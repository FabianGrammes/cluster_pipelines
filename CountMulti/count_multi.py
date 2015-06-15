import sys, itertools, warnings, traceback, os.path, collections, HTSeq
from optparse import OptionParser

# basic idea:
## - use the code from HTSeq count
## - modify that ONLY multimapping reads are counted
## Basic question afterwards is:
##     Should there be any way to preserve read information : Read pair A mapps to all the following loactions,
##     this could potentially lead to a hughe file

##     To ask a more specific question like: Are there many reads that are multimapping
##     for homoeolog regions might in practice be easier...

## Maybe just restrict the whole thing to reads mapping to exactly 2 genomic positions

parser = OptionParser()
parser.add_option("-g", "--gff", dest="gff_filename", type="string",
                  help="Input GTF/GFF file", metavar="GTF/GFF FILE")

parser.add_option("-o", "--out", dest="out_filename", type="string",
                  default = "multi_count.count", metavar = "OUTPUT FILE",
                  help="Name of the output file")

parser.add_option( "-s", "--stranded", type="choice", dest="stranded",
                  choices = ( "yes", "no", "reverse" ), default = "yes",
                  help = "whether the data is from a strand-specific assay." +
                  "'yes' is default and currently the only possible option.")

parser.add_option( "-t", "--type", type="string", dest="featuretype",
                  default = "exon",
                  help = "feature type (3rd column in GFF file) to be used, " +
                  "all features of other type are ignored (default, suitable for Ensembl " +
                  "GTF files: exon). Currently 'exon' is the only option" )

parser.add_option( "-a", "--idattr", type="string", dest="idattr",
                  default = "gene_id",
                  help = "GFF attribute to be used as feature ID (default, " +
                  "suitable for Ensembl GTF files: gene_id). Currently 'gene_id' is the only option" )

parser.add_option( "-q", "--quiet", action="store_true", dest="quiet",
                  help = "suppress progress report" ) # and warnings" )

(options, args) = parser.parse_args()

if len( args ) != 1:
    sys.stderr.write( sys.argv[0] + ": Error: Please provide the input .sam.\n" )
    sys.stderr.write( "  Call with '-h' to get usage information.\n" )
    sys.exit( 1 )

quiet = options.quiet

#-------------------------------------------------------------------------------
# Define the FUN to read the gff/gtf file
def gff_read( gff_filename, stranded, id_attribute, feature_type ):
    features = HTSeq.GenomicArrayOfSets( "auto", stranded )
    # read in the gff/gtf file 
    gff = HTSeq.GFF_Reader( gff_filename )   
    i = 0
    try:
        for f in gff:
            if f.type == feature_type:
                try:
                    feature_id = f.attr[ id_attribute ]
                except KeyError:
                    raise ValueError, ( "Feature %s does not contain a '%s' attribute" % ( f.name, id_attribute ) )
                if stranded != "no" and f.iv.strand == ".":
                    raise ValueError, ( "Feature %s at %s does not have strand information but you are " "running htseq-count in stranded mode. Use '--stranded=no'." % ( f.name, f.iv ) )
                features[ f.iv ] += feature_id
            i += 1
            if i % 200000 == 0 and not quiet:
                sys.stderr.write( "%d GFF lines processed.\n" % i )
    except KeyError:
        raise ValueError, ( "Something's fishy with the file")
    return( features )

# READ the GTF/GFF file
features = gff_read( options.gff_filename, options.stranded, options.idattr, options.featuretype )

#-------------------------------------------------------------------------------
# Function to referse strand used in the run_through_sam function
def invert_strand( iv ):
   iv2 = iv.copy()
   if iv2.strand == "+":
      iv2.strand = "-"
   elif iv2.strand == "-":
      iv2.strand = "+"
   else:
      raise ValueError, "Illegal strand"
   return iv2

'''
Function to read a NAME SORTED SAM file line by line.
For reads that map to 2 places in the genome it reports these places, more
percise it reports the gene ids of of the mapping hits.
So far the function considers strandedness, it expects:
- first read = reverse
- second read = forward
Which is the standard for stranded illumina reads...
'''
def run_through_sam( sam_filename ): 
    try:
        almnt_file = HTSeq.SAM_Reader( sam_filename )
    except KeyError:
        raise ValueError, ( "Can't find file %s" % (sam_filename))
    count_reads = collections.Counter()
    i = 0
    for bundle in HTSeq.pair_SAM_alignments( almnt_file, bundle=True ):
        if len(bundle) != 0:
            i += 1
            if i > 0 and i % 200000 == 0 and not quiet:
                sys.stderr.write( "%d SAM alignment records processed.\n" % ( i ) )
            rs = set()
            # Loop for multimapping reads: Reads that map to more than 3 positions
            if len(bundle) > 2:
                 count_reads[ '__Ambigios_read' ] += 1
                 continue
            # Loop for Singles: Reads that map to 1 genomic postion 
            elif len(bundle) == 1:
                for r1,r2 in bundle:
                    if r1 is None or r2 is None:
                        count_reads[ '__Single_hit:Not_aligned' ] += 1
                        continue
                    else:
                        try:
                            iv_seq1 = ( co.ref_iv for co in r1.cigar if co.type == "M" and co.size > 0 )
                            iv_seq2 = ( co.ref_iv for co in r2.cigar if co.type == "M" and co.size > 0 )
                        except AttributeError:
                            raise ValueError, ( "Single:Someting wrong with read %s" % (r1))
                            continue
                        for iv in iv_seq1:
                            for iv2, fs2 in features[ invert_strand(iv) ].steps():
                                rs = rs.union( fs2 )
                        for iv in iv_seq2:
                            for iv2, fs2 in features[ iv ].steps():
                                rs = rs.union( fs2 )
                                # Parsing through the set  
                if len(rs) == 0:
                    count_reads[ '__Single_hit:No_feature' ] += 1
                elif len(rs) == 1:
                    count_reads[ '__Single_hit:Feature_found' ] += 1
                elif len(rs) > 1:
                    count_reads[ '__Single_hit:Ambigous_features' ] += 1
            # Loop for Doubles: Reads that map to 2 genomic postion 
            elif len(bundle) == 2:
                found = []
                for r1,r2 in bundle:
                    if r1 is None or r2 is None:
                        found.append(False)
                        continue
                    else:
                        found.append(True)  
                        try:
                            iv_seq1 = ( co.ref_iv for co in r1.cigar if co.type == "M" and co.size > 0 )
                            iv_seq2 = ( co.ref_iv for co in r2.cigar if co.type == "M" and co.size > 0 )
                        except AttributeError:
                            raise ValueError, ( "Double:Someting wrong with read %s" % (r1))
                            continue
                        for iv in iv_seq1:
                            for iv2, fs2 in features[ invert_strand(iv) ].steps():
                                rs = rs.union( fs2 )
                        for iv in iv_seq2:
                            for iv2, fs2 in features[ iv ].steps():
                                rs = rs.union( fs2 )
                if all(found) == False:
                    count_reads[ '__Double_hit:Not_aligned' ] += 1
                    continue
                if any(found):
                    if len(rs) == 0:
                        count_reads[ '__Double_hit:No_feature' ] += 1
                    elif len(rs) == 1:
                        count_reads[ '_'.join(rs) ] += 1
                        count_reads[ '__Double_hit:Single_feature_found' ] += 1
                    elif len(rs) == 2:
                        count_reads[ '_'.join(rs) ] += 1
                        count_reads[ '__Double_hit:Double_feature_found' ] += 1
                    elif len(rs) > 2:
                        count_reads[ '__Double_hit:Ambigous_features' ] += 1
        else:
            continue
    # this sorts the collections.counter
    count_reads['__Total_reads' ] = i
    com_coll = sorted(count_reads.items(), key=lambda pair: pair[0], reverse=False)
    return( com_coll )

multi_counts = run_through_sam( args[0] )

#-------------------------------------------------------------------------------
# Write the results to file

handle = open(options.out_filename, 'w')
for line in multi_counts:
    handle.write("%s\t%d\n" % (line[0], line[1]))

handle.close()
