import sys, optparse, itertools, warnings, traceback, os.path, collections

import HTSeq

# basic idea:
## - use the code from HTSeq count
## - modify that ONLY multimapping reads are counted
## Basic question afterwards is:
##     Should there be any way to preserve read information : Read pair A mapps to all the following loactions,
##     this could potentially lead to a hughe file

##     To ask a more specific question like: Are there many reads that are multimapping
##     for homoeolog regions might in practice be easier...

## Maybe just restrict the whole thing to reads mapping to exactly 2 genomic positions


#==============================================================================
# define input files
gff_filename = 'Salmon_3p6_Chr_051214.gtf'
sam_filename = 's10_test2.sam'
stranded = 'yes'
feature_type = 'exon'
id_attribute = 'gene_id'
quiet = False

# Function to referse strand 
def invert_strand( iv ):
   iv2 = iv.copy()
   if iv2.strand == "+":
      iv2.strand = "-"
   elif iv2.strand == "-":
      iv2.strand = "+"
   else:
      raise ValueError, "Illegal strand"
   return iv2

# Define the FUN to read the gff/gtf file
def gff_read( gff_filename, stranded ):
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
            if i % 100000 == 0 and not quiet:
                sys.stderr.write( "%d GFF lines processed.\n" % i )
    except KeyError:
        raise ValueError, ( "Something's fishy with the file")
    return( features )

# READ the GTF/GFF file
features = gff_read( gff_filename, stranded )

# quick check if names are there
for iv, val in features[ HTSeq.GenomicInterval( "ssa01", 21874065,  21874200, '+')].steps():
    print iv, val

for iv, val in features[ invert_strand(HTSeq.GenomicInterval( "ssa01", 21874065,  21874200, '+'))].steps():
    print iv, val

# INPROGRESS:
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
    count_double = collections.Counter()
    count_single = collections.Counter()
    for bundle in HTSeq.pair_SAM_alignments( almnt_file, bundle=True ):
        rs = set()
        if len(bundle) == 2:
            for r1,r2 in bundle:
                if r1 is None or r2 is None:
                    count_double[ '__not_aligned' ] += 1
                    continue
                else:
                    try:
                        iv_seq1 = ( co.ref_iv for co in r1.cigar if co.type == "M" and co.size > 0 )
                        iv_seq2 = ( co.ref_iv for co in r2.cigar if co.type == "M" and co.size > 0 )
                    except AttributeError:
                        raise ValueError, ( "Someting wrong with read %s" % (r1))
                        continue
                    for iv in iv_seq1:
                        for iv2, fs2 in features[ invert_strand(iv) ].steps():
                            rs = rs.union( fs2 )
                    for iv in iv_seq2:
                        for iv2, fs2 in features[ iv ].steps():
                            rs = rs.union( fs2 )
            if len(rs) == 0:
                count_double[ '__no_feature' ] += 1
            elif len(rs) == 1:
                count_single[ '_'.join(rs) ] += 1
            elif len(rs) == 2:
                count_double[ '-'.join(rs) ] += 1
            else:
                count_double[ '__too_many_features' ] += 1
        elif len(bundle) == 1:
            count_double[ '__single_hit' ] += 1
        elif len(bundle) > 2:
            count_double[ '__ambigous' ] += 1
    # join the collections
    com_coll = count_single+count_double
    # this sorts the collections.counter 
    com_coll = sorted(com_coll.items(), key=lambda pair: pair[0], reverse=False)
    return( com_coll )

cd = run_through_sam( 'test3.sam' )
for i in cd:
    print i
