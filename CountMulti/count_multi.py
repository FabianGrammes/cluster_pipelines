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


for i in algn_counter:
    print i, algn_counter[i]

# this sorts the collections.counter 
a = sorted(count_double.items(), key=lambda pair: pair[0], reverse=False)
for i in a:
    print i
        
    
cig_out = []
if (type(rs) == str):
    out = rs
else:
    for f in rs:
        if ':' in f.name:
            cig_out.append( f.name)
    out = '-'.join(sorted(cig_out))
print out


    print rs, len(rs)
    print '---\n'


for iv in iv_seq2:
    for iv2, fs2 in features[ iv ].steps():
        print fs2

        
        rs = rs.union( fs2 )
print rs

            
            
            
            iv_seq1 = itertools.chain( iv_seq1, ( co1.ref_iv for co1 in r1.cigar if co1.type == "M" and co1.size > 0 ) )



        else:
            [ co1.ref_iv for co1 in r1.cigar if co1.type == "M" and co1.size > 0 ] # CHECK
            [ co2.ref_iv for co2 in r2.cigar if co2.type == "M" and co2.size > 0 ]
            rs = set()
            for iv2, fs2 in gff[ co2.ref_iv ].steps():
                rs = rs.union( fs2 )
                rs
                 
                rs = rs.union(s)
            
    #if report_bundle == True:  # unique paired hits eventually may be that one read is aligned but the other not ??
    #    # target bundle length
    #   unique += 1
    #   continue
    for r1,r2 in bundle:
        r1
        r2
        [ co1.ref_iv for co1 in r1.cigar if co1.type == "M" and co1.size > 0 ] # CHECK
        [ co2.ref_iv for co2 in r2.cigar if co2.type == "M" and co2.size > 0 ] # CHECK

        elif len(bundle) > 2:  # hitting too many positions
            too-ambiguous += 1
            continue
        else:
            empty += 1
            


        else 
        a.append( len(bundle) )


            
#==============================================================================
# Chunk to match read positions to the Genome...
# has to refined, maybe make less fancy...
     try:
            if overlap_mode == "union":
               fs = set()
               for iv in iv_seq:
                  if iv.chrom not in features.chrom_vectors:
                     raise UnknownChrom
                  for iv2, fs2 in features[ iv ].steps():
                     fs = fs.union( fs2 )
            elif overlap_mode == "intersection-strict" or overlap_mode == "intersection-nonempty":
               fs = None
               for iv in iv_seq:
                  if iv.chrom not in features.chrom_vectors:
                     raise UnknownChrom
                  for iv2, fs2 in features[ iv ].steps():
                     if len(fs2) > 0 or overlap_mode == "intersection-strict":
                        if fs is None:
                           fs = fs2.copy()
                        else:
                           fs = fs.intersection( fs2 )
            else:
               sys.exit( "Illegal overlap mode." )
            if fs is None or len( fs ) == 0:
               write_to_samout( r, "__no_feature" )
               empty += 1
            elif len( fs ) > 1:
               write_to_samout( r, "__ambiguous[" + '+'.join( fs ) + "]" )
               ambiguous += 1
            else:


#------------------------------------------------------------------------------
almnt_file = HTSeq.SAM_Reader( "s10_test2.sam" )
#counts = collections.Counter( )
c_mgood = 0
c_mbad = 0
a = []

for bundle in HTSeq.pair_SAM_alignments( almnt_file, bundle=True ):
    if (len(bundle) > 1 and len(bundle) < 3):
        a.append( len(bundle) )
    
    
    if len(bundle) > 1 & len(bundle) < 3:
        c_mgood += 1
    elif len(bundle) > 3:
        c_mbad +=1

print "OK multimappers %d; Too heavy multimappers %d" % (c_mgood, c_mbad) 


      continue  # Skip multiple alignments
   first_almnt, second_almnt = bundle[0]  # extract pair
   if not first_almnt.aligned and second_almnt.aligned:
      count[ "_unmapped" ] += 1
      continue
   gene_ids = set()
   for iv, val in features[ left_almnt.iv ].steps():
      gene_ids |= val
   for iv, val in features[ right_almnt.iv ].steps():
      gene_ids |= val
   if len(gene_ids) == 1:
      gene_id = list(gene_ids)[0]
      counts[ gene_id ] += 1
   elif len(gene_ids) == 0:
      counts[ "_no_feature" ] += 1
   else:
      counts[ "_ambiguous" ] += 1




#------------------------------------------------------------------------------







read_seq = HTSeq.SAM_Reader( 's10_test.sam' )
read_seq2 = HTSeq.pair_SAM_alignments( read_seq )
for r in read_seq:
    # identify non unique reads by the NH field
    if r.optional_field( "NH" ) > 1:
        print r

try:
    if pe_mode:
        if order == "name":
            read_seq = HTSeq.pair_SAM_alignments( read_seq )
        elif order == "pos":
            read_seq = HTSeq.pair_SAM_alignments_with_buffer( read_seq )
        else:
            raise ValueError, "Illegal order specified."
    empty = 0
    ambiguous = 0
    notaligned = 0
    lowqual = 0
    nonunique = 0
    i = 0   
    for r in read_seq:
        if i > 0 and i % 100000 == 0 and not quiet:
            sys.stderr.write( "%d SAM alignment record%s processed.\n" % ( i, "s" if not pe_mode else " pairs" ) )

        i += 1
        if not pe_mode:
            if not r.aligned:
                notaligned += 1
                write_to_samout( r, "__not_aligned" )
                continue
            try:
                if r.optional_field( "NH" ) > 1:
                    nonunique += 1
                    write_to_samout( r, "__alignment_not_unique" )
                    continue
            except KeyError:
                pass
        if r.aQual < minaqual:
            lowqual += 1
            write_to_samout( r, "__too_low_aQual" )
            continue
        if stranded != "reverse":
            iv_seq = ( co.ref_iv for co in r.cigar if co.type == "M" and co.size > 0 )
        else:
            iv_seq = ( invert_strand( co.ref_iv ) for co in r.cigar if co.type == "M" and co.size > 0 )            
        else:
            if r[0] is not None and r[0].aligned:
                if stranded != "reverse":
                    iv_seq = ( co.ref_iv for co in r[0].cigar if co.type == "M" and co.size > 0 )
                else:
                    iv_seq = ( invert_strand( co.ref_iv ) for co in r[0].cigar if co.type == "M" and co.size > 0 )
            else:
                iv_seq = tuple()
            if r[1] is not None and r[1].aligned:            
                if stranded != "reverse":
                    iv_seq = itertools.chain( iv_seq,
                                            ( invert_strand( co.ref_iv ) for co in r[1].cigar if co.type == "M" and co.size > 0 ) )
               else:
                  iv_seq = itertools.chain( iv_seq, 
                     ( co.ref_iv for co in r[1].cigar if co.type == "M" and co.size > 0 ) )
            else:
               if ( r[0] is None ) or not ( r[0].aligned ):
                  write_to_samout( r, "__not_aligned" )
                  notaligned += 1
                  continue         
            try:
               if ( r[0] is not None and r[0].optional_field( "NH" ) > 1 ) or \
                     ( r[1] is not None and r[1].optional_field( "NH" ) > 1 ):
                  nonunique += 1
                  write_to_samout( r, "__alignment_not_unique" )
                  continue
            except KeyError:
               pass
            if ( r[0] and r[0].aQual < minaqual ) or ( r[1] and r[1].aQual < minaqual ):
               lowqual += 1
               write_to_samout( r, "__too_low_aQual" )
               continue         
         
         try:
            if overlap_mode == "union":
               fs = set()
               for iv in iv_seq:
                  if iv.chrom not in features.chrom_vectors:
                     raise UnknownChrom
                  for iv2, fs2 in features[ iv ].steps():
                     fs = fs.union( fs2 )
            elif overlap_mode == "intersection-strict" or overlap_mode == "intersection-nonempty":
               fs = None
               for iv in iv_seq:
                  if iv.chrom not in features.chrom_vectors:
                     raise UnknownChrom
                  for iv2, fs2 in features[ iv ].steps():
                     if len(fs2) > 0 or overlap_mode == "intersection-strict":
                        if fs is None:
                           fs = fs2.copy()
                        else:
                           fs = fs.intersection( fs2 )
            else:
               sys.exit( "Illegal overlap mode." )
            if fs is None or len( fs ) == 0:
               write_to_samout( r, "__no_feature" )
               empty += 1
            elif len( fs ) > 1:
               write_to_samout( r, "__ambiguous[" + '+'.join( fs ) + "]" )
               ambiguous += 1
            else:
               write_to_samout( r, list(fs)[0] )
               counts[ list(fs)[0] ] += 1
         except UnknownChrom:
            write_to_samout( r, "__no_feature" )
            empty += 1

   except:
      sys.stderr.write( "Error occured when processing SAM input (%s):\n" % read_seq_file.get_line_number_string() )
      raise

   if not quiet:
      sys.stderr.write( "%d SAM %s processed.\n" % ( i, "alignments " if not pe_mode else "alignment pairs" ) )
         
   if samoutfile is not None:
      samoutfile.close()

   for fn in sorted( counts.keys() ):
      print "%s\t%d" % ( fn, counts[fn] )
   print "__no_feature\t%d" % empty
   print "__ambiguous\t%d" % ambiguous
   print "__too_low_aQual\t%d" % lowqual
   print "__not_aligned\t%d" % notaligned
   print "__alignment_not_unique\t%d" % nonunique
