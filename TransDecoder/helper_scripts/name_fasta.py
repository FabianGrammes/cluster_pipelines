from Bio import SeqIO
import argparse

"""
Short script to append fasta ids by a unique species id.
Script removes all fasta descriptions. 
"""
#=====================================================================

# create pareser object
parser = argparse.ArgumentParser()

# add command line options
parser.add_argument('-i', '-infile', type=str, help='input .fa file')
parser.add_argument('-o', '-outfile', type=str, default='out.gff', help='output .fa file')
parser.add_argument('-p', '-prefix', type=str, default='empty', help='fasta id prefix')

# read the command line inputs
args = parser.parse_args()

#=====================================================================

def append_ids( in_fa, out_fa, prefix ):
    clean_fa = []
    # parse fasta
    for record in SeqIO.parse( in_fa, "fasta" ):
        record.id = prefix+"|"+record.id
        record.description = ''
        clean_fa.append( record )
    # write out
    output_handle = open( out_fa, "w")
    SeqIO.write( clean_fa, output_handle, "fasta")
    output_handle.close()


# call the function
append_ids( args.i, args.o, args.p )


