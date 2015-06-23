import sys, random, itertools, argparse, gzip
from Bio import SeqIO


"""
Script to remove reads that are < cutoff [bp] form PE fastq files.
The script expects trimmed reads
"""


# create pareser object
parser = argparse.ArgumentParser()

# add command line options
parser.add_argument('--f1', '--fastq_1', type=str, help='Read 1')
parser.add_argument('--f2', '--fastq_2', type=str, help='Read 2')
parser.add_argument('--c', '--cut', type=int, default=500, help='sequence length cutoff')
parser.add_argument('--p', '--path', type=str, help='output path')

# read the command line inputs
args = parser.parse_args()


#---------------------------------------------------------------------

# input
def make_out( infile, out_path ):
    inl = infile.split( '/' )
    if len( inl ) == 1:
        inll = inl[0]
    else:
        inll = inl[-1]
    out = out_path+'/'+inll
    return out

infile1 = args.f1

# open files depending on gzip or not
if infile1.split('.')[-1] == 'gz':
    in1 = SeqIO.parse( gzip.open(args.f1, 'r'), 'fastq')
    in2 = SeqIO.parse( gzip.open(args.f2, 'r'), 'fastq')
    out1 = gzip.open(make_out( args.f1, args.p), 'wb', compresslevel = 6)
    out2 = gzip.open(make_out( args.f2, args.p), 'wb', compresslevel = 6)
elif infile1.split('.')[-1] == 'fastq':
    print 'fastq'
    in1 = SeqIO.parse( args.f1, 'fastq')
    in2 = SeqIO.parse( args.f2, 'fastq')
    out1 = open(make_out( args.f1, args.p), 'w')
    out2 = open(make_out( args.f2, args.p), 'w')

#---------------------------------------------------------------------

cutoff = args.c

counter_R = 0
counter_M = 0
counter = 0
for read1, read2 in itertools.izip( in1, in2 ):
    counter += 1
    if (len(read1.seq) < cutoff) or (len(read2.seq) < cutoff):
        counter_R += 1
        continue
    r1_name = read1.description.split( ' ' )[0]
    r2_name = read2.description.split( ' ' )[0]
    if r1_name != r2_name:
        counter_M += 1
        continue
    out1.write( read1.format("fastq") )
    out2.write( read2.format("fastq") )
print('%i Reads of totally %i Reads are shorter than %i bp' % (counter_R, counter, cutoff ))
print('%i Reads are NOT properly paired in %s and %s' % (counter_M, args.f1, args.f2 ))
out1.close()
out2.close()
