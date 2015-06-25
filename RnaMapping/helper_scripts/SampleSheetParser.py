#!/usr/bin/python
import argparse, sys

# create pareser object
parser = argparse.ArgumentParser()

# add command line options
parser.add_argument('-s', '-sheet', type=str, help='input sample sheet')
parser.add_argument('-o', '-output', type=str, help='output file')

# read the command line inputs
args = parser.parse_args()


'''
Function to parse the relevant information from an Illumina MiSeq sample sheet (.csv)
into a more usable file

'''



#-------------------------------------------------------------------------------
# Functions

def filter_data( infile ):   
    read_file = open( infile, 'r' )
    keeper = False
    sample_data = []
    for line in read_file:
        if line.startswith('[Data]'):
            keeper = True
            c = -2
        if keeper:
            c += 1
            if c > 0:
                m_entry = [i for i in line.split(',') ]
                m_entry.append(str(c))
                sample_data.append( m_entry )
    return sample_data

data = filter_data( args.s )

# Write output 
out = open( args.o, 'w')
for i in data:
    if i[0] != '' and i[1] == '':
        out.write('%s\t%s\t%s\t%s\n' % ( i[0], i[0].replace('_','-')+'_S'+i[8]+'_L001' , i[4], i[5] ))
    elif i[0] == '' and i[1] != '':
        out.write('%s\t%s\t%s\t%s\n' % ( i[1], i[1].replace('_','-')+'_S'+i[8]+'_L001' , i[4], i[5] ))
    elif i[0] != '' and i[1] != '':
        out.write('%s\t%s\t%s\t%s\n' % ( i[1], i[1].replace('_','-')+'_S'+i[8]+'_L001' , i[4], i[5] ))
        
out.close()

