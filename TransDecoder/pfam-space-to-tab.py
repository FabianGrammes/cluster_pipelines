#!/usr/bin/env python

# (c) CIGENE, Torfinn Nome, 2014-01

# Version 0.1: Initial version.

''' 
Input: pfam dat file
Output: pfam dat file with spaces replaced with tabulators, except in last column.
'''

import sys

with open(sys.argv[1]) as file:
    for lineIn in file:
        lineOut = lineIn.strip().split()
        print "\t".join(lineOut[0:22]) + "\t" + " ".join(lineOut[22:])
        