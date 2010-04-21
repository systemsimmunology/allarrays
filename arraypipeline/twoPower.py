#!/usr/bin/env python
#
# twoPower.py 
# Evaluate 2 to the power of data values
# ( Inverse of the log base 2 function )
#
# Assumes tab-separation
# Data columns are to be specified in dataColumnHeader file
# One line per data column. Needs exact match, including spaces.
# Very few format checks.
 
import sys
import re
 
if (len (sys.argv) != 3):
  print 'error!  usage:  twoPower.py <input file> <dataColumnHeaders>\n'
  sys.exit ()

inputFile = sys.argv[1]
dataHeadingFile = sys.argv[2]

## Read in data file
lines = open(inputFile).read().split('\n')
lines = lines[0:-1]
nlines = len(lines)
sys.stderr.write ('Found %d lines\n' % nlines )

## Read in data header files
dataHeaders = open(dataHeadingFile).read().split('\n')
dataHeaders = dataHeaders[0:-1]
ndataHeaders = len(dataHeaders)
sys.stderr.write ('Found %d dataHeaders\n' % ndataHeaders )

## Determine indices for which columns are to be taken to power two
fullHeader = lines[0]
fields = fullHeader.split('\t')
colsToTransform = []
for header in dataHeaders:
    ind=fields.index(header)
    colsToTransform.append(ind)

## Print out data lines
print lines[0]
for line in lines[1:]:
    fields = line.split('\t')
    nfields = len(fields)
    outString = ''
    for i in range(nfields):
      if ( i in colsToTransform ):
        outString += str(pow(2,float(fields[i])))
      else:
        outString += fields[i]
      if ( not(i== (nfields-1)) ):
        outString += '\t'
      else:
        outString += '\n'
    print outString,
