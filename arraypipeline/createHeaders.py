#!/usr/bin/env python
#
# Create specialized header files for VERA and SAM
# 
# Input
# 1. Group file: Column1: individual Column2: group

import sys

if (len (sys.argv) != 2):
    print 'error!  usage: createHeaders.py <Group file>\n'
    sys.exit ()
groupfile = sys.argv[1]

## Read in group file
lines = open(groupfile).read().split('\n')
lines = lines[0:-1]
nlines = len(lines)
sys.stderr.write ('Found %d lines in group file\n' % nlines )

# http://www.peterbe.com/plog/uniqifiers-benchmark
def unique(seq):
   # order preserving 
    checked = [] 
    for e in seq: 
        if e not in checked: 
            checked.append(e) 
    return checked 

## Get unique groups
grs = []
for line in lines:
    gr = line.split('\t')[1]
    grs.append(gr)
groups = unique(grs) ## Contains order, while dictionary will not!

## Populate group dictionary
groupdict = {}
for gr in groups:
    groupdict[gr] = []
for line in lines:
    toks = line.split('\t')
    gr = toks[1]
    indiv = toks[0]
    groupdict[gr].append(indiv)

### Create header files
for gr in groups:
    membs = groupdict[gr]
    fname = gr + '.headers'
    fp = open(fname,'w')
    fp.write('ProbesetID\tProbesetID\n')
    nm = len(membs)   
    fname2 = gr + '.dataColumnHeaders'
    fp2 = open(fname2,'w')
    fname3 = gr + '.XtoYheaders'
    fp3 = open(fname3,'w')
    for i in range(nm):
        mem = membs[i]
        ostring = mem + '\t' + 'X' + str(i) + '\n'
        fp.write(ostring)
        ostring2 = 'X' + str(i) + '\n'
        fp2.write(ostring2)
        ostring3 = 'X' + str(i) + '\tY' + str(i) + '\n'
        fp3.write(ostring3)
    fp.close()
    fp2.close()
    fp3.close()
    
