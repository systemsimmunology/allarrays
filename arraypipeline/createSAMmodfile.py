#!/usr/bin/env python
#
# creatSAMmodfile.py
#
# From two VERA_OC model files, create two-channel model file
# Correlation is set to zero.
import sys
import re

if (len (sys.argv) != 3):
  print 'error!  usage:  createSAMmodfile.py <VERA_OC modfile 1> <VERA_OC modfile 2>\n'
  sys.exit ()
 
vocfile1 = sys.argv[1]
vocfile2 = sys.argv[2]

# Read in VERA_OC model files
lines = open(vocfile1).read().split('\n')
vals = lines[1]
tokens  = vals.split()
sig_eps_x = tokens[0]
sig_del_x = tokens[1]

lines = open(vocfile2).read().split('\n')
vals = lines[1]
tokens  = vals.split()
sig_eps_y = tokens[0]
sig_del_y = tokens[1]

rho_eps = str(0.)

header = 'sig_eps_x sig_eps_y  rho_eps  sig_del_x  sig_del_y'
print header
print sig_eps_x+' '+sig_eps_y+' '+rho_eps+' '+ sig_del_x+' '+ sig_del_y
