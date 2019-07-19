#!/usr/bin/python
import glob
import os
import numpy as np
import sys
import argparse
parser = argparse.ArgumentParser(description='construct mainchan pdb using bbq algorithm')
parser.add_argument('--capdbfn', type=str, required=True)
parser.add_argument('--bbpdbfn', type=str, required=True)
args = parser.parse_args()

with open(args.capdbfn) as f:
	lines = f.readlines()
for i in range(len(lines)):
	lines[i] = lines[i].replace(' CA   ', '  CA  ')
prebbqfn = args.capdbfn+'.prebbq.pdb'
with open(prebbqfn, 'w') as f:
	f.writelines(lines)

bbqdir = os.path.dirname(os.path.realpath(__file__))+'/BBQ'
os.system( '/usr/bin/java -classpath '+bbqdir+':'+bbqdir+'/jbcl.jar BBQ -d='+bbqdir+'/q_50_xyz.dat -r='+prebbqfn+' 1> '+args.bbpdbfn+' 2> /dev/null' )
