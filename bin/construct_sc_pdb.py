#!/usr/bin/python
import os
import argparse
parser = argparse.ArgumentParser(description='construct mainchan pdb using scwrl')
parser.add_argument('--bbpdbfn', type=str, required=True)
parser.add_argument('--scpdbfn', type=str, required=True)
parser.add_argument('--scwrlpath', type=str, required=True)
args = parser.parse_args()

os.system(args.scwrlpath+' -h -i '+args.bbpdbfn+' -o '+args.scpdbfn+' 1>/dev/null' )
