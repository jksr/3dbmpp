#!/usr/bin/python

''' this script is to calculate the reg-ec scores
'''

import argparse
parser = argparse.ArgumentParser(description='compute register scores from sequence covariation info')
parser.add_argument('--weights', type=str, nargs=3, default=['110','030','080'])#
parser.add_argument('--logp', type=int, default=0, choices=[0,1,2])#
parser.add_argument('--theta', type=float, default=0)#

parser.add_argument('--ec_file', type=str)#
parser.add_argument('--strands_file', type=argparse.FileType('r'), required=True)#
parser.add_argument('--peris_file', type=argparse.FileType('r'))#
parser.add_argument('--output_file', type=argparse.FileType('w'), required=True)#
args = parser.parse_args()


lines = args.strands_file.readlines()
for line in lines:
	if len(line.strip())!=0:
		args.output_file.write('0\n')
