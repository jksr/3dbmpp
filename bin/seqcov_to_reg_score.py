#!/usr/bin/python

''' this script is to calculate the reg-ec scores
'''

import argparse
parser = argparse.ArgumentParser(description='compute register scores from sequence covariation info')
parser.add_argument('--weights', type=str, nargs=3, default=['110','030','080'])#
parser.add_argument('--logp', type=int, default=0, choices=[0,1,2])#
parser.add_argument('--theta', type=float, default=0)#

parser.add_argument('--ec_file', type=argparse.FileType('r'), required=True)#
parser.add_argument('--strands_file', type=argparse.FileType('r'), required=True)#
parser.add_argument('--peris_file', type=argparse.FileType('r'))#
parser.add_argument('--output_file', type=argparse.FileType('w'), required=True)#
args = parser.parse_args()


import glob
import sys
import math
import parse_ec
import regpy

def get_ecs(strandfn, ecfn, weights, theta, log, peris):
	## get ec dict
	ec_dict = parse_ec.get_neighbor_ec_dict(strandfn,ecfn)
	## construct circular barrel
	barrel = parse_ec.get_circular_barrel(strandfn)
	strand_ranges = parse_ec.get_strand_ranges(strandfn)
	barrel = []
	for i in range(len(strand_ranges)):
		if i%2==0:
			strand = list(reversed(range(strand_ranges[i][0], strand_ranges[i][1])))
		else:
			strand = range(strand_ranges[i][0], strand_ranges[i][1])
		barrel.append(strand)
	rtn= {}
	## scan all the strands
	for i in range(1,len(barrel)-1):
		rtn[i] = []
		offset1 = barrel[i][0]-peris[i-1]
		try:
			offset2 = barrel[i+1][0]-peris[i]
		except IndexError:
			offset2 = barrel[i+1][0]-peris[0]
		## enumerate registration
		for reg in range( -len(barrel[i])+1, len(barrel[i+1]) ):
			strand1, strand2 = regpy.construct_hairpin( barrel[i], barrel[i+1], reg )
			## adjustment of reg to interface
			if i%2==0:
				reg = reg+offset1+offset2
			else:
				reg = reg-offset1-offset2
			## get distant contact list
			dist_cont_list = []
			dist_cont_list.append( regpy.get_dist_contacts(strand1, strand2, 0) )
			dist_cont_list.append( regpy.get_dist_contacts(strand1, strand2, 1) )
			dist_cont_list.append( regpy.get_dist_contacts(strand1, strand2, 2) )
		
			score = 0
			contnum = 0
			for j in range(len(dist_cont_list)):
				conts = dist_cont_list[j]
				for cont in conts:
					## strand indexes are i-1, i instead of i, i+1
					try:
						prob = ec_dict[i-1,i%(len(barrel)-2)][cont]
						if prob < theta:
							continue
						if log == 1:
							score -= weights[j] * math.log(prob)
						else:
							score += weights[j] * prob
						contnum += 1
					except KeyError:
						pass

			rtn[i].append([reg,score,contnum])
	if log==2:
		minp = 1e-10
		for i in range(1,len(barrel)-1):
			tot = 0.0
			for dummy,score,dummy in rtn[i]:
				tot+=score
			if tot==0.0:
				continue
			for j in range(len(rtn[i])):
				if rtn[i][j][1]==0.0:
					rtn[i][j][1]=math.log(minp/tot)
				else:
					rtn[i][j][1]=math.log(rtn[i][j][1]/tot)
	return rtn


##################################################################
############################ MAIN ################################
##################################################################

#~~~~~~~~~~~~~~~~~~~ load inputs ~~~~~~~~~~~~~~~~~~~~
weights = [float(w) for w in args.weights]
dummy = sum(weights)
weights = [ w/dummy for w in weights ]
wtag = '_' + args.weights[0].zfill(3)+ args.weights[1].zfill(3)+ args.weights[2].zfill(3)

log = args.logp
if log==1:
	ltag = '_lp'
elif log==2:
	ltag = '_lnp'
else:
	ltag = '_p'

theta = args.theta
ttag = '_%1.2E'%theta


if args.peris_file is not None:
	lines = args.peris_file.readlines()
	peris = [int(line.strip()) for line in lines]
else:
	peris = []
	lines = args.strands_file.readlines()
	for i in range(len(lines)):
		if i%2==0:
			peris.append( int(lines[i].split()[0]) )
		else:
			peris.append( int(lines[i].split()[1]) )



fout = args.output_file

ecs = get_ecs(args.strands_file.name, args.ec_file.name, weights, theta, log, peris)
for i in sorted(ecs.keys()):
	fout.write(str(len(ecs[i])))
	for itm in ecs[i]:
		fout.write('\t'+str(itm[0])+'\t'+str(itm[1]))
	fout.write('\n')

fout.close()


