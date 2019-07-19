#!/usr/bin/python

import sys
import os
import numpy as np
import glob
import argparse
parser = argparse.ArgumentParser(description='Adjust local register prediction using global information')
parser.add_argument('--opt', type=int, default=-1, choices=[0,1,2,3,4,5,6,7], required=True)
parser.add_argument('--folder', type=str, required=True)
parser.add_argument('--oddsdir', type=str, required=True)
args = parser.parse_args()


oddsei = np.loadtxt(glob.glob(args.oddsdir+'/ExtraIn.odds' )[0])
oddseo = np.loadtxt(glob.glob(args.oddsdir+'/ExtraOut.odds')[0])
oddspi = np.loadtxt(glob.glob(args.oddsdir+'/PeriIn.odds'  )[0])
oddspo = np.loadtxt(glob.glob(args.oddsdir+'/PeriOut.odds' )[0])
oddsci = np.loadtxt(glob.glob(args.oddsdir+'/CoreIn.odds'  )[0])
oddsco = np.loadtxt(glob.glob(args.oddsdir+'/CoreOut.odds' )[0])


def determine_facings_odds(periids, extraids, res):
	predicted_peri_facings = []
	for i in range(len(periids)):
		if periids[i]<extraids[i]:
			strand = range(periids[i], extraids[i]+1, 1)
		else:
			strand = range(periids[i], extraids[i]-1, -1)

		for j in range(len(strand)):
			strand[j] = res[ strand[j]-1 ]

		bestscore = -10000000000000
		for peri_num,core_num,extra_num in [ (2,5,2), (2,5,3), (2,6,2), (3,5,2), (2,6,3), (3,6,2), (3,5,3), (3,6,3), (2,5,1), (2,4,2), (1,5,2), (2,4,1), (1,5,1), (1,4,2) ]:
			for peri_facing_out in [0,1]:
				score = 0.0
				for j in range(peri_num):
					try:
						strand[j]
					except IndexError:
						continue
					if j%2==peri_facing_out:
						score += oddspi[strand[j]]
					else:
						score += oddspo[strand[j]]
				for j in range(peri_num, peri_num+core_num):
					try:
						strand[j]
					except IndexError:
						continue
					if j%2==peri_facing_out:
						score += oddsci[strand[j]]
					else:
						score += oddsco[strand[j]]
				for j in range(peri_num+core_num, peri_num+core_num+extra_num):
					try:
						strand[j]
					except IndexError:
						continue
					if j%2==peri_facing_out:
						score += oddsei[strand[j]]
					else:
						score += oddseo[strand[j]]
				if score > bestscore:
					bestscore = score
					best_peri_facing_out = peri_facing_out
		if best_peri_facing_out:
			predicted_peri_facings.append('OUT')
		else:
			predicted_peri_facings.append('IN')
	return predicted_peri_facings


getfep_score_dict = { 0: 0, 1:-1, 2:-1, 3:-1, 4:-1, 5:-1, 6:-1, 7:-1, 8:-1, 9: 0, 10: 0, 11:-1, 12: 0, 13: 0, 14:-1, 15:-1, 16:-1, 17: 0, 18: 0, 19: 0 }
def determine_facings_GeTFEP(periids, extraids,res):
	predicted_peri_facings=[]

	for i in range(len(periids)):
		if periids[i]<extraids[i]:
			strand = range(periids[i], extraids[i]+1, 1)
		else:
			strand = range(periids[i], extraids[i]-1, -1)

		for j in range(len(strand)):
			strand[j] = res[ strand[j]-1 ]

		periout_score = 0
		for j in [0,2,4]:
			try:
				periout_score += getfep_score_dict[ strand[j] ]
			except:
				pass
		for j in [1,3,5]:
			try:
				periout_score -= getfep_score_dict[ strand[j] ]
			except:
				pass

		if periout_score>=0:
			predicted_peri_facings.append('OUT')
		else:
			predicted_peri_facings.append('IN')

	return predicted_peri_facings



def get_possible_targetshears(pdbdata):
	N = len(pdbdata)
	if N <= 8:
		targetshears = [N,N+2]
	elif N <= 10:
		targetshears = [N+2]
	elif N <= 12:
		targetshears = [N+2,N+4]
	elif N <= 14:
		targetshears = [N,N+2]
	elif N <= 16:
		targetshears = [N+4]
	elif N <= 18:
		targetshears = [N+2,N+4]
	elif N <= 24:
		targetshears = [N+2]
	elif N <= 26:
		targetshears = [N+4]

def get_common_targetshears(pdbdata):
	N = len(pdbdata)
	if N <= 12:
		targetshears = [N+2]
	elif N <= 14:
		targetshears = [N]
	elif N <= 18:
		targetshears = [N+4]
	elif N <= 24:
		targetshears = [N+2]
	elif N <= 26:
		targetshears = [N+2]
	#	targetshears = [N+4]
	else:
		pass #TODO
	return targetshears



def parity_first_shear_adjustment(pdbdata, targetshears):
	N = len(pdbdata)
	currshear = sum(pdbdata[:,3].astype(int))

	newdata = []
	for i in range(N):
		score_diff, pred0_score, pred1_score, pred0_reg, pred1_reg, predface, strandid = pdbdata[i]
		score_diff = float(score_diff)
		pred0_reg = int(pred0_reg)
		pred1_reg = int(pred1_reg)
		regmod01 = (pred0_reg-pred1_reg)%2

		if currshear < targetshears[0]:
			if pred0_reg < pred1_reg:
				## adjust shear parity first unless top 2 regs can make the shear larger
				if i<2 or (currshear%2==regmod01 and targetshears[-1]-currshear>=pred1_reg-pred0_reg):
				# #adjust shear parity first
				#if (currshear%2==regmod01 and targetshears[-1]-currshear>=pred1_reg-pred0_reg):
					currshear += pred1_reg-pred0_reg
					pred0_reg = pred1_reg
					pred0_score = pred1_score
					score_diff = 0
		elif currshear > targetshears[-1]:
			if pred0_reg > pred1_reg:
				# #adjust shear parity first unless top 2 regs can make the shear smaller
				if i<2 or (currshear%2==regmod01 and targetshears[0]-currshear<=pred1_reg-pred0_reg):
				# #adjust shear parity first
				#if (currshear%2==regmod01 and targetshears[0]-currshear<=pred1_reg-pred0_reg):
					currshear += pred1_reg-pred0_reg
					pred0_reg = pred1_reg
					pred0_score = pred1_score
					score_diff = 0

		newdata.append( (score_diff, pred0_score, pred1_score, pred0_reg, pred1_reg, predface, strandid) )

	newdata = np.array(sorted(newdata))
	return newdata



def brute_shear_adjustment(alldata, targetshears):
	N = len(alldata)
	currshear = sum(alldata[:,3].astype(int))

	newdata = []
	for i in range(N):
		score_diff, pred0_score, pred1_score, pred0_reg, pred1_reg, predface, strandid = alldata[i]
		score_diff = float(score_diff)
		pred0_reg = int(pred0_reg)
		pred1_reg = int(pred1_reg)
		peri = int(peri)

		if currshear < targetshears[0]:
			if pred0_reg < pred1_reg:
				if i<2 or targetshears[-1]-currshear>=pred1_reg-pred0_reg:
					currshear += pred1_reg-pred0_reg
					pred0_reg = pred1_reg
					pred0_score = pred1_score
					score_diff = 0
		elif currshear > targetshears[-1]:
			if pred0_reg > pred1_reg:
				if i<2 or targetshears[0]-currshear<=pred1_reg-pred0_reg:
					currshear += pred1_reg-pred0_reg
					pred0_reg = pred1_reg
					pred0_score = pred1_score
					score_diff = 0

		newdata.append( (score_diff, pred0_score, pred1_score, pred0_reg, pred1_reg, predface, strandid) )

	newdata = np.array(sorted(newdata))
	return newdata


def get_final(alldata):
	N = len(alldata)
	newdata = []
	for i in range(N):
		score_diff, pred0_score, pred1_score, pred0_reg, pred1_reg, predface, strandid = alldata[i]
		newdata.append((int(strandid), pred0_reg, predface))
	newdata = sorted(newdata)
	return np.array(newdata)[:,1:]

## correct facing according to reg
##    what it actually does:
##    if 2 consuctive pred_reg0 needs to be corrected by pred_reg1 according to the pred_fac, it is probably that the pred_fac is wrong
def correct_fac(pdbdata):
	newdata = []
	N = len(pdbdata)
	facingcorrecting_cond1_strandids = []
	facingcorrecting_cond2_strandids = []

	for i in range(N):
		ip1 = (i+1)%N
		score_diff, pred0_score, pred1_score, pred0_reg, pred1_reg, predface, strandid = pdbdata[i]
		dummy,      dummy,       dummy,       dummy,     dummy, predface_ip1, dummy    = pdbdata[ip1]
		pred0_reg = int(pred0_reg)
		pred1_reg = int(pred1_reg)
		#truereg = int(truereg)
		if predface==predface_ip1:
			if pred0_reg%2!=0 and pred1_reg%2==0:
				###------------------------------------------
				### facing correction condition 1
				#print '*wew', i, predface, trueface
				#print ' ', (i+1)%N, pdbdata[(i+1)%N][5], pdbdata[(i+1)%N][7]
				facingcorrecting_cond1_strandids.append(i)
				###------------------------------------------
			elif pred0_reg%2==0 and pred1_reg%2!=0:
				###------------------------------------------
				pass
				###------------------------------------------
			elif pred0_reg%2!=0 and pred1_reg%2!=0:
				###------------------------------------------
				### facing correction condition 2
				#print '*hew', i, predface, trueface
				#print ' ', (i+1)%N, pdbdata[(i+1)%N][5], pdbdata[(i+1)%N][7]
				facingcorrecting_cond2_strandids.append(i)
				###------------------------------------------
		else: ## predface != predface_ip1
			if pred0_reg%2==0 and pred1_reg%2==0:
				###------------------------------------------
				### facing correction condition 2
				#print '*how', i, predface, trueface
				#print ' ', (i+1)%N, pdbdata[(i+1)%N][5], pdbdata[(i+1)%N][7]
				facingcorrecting_cond2_strandids.append(i)
				###------------------------------------------
			elif pred0_reg%2!=0 and pred1_reg%2==0:
				###------------------------------------------
				pass
				###------------------------------------------
			elif pred0_reg%2==0 and pred1_reg%2!=0:
				###------------------------------------------
				### facing correction condition 1
				#print '*wow', i, predface, trueface
				#print ' ', (i+1)%N, pdbdata[(i+1)%N][5], pdbdata[(i+1)%N][7]
				facingcorrecting_cond1_strandids.append(i)
				###------------------------------------------
		newdata.append( (score_diff, pred0_score, pred1_score, pred0_reg, pred1_reg, predface, strandid) )

	## correct facings here
	for i in range(len(newdata)):
		im1 = (i-1)%N
		if i in facingcorrecting_cond1_strandids and im1 in facingcorrecting_cond1_strandids:
			if newdata[i][5] == 'OUT':
				newdata[i][5] == 'IN'
			else:
				newdata[i][5] == 'OUT'
	#	if i in facingcorrecting_cond2_strandids and im1 in facingcorrecting_cond2_strandids:
	#		if newdata[i][5] == 'OUT':
	#			newdata[i][5] == 'IN'
	#		else:
	#			newdata[i][5] == 'OUT'

	return np.array(newdata)


def filter_by_fac(pdbdata):
	newdata = []
	N = len(pdbdata)

	for i in range(N):
		ip1 = (i+1)%N
		score_diff, pred0_score, pred1_score, pred0_reg, pred1_reg, predface, strandid = pdbdata[i]
		dummy,      dummy,       dummy,       dummy,     dummy, predface_ip1, dummy    = pdbdata[ip1]
		pred0_reg = int(pred0_reg)
		pred1_reg = int(pred1_reg)
		#truereg = int(truereg)
		if predface==predface_ip1:
			if pred0_reg%2!=0 and pred1_reg%2==0:
				###------------------------------------------
				### filter regs using facing
				score_diff = 0
				pred0_score = pred1_score
				pred0_reg = pred1_reg
				#print pdb,i,'wew'
				pass
				###------------------------------------------
			elif pred0_reg%2==0 and pred1_reg%2!=0:
				###------------------------------------------
				### filter regs using facing
				score_diff = 0
				pred1_score = pred0_score
				pred1_reg = pred0_reg
				pass
				###------------------------------------------
			elif pred0_reg%2!=0 and pred1_reg%2!=0:
				###------------------------------------------
				pass
				###------------------------------------------
		else: ## predface != predface_ip1
			if pred0_reg%2==0 and pred1_reg%2==0:
				###------------------------------------------
				pass
				###------------------------------------------
			elif pred0_reg%2!=0 and pred1_reg%2==0:
				###------------------------------------------
				### filter regs using facing
				score_diff = 0
				pred1_score = pred0_score
				pred1_reg = pred0_reg
				#print pdb,i,'how'
				pass
				###------------------------------------------
			elif pred0_reg%2==0 and pred1_reg%2!=0:
				###------------------------------------------
				### in this case, correcting pred0 does more bad than good, according to manual counted results
				#score_diff = 0
				#pred0_score = pred1_score
				#pred0_reg = pred1_reg
				#print pdb,i,'wow'
				pass
				###------------------------------------------
		newdata.append( (score_diff, pred0_score, pred1_score, pred0_reg, pred1_reg, predface, strandid) )

	return np.array(newdata)


if __name__ == '__main__':
	#print 'Optimization #:'
	#print '0             : fac_crct'
	#print '1             : fac_crct + fac_fltr'
	#print '2             : fac_crct + fac_fltr + parity'
	#print '3             : fac_crct + parity'
	#print '4             : fac_crct + fac_fltr + brute'
	#print '5             : fac_crct + brute'
	#print '6             : fac_crct + fac_fltr + parity + brute'
	#print '7             : fac_crct + parity + brute'
	#print 'o/w (default) : no optimization'

	scorefn = glob.glob(args.folder+'/register.score')[0]
	fout = open(args.folder+'/register.pred', 'w')
	res = np.loadtxt(glob.glob(args.folder+'/*.res')[0]).astype(int)
	strandsfn = glob.glob(args.folder+'/*.strands')[0]

	peris = []
	extras = []
	with open(strandsfn) as f:
		lines = f.readlines()
	for i in range(len(lines)):
		end1, end2 = lines[i].split()
		if i%2==0:
			peris.append(int(end1))
			extras.append(int(end2))
		else:
			peris.append(int(end2))
			extras.append(int(end1))
	pred_facing = determine_facings_odds(peris, extras, res) # accuracy 723/736
	perifacing_dict = pred_facing

	###################################
	### adjustment starts from here ###
	###################################


	with open(scorefn) as f:
		lines = f.readlines()
	pdbdata_ori = []
	for strandid in range(len(lines)):
		split = lines[strandid].split()
		strandpairs = []
		for pair in split:
			perireg,score = pair.split(':')
			perireg = int(perireg)
			score = float(score)
			strandpairs.append((score,perireg))
		cand1, cand0 = sorted(strandpairs)[-2:]


		perifacing = perifacing_dict[strandid]
		## original data without adjustment
		## 0,			1,		2,		3,			4,			5,			6,			
		## score_diff,	score0,	score1,	pred_reg0,	pred_reg1,	pred_fac,	strand_id	
		pdbdata_ori.append( (float('%.3f'%((cand0[0]-cand1[0])/abs(cand0[0]+cand1[0]))),float('%.3f'%(cand0[0])),float('%.3f'%(cand1[0])),cand0[1],cand1[1],perifacing,strandid) )

		# all these tried c0+c1, abs(c0+c1), c0, 1
		# best is abs(c0+c1)


	################################################
	### process info of this pdb start from here ###
	################################################

	## get target shear
	targetshears = get_common_targetshears(pdbdata_ori)

	## get statistics of original prediction
	pdbdata_ori = np.array(pdbdata_ori)

	if args.opt == 0:
		## opt 0 : no adjustment, just a little bit facing correction
		pdbdata0 = np.array(pdbdata_ori)
		tofinal = correct_fac(pdbdata0)
		
	elif  args.opt == 1:
		## opt 1 : facing filtering
		pdbdata1 = np.array(pdbdata_ori)
		pdbdata1 = correct_fac(pdbdata1)
		tofinal = filter_by_fac(pdbdata1)
        
	elif  args.opt == 2:
		## opt 2 : facing filtering + parity first adjustment
		pdbdata2 = np.array(pdbdata_ori)
		pdbdata2 = correct_fac(pdbdata2)
		pdbdata2 = filter_by_fac(pdbdata2)
		pdbdata2 = np.array(sorted(pdbdata2.tolist()))
		tofinal = parity_first_shear_adjustment(pdbdata2,targetshears)
        
	elif  args.opt == 3:
		## opt 3 : parity first adjustment
		pdbdata3 = np.array(pdbdata_ori)
		pdbdata3 = correct_fac(pdbdata3)
		pdbdata3 = np.array(sorted(pdbdata3.tolist()))
		tofinal = parity_first_shear_adjustment(pdbdata3,targetshears)
        
	elif  args.opt == 4:
		## opt 4 : facing filtering + brute adjustment
		pdbdata4 = np.array(pdbdata_ori)
		pdbdata4 = correct_fac(pdbdata4)
		pdbdata4 = filter_by_fac(pdbdata4)
		pdbdata4 = np.array(sorted(pdbdata4.tolist()))
		tofinal = brute_shear_adjustment(pdbdata4,targetshears)
        
	elif  args.opt == 5:
		## opt 5 : brute adjustment
		pdbdata5 = np.array(pdbdata_ori)
		pdbdata5 = correct_fac(pdbdata5)
		pdbdata5 = np.array(sorted(pdbdata5.tolist()))
		tofinal = brute_shear_adjustment(pdbdata5,targetshears)
        
	elif  args.opt == 6:
		## opt 6 : facing filtering + parity first adjustment + brute adjustment
		pdbdata6 = np.array(pdbdata_ori)
		pdbdata6 = correct_fac(pdbdata6)
		pdbdata6 = filter_by_fac(pdbdata6)
		pdbdata6 = np.array(sorted(pdbdata6.tolist()))
		pdbdata6 = parity_first_shear_adjustment(pdbdata6,targetshears)
		tofinal = brute_shear_adjustment(pdbdata6,targetshears)

	elif  args.opt == 7:
		## opt 7 : parity first adjustment + brute adjustment
		pdbdata7 = np.array(pdbdata_ori)
		pdbdata7 = correct_fac(pdbdata7)
		pdbdata7 = np.array(sorted(pdbdata7.tolist()))
		pdbdata7 = parity_first_shear_adjustment(pdbdata7,targetshears)
		tofinal = brute_shear_adjustment(pdbdata7,targetshears)

	else:
		## opt -1 : no adjustment
		pdbdatam1 = np.array(pdbdata_ori)
		tofinal = pdbdatam1


	finaldata = get_final(tofinal)

	if fout is not None:
		for i in range(len(finaldata)):
			for itm in finaldata[i]:
				fout.write(' '+str(itm))
			fout.write('\n')
	else:
		#print pdb, len(finaldata), sum(finaldata[:,1].astype(int))
		#labs = [ 'strdid', 'predreg', 'predfac', 'truereg', 'truefac', 'peri' ]
		#print pandas.DataFrame(finaldata, columns=labs)
		#print '='*40
		pass


	fout.close()


