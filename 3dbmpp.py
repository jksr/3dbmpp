#!/usr/bin/python

import argparse
parser = argparse.ArgumentParser(description='Predict 3D structures of beta barrel membrane proteins')
parser.add_argument('--group', type=int, default=0, choices=[0,1,2,3,4], required=True)
parser.add_argument('--folder', type=str, required=True)
parser.add_argument('--scwrlpath', type=str, required=True)
args = parser.parse_args()

import glob
import os 

fastafn = glob.glob(args.folder+'/*.fasta')[0]
resfn = fastafn+'.res'
os.system('python bin/fasta_to_res.py '+fastafn+' '+resfn)


strandsfn = glob.glob(args.folder+'/*.strands')[0]
#if len(glob.glob(args.folder+'/*.peris')) == 0:
#	perisfn = strandsfn+'.peris'
#	with open(strandsfn) as f:
#		lines = f.readlines()
#	with open(perisfn, 'w') as f:
#		for i in range(len(lines)):
#			f.write( lines[i].strip().split()[i%2]+'\n' )
#perisfn = glob.glob(args.folder+'/*.peris')[0]

if len(glob.glob(args.folder+'/*.psicov'))==0:
	os.system('touch '+args.folder+'/dummy.psicov')
psicovfn = glob.glob(args.folder+'/*.psicov')[0]
seqcovfn = psicovfn+'.seqcov'

print 'Loading sequence convariation information'
#os.system('python bin/seqcov_to_reg_score.py --strands_file '+strandsfn+' --ec_file '+psicovfn+' --peris_file '+perisfn+' --output_file '+seqcovfn)
os.system('python bin/seqcov_to_reg_score.py --strands_file '+strandsfn+' --ec_file '+psicovfn+' --output_file '+seqcovfn)

print 'Predicting registers'
scorefn = args.folder+'/register.score'
os.system('bin/pred_reg --group '+str(args.group)+' --res '+resfn+' --strand '+strandsfn+' --seqcov '+seqcovfn+' > '+scorefn)

print 'Adjusting registers using global information'
os.system('bin/shear_optimization.py --opt 2 --folder '+args.folder+' --oddsdir  odds')


tmpfolder = args.folder+'/tmp'
if not os.path.exists(tmpfolder):
	os.mkdir(tmpfolder)

print 'Constructing C-alpha trace'
os.system('bin/construct_ca_pdb.py --folder '+args.folder+' --tmpfolder '+tmpfolder)
print 'Constructing mainchain residues using BBQ algorithm'
os.system('bin/construct_bb_pdb.py --capdbfn '+tmpfolder+'/ca_reidx_ext.pdb --bbpdbfn '+tmpfolder+'/bb_reidx_ext.pdb')
print 'Constructing sidechain residues using Scwrl4'
os.system('bin/construct_sc_pdb.py --bbpdbfn '+tmpfolder+'/bb_reidx_ext.pdb --scpdbfn '+tmpfolder+'/sc_reidx_ext.pdb --scwrl '+args.scwrlpath)
print 'Finalizing the pdb structure'
os.system('bin/finalize_pdb.py --scpdbfn '+tmpfolder+'/sc_reidx_ext.pdb --pdbfn '+args.folder+'/result.pdb --indexmap '+tmpfolder+'/reidx.map' )
print 'Done'
print 'The predicted structure is saved in '+args.folder+'/result.pdb'
