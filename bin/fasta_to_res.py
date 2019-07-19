#!/usr/bin/python
import sys


if len(sys.argv) != 3:
	print sys.argv[0], 'fasta_file', 'output_res_file'
	sys.exit(0)

fastafn = sys.argv[1]
resfn = sys.argv[2]

with open(fastafn) as f:
	readstart = False
	seq = ''
	while True:
		line = f.readline()
		if not line:
			break

		line = line.strip()
		if len(line)==0:
			continue
		if line[0] == '>':
			if readstart:
				break
			else:
				readstart = True
		else:
			if not readstart:
				readstart = True
			seq += line

seq = seq.upper()

one_list = [ 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' ]
index_list = range(20)
one_to_index_dict = dict( zip(one_list, index_list) )


with open(resfn, 'w') as f:
	for aa in seq:
		f.write(str(one_to_index_dict[aa])+'\n')
