import numpy as np 
import os
import sys
import pdb


def get_af_bin(maf):
	if maf < .05:
		print('assumption error')
		pdb.set_trace()
	elif maf >= .05 and maf < .1:
		maf_bin = 0
	elif maf >= .1 and maf < .15:
		maf_bin = 1
	elif maf >= .15 and maf < .2:
		maf_bin = 2
	elif maf >= .2 and maf < .25:
		maf_bin = 3
	elif maf >= .25 and maf < .3:
		maf_bin = 4
	elif maf >= .3 and maf < .35:
		maf_bin = 5
	elif maf >= .35 and maf < .4:
		maf_bin = 6
	elif maf >= .4 and maf < .45:
		maf_bin = 7
	elif maf >= .45 and maf <= .5:
		maf_bin = 8
	else:
		print('assumption error')
		pdb.set_trace()
	return maf_bin


uldsc_variant_file = sys.argv[1]
allele_frequency_output_file = sys.argv[2]


f = open(uldsc_variant_file)
t = open(allele_frequency_output_file,'w')
t.write('variant_id\tmaf\tmaf_bin\n')

head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	if len(data) != 6:
		print('assumptoin eororr')
		pdb.set_trace()
	variant_id = data[0] + ':' + data[1] + ':' + data[4].split(':')[0] + ':' + data[5].split(':')[0]
	af = float(data[4].split(':')[1])
	if af > .5:
		maf = 1.0 -af 
	else:
		maf = af
	maf_bin = get_af_bin(maf)
	t.write(variant_id + '\t' + str(maf) + '\t' + str(maf_bin) + '\n')
f.close()
t.close()



