import numpy as np 
import os
import sys
import pdb





one_k_genomes_sample_annotation_file = sys.argv[1]
eur_1k_genomes_samples_file = sys.argv[2]



f = open(one_k_genomes_sample_annotation_file)
t = open(eur_1k_genomes_samples_file, 'w')

head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split()
	if len(data) != 4:
		print('assumption error')
		pdb.set_trace()
	if head_count == 0:
		head_count = head_count + 1
		continue
	if data[2] != 'EUR':
		continue
	t.write(data[0] + '\n')

f.close()
t.close()