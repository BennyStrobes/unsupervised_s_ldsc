import numpy as np 
import os
import sys
import pdb




def get_study_sample_size(study_file):
	g = open(study_file)
	head_count = 0
	for line in g:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		if head_count == 1:
			head_count = head_count + 1
			sample_size = data[4]
			continue
		if data[4] != sample_size:
			print('assumption error')
			pdb.set_trace()
	g.close()
	return sample_size







old_ukbb_studies_file = sys.argv[1]
new_ukbb_studies_file = sys.argv[2]


f = open(old_ukbb_studies_file)
t = open(new_ukbb_studies_file,'w')

head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\t' + 'sample_size' + '\n')
		continue
	study_name = data[0]
	study_file = data[1]
	print(study_name)
	study_sample_size = get_study_sample_size(study_file)
	t.write(study_name + '\t' + study_file + '\t' + study_sample_size + '\n')
f.close()
t.close()
