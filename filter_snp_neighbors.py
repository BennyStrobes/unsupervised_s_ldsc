import numpy as np 
import os
import sys
import pdb
import time









# Command line args
r_squared_threshold = float(sys.argv[1])
input_neighbors_file = sys.argv[2]
output_neighbors_file = sys.argv[3]


f = open(input_neighbors_file)
t = open(output_neighbors_file,'w')
counter = 0
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\n')
		continue
	counter = counter + 1
	ld_score = float(data[3])
	neighbor_indices = np.asarray(data[4].split(','))
	neighbor_r_squared_string = np.asarray(data[5].split(','))
	neighbor_r_squared = neighbor_r_squared_string.astype(float)
	new_indices = np.where(neighbor_r_squared > r_squared_threshold)[0]
	new_neighbor_r_squared = neighbor_r_squared[new_indices]
	new_neighbor_r_squared_string = neighbor_r_squared_string[new_indices]
	new_neighbor_indices = neighbor_indices[new_indices]
	new_ld_score = np.sum(new_neighbor_r_squared) 
	t.write(data[0] + '\t' + data[1] + '\t' + data[2] + '\t' + str(new_ld_score) + '\t')
	t.write(','.join(new_neighbor_indices) + '\t' + ','.join(new_neighbor_r_squared_string) + '\n')
f.close()
t.close()