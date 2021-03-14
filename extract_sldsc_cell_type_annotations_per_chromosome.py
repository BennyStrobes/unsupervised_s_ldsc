import numpy as np 
import os
import sys
import pdb






def add_cell_type_annotations_to_dictionary(input_file, variant_to_anno, column_index):
	f = open(input_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		variant_id = data[0] + ':' + data[1]
		if variant_id not in variant_to_anno:
			continue
		variant_to_anno[variant_id][column_index] = float(data[4])
	f.close()
	return variant_to_anno



sldsc_cell_type_annotation_dir = sys.argv[1]
chrom_num = sys.argv[2]
cell_type_sldsc_output_file = sys.argv[3]
uldsc_variant_file = sys.argv[4]


# Create mapping from number to names
mapping = {1:'Adrenal_Pancreas', 2: 'Cardiovascular', 3: 'CNS', 4: 'Connective_Bone', 5: 'GI', 6: 'Hematopoietic', 7: 'Kidney', 8: 'Liver', 9: 'Other', 10: 'SkeletalMuscle'}
int_names = range(1,11)


# Open file handle for list of variants used in our analysis
# Grab $window_size variants at a time
f = open(uldsc_variant_file)
head_count = 0

# Keep track of number of variants currently in window
window_counter = 0
# Create array containing variant ids (in order)
variant_arr = []
# Create dictionary mapping variant ids to string of genomic annotations describing the variant
variant_to_anno = {}

counter = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	if len(data) != 6:
		print('assumption error')
		pdb.set_trace()
	variant_id_long = data[0] + ':' + data[1] + ':' + data[4].split(':')[0] + ':' + data[5].split(':')[0]
	variant_id_short = data[0] + ':' + data[1]
	variant_arr.append(variant_id_long)
	variant_to_anno[variant_id_short] = np.zeros(10) - 2.2
	window_counter = window_counter + 1
f.close()

cell_types = []
for int_name in int_names:
	input_file = sldsc_cell_type_annotation_dir + 'cell_type_group.' + str(int_name) + '.' + chrom_num + '.annot'
	cell_type = mapping[int_name]
	cell_types.append(cell_type)
	variant_to_anno = add_cell_type_annotations_to_dictionary(input_file, variant_to_anno, int_name - 1)

cell_types = np.asarray(cell_types)

missing_string = '\t'.join(['NA']*len(cell_types))


t = open(cell_type_sldsc_output_file, 'w')
t.write('variant_id\t' + '\t'.join(cell_types) + '\n')
for variant_id_long in variant_arr:
	variant_id = variant_id_long.split(':')[0] + ':' + variant_id_long.split(':')[1]
	output_vec = variant_to_anno[variant_id]
	if len(np.where(output_vec == -2.2)[0]) == 0:
		t.write(variant_id_long + '\t' + '\t'.join(output_vec.astype(str)) + '\n')
	else:
		if len(np.where(output_vec == -2.2)[0]) != 10:
			print('assumption erorr')
			pdb.set_trace()
		t.write(variant_id_long + '\t' + missing_string + '\n')
t.close()


