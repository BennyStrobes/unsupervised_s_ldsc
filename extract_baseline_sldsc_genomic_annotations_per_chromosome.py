import numpy as np 
import os
import sys
import pdb




def extract_genomic_annotations_for_set_of_variants(variant_to_anno, baseline_ldsc_input_file):
	g = open(baseline_ldsc_input_file)
	head_county = 0
	for line in g:
		line = line.rstrip()
		data = line.split('\t')
		if len(data) != 101:
			print('assumption eeorror')
			pdb.set_trace()
		if head_county == 0:
			head_county = head_county + 1
			continue
		variant_id = data[0] + ':' + data[1]
		if variant_id not in variant_to_anno:
			continue
		variant_to_anno[variant_id] = '\t'.join(data[5:])
	g.close()
	return variant_to_anno



baseline_ldsc_input_file = sys.argv[1]
baseline_sldsc_output_file = sys.argv[2]
uldsc_variant_file = sys.argv[3]

# Extract annotations for $window_size variants at one time (for memory purposes)
window_size = 8000

# Open output file handle
t = open(baseline_sldsc_output_file, 'w')
# print header
f = open(baseline_ldsc_input_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		header_info = data[5:]
		continue
	break
f.close()

t.write('variant_id\t' + '\t'.join(header_info) + '\n')
missing_string = '\t'.join(['NA']*len(header_info))

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
	variant_to_anno[variant_id_short] = 0
	window_counter = window_counter + 1
	if window_counter == window_size:
		print(counter)
		misses = 0
		counter = counter + 1
		variant_to_anno = extract_genomic_annotations_for_set_of_variants(variant_to_anno, baseline_ldsc_input_file)
		if len(variant_to_anno) != window_size or len(variant_arr) != window_size:
			print('assumptoin eroror')
			pdb.set_trace()
		for variant_id in variant_arr:
			short_variant_id = variant_id.split(':')[0] + ':' + variant_id.split(':')[1]
			if variant_to_anno[short_variant_id] == 0:
				t.write(variant_id + '\t' + missing_string + '\n')
				misses = misses + 1
			else:
				t.write(variant_id + '\t' + variant_to_anno[short_variant_id] + '\n')
		print(misses)
		# print line
		window_counter = 0
		variant_arr = []
		variant_to_anno = {}
f.close()

for variant_id in variant_arr:
	short_variant_id = variant_id.split(':')[0] + ':' + variant_id.split(':')[1]
	if variant_to_anno[short_variant_id] == 0:
		t.write(variant_id + '\t' + missing_string +'\n')
	else:
		t.write(variant_id + '\t' + variant_to_anno[short_variant_id] + '\n')
t.close()