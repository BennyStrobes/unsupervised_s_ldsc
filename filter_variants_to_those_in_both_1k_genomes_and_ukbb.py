import numpy as np 
import os
import sys
import pdb
from scipy.stats import chi2


def extract_1k_genomes_variants(genotype_1k_genomes_file):
	f = open(genotype_1k_genomes_file)
	arr = []
	dicti = {}
	head_count = 0
	position_used = {}
	for line in f:
		if head_count == 0:
			head_count = head_count + 1
			continue
		line = line.rstrip()
		data = line.split()
		if len(data) != 6:
			print('assumtpoin eororor')
			pdb.set_trace()
		if data[2] != '2':
			print('assumption error')
			pdb.set_trace()
		chrom_num = data[0]
		variant_position = data[1]
		variant_position_id = chrom_num + ':' + variant_position
		ref_allele = data[4].split(':')[0]
		alt_allele = data[5].split(':')[0]
		# throw out indels
		#if len(ref_allele) != 1 or len(alt_allele) != 1:
		#	continue
		if variant_position_id in position_used:
			print('assumptin error')
			pdb.set_trace()
		position_used[variant_position_id] = 1
		variant_id = chrom_num + ':' + variant_position + ':' + ref_allele + ':' + alt_allele
		if variant_id in dicti:
			print('assumption errror')
			pdb.set_trace()
		arr.append(variant_id)
		dicti[variant_id] = 0
	f.close()
	return np.asarray(arr), dicti

def extract_ukbb_studies(ukbb_studies_file):
	f = open(ukbb_studies_file)
	ukbb_study_names = []
	ukbb_study_files = []
	ukbb_sample_sizes = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		study_name = data[0]
		study_file = data[1]
		num_samples = int(data[2])
		ukbb_study_names.append(study_name)
		ukbb_study_files.append(study_file)
		ukbb_sample_sizes.append(num_samples)
	f.close()
	return np.asarray(ukbb_study_names), np.asarray(ukbb_study_files), np.asarray(ukbb_sample_sizes)

def pvalue_to_z_score(pvalue):
	return np.sqrt(chi2.isf(pvalue, 1))

def update_variant_dicti_with_ukbb_data_from_single_study(variant_dicti, chrom_num, ukbb_study_file, ukbb_sample_size):
	chrom_string = chrom_num + ':'
	max_chi_squared = np.max([80, 0.001*ukbb_sample_size])
	f = open(ukbb_study_file)
	head_count = 0
	for line in f:
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Limit to variants from desired chromosome
		if line.startswith(chrom_string) == False:
			continue
		# Parse line
		line = line.rstrip()
		data = line.split()
		# error checking to make sure we have the apprpriate number of lines
		if len(data) != 11:
			print('assumption error')
			pdb.set_trace()
		# Throw out variants not found in 1K genomes
		if data[0] not in variant_dicti:
			continue
		# Check if the variant passes filters in UKBB
		maf = float(data[2])
		if maf < .05:
			continue
		if maf > .5:
			print('assumption error')
			pdb.set_trace()
		low_confidence_boolean = data[3]
		if low_confidence_boolean == "true":
			continue
		#pvalue = float(data[10])
		#z_score = pvalue_to_z_score(pvalue)
		z_score = float(data[9])
		chi_sq = np.square(z_score)
		if chi_sq > max_chi_squared:
			continue
		# If we made it here, we have passed all of our filters
		# And thus we update the variant_dicti
		variant_dicti[data[0]] = variant_dicti[data[0]] + 1
	f.close()
	return variant_dicti

def update_variant_dicti_with_ukbb_data(variant_dicti, ukbb_studies_file, chrom_num):
	ukbb_study_names, ukbb_study_files, ukbb_sample_sizes = extract_ukbb_studies(ukbb_studies_file)
	num_studies = len(ukbb_study_names)
	# Loop through UKBB studies
	for study_index, study_name in enumerate(ukbb_study_names):
		print(str(study_index) + ' : ' + study_name)
		ukbb_study_file = ukbb_study_files[study_index]
		ukbb_sample_size = ukbb_sample_sizes[study_index]
		variant_dicti = update_variant_dicti_with_ukbb_data_from_single_study(variant_dicti, chrom_num, ukbb_study_file, ukbb_sample_size)

	return variant_dicti, num_studies




chrom_num = sys.argv[1]  # Chromosome to run analysis on 
processed_1k_genomes_genotype_dir = sys.argv[2]  #1kG genotype dir
ukbb_studies_file = sys.argv[3]  # file with line for each ukbb study
variant_file = sys.argv[4]  # output file
variant_position_file = sys.argv[5] # output file 2

# Genotype file from UKBB
genotype_1k_genomes_file = processed_1k_genomes_genotype_dir + 'chr_' + chrom_num + '_1k_genomes_european_only_maf_thresh_05.frq'

# Create both:
## 1. ordered array of 1KG variants
## 2. Dictionary of 1KG variants
ordered_1k_genomes_variants, variant_dicti = extract_1k_genomes_variants(genotype_1k_genomes_file)

# add info to variant dicti on whether each variant was found in all UKBB studies
variant_dicti, num_studies = update_variant_dicti_with_ukbb_data(variant_dicti, ukbb_studies_file, chrom_num)

# Print variants to output file
t = open(variant_file, 'w')
for variant in ordered_1k_genomes_variants:
	if variant_dicti[variant] == num_studies:
		t.write(variant + '\n')
t.close()


# Print output file 2
f = open(variant_file)
t = open(variant_position_file,'w')
used_positions = {}
for line in f:
	line = line.rstrip()
	data = line.split(':')
	position_id = data[0] + ':' + data[1]
	if position_id in used_positions:
		print('assumption error')
		pdb.set_trace()
	used_positions[position_id] = 1
	t.write(data[0] + '\t' + data[1] + '\n')
f.close()
t.close()
