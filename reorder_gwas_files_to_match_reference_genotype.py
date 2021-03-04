import numpy as np 
import os
import sys
import pdb
import gzip



def extract_reference_genotype_snps(ordered_snp_file):
	f = open(ordered_snp_file)
	snp_arr = []
	snp_dicti = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		snp_name = data[0] + ':' + data[1] + ':' + data[4].split(':')[0] + ':' + data[5].split(':')[0]
		snp_arr.append(snp_name)
		snp_dicti[snp_name] = (0.0, 0.0)
	f.close()
	return np.asarray(snp_arr), snp_dicti





def extract_ukbb_study_names(gwas_studies_file):
	f = open(gwas_studies_file)
	study_names = []
	study_files = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		study_name = data[0]
		study_file = data[1]
		study_names.append(study_name)
		study_files.append(study_file)
	f.close()
	return np.asarray(study_names), np.asarray(study_files)



def filter_and_print_gwas(study_file, ordered_snp_arr, snp_dicti, chrom_num, output_file):
	chrom_string = chrom_num + ':'
	f = gzip.open(study_file)
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
		if len(data) == 11:
			adder = 0
		elif len(data) == 12:
			adder = 1
		else:
			print('assumption error')
			pdb.set_trace()
		# Throw out variants not found in 1K genomes
		if data[0] not in snp_dicti:
			continue
		z_score = float(data[(9+adder)])
		chi_sq = np.square(z_score)
		# If we made it here, we have passed all of our filters
		# And thus we update the variant_dicti
		snp_dicti[data[0]] = (1.0, chi_sq)
		#variant_dicti[data[0]] = variant_dicti[data[0]] + 1
	f.close()
	t = open(output_file, 'w')
	t.write('snp_id\tchi_squared\n')
	arr = []
	missing = 0
	for snp in ordered_snp_arr:
		tupler = snp_dicti[snp]
		t.write(snp + '\t' + str(tupler[1]) + '\n')
		arr.append(tupler[1])
		if tupler[0] == 0.0:
			missing = missing + 1
			print('assumption error missing snps')
			pdb.set_trace()
		snp_dicti[snp] = (0.0, 0.0)
	t.close()
	print('mean: ' + str(np.mean(arr)))
	return 



######################
# Command line args
######################
chrom_num = sys.argv[1]
ordered_snp_file = sys.argv[2]
gwas_studies_file = sys.argv[3]
processed_ukbb_dir = sys.argv[4]



ordered_snp_arr, snp_dicti = extract_reference_genotype_snps(ordered_snp_file)


study_names, study_files = extract_ukbb_study_names(gwas_studies_file)


for study_iter in range(len(study_names)):
	study_name = study_names[study_iter]
	study_file = study_files[study_iter]
	print(study_name)
	output_file = processed_ukbb_dir + 'chr' + chrom_num + '_' + study_name + '_chi_squared.txt'
	filter_and_print_gwas(study_file, ordered_snp_arr, snp_dicti, chrom_num, output_file)




