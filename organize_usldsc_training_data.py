import numpy as np 
import os
import sys
import pdb
import pickle

def extract_array_of_chromosomes(chromosomes_for_training_file):
	arr = []
	f = open(chromosomes_for_training_file)
	for line in f:
		line = line.rstrip()
		arr.append(line)
	f.close()
	return np.asarray(arr)


def save_ukbb_study_data_to_npy(study, chromosome_arr, processed_ukbb_dir, npy_output_file, variant_ids):
	arr = []
	counter = 0
	for chrom_num in chromosome_arr:
		file_name = processed_ukbb_dir + 'chr' + chrom_num + '_' + study + '_chi_squared.txt'
		g = open(file_name)
		head_count = 0
		for line in g:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			if variant_ids[counter] != data[0]:
				print('assumption errro')
				pdb.set_trace()
			counter = counter + 1
			arr.append(float(data[1]))
		g.close()
	arr = np.asarray(arr)
	np.save(npy_output_file, arr)

def extract_variant_ids(input_study_file, chromosome_arr, processed_ukbb_dir):
	f = open(input_study_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if head_count == 1:
			study = data[0]
			head_count = head_count + 1
			continue
	f.close()
	arr = []
	for chrom_num in chromosome_arr:
		file_name = processed_ukbb_dir + 'chr' + chrom_num + '_' + study + '_chi_squared.txt'
		g = open(file_name)
		head_count = 0
		for line in g:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			arr.append(data[0])
		g.close()
	return np.asarray(arr)


def pickle_ukbb_chi_squared_stats(chromosome_arr, processed_ukbb_dir, organized_training_data_dir):
	input_study_file = processed_ukbb_dir + 'updated_ukbb_studies.txt'
	output_study_file = organized_training_data_dir + 'usldsc_training_studies.txt'
	variant_ids = extract_variant_ids(input_study_file, chromosome_arr, processed_ukbb_dir)
	npy_variant_output_file = organized_training_data_dir + 'usldsc_training_ukbb_variant_ids.npy'
	np.save(npy_variant_output_file, variant_ids)

	f = open(input_study_file)
	t = open(output_study_file,'w')
	t.write('study\tfile_name\tsample_size\n')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		study = data[0]
		sample_size = data[2]
		npy_output_file = organized_training_data_dir + 'usldsc_training_' + study + '_chi_squared_stats.npy'
		t.write(study + '\t' + npy_output_file + '\t' + sample_size + '\n')
		save_ukbb_study_data_to_npy(study, chromosome_arr, processed_ukbb_dir, npy_output_file, variant_ids)
	f.close()
	t.close()

def pickle_pairwise_ld_info(chromosome_arr, processed_ld_score_dir, organized_training_data_dir):
	aggregate_ld_file = organized_training_data_dir + 'usldsc_training_pairwise_ld_files.txt'
	t = open(aggregate_ld_file,'w')
	t.write('chromosome\tpairwise_indices_file\tpairwise_ld_file\n')
	snp_counter = 0
	for chromosome in chromosome_arr:
		print(chromosome)
		# Pairwise ld file for this chromosome
		chromosome_pairwise_ld_file = processed_ld_score_dir + 'chr_' + chromosome + '_ld_scores_pairwise_l2_neighbors_filtered.2_clustered.txt'
		# Initialize array where each position is a variant and it maps to neighboring indices of variant
		neighboring_indices_arr = []
		# Initialize array where each position is a variant and it maps to pairwise lds of each of its neighbors
		neighboring_pairwise_lds_arr = []
		head_count = 0
		f = open(chromosome_pairwise_ld_file)
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				starting_index = snp_counter
				continue
			# Extract relevent fields
			variant_name = data[0] + ':' + data[2]
			cluster_id = 'chr' + chromosome + '_' + data[6]
			neighboring_indices = np.asarray(data[4].split(',')).astype(int) + starting_index
			pairwise_lds = np.asarray(data[5].split(',')).astype(float)
			variant_index = snp_counter
			neighboring_indices_arr.append(neighboring_indices)
			neighboring_pairwise_lds_arr.append(pairwise_lds)
			snp_counter = snp_counter + 1
		f.close()
		
		indices_output_file = organized_training_data_dir + 'usldsc_training_chr_' + chromosome + '_pairwise_ld_indices.npy'
		ld_output_file = organized_training_data_dir + 'usldsc_training_chr_' + chromosome + '_pairwise_lds.npy'

		np.save(indices_output_file, neighboring_indices_arr)
		np.save(ld_output_file, neighboring_pairwise_lds_arr)

		t.write(chromosome + '\t' + indices_output_file + '\t' + ld_output_file + '\n')
	t.close()

def extract_chromosome_pairwise_ld_info(raw_pairwise_ld_file):
	# Initialize array where each position is a variant and it maps to neighboring indices of variant
	neighboring_indices_arr = []
	# Initialize array where each position is a variant and it maps to pairwise lds of each of its neighbors
	neighboring_pairwise_lds_arr = []
	head_count = 0
	f = open(raw_pairwise_ld_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		neighboring_indices = np.asarray(data[4].split(',')).astype(int)
		pairwise_lds = np.asarray(data[5].split(',')).astype(float)
		neighboring_indices_arr.append(neighboring_indices)
		neighboring_pairwise_lds_arr.append(pairwise_lds)
	return neighboring_indices_arr, neighboring_pairwise_lds_arr


def pickle_cluster_level_pairwise_ld_data(chromosome_arr, processed_ld_score_dir, organized_training_data_dir):
	for chromosome in chromosome_arr:
		# Extract pairwise ld info on this chromosome
		raw_pairwise_ld_file = processed_ld_score_dir + 'chr_' + chromosome + '_ld_scores_pairwise_l2_neighbors_filtered.2_clustered.txt'
		neighboring_indices_arr, neighboring_pairwise_lds_arr = extract_chromosome_pairwise_ld_info(raw_pairwise_ld_file)
		# Now on this chromosome loop through snp clusters
		# And save npy matrix corresponding to LD-matrix across all snps in the cluster
		cluster_mapping_file = processed_ld_score_dir + 'chr_' + str(chromosome) + '_ld_scores_pairwise_l2_neighbors_filtered.2_cluster_mapping.txt'
		f = open(cluster_mapping_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			# Skip header
			if head_count == 0:
				head_count = head_count + 1
				continue
			# Name of cluster
			cluster_id = data[0]
			# Number of snps assigned to cluster
			num_snps = int(data[1])
			# Get indices of snps assigned to this cluster
			snp_indices = np.asarray(data[2].split(',')).astype(int)
			# Simple error checking
			if len(snp_indices) != num_snps:
				print('assumption error')
				pdb.set_trace()
			pairwise_correlation_mat = np.zeros((num_snps, num_snps))
			name_to_position = {}
			for position, name in enumerate(snp_indices):
				name_to_position[name] = position
			neighbor_positions_arr = []
			for snp_iter in range(num_snps):
				snp_neighbor_positions = []
				# wish to fill in pairwise_correlation_mat[snp_iter, :]
				snp_name = snp_indices[snp_iter]
				snp_neighbor_names = neighboring_indices_arr[snp_name]
				snp_neighbor_lds = neighboring_pairwise_lds_arr[snp_name]
				for neighbor_iter in range(len(snp_neighbor_names)):
					neighbor_name = snp_neighbor_names[neighbor_iter]
					neighbor_ld = snp_neighbor_lds[neighbor_iter]
					pairwise_correlation_mat[snp_iter, name_to_position[neighbor_name]] = neighbor_ld
					snp_neighbor_positions.append(name_to_position[neighbor_name])
				neighbor_positions_arr.append(np.asarray(snp_neighbor_positions))
			npy_output_file = organized_training_data_dir + 'usldsc_training_chr_' + chromosome + '_' + cluster_id + '_pairwise_ld_correlation_mat.npy'
			np.save(npy_output_file, pairwise_correlation_mat)
			npy_output_file = organized_training_data_dir + 'usldsc_training_chr_' + chromosome + '_' + cluster_id + '_pairwise_ld_neighbor_positions.npy'
			np.save(npy_output_file, neighbor_positions_arr)
		f.close()

def pickle_cluster_level_ukbb_data(chromosome_arr, study_names, processed_ukbb_dir, processed_ld_score_dir, organized_training_data_dir):
	# Loop through chromosomes
	for chromosome in chromosome_arr:
		# For this chromosome first create matrix of dimension (number_of_snps_on_chromosomeXnumber_of_ukbb_studies)
		# Where each element fo this matrix the chi squared statistic for that test in that study
		chromosome_ukbb = []
		for study in study_names:
			file_name = processed_ukbb_dir + 'chr' + chromosome + '_' + study + '_chi_squared.txt'
			g = open(file_name)
			arr = []
			head_count = 0
			for line in g:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				arr.append(float(data[1]))
			g.close()
			chromosome_ukbb.append(np.asarray(arr))
		chromosome_ukbb = np.transpose(np.asarray(chromosome_ukbb))
		num_snps = chromosome_ukbb.shape[0]
		used_snps = np.zeros(num_snps)

		# Now on this chromosome loop through snp clusters
		# And save npy array corresponding to chi-squared stats for that cluster
		cluster_mapping_file = processed_ld_score_dir + 'chr_' + str(chromosome) + '_ld_scores_pairwise_l2_neighbors_filtered.2_cluster_mapping.txt'
		f = open(cluster_mapping_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			# Skip header
			if head_count == 0:
				head_count = head_count + 1
				continue
			# Name of cluster
			cluster_id = data[0]
			# Number of snps assigned to cluster
			num_snps = int(data[1])
			# Get indices of snps assigned to this cluster
			snp_indices = np.asarray(data[2].split(',')).astype(int)
			# Save J*S length array to npy object
			npy_output_file = organized_training_data_dir + 'usldsc_training_chr_' + chromosome + '_' + cluster_id + '_chi_squared_stats.npy'
			np.save(npy_output_file, chromosome_ukbb[snp_indices,:].flatten(order='F'))
			# Fill in used snp array for error checking
			used_snps[snp_indices] = used_snps[snp_indices] + 1.0
		f.close()
		# need to make sure all snps are seen once and only once on this chromosome(error checking)
		if sum(used_snps != 1) != 0:
			print('assumption error')
			pdb.set_trace()

def get_number_of_snps_in_this_file(file_name):
	f = open(file_name)
	snp_counter = 0
	head_count = 0
	for line in f:
		line = line.rstrip()
		if head_count == 0:
			head_count = head_count + 1
			continue
		snp_counter = snp_counter + 1
	f.close()
	return snp_counter


def pickle_cluster_level_variant_names(chromosome_arr, processed_ld_score_dir, organized_training_data_dir):
	chromosome_starting_index = 0
	#used_snps = np.zeros(2070780)
	for chromosome in chromosome_arr:
		#print(chromosome)
		# Now on this chromosome loop through snp clusters
		# And save snp indices that belong tot that cluster
		cluster_mapping_file = processed_ld_score_dir + 'chr_' + str(chromosome) + '_ld_scores_pairwise_l2_neighbors_filtered.2_cluster_mapping.txt'
		f = open(cluster_mapping_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			# Skip header
			if head_count == 0:
				head_count = head_count + 1
				continue
			# Name of cluster
			cluster_id = data[0]
			# Get indices of snps assigned to this cluster
			snp_indices = np.asarray(data[2].split(',')).astype(int) + chromosome_starting_index
			# Save snp indices to output file
			npy_output_file = organized_training_data_dir + 'usldsc_training_chr_' + chromosome + '_' + cluster_id + '_snp_indices.npy'
			np.save(npy_output_file, snp_indices)
			# Error checking
			#for snp_index in snp_indices:
			#	if used_snps[snp_index] == 1.0:
			#		print('assumption errror')
			#		pdb.set_trace()
			#	used_snps[snp_index] = 1.0
		f.close()
		chromosome_starting_index = chromosome_starting_index + get_number_of_snps_in_this_file(processed_ld_score_dir + 'chr_' + chromosome + '_ld_scores_pairwise_l2_neighbors_filtered.2_clustered.txt')
	#pdb.set_trace()




def pickle_cluster_level_data(chromosome_arr, processed_ukbb_dir, processed_ld_score_dir, organized_training_data_dir):
	input_study_file = processed_ukbb_dir + 'updated_ukbb_studies.txt'
	temp_data = np.loadtxt(input_study_file,dtype=str,delimiter='\t')
	study_names = temp_data[1:,0]
	sample_sizes = temp_data[1:,2]
	pickle_cluster_level_ukbb_data(chromosome_arr, study_names, processed_ukbb_dir, processed_ld_score_dir, organized_training_data_dir)
	pickle_cluster_level_pairwise_ld_data(chromosome_arr, processed_ld_score_dir, organized_training_data_dir)
	pickle_cluster_level_variant_names(chromosome_arr, processed_ld_score_dir, organized_training_data_dir)
	aggregate_cluster_file = organized_training_data_dir + 'usldsc_training_snp_cluster_files.txt'
	t = open(aggregate_cluster_file,'w')
	t.write('cluster_id\tcluster_ukbb_file\tcluster_pairwise_ld_matrix_file\tcluster_variant_names_file\tcluster_pairwise_ld_neighbor_positions_file\n')
	for chromosome in chromosome_arr:
		#print(chromosome)
		# Now on this chromosome loop through snp clusters
		# And save snp indices that belong tot that cluster
		cluster_mapping_file = processed_ld_score_dir + 'chr_' + str(chromosome) + '_ld_scores_pairwise_l2_neighbors_filtered.2_cluster_mapping.txt'
		f = open(cluster_mapping_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			# Skip header
			if head_count == 0:
				head_count = head_count + 1
				continue
			# Name of cluster
			cluster_id = data[0]
			full_cluster_id = 'chr' + chromosome + ':' + cluster_id
			cluster_ukbb_file = organized_training_data_dir + 'usldsc_training_chr_' + chromosome + '_' + cluster_id + '_chi_squared_stats.npy'
			cluster_pairwise_ld_matrix_file = organized_training_data_dir + 'usldsc_training_chr_' + chromosome + '_' + cluster_id + '_pairwise_ld_correlation_mat.npy'
			cluster_variant_names_file = organized_training_data_dir + 'usldsc_training_chr_' + chromosome + '_' + cluster_id + '_snp_indices.npy'
			cluster_pairwise_ld_neighbor_positions_file = organized_training_data_dir + 'usldsc_training_chr_' + chromosome + '_' + cluster_id + '_pairwise_ld_neighbor_positions.npy'
			t.write(full_cluster_id + '\t' + cluster_ukbb_file + '\t' + cluster_pairwise_ld_matrix_file + '\t' + cluster_variant_names_file + '\t' + cluster_pairwise_ld_neighbor_positions_file + '\n')
		f.close()
	t.close()

######################
# command line args
########################
chromosomes_for_training_file = sys.argv[1]  # File containing list of chromosomes to train model on
processed_ukbb_dir = sys.argv[2]  # directory containing ukbb data
processed_ld_score_dir = sys.argv[3]  # Directory containing info on pairwise ld relationships
organized_training_data_dir = sys.argv[4]  # Output dir


# Extract array of chromosomes to run analysis on 
chromosome_arr = extract_array_of_chromosomes(chromosomes_for_training_file)
pickle_ukbb_chi_squared_stats(chromosome_arr, processed_ukbb_dir, organized_training_data_dir)

pickle_pairwise_ld_info(chromosome_arr, processed_ld_score_dir, organized_training_data_dir)

# Generate cluster chi squared stats
pickle_cluster_level_data(chromosome_arr, processed_ukbb_dir, processed_ld_score_dir, organized_training_data_dir)
