import numpy as np 
import os
import sys
import pdb
import usldsc_als
import usldsc_vi
import usldsc_shared_factor_vi
import usldsc_study_variance_vi


def extract_study_info(studies_file):
	f = open(studies_file)
	chi_squared_files = []
	study_sample_sizes = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		chi_squared_files.append(data[1])
		study_sample_sizes.append(int(data[2]))
	f.close()
	return np.asarray(chi_squared_files), np.asarray(study_sample_sizes)

def extract_pairwise_ld_files(pairwise_ld_summary_file):
	f = open(pairwise_ld_summary_file)
	pairwise_ld_indices_files = []
	pairwise_ld_files = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		pairwise_ld_indices_files.append(data[1])
		pairwise_ld_files.append(data[2])
	f.close()
	return np.asarray(pairwise_ld_indices_files), np.asarray(pairwise_ld_files)

def extract_number_of_snps(chi_squared_file):
	return len(np.load(chi_squared_file))

def extract_cluster_info_files(training_data_cluster_info_file):
	f = open(training_data_cluster_info_file)
	head_count = 0
	cluster_ids = []
	cluster_ukbb_files = []
	cluster_pairwise_ld_matrix_files = []
	cluster_variant_names_files = []
	cluster_variant_neighbor_positions_files = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		cluster_ids.append(data[0])
		cluster_ukbb_files.append(data[1])
		cluster_pairwise_ld_matrix_files.append(data[2])
		cluster_variant_names_files.append(data[3])
		cluster_variant_neighbor_positions_files.append(data[4])
	f.close()
	return np.asarray(cluster_ids), np.asarray(cluster_ukbb_files), np.asarray(cluster_pairwise_ld_matrix_files), np.asarray(cluster_variant_names_files), np.asarray(cluster_variant_neighbor_positions_files)

def run_model(studies_file, pairwise_ld_summary_file, training_data_cluster_info_file, k, b_v, model_version, output_root):
	chi_squared_files, study_sample_sizes = extract_study_info(studies_file)
	pairwise_ld_indices_files, pairwise_ld_files = extract_pairwise_ld_files(pairwise_ld_summary_file)
	cluster_ids, cluster_ukbb_files, cluster_pairwise_ld_matrix_files, cluster_variant_names_files, cluster_variant_neighbor_positions_files = extract_cluster_info_files(training_data_cluster_info_file)
	num_snps = extract_number_of_snps(chi_squared_files[0])
	if model_version == 'als':
		usldsc = usldsc_als.USLDSC(K=k)
		usldsc.fit(chi_squared_files=chi_squared_files, study_sample_sizes=study_sample_sizes, pairwise_ld_files=pairwise_ld_files, pairwise_ld_indices_files=pairwise_ld_indices_files, num_snps=num_snps, cluster_ukbb_files=cluster_ukbb_files, cluster_pairwise_ld_matrix_files=cluster_pairwise_ld_matrix_files, cluster_variant_names_files=cluster_variant_names_files, cluster_variant_neighbor_positions_files=cluster_variant_neighbor_positions_files)
	elif model_version == 'vi':
		usldsc = usldsc_vi.USLDSC(K=k, b_v=b_v)
		usldsc.fit(chi_squared_files=chi_squared_files, study_sample_sizes=study_sample_sizes, pairwise_ld_files=pairwise_ld_files, pairwise_ld_indices_files=pairwise_ld_indices_files, num_snps=num_snps, cluster_ukbb_files=cluster_ukbb_files, cluster_pairwise_ld_matrix_files=cluster_pairwise_ld_matrix_files, cluster_variant_names_files=cluster_variant_names_files, cluster_variant_neighbor_positions_files=cluster_variant_neighbor_positions_files, output_root=output_root)
	elif model_version == 'shared_factor_vi':
		usldsc = usldsc_shared_factor_vi.USLDSC(K=k, b_v=b_v)
		usldsc.fit(chi_squared_files=chi_squared_files, study_sample_sizes=study_sample_sizes, pairwise_ld_files=pairwise_ld_files, pairwise_ld_indices_files=pairwise_ld_indices_files, num_snps=num_snps, cluster_ukbb_files=cluster_ukbb_files, cluster_pairwise_ld_matrix_files=cluster_pairwise_ld_matrix_files, cluster_variant_names_files=cluster_variant_names_files, cluster_variant_neighbor_positions_files=cluster_variant_neighbor_positions_files, output_root=output_root)	
	elif model_version == 'study_variance_vi':
		usldsc = usldsc_study_variance_vi.USLDSC(K=k, b_v=b_v)
		usldsc.fit(chi_squared_files=chi_squared_files, study_sample_sizes=study_sample_sizes, pairwise_ld_files=pairwise_ld_files, pairwise_ld_indices_files=pairwise_ld_indices_files, num_snps=num_snps, cluster_ukbb_files=cluster_ukbb_files, cluster_pairwise_ld_matrix_files=cluster_pairwise_ld_matrix_files, cluster_variant_names_files=cluster_variant_names_files, cluster_variant_neighbor_positions_files=cluster_variant_neighbor_positions_files, output_root=output_root)

	pdb.set_trace()









###################
# Command line args
###################
studies_file = sys.argv[1]
pairwise_ld_summary_file = sys.argv[2]
training_data_cluster_info_file = sys.argv[3]
k = int(sys.argv[4])
model_version = sys.argv[5]
output_root = sys.argv[6]
b_v = float(sys.argv[7])



run_model(studies_file, pairwise_ld_summary_file, training_data_cluster_info_file, k, b_v, model_version, output_root)


