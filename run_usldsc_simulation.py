import numpy as np 
import os
import sys
import pdb
import usldsc_debug_vi
import usldsc_debug_fixed_U_vi
import usldsc_debug_ard_vi
import usldsc_debug2_vi
import usldsc_debug_vi3

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

def simulate_data(training_data_study_file, training_data_pairwise_ld_file, training_data_cluster_info_file, simulation_k, simulation_output_root):
	# Load in study sample sizes
	temp_study_data = np.loadtxt(training_data_study_file,dtype=str, delimiter='\t')
	study_names = temp_study_data[1:,0]  # Will remain same
	study_sample_sizes = temp_study_data[1:,2].astype(float)  # Will remain same
	normalized_study_sample_sizes = study_sample_sizes/np.mean(study_sample_sizes)
	# Extract single cluster with  > 2000 snps in it
	f = open(training_data_cluster_info_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		cluster_indices = np.load(data[3])
		if len(cluster_indices) > 2000 and len(cluster_indices) < 3000:
			cluster_id_real = data[0]  # Will remain same
			cluster_chi_squared_npy_file_real = data[1]
			cluster_correlation_mat_npy_file_real = data[2] # Will remain same
			cluster_indices_npy_file_real = data[3] 
			cluster_neighbor_pos_npy_file_real = data[4] # will remain same
			break
	f.close()
	# Load in some data
	corr_mat = np.load(cluster_correlation_mat_npy_file_real)

	# Simulate data
	num_snps = corr_mat.shape[0]
	num_studies = len(study_names)
	U_sim = np.random.randn(num_snps, simulation_k)
	S_U_sim = np.random.randint(2, size=(num_snps, simulation_k))
	V_sim = np.random.randn(simulation_k, num_studies)
	S_V_sim = np.random.randint(2, size=(simulation_k, num_studies))
	intercept_sim = np.random.randn(num_studies) 
	intercept_sim = intercept_sim - np.min(intercept_sim)
	tau_sim = 1.5
	ld_scores_sim = np.dot(corr_mat, U_sim*S_U_sim)
	factor_predicted_mean_sim = np.dot(ld_scores_sim, V_sim*S_V_sim)
	chi_squared_sim = np.zeros((num_snps, num_studies))

	for study_num in range(num_studies):
		study_sample_size = normalized_study_sample_sizes[study_num]

		study_means = 1.0 + study_sample_size*intercept_sim[study_num] + study_sample_size*factor_predicted_mean_sim[:, study_num]
		chi_squared_sim[:, study_num] = np.random.normal(study_means, scale=np.sqrt(1.0/tau_sim))

	simulation_data = {}
	simulation_data['U_sim'] = U_sim
	simulation_data['S_U_sim'] = S_U_sim
	simulation_data['V_sim'] = V_sim
	simulation_data['S_V_sim'] = S_V_sim
	simulation_data['intercept_sim'] = intercept_sim
	simulation_data['tau_sim'] = tau_sim
	simulation_data['corr_mat'] = corr_mat

	# print training data study file
	f = open(training_data_study_file)
	simulated_training_data_study_file = simulation_output_root + 'training_data_study_file.txt'
	t = open(simulated_training_data_study_file, 'w')
	study_num = 0
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		study_name = data[0]
		study_npy_file = simulation_output_root + study_name + '_chi_squared_stats.npy'
		np.save(study_npy_file, chi_squared_sim[:, study_num])
		study_num = study_num + 1
		t.write(data[0] + '\t' + study_npy_file + '\t' + data[2] + '\n')
	f.close()
	t.close()
	# Print training data cluster info file
	simulated_training_data_cluster_info_file = simulation_output_root + 'training_data_cluster_info_file.txt'
	t = open(simulated_training_data_cluster_info_file, 'w')
	cluster_chi_squared_npy_file_sim = simulation_output_root + cluster_id_real + '_cluster_chi_squared.npy'
	np.save(cluster_chi_squared_npy_file_sim, chi_squared_sim.flatten(order='F'))
	cluster_variant_names_npy_file_sim = simulation_output_root + cluster_id_real + '_cluster_variant_names.npy'
	np.save(cluster_variant_names_npy_file_sim, np.arange(num_snps))
	t.write('cluster_id\tcluster_ukbb_file\tcluster_pairwise_ld_matrix_file\tcluster_variant_names_file\tcluster_pairwise_ld_neighbor_positions_file\n')
	t.write(cluster_id_real + '\t' + cluster_chi_squared_npy_file_sim + '\t' + cluster_correlation_mat_npy_file_real + '\t' + cluster_variant_names_npy_file_sim + '\t' + cluster_neighbor_pos_npy_file_real + '\n')
	t.close()

	# Extract relevent files from pairwise ld file
	temp_data = np.loadtxt(training_data_pairwise_ld_file, dtype=str,delimiter='\t')
	pairwise_indices_npy_file = temp_data[1,1]
	pairwise_ld_npy_file = temp_data[1,2]

	# Print training data pairwise ld file
	sim_pairwise_indices = []
	sim_pairwise_ld = []
	for snp_num in range(num_snps):
		neighbor_indices = np.where(corr_mat[snp_num,:] > 0.0)[0]
		neighbor_ld = corr_mat[snp_num, neighbor_indices]
		sim_pairwise_indices.append(neighbor_indices)
		sim_pairwise_ld.append(neighbor_ld)
	sim_pairwise_indices_npy_file = simulation_output_root + 'pairwise_indices.npy'
	np.save(sim_pairwise_indices_npy_file, sim_pairwise_indices)
	sim_pairwise_ld_npy_file = simulation_output_root + 'pairwise_ld.npy'
	np.save(sim_pairwise_ld_npy_file, sim_pairwise_ld)

	simulated_training_data_pairwise_ld_file = simulation_output_root + 'training_data_pairwise_ld.txt'
	t = open(simulated_training_data_pairwise_ld_file,'w')
	t.write('\t'.join(temp_data[0,:]) + '\n')
	t.write('4\t' + sim_pairwise_indices_npy_file + '\t' + sim_pairwise_ld_npy_file + '\n')
	t.close()

	return simulation_data, simulated_training_data_study_file, simulated_training_data_pairwise_ld_file, simulated_training_data_cluster_info_file



def simulate_independent_data(training_data_study_file, training_data_pairwise_ld_file, training_data_cluster_info_file, simulation_k, simulation_output_root):
	# Load in study sample sizes
	temp_study_data = np.loadtxt(training_data_study_file,dtype=str, delimiter='\t')
	study_names = temp_study_data[1:,0]  # Will remain same
	study_sample_sizes = temp_study_data[1:,2].astype(float)  # Will remain same
	normalized_study_sample_sizes = study_sample_sizes/np.mean(study_sample_sizes)
	# Extract single cluster with  > 2000 snps in it
	f = open(training_data_cluster_info_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		cluster_indices = np.load(data[3])
		if len(cluster_indices) > 2000 and len(cluster_indices) < 3000:
			cluster_id_real = data[0]  # Will remain same
			cluster_chi_squared_npy_file_real = data[1]
			cluster_correlation_mat_npy_file_real = data[2] # Will remain same
			cluster_indices_npy_file_real = data[3] 
			cluster_neighbor_pos_npy_file_real = data[4] # will remain same
			break
	f.close()
	np.random.seed(1)
	# Load in some data
	corr_mat_orig = np.load(cluster_correlation_mat_npy_file_real)
	corr_mat = np.eye(corr_mat_orig.shape[0])

	# Simulate data
	num_snps = corr_mat.shape[0]
	num_studies = len(study_names)
	U_sim = np.random.randn(num_snps, simulation_k)
	S_U_sim = np.random.randint(2, size=(num_snps, simulation_k))
	theta_U = []
	for k in range(simulation_k):
		p_k = np.random.beta(2, 2)
		theta_U.append(p_k)
		S_U_sim[:, k] = np.random.binomial(1, p=p_k, size=S_U_sim.shape[0])

	#S_U_sim = np.ones((num_snps, simulation_k))
	V_sim = np.random.normal(0, .85, size=(simulation_k, num_studies))
	#V_sim = np.ones((simulation_k, num_studies))*3.0
	S_V_sim = np.random.randint(2, size=(simulation_k, num_studies))

	theta_V = []
	for k in range(simulation_k):
		p_k = np.random.beta(2,2)
		theta_V.append(p_k)
		S_V_sim[k, :] = np.random.binomial(1, p=p_k, size=S_V_sim.shape[1])


	intercept_sim = np.random.randn(num_studies) 
	intercept_sim = intercept_sim - np.min(intercept_sim)
	tau_sim = 1/10.0
	ld_scores_sim = np.dot(corr_mat, U_sim*S_U_sim)
	factor_predicted_mean_sim = np.dot(ld_scores_sim, S_V_sim)
	chi_squared_sim = np.zeros((num_snps, num_studies))

	for study_num in range(num_studies):
		study_sample_size = normalized_study_sample_sizes[study_num]

		study_means = 1.0 + 1.0*factor_predicted_mean_sim[:, study_num]
		chi_squared_sim[:, study_num] = np.random.normal(study_means, scale=np.sqrt(1.0/tau_sim))


	simulation_data = {}
	simulation_data['U_sim'] = U_sim
	simulation_data['S_U_sim'] = S_U_sim
	simulation_data['V_sim'] = V_sim
	simulation_data['S_V_sim'] = S_V_sim
	simulation_data['intercept_sim'] = intercept_sim
	simulation_data['tau_sim'] = tau_sim
	simulation_data['corr_mat'] = corr_mat
	simulation_data['Y'] = chi_squared_sim - 1.0

	# print training data study file
	f = open(training_data_study_file)
	simulated_training_data_study_file = simulation_output_root + 'training_data_study_file.txt'
	t = open(simulated_training_data_study_file, 'w')
	study_num = 0
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		study_name = data[0]
		study_npy_file = simulation_output_root + study_name + '_chi_squared_stats.npy'
		np.save(study_npy_file, chi_squared_sim[:, study_num])
		study_num = study_num + 1
		t.write(data[0] + '\t' + study_npy_file + '\t' + '1' + '\n')
	f.close()
	t.close()
	# Print training data cluster info file
	simulated_training_data_cluster_info_file = simulation_output_root + 'training_data_cluster_info_file.txt'
	t = open(simulated_training_data_cluster_info_file, 'w')
	cluster_chi_squared_npy_file_sim = simulation_output_root + cluster_id_real + '_cluster_chi_squared.npy'
	np.save(cluster_chi_squared_npy_file_sim, chi_squared_sim.flatten(order='F'))
	cluster_variant_names_npy_file_sim = simulation_output_root + cluster_id_real + '_cluster_variant_names.npy'
	np.save(cluster_variant_names_npy_file_sim, np.arange(num_snps))
	cluster_corr_mat_npy_file_sim = simulation_output_root + cluster_id_real + '_cluster_corr_mat.npy'
	np.save(cluster_corr_mat_npy_file_sim, corr_mat)
	temp_arr = []
	for i in range(num_snps):
		temp_arr.append(np.asarray([i]))
	cluster_neighbor_pos_npy_file_sim = simulation_output_root + cluster_id_real + '_cluster_neighbor_position.npy'
	np.save(cluster_neighbor_pos_npy_file_sim, temp_arr)

	t.write('cluster_id\tcluster_ukbb_file\tcluster_pairwise_ld_matrix_file\tcluster_variant_names_file\tcluster_pairwise_ld_neighbor_positions_file\n')
	t.write(cluster_id_real + '\t' + cluster_chi_squared_npy_file_sim + '\t' + cluster_corr_mat_npy_file_sim + '\t' + cluster_variant_names_npy_file_sim + '\t' + cluster_neighbor_pos_npy_file_sim + '\n')
	t.close()

	# Extract relevent files from pairwise ld file
	temp_data = np.loadtxt(training_data_pairwise_ld_file, dtype=str,delimiter='\t')
	pairwise_indices_npy_file = temp_data[1,1]
	pairwise_ld_npy_file = temp_data[1,2]

	# Print training data pairwise ld file
	sim_pairwise_indices = []
	sim_pairwise_ld = []
	for snp_num in range(num_snps):
		neighbor_indices = np.where(corr_mat[snp_num,:] > 0.0)[0]
		neighbor_ld = corr_mat[snp_num, neighbor_indices]
		sim_pairwise_indices.append(neighbor_indices)
		sim_pairwise_ld.append(neighbor_ld)
	sim_pairwise_indices_npy_file = simulation_output_root + 'pairwise_indices.npy'
	np.save(sim_pairwise_indices_npy_file, sim_pairwise_indices)
	sim_pairwise_ld_npy_file = simulation_output_root + 'pairwise_ld.npy'
	np.save(sim_pairwise_ld_npy_file, sim_pairwise_ld)

	simulated_training_data_pairwise_ld_file = simulation_output_root + 'training_data_pairwise_ld.txt'
	t = open(simulated_training_data_pairwise_ld_file,'w')
	t.write('\t'.join(temp_data[0,:]) + '\n')
	t.write('4\t' + sim_pairwise_indices_npy_file + '\t' + sim_pairwise_ld_npy_file + '\n')
	t.close()

	return simulation_data, simulated_training_data_study_file, simulated_training_data_pairwise_ld_file, simulated_training_data_cluster_info_file



def simulate_non_independent_data(training_data_study_file, training_data_pairwise_ld_file, training_data_cluster_info_file, simulation_k, simulation_output_root):
	# Load in study sample sizes
	temp_study_data = np.loadtxt(training_data_study_file,dtype=str, delimiter='\t')
	study_names = temp_study_data[1:,0]  # Will remain same
	study_sample_sizes = temp_study_data[1:,2].astype(float)  # Will remain same
	normalized_study_sample_sizes = study_sample_sizes/np.mean(study_sample_sizes)
	# Extract single cluster with  > 2000 snps in it
	f = open(training_data_cluster_info_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		cluster_indices = np.load(data[3])
		if len(cluster_indices) > 2000 and len(cluster_indices) < 3000:
			cluster_id_real = data[0]  # Will remain same
			cluster_chi_squared_npy_file_real = data[1]
			cluster_correlation_mat_npy_file_real = data[2] # Will remain same
			cluster_indices_npy_file_real = data[3] 
			cluster_neighbor_pos_npy_file_real = data[4] # will remain same
			break
	f.close()
	# Load in some data
	corr_mat_orig = np.load(cluster_correlation_mat_npy_file_real)
	corr_mat = np.copy(corr_mat_orig)

	# Simulate data
	num_snps = corr_mat.shape[0]
	num_studies = len(study_names)
	U_sim = np.random.randn(num_snps, simulation_k)
	S_U_sim = np.random.randint(2, size=(num_snps, simulation_k))
	theta_U = []
	for k in range(simulation_k):
		p_k = np.random.beta(2, 2)
		theta_U.append(p_k)
		S_U_sim[:, k] = np.random.binomial(1, p=p_k, size=S_U_sim.shape[0])

	#S_U_sim = np.ones((num_snps, simulation_k))
	V_sim = np.random.randn(simulation_k, num_studies)
	S_V_sim = np.random.randint(2, size=(simulation_k, num_studies))

	theta_V = []
	for k in range(simulation_k):
		p_k = np.random.beta(2,2)
		theta_V.append(p_k)
		S_V_sim[k, :] = np.random.binomial(1, p=p_k, size=S_V_sim.shape[1])


	intercept_sim = np.random.randn(num_studies) 
	intercept_sim = intercept_sim - np.min(intercept_sim)
	tau_sim = 1.0
	ld_scores_sim = np.dot(corr_mat, U_sim*S_U_sim)
	factor_predicted_mean_sim = np.dot(ld_scores_sim, V_sim*S_V_sim)
	chi_squared_sim = np.zeros((num_snps, num_studies))

	for study_num in range(num_studies):
		study_sample_size = normalized_study_sample_sizes[study_num]

		study_means = 1.0 + 1.0*factor_predicted_mean_sim[:, study_num]
		chi_squared_sim[:, study_num] = np.random.normal(study_means, scale=np.sqrt(1.0/tau_sim))


	simulation_data = {}
	simulation_data['U_sim'] = U_sim
	simulation_data['S_U_sim'] = S_U_sim
	simulation_data['V_sim'] = V_sim
	simulation_data['S_V_sim'] = S_V_sim
	simulation_data['intercept_sim'] = intercept_sim
	simulation_data['tau_sim'] = tau_sim
	simulation_data['corr_mat'] = corr_mat
	simulation_data['Y'] = chi_squared_sim - 1.0

	# print training data study file
	f = open(training_data_study_file)
	simulated_training_data_study_file = simulation_output_root + 'training_data_study_file.txt'
	t = open(simulated_training_data_study_file, 'w')
	study_num = 0
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		study_name = data[0]
		study_npy_file = simulation_output_root + study_name + '_chi_squared_stats.npy'
		np.save(study_npy_file, chi_squared_sim[:, study_num])
		study_num = study_num + 1
		t.write(data[0] + '\t' + study_npy_file + '\t' + '1' + '\n')
	f.close()
	t.close()
	# Print training data cluster info file
	simulated_training_data_cluster_info_file = simulation_output_root + 'training_data_cluster_info_file.txt'
	t = open(simulated_training_data_cluster_info_file, 'w')
	cluster_chi_squared_npy_file_sim = simulation_output_root + cluster_id_real + '_cluster_chi_squared.npy'
	np.save(cluster_chi_squared_npy_file_sim, chi_squared_sim.flatten(order='F'))
	cluster_variant_names_npy_file_sim = simulation_output_root + cluster_id_real + '_cluster_variant_names.npy'
	np.save(cluster_variant_names_npy_file_sim, np.arange(num_snps))
	cluster_corr_mat_npy_file_sim = simulation_output_root + cluster_id_real + '_cluster_corr_mat.npy'
	np.save(cluster_corr_mat_npy_file_sim, corr_mat)
	temp_arr = []
	for i in range(num_snps):
		temp_arr.append(np.where(corr_mat[i,:] > 0.0)[0])
	cluster_neighbor_pos_npy_file_sim = simulation_output_root + cluster_id_real + '_cluster_neighbor_position.npy'
	np.save(cluster_neighbor_pos_npy_file_sim, temp_arr)

	t.write('cluster_id\tcluster_ukbb_file\tcluster_pairwise_ld_matrix_file\tcluster_variant_names_file\tcluster_pairwise_ld_neighbor_positions_file\n')
	t.write(cluster_id_real + '\t' + cluster_chi_squared_npy_file_sim + '\t' + cluster_corr_mat_npy_file_sim + '\t' + cluster_variant_names_npy_file_sim + '\t' + cluster_neighbor_pos_npy_file_sim + '\n')
	t.close()

	# Extract relevent files from pairwise ld file
	temp_data = np.loadtxt(training_data_pairwise_ld_file, dtype=str,delimiter='\t')
	pairwise_indices_npy_file = temp_data[1,1]
	pairwise_ld_npy_file = temp_data[1,2]

	# Print training data pairwise ld file
	sim_pairwise_indices = []
	sim_pairwise_ld = []
	for snp_num in range(num_snps):
		neighbor_indices = np.where(corr_mat[snp_num,:] > 0.0)[0]
		neighbor_ld = corr_mat[snp_num, neighbor_indices]
		sim_pairwise_indices.append(neighbor_indices)
		sim_pairwise_ld.append(neighbor_ld)
	sim_pairwise_indices_npy_file = simulation_output_root + 'pairwise_indices.npy'
	np.save(sim_pairwise_indices_npy_file, sim_pairwise_indices)
	sim_pairwise_ld_npy_file = simulation_output_root + 'pairwise_ld.npy'
	np.save(sim_pairwise_ld_npy_file, sim_pairwise_ld)

	simulated_training_data_pairwise_ld_file = simulation_output_root + 'training_data_pairwise_ld.txt'
	t = open(simulated_training_data_pairwise_ld_file,'w')
	t.write('\t'.join(temp_data[0,:]) + '\n')
	t.write('4\t' + sim_pairwise_indices_npy_file + '\t' + sim_pairwise_ld_npy_file + '\n')
	t.close()

	return simulation_data, simulated_training_data_study_file, simulated_training_data_pairwise_ld_file, simulated_training_data_cluster_info_file





training_data_study_file = sys.argv[1]
training_data_pairwise_ld_file = sys.argv[2]
training_data_cluster_info_file = sys.argv[3]
simulation_output_root = sys.argv[4]

simulation_k = 4
#simulation_data, simulated_training_data_study_file, simulated_training_data_pairwise_ld_file, simulated_training_data_cluster_info_file = simulate_data(training_data_study_file, training_data_pairwise_ld_file, training_data_cluster_info_file, simulation_k, simulation_output_root + 'simulation_k_' + str(simulation_k) + '_')
#simulation_data, simulated_training_data_study_file, simulated_training_data_pairwise_ld_file, simulated_training_data_cluster_info_file = simulate_independent_data(training_data_study_file, training_data_pairwise_ld_file, training_data_cluster_info_file, simulation_k, simulation_output_root + 'simulation_k_' + str(simulation_k) + '_')
simulation_data, simulated_training_data_study_file, simulated_training_data_pairwise_ld_file, simulated_training_data_cluster_info_file = simulate_independent_data(training_data_study_file, training_data_pairwise_ld_file, training_data_cluster_info_file, simulation_k, simulation_output_root + 'simulation_k_' + str(simulation_k) + '_')


np.save(simulation_output_root + 'sim_data_Y.npy', simulation_data['Y'])
np.save(simulation_output_root + 'sim_data_S_U.npy', simulation_data['S_U_sim'])
np.save(simulation_output_root + 'sim_data_S_V.npy', simulation_data['S_V_sim'])
np.save(simulation_output_root + 'sim_data_V.npy', simulation_data['V_sim'])
np.save(simulation_output_root + 'sim_data_U.npy', simulation_data['U_sim'])



chi_squared_files, study_sample_sizes = extract_study_info(simulated_training_data_study_file)
pairwise_ld_indices_files, pairwise_ld_files = extract_pairwise_ld_files(simulated_training_data_pairwise_ld_file)
cluster_ids, cluster_ukbb_files, cluster_pairwise_ld_matrix_files, cluster_variant_names_files, cluster_variant_neighbor_positions_files = extract_cluster_info_files(simulated_training_data_cluster_info_file)
num_snps = extract_number_of_snps(chi_squared_files[0])



usldsc = usldsc_debug_vi3.USLDSC(K=4, b_v=1.0)
usldsc.fit(chi_squared_files=chi_squared_files, study_sample_sizes=study_sample_sizes, pairwise_ld_files=pairwise_ld_files, pairwise_ld_indices_files=pairwise_ld_indices_files, num_snps=num_snps, cluster_ukbb_files=cluster_ukbb_files, cluster_pairwise_ld_matrix_files=cluster_pairwise_ld_matrix_files, cluster_variant_names_files=cluster_variant_names_files, cluster_variant_neighbor_positions_files=cluster_variant_neighbor_positions_files, output_root=simulation_output_root)


sim_U = simulation_data['U_sim']*simulation_data['S_U_sim']
pred_U = usldsc.U_mu*usldsc.S_U
sim_ld_scores = np.dot(simulation_data['corr_mat'], sim_U)
pred_ld_scores = np.dot(simulation_data['corr_mat'], pred_U)


print(np.corrcoef(np.dot(sim_U, simulation_data['S_V_sim'])[:,5], np.dot(pred_U, usldsc.S_V)[:,5]))


sim_U = simulation_data['U_sim']*simulation_data['S_U_sim']
pred_U = usldsc.U_mu*usldsc.S_U
sim_ld_scores = np.dot(simulation_data['corr_mat'], sim_U)
pred_ld_scores = np.dot(simulation_data['corr_mat'], pred_U)
sim_V = simulation_data['V_sim']*simulation_data['S_V_sim']
sim_factor_predicted = np.dot(sim_ld_scores, sim_V)
pred_V = usldsc.V_mu
pred_factor_predicted = np.dot(pred_ld_scores, pred_V)
theta_U = usldsc.theta_U_a/(usldsc.theta_U_a + usldsc.theta_U_b)
gamma_V = usldsc.gamma_V_alpha/usldsc.gamma_V_beta
gamma_U = usldsc.gamma_U_alpha/usldsc.gamma_U_beta

print(np.corrcoef(np.dot(sim_U, sim_V)[:,1], np.dot(pred_U, pred_V)[:,1]))
print(np.corrcoef(np.dot(sim_U, sim_V)[:,2], np.dot(pred_U, pred_V)[:,2]))
print(np.corrcoef(np.dot(sim_U, sim_V)[:,3], np.dot(pred_U, pred_V)[:,3]))

pdb.set_trace()

