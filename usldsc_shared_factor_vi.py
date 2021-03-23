import numpy as np 
import os
import sys
import pdb
import sklearn.decomposition
from sklearn.linear_model import Ridge
from sklearn.linear_model import Lasso
from sklearn.linear_model import ElasticNet
import scipy.special as special
import time
from sklearn.linear_model import LinearRegression


def sigmoid_function(x):
	return 1.0/(1.0 + np.exp(-x))

def generate_ld_scores(U, pairwise_ld_files, pairwise_ld_indices_files):
	ld_scores = np.zeros(U.shape)
	global_snp_counter = 0
	# We seperated ld files by chromosome to reduce memory burden
	num_files = len(pairwise_ld_files)
	# Loop through ld files
	for file_num in range(num_files):
		# Load in data
		pairwise_ld_file = pairwise_ld_files[file_num]
		pairwise_ld_indices_file = pairwise_ld_indices_files[file_num]
		pairwise_ld = np.load(pairwise_ld_file)
		pairwise_ld_indices = np.load(pairwise_ld_indices_file)
		# Simple error checking to make sure reasonably sized
		if len(pairwise_ld) != len(pairwise_ld_indices):
			print('assumption error')
			pdb.set_trace()
		# Loop through variants on this chromosome
		for variant_iter in range(len(pairwise_ld)):
			variant_pairwise_ld = pairwise_ld[variant_iter]
			variant_pairwise_ld_indices = pairwise_ld_indices[variant_iter]
			# Compute annotation (U) weighted ld scores for each variant
			ld_scores[global_snp_counter, :] = np.dot(variant_pairwise_ld, U[variant_pairwise_ld_indices,:])
			global_snp_counter = global_snp_counter + 1
	return ld_scores


def generate_ld_scores_and_squared_ld_scores(U, U_squared, pairwise_ld_files, pairwise_ld_indices_files):
	ld_scores = np.zeros(U.shape)
	squared_ld_scores = np.zeros(U.shape)
	global_snp_counter = 0
	# We seperated ld files by chromosome to reduce memory burden
	num_files = len(pairwise_ld_files)
	# Loop through ld files
	for file_num in range(num_files):
		# Load in data
		pairwise_ld_file = pairwise_ld_files[file_num]
		pairwise_ld_indices_file = pairwise_ld_indices_files[file_num]
		pairwise_ld = np.load(pairwise_ld_file)
		pairwise_ld_indices = np.load(pairwise_ld_indices_file)
		# Simple error checking to make sure reasonably sized
		if len(pairwise_ld) != len(pairwise_ld_indices):
			print('assumption error')
			pdb.set_trace()
		# Loop through variants on this chromosome
		for variant_iter in range(len(pairwise_ld)):
			variant_pairwise_ld = pairwise_ld[variant_iter]
			variant_pairwise_ld_indices = pairwise_ld_indices[variant_iter]
			# Compute annotation (U) weighted ld scores for each variant
			ld_scores[global_snp_counter, :] = np.dot(variant_pairwise_ld, U[variant_pairwise_ld_indices,:])
			squared_ld_scores[global_snp_counter, :] = np.square(ld_scores[global_snp_counter, :]) - np.dot(np.square(variant_pairwise_ld), np.square(U[variant_pairwise_ld_indices,:])) + np.dot(np.square(variant_pairwise_ld), U_squared[variant_pairwise_ld_indices,:])
			global_snp_counter = global_snp_counter + 1
	return ld_scores, squared_ld_scores



class USLDSC(object):
	def __init__(self, K=5, a_u=1.0, b_u=1.0, a_v=1.0, b_v=1.0, alpha=1e-3, beta=1e-3, max_iter=1000):
		self.K = K
		self.max_iter = max_iter
		self.a_u_prior = a_u
		self.b_u_prior = b_u
		self.a_v_prior = a_v
		self.b_v_prior = b_v
		self.alpha_prior = alpha
		self.beta_prior = beta
	def fit(self, chi_squared_files, study_sample_sizes, pairwise_ld_files, pairwise_ld_indices_files, num_snps, cluster_ukbb_files, cluster_pairwise_ld_matrix_files, cluster_variant_names_files, cluster_variant_neighbor_positions_files, output_root):
		# Load in data
		self.chi_squared_files = chi_squared_files
		self.study_sample_sizes = study_sample_sizes
		self.study_sample_sizes = self.study_sample_sizes/np.mean(self.study_sample_sizes)
		self.num_studies = len(study_sample_sizes)
		self.pairwise_ld_files = pairwise_ld_files
		self.pairwise_ld_indices_files = pairwise_ld_indices_files
		self.num_snps = num_snps
		self.cluster_ukbb_files = cluster_ukbb_files
		self.cluster_pairwise_ld_matrix_files = cluster_pairwise_ld_matrix_files
		self.cluster_variant_names_files = cluster_variant_names_files
		self.cluster_variant_neighbor_positions_files = cluster_variant_neighbor_positions_files
		self.num_snp_clusters = len(cluster_variant_names_files) 
		self.output_root = output_root
		# Initialize variables
		self.initialize_variables()
		# Run iterative optimization
		for vi_iter in range(self.max_iter):
			print('ITERATION: ' + str(vi_iter) + '  ' + str(time.time()))
			# Perform Updates
			self.update_intercept()
			self.update_V()
			print('update F here')
			self.update_F()
			self.update_U()
			self.update_theta_V()
			self.update_theta_U()
			self.update_gamma_U()
			self.update_gamma_V()
			self.update_tau()
			# Save results
			np.save(self.output_root + 'U.npy', self.U_mu)
			np.save(self.output_root + 'U_S.npy', self.U_mu*self.S_U)
			np.save(self.output_root + 'V.npy', self.V_mu)
			np.save(self.output_root + 'V_S.npy', self.V_mu*self.S_V)
			np.save(self.output_root + 'theta_U.npy', self.theta_U_a/(self.theta_U_a + self.theta_U_b))
			np.save(self.output_root + 'theta_V.npy', self.theta_V_a/(self.theta_V_a + self.theta_V_b))
			np.save(self.output_root + 'intercept_mu.npy', self.intercept_mu)
			np.save(self.output_root + 'tau.npy', self.tau_alpha/self.tau_beta)
			np.save(self.output_root + 'gamma_U.npy', self.gamma_U_alpha/self.gamma_U_beta)
			np.save(self.output_root + 'gamma_V.npy', self.gamma_V_alpha/self.gamma_V_beta)
			np.savetxt(self.output_root + 'iter.txt', np.asmatrix(vi_iter), fmt="%s", delimiter='\t')

	def update_gamma_U(self):
		U_var_s_0 = self.gamma_U_beta/self.gamma_U_alpha
		U_squared_expected_val = ((np.square(self.U_mu) + self.U_var)*self.S_U) + (1.0-self.S_U)*(U_var_s_0)
		self.gamma_U_alpha = self.alpha_prior + ((self.num_snps*self.K)/2.0)
		self.gamma_U_beta = self.beta_prior + (np.sum(U_squared_expected_val)/2.0)
	def update_gamma_V(self):
		V_var_s_0 = self.gamma_V_beta/self.gamma_V_alpha
		V_squared_expected_val = ((np.square(self.V_mu) + self.V_var)*self.S_V) + (1.0-self.S_V)*(V_var_s_0)
		self.gamma_V_alpha = self.alpha_prior + ((self.num_studies*self.K)/2.0)
		self.gamma_V_beta = self.beta_prior + (np.sum(V_squared_expected_val)/2.0)
	def update_theta_V(self):
		# Loop through factors
		for k in range(self.K):
			self.theta_V_a[k] = self.a_v_prior + np.sum(self.S_V[k, :])
			self.theta_V_b[k] = self.b_v_prior + self.num_studies - np.sum(self.S_V[k, :])

	def update_theta_U(self):
		# Loop through factors
		for k in range(self.K):
			self.theta_U_a[k] = self.a_u_prior + np.sum(self.S_U[:, k])
			self.theta_U_b[k] = self.b_u_prior + self.num_snps - np.sum(self.S_U[:, k])
	def update_tau_study_specific(self):
		# Compute other useful expectations
		V_expected_val = self.V_mu*self.S_V
		V_squared_expected_val = (np.square(self.V_mu) + self.V_var)*self.S_V
		U_expected_val = self.U_mu*self.S_U
		U_squared_expected_val = (np.square(self.U_mu) + self.U_var)*self.S_U
		intercept_expected_val = self.intercept_mu
		intercept_squared_expected_val = np.square(self.intercept_mu) + self.intercept_var
		squared_study_sample_size = np.square(self.study_sample_sizes)

		# Initizlize variables to keep track of variance info
		a_temp = np.zeros(self.num_studies)
		b_temp = np.zeros(self.num_studies)
		# loop through clusters
		for cluster_iter in range(self.num_snp_clusters):
			# For this cluster load in relevent data
			cluster_pairwise_ld_matrix = np.load(self.cluster_pairwise_ld_matrix_files[cluster_iter])
			cluster_variant_names = np.load(self.cluster_variant_names_files[cluster_iter])
			cluster_variant_neighbor_positions = np.load(self.cluster_variant_neighbor_positions_files[cluster_iter])
			# Number of snps assigned to this snp cluster
			num_cluster_snps = len(cluster_variant_names)
			cluster_chi_squared = np.load(self.cluster_ukbb_files[cluster_iter]).reshape((num_cluster_snps, self.num_studies), order='F')

			a_temp = a_temp + ((num_cluster_snps)/2.0) #

			b_temp = b_temp + np.sum(np.square(cluster_chi_squared), axis=0)/2.0 #
			b_temp = b_temp + (num_cluster_snps)/2.0 #
			b_temp = b_temp + (num_cluster_snps*intercept_squared_expected_val*squared_study_sample_size)/2.0 #
			b_temp = b_temp - np.sum(cluster_chi_squared, axis=0) #
			b_temp = b_temp - np.sum(cluster_chi_squared*((self.study_sample_sizes)*intercept_expected_val), axis=0) #

			U_expected_cluster = U_expected_val[cluster_variant_names, :]
			U_squared_expected_cluster = U_squared_expected_val[cluster_variant_names, :]
			expected_ld_scores = np.dot(cluster_pairwise_ld_matrix, U_expected_cluster)
			factor_predictions = np.dot(expected_ld_scores, V_expected_val)

			b_temp = b_temp - np.sum((factor_predictions*cluster_chi_squared)*(self.study_sample_sizes), axis=0) #
			b_temp = b_temp + (self.study_sample_sizes)*intercept_expected_val*num_cluster_snps #
			b_temp = b_temp + np.sum(factor_predictions*self.study_sample_sizes, axis=0)  #
			b_temp = b_temp + np.sum(factor_predictions*(squared_study_sample_size*intercept_expected_val), axis=0)

			b_temp = b_temp + np.sum(np.square(factor_predictions*self.study_sample_sizes), axis=0)/2.0
			b_temp = b_temp - np.sum(np.dot(np.dot(np.square(cluster_pairwise_ld_matrix), np.square(U_expected_cluster)), np.square(V_expected_val))*squared_study_sample_size, axis=0)/2.0
			b_temp = b_temp + np.sum(np.dot(np.dot(np.square(cluster_pairwise_ld_matrix), U_squared_expected_cluster), V_squared_expected_val)*squared_study_sample_size, axis=0)/2.0
		self.tau_alpha = self.alpha_prior + a_temp
		self.tau_beta = self.beta_prior + b_temp

	def update_tau(self):
		# Compute other useful expectations
		V_expected_val = self.V_mu*self.S_V
		V_squared_expected_val = (np.square(self.V_mu) + self.V_var)*self.S_V
		U_expected_val = self.U_mu*self.S_U
		U_squared_expected_val = (np.square(self.U_mu) + self.U_var)*self.S_U
		intercept_expected_val = self.intercept_mu
		intercept_squared_expected_val = np.square(self.intercept_mu) + self.intercept_var
		squared_study_sample_size = np.square(self.study_sample_sizes)

		# Initizlize variables to keep track of variance info
		a_temp = 0
		b_temp = 0
		# loop through clusters
		for cluster_iter in range(self.num_snp_clusters):
			# For this cluster load in relevent data
			cluster_pairwise_ld_matrix = np.load(self.cluster_pairwise_ld_matrix_files[cluster_iter])
			cluster_variant_names = np.load(self.cluster_variant_names_files[cluster_iter])
			cluster_variant_neighbor_positions = np.load(self.cluster_variant_neighbor_positions_files[cluster_iter])
			# Number of snps assigned to this snp cluster
			num_cluster_snps = len(cluster_variant_names)
			cluster_chi_squared = np.load(self.cluster_ukbb_files[cluster_iter]).reshape((num_cluster_snps, self.num_studies), order='F')

			a_temp = a_temp + ((num_cluster_snps*self.num_studies)/2.0)

			b_temp = b_temp + np.sum(np.square(cluster_chi_squared))/2.0
			b_temp = b_temp + (num_cluster_snps*self.num_studies)/2.0
			b_temp = b_temp + (num_cluster_snps*np.sum(intercept_squared_expected_val*squared_study_sample_size))/2.0
			b_temp = b_temp - np.sum(cluster_chi_squared)
			b_temp = b_temp - np.sum(cluster_chi_squared*((self.study_sample_sizes)*intercept_expected_val))

			U_expected_cluster = U_expected_val[cluster_variant_names, :]
			U_squared_expected_cluster = U_squared_expected_val[cluster_variant_names, :]
			expected_ld_scores = np.dot(cluster_pairwise_ld_matrix, U_expected_cluster)
			factor_predictions = np.dot(expected_ld_scores, V_expected_val)

			b_temp = b_temp - np.sum((factor_predictions*cluster_chi_squared)*(self.study_sample_sizes))
			b_temp = b_temp + np.sum((self.study_sample_sizes)*intercept_expected_val)*num_cluster_snps
			b_temp = b_temp + np.sum(factor_predictions*self.study_sample_sizes)
			b_temp = b_temp + np.sum(factor_predictions*(squared_study_sample_size*intercept_expected_val))

			b_temp = b_temp + np.sum(np.square(factor_predictions*self.study_sample_sizes))/2.0
			b_temp = b_temp - np.sum(np.dot(np.dot(np.square(cluster_pairwise_ld_matrix), np.square(U_expected_cluster)), np.square(V_expected_val))*squared_study_sample_size)/2.0
			b_temp = b_temp + np.sum(np.dot(np.dot(np.square(cluster_pairwise_ld_matrix), U_squared_expected_cluster), V_squared_expected_val)*squared_study_sample_size)/2.0
		self.tau_alpha = self.alpha_prior + a_temp
		self.tau_beta = self.beta_prior + b_temp
	def update_F(self):
		# Compute other useful expectations
		tau_expected_val = self.tau_alpha/self.tau_beta
		V_expected_val = self.V_mu*self.S_V
		U_expected_val = self.U_mu*self.S_U
		squared_study_sample_size = np.square(self.study_sample_sizes)
		gamma_F_expected_val = self.gamma_F_alpha/self.gamma_F_beta

		# loop through clusters
		for cluster_iter in range(self.num_snp_clusters):
			# For this cluster load in relevent data
			cluster_pairwise_ld_matrix = np.load(self.cluster_pairwise_ld_matrix_files[cluster_iter])
			cluster_variant_names = np.load(self.cluster_variant_names_files[cluster_iter])
			cluster_variant_neighbor_positions = np.load(self.cluster_variant_neighbor_positions_files[cluster_iter])
			# Number of snps assigned to this snp cluster
			num_cluster_snps = len(cluster_variant_names)
			cluster_chi_squared = np.load(self.cluster_ukbb_files[cluster_iter]).reshape((num_cluster_snps, self.num_studies), order='F')

			predicted_U = (self.U_mu[cluster_variant_names,:]*self.S_U[cluster_variant_names,:])
			ld_scores = np.dot(cluster_pairwise_ld_matrix, predicted_U)

			predicted_F = self.F_mu[cluster_variant_names,:]*self.S_F[cluster_variant_names,:]
			F_ld_scores = np.dot(cluster_pairwise_ld_matrix, predicted_F)
			# Loop through snps
			for snp_iter in np.random.permutation(range(num_cluster_snps)):
				# Extract relevent info from snp
				snp_name = cluster_variant_names[snp_iter]
				snp_neighbor_indices = cluster_variant_neighbor_positions[snp_iter]
				num_neighbor_snps = len(snp_neighbor_indices)
				pdb.set_trace()
			'''	
				sum_neighbor_ld = np.sum(cluster_pairwise_ld_matrix[snp_iter,:])
				sum_neighbor_ld_squared = np.sum(np.square(cluster_pairwise_ld_matrix[snp_iter,:]))
				ld_weighted_chi_squared = np.dot(cluster_pairwise_ld_matrix[snp_iter,:], cluster_chi_squared)
				ld_weighted_F = np.sum((F_ld_scores[:,0])*(cluster_pairwise_ld_matrix[snp_iter,:]))
				
				current_U = (self.U_mu[snp_name,:]*self.S_U[snp_name,:])
				other_snps = other_snps - np.dot(cluster_pairwise_ld_matrix[snp_iter,:, np.newaxis], current_U[np.newaxis,:])


				for kk in range(self.K):
					# This is the update point for U_mk
					a_term = (-gamma_U_expected_val/2.0) - (1.0/2.0)*sum_neighbor_ld_squared*sum_squared_Vs[kk] #
					b_term = np.sum(ld_weighted_chi_squared*V_expected_val[kk,:]*self.study_sample_sizes*tau_expected_val) #
					b_term = b_term - sum_neighbor_ld*sum_Vs[kk] #
					b_term = b_term - ld_weighted_F*np.sum(tau_expected_val*squared_study_sample_size*V_expected_val[kk,:])
					b_term = b_term - sum_neighbor_ld*np.sum(tau_expected_val*squared_study_sample_size*self.intercept_mu*V_expected_val[kk,:]) #
					temp_U = self.U_mu[snp_name, :]*self.S_U[snp_name,:]
					other_components = np.dot(temp_U, V_expected_val) - temp_U[kk]*V_expected_val[kk,:]
					b_term = b_term - sum_neighbor_ld_squared*np.sum(tau_expected_val*squared_study_sample_size*other_components*V_expected_val[kk,:]) #

					b_term = b_term - np.sum(np.dot(other_snps, interaction_Vs[kk])*cluster_pairwise_ld_matrix[snp_iter,:])

					self.U_mu[snp_name, kk] = (-b_term)/(2.0*a_term)
					self.U_var[snp_name, kk] = (-1.0)/(2.0*a_term)
					# Sparsity
					ln_theta_U_expected_val = special.digamma(self.theta_U_a[kk]) - special.digamma(self.theta_U_a[kk]+self.theta_U_b[kk])
					ln_1_minus_theta_U_expected_val = special.digamma(self.theta_U_b[kk]) - special.digamma(self.theta_U_a[kk]+self.theta_U_b[kk])
					z_term = ln_theta_U_expected_val - ln_1_minus_theta_U_expected_val + .5*np.log(gamma_U_expected_val) - .5*np.log(-2.0*a_term) + (np.square(b_term)/(-4.0*a_term))
					self.S_U[snp_name, kk] = sigmoid_function(z_term)
				current_U = (self.U_mu[snp_name,:]*self.S_U[snp_name,:])
				other_snps = other_snps + np.dot(cluster_pairwise_ld_matrix[snp_iter,:, np.newaxis], current_U[np.newaxis,:])
			'''

	def update_U(self):
		# Compute other useful expectations
		tau_expected_val = self.tau_alpha/self.tau_beta
		V_expected_val = self.V_mu*self.S_V
		V_squared_expected_val = (np.square(self.V_mu) + self.V_var)*self.S_V
		squared_study_sample_size = np.square(self.study_sample_sizes)
		gamma_U_expected_val = self.gamma_U_alpha/self.gamma_U_beta
		# Precompute some quantities
		interaction_Vs = []
		sum_squared_Vs = []
		sum_Vs = []
		for kk in range(self.K):
			temp_V = (V_expected_val*V_expected_val[kk,:])
			temp_V[kk,:] = V_squared_expected_val[kk,:]
			interaction_Vs.append(np.sum(temp_V*(squared_study_sample_size*tau_expected_val),axis=1))
			sum_squared_Vs.append(np.sum(tau_expected_val*squared_study_sample_size*V_squared_expected_val[kk,:]))
			sum_Vs.append(np.sum(tau_expected_val*self.study_sample_sizes*V_expected_val[kk,:]))

		# loop through clusters
		for cluster_iter in range(self.num_snp_clusters):
			# For this cluster load in relevent data
			cluster_pairwise_ld_matrix = np.load(self.cluster_pairwise_ld_matrix_files[cluster_iter])
			cluster_variant_names = np.load(self.cluster_variant_names_files[cluster_iter])
			cluster_variant_neighbor_positions = np.load(self.cluster_variant_neighbor_positions_files[cluster_iter])
			# Number of snps assigned to this snp cluster
			num_cluster_snps = len(cluster_variant_names)
			cluster_chi_squared = np.load(self.cluster_ukbb_files[cluster_iter]).reshape((num_cluster_snps, self.num_studies), order='F')

			predicted_U = (self.U_mu[cluster_variant_names,:]*self.S_U[cluster_variant_names,:])
			other_snps = np.dot(cluster_pairwise_ld_matrix, predicted_U)

			predicted_F = self.F_mu[cluster_variant_names,:]*self.S_F[cluster_variant_names,:]
			F_ld_scores = np.dot(cluster_pairwise_ld_matrix, predicted_F)
			# Loop through snps
			for snp_iter in np.random.permutation(range(num_cluster_snps)):
				# Extract relevent info from snp
				snp_name = cluster_variant_names[snp_iter]
				snp_neighbor_indices = cluster_variant_neighbor_positions[snp_iter]
				num_neighbor_snps = len(snp_neighbor_indices)

				sum_neighbor_ld = np.sum(cluster_pairwise_ld_matrix[snp_iter,:])
				sum_neighbor_ld_squared = np.sum(np.square(cluster_pairwise_ld_matrix[snp_iter,:]))
				ld_weighted_chi_squared = np.dot(cluster_pairwise_ld_matrix[snp_iter,:], cluster_chi_squared)
				ld_weighted_F = np.sum((F_ld_scores[:,0])*(cluster_pairwise_ld_matrix[snp_iter,:]))
				
				current_U = (self.U_mu[snp_name,:]*self.S_U[snp_name,:])
				other_snps = other_snps - np.dot(cluster_pairwise_ld_matrix[snp_iter,:, np.newaxis], current_U[np.newaxis,:])


				for kk in range(self.K):
					# This is the update point for U_mk
					a_term = (-gamma_U_expected_val/2.0) - (1.0/2.0)*sum_neighbor_ld_squared*sum_squared_Vs[kk] #
					b_term = np.sum(ld_weighted_chi_squared*V_expected_val[kk,:]*self.study_sample_sizes*tau_expected_val) #
					b_term = b_term - sum_neighbor_ld*sum_Vs[kk] #
					b_term = b_term - ld_weighted_F*np.sum(tau_expected_val*squared_study_sample_size*V_expected_val[kk,:])
					b_term = b_term - sum_neighbor_ld*np.sum(tau_expected_val*squared_study_sample_size*self.intercept_mu*V_expected_val[kk,:]) #
					temp_U = self.U_mu[snp_name, :]*self.S_U[snp_name,:]
					other_components = np.dot(temp_U, V_expected_val) - temp_U[kk]*V_expected_val[kk,:]
					b_term = b_term - sum_neighbor_ld_squared*np.sum(tau_expected_val*squared_study_sample_size*other_components*V_expected_val[kk,:]) #

					b_term = b_term - np.sum(np.dot(other_snps, interaction_Vs[kk])*cluster_pairwise_ld_matrix[snp_iter,:])

					self.U_mu[snp_name, kk] = (-b_term)/(2.0*a_term)
					self.U_var[snp_name, kk] = (-1.0)/(2.0*a_term)
					# Sparsity
					ln_theta_U_expected_val = special.digamma(self.theta_U_a[kk]) - special.digamma(self.theta_U_a[kk]+self.theta_U_b[kk])
					ln_1_minus_theta_U_expected_val = special.digamma(self.theta_U_b[kk]) - special.digamma(self.theta_U_a[kk]+self.theta_U_b[kk])
					z_term = ln_theta_U_expected_val - ln_1_minus_theta_U_expected_val + .5*np.log(gamma_U_expected_val) - .5*np.log(-2.0*a_term) + (np.square(b_term)/(-4.0*a_term))
					self.S_U[snp_name, kk] = sigmoid_function(z_term)
				current_U = (self.U_mu[snp_name,:]*self.S_U[snp_name,:])
				other_snps = other_snps + np.dot(cluster_pairwise_ld_matrix[snp_iter,:, np.newaxis], current_U[np.newaxis,:])
	def update_intercept(self):
		# Compute other useful expectations
		tau_expected_val = self.tau_alpha/self.tau_beta
		V_expected_val = self.V_mu*self.S_V
		# Calculate LD SCORES
		#ld_scores = generate_ld_scores(self.U_mu*self.S_U, self.pairwise_ld_files, self.pairwise_ld_indices_files)
		#F_ld_scores = generate_ld_scores(self.F_mu*self.S_F, self.pairwise_ld_files, self.pairwise_ld_indices_files)
		shared_ld_scores = generate_ld_scores(np.hstack((self.F_mu*self.S_F, self.U_mu*self.S_U)), self.pairwise_ld_files, self.pairwise_ld_indices_files)
		F_ld_scores = shared_ld_scores[:,0, np.newaxis]
		ld_scores = shared_ld_scores[:,1:]

		# Perform VI updates independently for each study
		for study_num in range(self.num_studies):
			# Extract relevent data for this study
			study_chi_sq = np.load(self.chi_squared_files[study_num])
			study_sample_size = self.study_sample_sizes[study_num]
			# Simple error checking
			if len(study_chi_sq) != ld_scores.shape[0]:
				print('assumption error')
				pdb.set_trace()
			a_term = (-0.0/2.0) - (tau_expected_val/2.0)*np.square(study_sample_size)*self.num_snps
			b_term = (tau_expected_val*study_sample_size*np.sum(study_chi_sq)) - (tau_expected_val*study_sample_size*self.num_snps) - (tau_expected_val*np.square(study_sample_size)*np.sum(np.dot(ld_scores, V_expected_val[:, study_num])))
			b_term = b_term - (tau_expected_val*np.square(study_sample_size)*np.sum(F_ld_scores))
			self.intercept_mu[study_num] = (-b_term)/(2.0*a_term)
			self.intercept_var[study_num] = (-1.0)/(2.0*a_term)


	def update_V(self):
		# Calculate expectation of ld scores and squared ld scores
		ld_scores, squared_ld_scores = generate_ld_scores_and_squared_ld_scores(self.U_mu*self.S_U, (np.square(self.U_mu) + self.U_var)*self.S_U, self.pairwise_ld_files, self.pairwise_ld_indices_files)
		F_ld_scores = generate_ld_scores(self.F_mu*self.S_F, self.pairwise_ld_files, self.pairwise_ld_indices_files)
		# Compute other useful expectations
		tau_expected_val = self.tau_alpha/self.tau_beta
		gamma_V_expected_val = self.gamma_V_alpha/self.gamma_V_beta
		# Perform VI updates independently for each study
		for study_num in range(self.num_studies):
			# Extract relevent data for this study
			study_chi_sq = np.load(self.chi_squared_files[study_num])
			study_sample_size = self.study_sample_sizes[study_num]
			# Simple error checking
			if len(study_chi_sq) != ld_scores.shape[0]:
				print('assumption error')
			# Loop through latent fractors
			for kk in range(self.K):
				# Compute VI updates for this (study, latent factor) pair
				a_term = (-gamma_V_expected_val/2.0) - (tau_expected_val/2.0)*np.square(study_sample_size)*np.sum(squared_ld_scores[:, kk])

				other_components = np.dot(ld_scores, (self.V_mu[:,study_num]*self.S_V[:, study_num])) - ld_scores[:,kk]*(self.V_mu[kk,study_num]*self.S_V[kk, study_num])

				b_term = tau_expected_val*study_sample_size*np.sum(study_chi_sq*ld_scores[:,kk] - ld_scores[:,kk] - study_sample_size*self.intercept_mu[study_num]*ld_scores[:,kk] - study_sample_size*ld_scores[:,kk]*other_components - study_sample_size*ld_scores[:,kk]*F_ld_scores[:,0])
				# Update global model parameters
				self.V_mu[kk, study_num] = (-b_term)/(2.0*a_term)
				self.V_var[kk, study_num] = (-1.0)/(2.0*a_term)
				# Sparsity
				ln_theta_V_expected_val = special.digamma(self.theta_V_a[kk]) - special.digamma(self.theta_V_a[kk]+self.theta_V_b[kk])  # expectation of ln(1-X)
				ln_1_minus_theta_V_expected_val = special.digamma(self.theta_V_b[kk]) - special.digamma(self.theta_V_a[kk]+self.theta_V_b[kk])
				
				z_term = ln_theta_V_expected_val - ln_1_minus_theta_V_expected_val + .5*np.log(gamma_V_expected_val) - .5*np.log(-2.0*a_term) + (np.square(b_term)/(-4.0*a_term))
				self.S_V[kk, study_num] = sigmoid_function(z_term)
	def initialize_variables(self):
		# initialization of V doesn't matter as we learn V on the first step conditioned on U
		self.V_mu = np.zeros((self.K, self.num_studies))
		self.V_mu = np.random.randn(self.K, self.num_studies)
		self.V_var = np.ones((self.K, self.num_studies))
		self.S_V = np.ones((self.K, self.num_studies))
		# Initialization of intercept doesn't matter as we learn intercept in the first step, conditioned on U
		self.intercept_mu = np.zeros(self.num_studies)
		self.intercept_var = np.ones(self.num_studies)
		# Initialization of U DOES matter
		self.U_mu = np.random.randn(self.num_snps, self.K)
		self.U_var = np.ones((self.num_snps, self.K))
		for k in range(self.K):
			self.U_mu[:,k] = ((self.U_mu[:,k]-np.mean(self.U_mu[:,k]))/(10.0*np.std(self.U_mu[:,k])))
		self.S_U = np.ones((self.num_snps, self.K))
		
		self.F_mu = np.random.randn(self.num_snps, 1)
		self.F_var = np.ones((self.num_snps, 1))
		self.F_mu[:,0] = ((self.F_mu[:,0]-np.mean(self.F_mu[:,0]))/(10.0*np.std(self.F_mu[:,0])))
		self.S_F = np.ones((self.num_snps, 1))

		# Smart init for Residual variance
		ld_score = generate_ld_scores(np.ones((self.U_mu.shape[0],1)), self.pairwise_ld_files, self.pairwise_ld_indices_files)
		# Perform regression analysis independently for each study
		resid_varz = []
		for study_num in range(self.num_studies):
			# Extract relevent data for this study
			study_chi_sq = np.load(self.chi_squared_files[study_num])
			study_sample_size = self.study_sample_sizes[study_num]
			# Simple error checking
			if len(study_chi_sq) != ld_score.shape[0]:
				print('assumption error')
			# Fit regression model
			#clf = Lasso(alpha=.1)
			reg = LinearRegression().fit(ld_score*study_sample_size, study_chi_sq-1.0)
			residual_var = np.sum(np.square(reg.predict(ld_score*study_sample_size) - (study_chi_sq-1.0)))/(len(study_chi_sq) -2)
			resid_varz.append(residual_var)
		mean_resid_var = np.mean(resid_varz)
		# Variance params
		# resid var
		#self.tau_alpha = np.ones(self.num_studies)
		#self.tau_beta = np.ones(self.num_studies)
		#self.tau_beta = np.asarray(resid_varz)
		self.tau_alpha = 1.0
		self.tau_beta = mean_resid_var
		# U var
		self.gamma_U_alpha = 1.0
		self.gamma_U_beta = 1.0
		# V var
		self.gamma_V_alpha = 1.0
		self.gamma_V_beta = 1.0
		# F var
		self.gamma_F_alpha = 1.0
		self.gamma_F_beta = 1.0

		# Sparsity parameters
		self.theta_U_a = np.ones(self.K)*(.2)
		self.theta_U_b = np.ones(self.K)
		self.theta_V_a = np.ones(self.K)*(.2)
		self.theta_V_b = np.ones(self.K)
		self.theta_F_a = 1.0*.2
		self.theta_F_b = 1.0

