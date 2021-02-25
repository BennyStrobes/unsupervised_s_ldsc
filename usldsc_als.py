import numpy as np 
import os
import sys
import pdb
import sklearn.decomposition
from sklearn.linear_model import Ridge
from sklearn.linear_model import Lasso
from sklearn.linear_model import ElasticNet




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




class USLDSC(object):
	def __init__(self, K=5, max_iter=1000):
		self.K = K
		self.max_iter = max_iter
	def fit(self, chi_squared_files, study_sample_sizes, pairwise_ld_files, pairwise_ld_indices_files, num_snps, cluster_ukbb_files, cluster_pairwise_ld_matrix_files, cluster_variant_names_files, cluster_variant_neighbor_positions_files):
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
		# Initialize variables
		self.initialize_variables()
		# Run iterative optimization
		for als_iter in range(self.max_iter):
			print('ITERATION: ' + str(als_iter))
			self.update_V_and_intercept()
			self.update_U()
			np.save('temp_U.npy', self.U)
			np.save('temp_V.npy', self.V)
			np.save('temp_intercept.npy', self.intercept)
			np.save('temp_itera.npy', als_iter)
	def update_U(self):
		# loop through clusters
		for cluster_iter in range(self.num_snp_clusters):
			# For this cluster load in relevent data
			cluster_pairwise_ld_matrix = np.load(self.cluster_pairwise_ld_matrix_files[cluster_iter])
			cluster_variant_names = np.load(self.cluster_variant_names_files[cluster_iter])
			cluster_variant_neighbor_positions = np.load(self.cluster_variant_neighbor_positions_files[cluster_iter])
			# Number of snps assigned to this snp cluster
			num_cluster_snps = len(cluster_variant_names)
			cluster_chi_squared = np.load(self.cluster_ukbb_files[cluster_iter]).reshape((num_cluster_snps, self.num_studies), order='F')
			print(str(cluster_iter) + ': ' + str(num_cluster_snps))

			# Generate predicted chi squared values in this cluster
			predicted_cluster_chi_squared = np.matmul(np.matmul(cluster_pairwise_ld_matrix, self.U[cluster_variant_names,:]), self.V)*self.study_sample_sizes + 1.0 + self.intercept
			# Loop through snps
			for snp_iter in range(num_cluster_snps):
				# Extract relevent info from snp
				snp_name = cluster_variant_names[snp_iter]
				snp_neighbor_indices = cluster_variant_neighbor_positions[snp_iter]
				num_neighbor_snps = len(snp_neighbor_indices)
				# QUICK ERROR CHECK (TO BE REMOVED)
				if np.array_equal(sorted(np.where(cluster_pairwise_ld_matrix[snp_iter,:] != 0.0)[0]), sorted(snp_neighbor_indices)) == False:
					print('assumption error')
				# Remove this snps predicted effect from predicted cluster chi squared
				predicted_cluster_chi_squared = predicted_cluster_chi_squared - np.matmul(np.matmul(cluster_pairwise_ld_matrix[:, snp_iter, np.newaxis], self.U[np.newaxis,snp_name,:]), self.V)*self.study_sample_sizes
				residual_chi_squared = (cluster_chi_squared[snp_neighbor_indices, :] - predicted_cluster_chi_squared[snp_neighbor_indices,:]).flatten(order='F')
				# Estimate predicted effect from this snp
				X = np.tile(cluster_pairwise_ld_matrix[snp_iter,:][snp_neighbor_indices,np.newaxis], (self.num_studies,self.K))
				# Generate vector of study sample sizes
				repeated_study_sample_sizes = np.repeat(self.study_sample_sizes,num_neighbor_snps)
				# Add study sample sizes to X
				X = X*repeated_study_sample_sizes[:, np.newaxis]
				X = X*np.repeat(np.transpose(self.V), repeats=num_neighbor_snps,axis=0)
				clf = Ridge(alpha=(1e-2)*num_neighbor_snps, fit_intercept=False)
				clf.fit(X, residual_chi_squared)
				# Add predicted from this new estimated snp to predicted cluster chi squared
				self.U[snp_name,:] = clf.coef_
				predicted_cluster_chi_squared = predicted_cluster_chi_squared + np.matmul(np.matmul(cluster_pairwise_ld_matrix[:, snp_iter, np.newaxis], self.U[np.newaxis,snp_name,:]), self.V)*self.study_sample_sizes
			

	def update_U_snps_jointly(self):
		# loop through clusters
		print(self.num_snp_clusters)
		for cluster_iter in range(self.num_snp_clusters):
			# For this cluster load in relevent data
			cluster_pairwise_ld_matrix = np.load(self.cluster_pairwise_ld_matrix_files[cluster_iter])
			cluster_variant_names = np.load(self.cluster_variant_names_files[cluster_iter])
			# Number of snps assigned to this snp cluster
			num_cluster_snps = len(cluster_variant_names)
			cluster_chi_squared = np.load(self.cluster_ukbb_files[cluster_iter])
			print(str(cluster_iter) + ': ' + str(num_cluster_snps))
				
			# Generate full coeficent matrix
			X = np.tile(cluster_pairwise_ld_matrix, (self.num_studies, self.K))
			# Generate vector of study sample sizes
			repeated_study_sample_sizes = np.repeat(self.study_sample_sizes,num_cluster_snps)
			# Add study sample sizes to X
			X = X*repeated_study_sample_sizes[:, np.newaxis]
			# Add V to X
			X = X*np.repeat(np.repeat(np.transpose(self.V), repeats=num_cluster_snps,axis=1), repeats=num_cluster_snps,axis=0)
			# Generate vector of intercepts
			repeated_intercepts = np.repeat(self.intercept, num_cluster_snps)
			# Remove intercepts
			cluster_chi_squared = cluster_chi_squared - 1.0 - repeated_intercepts
			#clf = Lasso(alpha=.001, fit_intercept=False)
			clf = Ridge(alpha=(1e-4)*num_cluster_snps, fit_intercept=False)
			clf.fit(X, cluster_chi_squared)
			self.U[cluster_variant_names,:] = clf.coef_.reshape((num_cluster_snps,self.K), order='F')

	def update_V_and_intercept(self):
		ld_scores = generate_ld_scores(self.U, self.pairwise_ld_files, self.pairwise_ld_indices_files)
		# Perform regression analysis independently for each study
		for study_num in range(self.num_studies):
			# Extract relevent data for this study
			study_chi_sq = np.load(self.chi_squared_files[study_num])
			study_sample_size = self.study_sample_sizes[study_num]
			# Simple error checking
			if len(study_chi_sq) != ld_scores.shape[0]:
				print('assumption error')
			# Fit regression model
			clf = Lasso(alpha=.1)
			clf.fit(ld_scores*study_sample_size, study_chi_sq-1.0)
			# Update parameters with MAP estimates from regression
			self.V[:, study_num] = clf.coef_
			self.intercept[study_num] = clf.intercept_
			'''	
			clf = Ridge(alpha=1.0)
			clf.fit(ld_scores*study_sample_size, study_chi_sq-1.0)
			# Update parameters with MAP estimates from regression
			self.V[:, study_num] = clf.coef_
			self.intercept[study_num] = clf.intercept_
			'''
	def initialize_variables(self):
		# initialization of V doesn't matter as we learn V on the first step conditioned on U
		self.V = np.zeros((self.K, self.num_studies))
		# Initialization of intercept doesn't matter as we learn intercept in the first step, conditioned on U
		self.intercept = np.zeros(self.num_studies)
		# Initialization of U DOES matter
		self.U = np.random.randn(self.num_snps, self.K)
		for k in range(self.K):
			self.U[:,k] = ((self.U[:,k]-np.mean(self.U[:,k]))/(np.std(self.U[:,k])))
