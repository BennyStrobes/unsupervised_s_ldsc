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

def load_in_chi_squared_data(chi_squared_files):
	chi_squared = []
	for chi_squared_file in chi_squared_files:
		chi_squared.append(np.load(chi_squared_file))
	return np.transpose(np.asarray(chi_squared))

class USLDSC(object):
	def __init__(self, K=5, pruning_thresh=.2):
		self.K = K
		self.pruning_thresh = pruning_thresh
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

		self.chi_squared = load_in_chi_squared_data(chi_squared_files)
		self.learn_pruned_variant_indices()
		self.chi_squared_pruned = self.chi_squared[self.pruned_snps, :]
		self.standardized_chi_squared_pruned = (self.chi_squared_pruned/self.study_sample_sizes)
		for study_num in range(self.num_studies):
			self.standardized_chi_squared_pruned[:, study_num] = (self.standardized_chi_squared_pruned[:, study_num] - np.mean(self.standardized_chi_squared_pruned[:, study_num]))/np.std(self.standardized_chi_squared_pruned[:, study_num])

		uuu, sss, vh = np.linalg.svd(np.transpose(self.standardized_chi_squared_pruned), full_matrices=False)
		svd_loadings = np.transpose(vh)[:,:self.K]

		V1 = uuu[:self.K, :]
		V2 = np.transpose(uuu[:, :self.K])
		np.save(self.output_root + 'V.npy', V2)
	def learn_pruned_variant_indices(self):
		pruned_snps = []
		for cluster_iter in range(self.num_snp_clusters):
			# For this cluster load in relevent data
			cluster_pairwise_ld_matrix = np.load(self.cluster_pairwise_ld_matrix_files[cluster_iter])
			cluster_variant_names = np.load(self.cluster_variant_names_files[cluster_iter])
			cluster_variant_neighbor_positions = np.load(self.cluster_variant_neighbor_positions_files[cluster_iter])
			# Number of snps assigned to this snp cluster
			num_cluster_snps = len(cluster_variant_names)

			valid_cluster_snps = np.arange(num_cluster_snps)
			pruned_cluster_snps = []
			converged = False
			while converged == False:
				new_snp = np.random.choice(valid_cluster_snps)
				pruned_cluster_snps.append(new_snp)
				new_invalid_snps = cluster_pairwise_ld_matrix[new_snp,:] > self.pruning_thresh
				new_valid_cluster_snps = []
				for snp in valid_cluster_snps:
					if new_invalid_snps[snp] == False:
						new_valid_cluster_snps.append(snp)
				valid_cluster_snps = np.copy(new_valid_cluster_snps)
				if len(valid_cluster_snps) == 0:
					converged = True
			pruned_cluster_snps = np.asarray(pruned_cluster_snps)
			for pruned_cluster_snp in pruned_cluster_snps:
				pruned_snps.append(cluster_variant_names[pruned_cluster_snp])

		self.pruned_snps = np.asarray(pruned_snps)

