import numpy as np 
import os
import sys
import pdb


def get_fraction_of_loaded_variants_in_ld_with_loaded_variant(S_matrix, pairwise_ld_file, r_squared, component):
	head_count = 0
	S_slice = S_matrix[:, component]
	f = open(pairwise_ld_file)
	fracs = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		pairwise_indices_npy = data[1]
		pairwise_ld_npy = data[2]
		pairwise_indices = np.load(pairwise_indices_npy)
		pairwise_ld = np.load(pairwise_ld_npy)
		for snp_index in range(len(pairwise_ld)):
			if S_slice[snp_index] < .5:
				continue
			num_neighbors = 0
			num_loaded_neighbors = 0
			for i, neighbor_index in enumerate(pairwise_indices[snp_index]):
				if neighbor_index == snp_index:
					continue
				if pairwise_ld[snp_index][i] < r_squared:
					continue
				num_neighbors = num_neighbors + 1
				if S_slice[neighbor_index] > .5:
					num_loaded_neighbors = num_loaded_neighbors + 1
			if num_neighbors > 0:
				fracs.append(num_loaded_neighbors*1.0/num_neighbors)
			else:
				fracs.append(0.0)
	f.close()
	fracs = np.asarray(fracs)
	return fracs

def generate_ld_scores_of_variants(pairwise_ld_file):
	ld_scores = []
	head_count = 0
	f = open(pairwise_ld_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		pairwise_ld_npy = data[2]
		pairwise_ld = np.load(pairwise_ld_npy)
		for snp_index in range(len(pairwise_ld)):
			ld_scores.append(np.sum(pairwise_ld[snp_index]))
	f.close()
	return np.asarray(ld_scores)

model_dir = sys.argv[1]
model_name = sys.argv[2]
pairwise_ld_file = sys.argv[3]
usldsc_visualize_results_dir = sys.argv[4]


model_stem = model_dir + model_name


###########################
# Generate LD Scores of variants
###########################
ld_scores = generate_ld_scores_of_variants(pairwise_ld_file)
np.savetxt(usldsc_visualize_results_dir + model_name + 'ld_scores.txt', ld_scores, delimiter='\t', fmt="%s")


###########################
# Load and save V matrix
###########################
V_matrix = np.load(model_stem + 'V_S.npy')
np.savetxt(usldsc_visualize_results_dir + model_name + 'V_S.txt', V_matrix, delimiter='\t', fmt="%s")

###########################
# Load and save U matrix and S matrix
###########################
U_S_matrix = np.load(model_stem + 'U_S.npy')
U_matrix = np.load(model_stem + 'U.npy')
np.savetxt(usldsc_visualize_results_dir + model_name + 'U_S.txt', U_S_matrix, delimiter='\t', fmt="%s")

S_matrix = U_S_matrix/U_matrix
np.savetxt(usldsc_visualize_results_dir + model_name + 'S_U.txt', S_matrix, delimiter='\t', fmt="%s")


r_squared = .5
for component_num in range(U_S_matrix.shape[1]):
	fraction_variants = get_fraction_of_loaded_variants_in_ld_with_loaded_variant(S_matrix, pairwise_ld_file, r_squared, component_num)
	np.savetxt(usldsc_visualize_results_dir + model_name + 'fraction_loaded_neighbors_component_' + str(component_num) + '.txt', fraction_variants, delimiter='\t', fmt="%s")

