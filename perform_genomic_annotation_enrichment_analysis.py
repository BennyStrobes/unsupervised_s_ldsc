import numpy as np
import os
import sys
import pdb
import statsmodels.api as sm




def extract_af_bins(allele_frequency_file):
	arr = []
	f = open(allele_frequency_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		arr.append(int(data[2]))
	f.close()
	return np.asarray(arr)

def extract_genomic_annotation_names(genomic_annotation_file):
	names = []
	head_count = 0
	f = open(genomic_annotation_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			names = data[1:]
			continue
		break
	return np.asarray(names)

def extract_genomic_annotation(annotation_iter, genomic_annotation_file):
	f = open(genomic_annotation_file)
	head_count = 0
	anno = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		anno_string = data[(annotation_iter + 1)]
		if anno_string == 'NA':
			anno.append(np.nan)
		else:
			anno.append(float(anno_string))
	return np.asarray(anno)


def annotation_lv_logistic_regression(anno, S_U_binary, af_bins):
	multiplier = 1
	permute = False
	# Get indices where annotation is observed
	observed_indices = np.isnan(anno) == False

	observed_anno = anno[observed_indices]
	observed_S_U_binary = S_U_binary[observed_indices]
	observed_af_bins = af_bins[observed_indices]
	if permute == True:
		num_ones = np.sum(observed_S_U_binary)
		pp = num_ones/len(observed_S_U_binary)
		observed_S_U_binary = np.random.binomial(1, pp, size=len(observed_S_U_binary))

	if np.var(observed_anno) == 0:
		test_info = {'beta': np.nan, 'beta_ub': np.nan, 'beta_lb': np.nan, 'pvalue': np.nan}
		return test_info

	observed_anno = (observed_anno - np.mean(observed_anno))/np.std(observed_anno)

	ys = []
	gs = []
	loaded_snps = np.where(observed_S_U_binary == 1.0)[0]
	ys.append([1]*len(loaded_snps))
	gs.append(observed_anno[loaded_snps])

	bin_counts = {}
	for snp_index in loaded_snps:
		snp_af_bin = observed_af_bins[snp_index]
		if snp_af_bin not in bin_counts:
			bin_counts[snp_af_bin] = 0
		bin_counts[snp_af_bin] = bin_counts[snp_af_bin] + 1

	for bin_num in bin_counts.keys():
		num_in_bin = bin_counts[bin_num]*multiplier
		null_indices = np.where((observed_S_U_binary != 1.0) & (observed_af_bins == bin_num))[0]
		if len(null_indices) < num_in_bin:
			print('assumption eororr')
			pdb.set_trace()
		randomly_sampled_null_indices = np.random.choice(null_indices, size=num_in_bin,replace=False)
		ys.append([0]*len(randomly_sampled_null_indices))
		gs.append(observed_anno[randomly_sampled_null_indices])
	ys = np.hstack(ys)
	gs = np.hstack(gs)
	# Standardize annotations
	# gs = (gs - np.mean(gs))/np.std(gs)
	try:
		model = sm.Logit(ys, sm.add_constant(gs))
		res = model.fit()
		ci = res.conf_int()[1,:]
		beta_lb = ci[0]
		beta_ub = ci[1]
		beta = res.params[1]
		pvalue = res.pvalues[1]
		if res.mle_retvals['converged'] == False:
			beta_lb = np.nan
			beta_ub = np.nan
			beta = np.nan
			pvalue = np.nan
	except:
		pvalue = np.nan
		beta_lb = np.nan
		beta_ub = np.nan
		beta = np.nan
	test_info = {'beta': beta, 'beta_ub': beta_ub, 'beta_lb': beta_lb, 'pvalue': pvalue}
	return test_info


allele_frequency_file = sys.argv[1]
genomic_annotation_file = sys.argv[2]
U_S_npy_file = sys.argv[3]
U_npy_file = sys.argv[4]
enrichment_results_file_stem = sys.argv[5]


# Loaded in U-LDSC S_U variable (ie bernoulli prob that U is loaded)
S_U = np.load(U_S_npy_file)/np.load(U_npy_file)
# Binarize S_U
S_U_binary = 1.0*(S_U > .5)

K = S_U_binary.shape[1]


# Extract AF bins
af_bins = extract_af_bins(allele_frequency_file)

# Extract genomic annotation data
annotation_names = extract_genomic_annotation_names(genomic_annotation_file)


t = open(enrichment_results_file_stem + 'logstic_regression_results.txt', 'w')
t.write('annotation\tlatent_factor\tpvalue\tbeta\tbeta_lb\tbeta_ub\n')

for annotation_iter, annotation_name in enumerate(annotation_names):
	anno = extract_genomic_annotation(annotation_iter, genomic_annotation_file)
	for lv_num in range(K):
		test_results = annotation_lv_logistic_regression(anno, S_U_binary[:, lv_num], af_bins)
		t.write(annotation_name + '\t' + str(lv_num) + '\t' + str(test_results['pvalue']) + '\t' + str(test_results['beta']) + '\t' + str(test_results['beta_lb']) + '\t' + str(test_results['beta_ub']) + '\n')

t.close()
