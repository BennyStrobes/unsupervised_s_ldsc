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
	union_other_anno = []
	avg_other_anno = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		anno_string = data[(annotation_iter + 1)]
		annos_string = data[1:]
		other_anno = []
		for i, annos in enumerate(annos_string):
			if i != annotation_iter:
				other_anno.append(annos)
		if anno_string == 'NA':
			anno.append(np.nan)
			union_other_anno.append(np.nan)
			avg_other_anno.append(np.nan)
		else:
			anno.append(float(anno_string))
			other_anno = np.asarray(other_anno).astype(float)
			union_other_anno.append(np.max(other_anno))
			avg_other_anno.append(np.mean(other_anno))
	return np.asarray(anno), np.asarray(union_other_anno), np.asarray(avg_other_anno)

def extract_shared_genomic_annotation(genomic_annotation_file):
	f = open(genomic_annotation_file)
	head_count = 0
	anno = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		hit = 0.0
		nan_count = 0
		for ele in data[1:]:
			if ele == 'NA':
				nan_count = nan_count + 1
			else:
				if ele == '0.0':
					continue
				elif ele == '1.0':
					hit = 1.0
				else:
					print('assumption error')
					pdb.set_trace()
		if nan_count == 10:
			anno.append(np.nan)
		elif nan_count == 0:
			anno.append(hit)
		else:
			print('assumption error')
			pdb.set_trace()
	return np.asarray(anno)

def annotation_lv_logistic_regression(anno, shared_anno, avg_anno, S_U_binary, af_bins):
	multiplier = 1
	permute = False
	# Get indices where annotation is observed
	observed_indices = np.isnan(anno) == False

	observed_anno = anno[observed_indices]
	observed_shared_anno = shared_anno[observed_indices]
	observed_avg_anno = avg_anno[observed_indices]
	observed_S_U_binary = S_U_binary[observed_indices]
	observed_af_bins = af_bins[observed_indices]


	if np.var(observed_anno) == 0:
		test_info = {'beta': np.nan, 'beta_ub': np.nan, 'beta_lb': np.nan, 'pvalue': np.nan}
		return test_info

	observed_anno = (observed_anno - np.mean(observed_anno))/np.std(observed_anno)
	observed_shared_anno = (observed_shared_anno - np.mean(observed_shared_anno))/np.std(observed_shared_anno)
	observed_avg_anno = (observed_avg_anno - np.mean(observed_avg_anno))/np.std(observed_avg_anno)

	ys = []
	gs = []
	gs_shared = []
	gs_avg = []
	loaded_snps = np.where(observed_S_U_binary == 1.0)[0]
	ys.append([1]*len(loaded_snps))
	gs.append(observed_anno[loaded_snps])
	gs_shared.append(observed_shared_anno[loaded_snps])
	gs_avg.append(observed_avg_anno[loaded_snps])

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
		gs_shared.append(observed_shared_anno[randomly_sampled_null_indices])
		gs_avg.append(observed_avg_anno[randomly_sampled_null_indices])
	ys = np.hstack(ys)
	gs = np.hstack(gs)
	gs_shared = np.hstack(gs_shared)
	gs_avg = np.hstack(gs_avg)
	# Standardize annotations
	# gs = (gs - np.mean(gs))/np.std(gs)
	try:
		G = np.transpose(np.vstack((gs_shared, gs_avg, gs)))
		model = sm.Logit(ys, sm.add_constant(G))
		res = model.fit()
		ci = res.conf_int()[3,:]
		beta_lb = ci[0]
		beta_ub = ci[1]
		beta = res.params[3]
		pvalue = res.pvalues[3]
		#G = np.transpose(np.vstack((gs)))
		#model = sm.Logit(ys, sm.add_constant(gs))
		#res = model.fit()
		#ci = res.conf_int()[1,:]
		#beta_lb = ci[0]
		#beta_ub = ci[1]
		#beta = res.params[1]
		#pvalue = res.pvalues[1]
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

def extract_ct_specific_genomic_annotation(annotation_iter, genomic_annotation_file):
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
			valz = np.asarray(data[1:]).astype(float)
			if np.sum(valz) == 1.0 and valz[annotation_iter] == 1.0:
				anno.append(1.0)
			else:
				anno.append(0.0)
			#anno.append(float(anno_string))
	return np.asarray(anno)


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

shared_anno = extract_shared_genomic_annotation(genomic_annotation_file)
for annotation_iter, annotation_name in enumerate(annotation_names):
	anno, union_other_anno, avg_other_anno =extract_genomic_annotation(annotation_iter, genomic_annotation_file)
	#ct_spec_anno = extract_ct_specific_genomic_annotation(annotation_iter, genomic_annotation_file)
	if np.array_equal(np.isnan(anno), np.isnan(shared_anno)) == False:
		print('assumption error')
	for lv_num in range(K):
		test_results = annotation_lv_logistic_regression(anno, union_other_anno, avg_other_anno, S_U_binary[:, lv_num], af_bins)
		t.write(annotation_name + '\t' + str(lv_num + 1) + '\t' + str(test_results['pvalue']) + '\t' + str(test_results['beta']) + '\t' + str(test_results['beta_lb']) + '\t' + str(test_results['beta_ub']) + '\n')

t.close()
