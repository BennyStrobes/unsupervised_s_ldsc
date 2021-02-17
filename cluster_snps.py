import numpy as np 
import os
import pdb
import sys

def get_number_of_snps(input_file):
	f = open(input_file)
	counter = 0
	for line in f:
		counter = counter + 1
	f.close()
	return counter -1 

def add_cluster_labels(input_file, output_file):
	num_snps = get_number_of_snps(input_file)
	#num_snps = 124359
	# Dictiontary mapping from snp_id to cluster id
	snp_to_cluster = np.asarray(['null']*num_snps).astype('U25')
	# Dictionary mapping from cluster id to set of snps
	cluster_to_snp_arr = {}
	# First pass
	f = open(input_file)
	head_count = 0
	snp_counter = 0
	cluster_counter = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		counter = counter +1
		snp_ids = np.asarray(data[4].split(',')).astype(int)
		cluster_ids = snp_to_cluster[snp_ids]
		unique_cluster_ids = np.unique(cluster_ids)
		num_current_clusters = np.sum(unique_cluster_ids != 'null')
		if num_current_clusters == 0:
			temp_cluster_id = 'cluster' + str(cluster_counter)
			cluster_counter = cluster_counter + 1
		elif num_current_clusters == 1:
			for cluster_id_temper in unique_cluster_ids:
				if cluster_id_temper != 'null':
					temp_cluster_id = cluster_id_temper
		elif num_current_clusters > 1:
			temp_cluster_id = 'cluster' + str(cluster_counter)
			cluster_counter = cluster_counter + 1
		else:
			print('assuption errro: should not be here')
			pdb.set_trace()
		for snp_iter in range(len(snp_ids)):
			snp_id = snp_ids[snp_iter]
			cluster_id = snp_to_cluster[snp_id]
			if cluster_id == 'null':
				snp_to_cluster[snp_id] = temp_cluster_id
				if temp_cluster_id not in cluster_to_snp_arr:
					cluster_to_snp_arr[temp_cluster_id] = []
				cluster_to_snp_arr[temp_cluster_id].append(snp_id)
			elif cluster_id != temp_cluster_id:
				# need to deal with cluster merge
				old_cluster_snps = cluster_to_snp_arr[cluster_id]
				for temp_snp_id in old_cluster_snps:
					snp_to_cluster[temp_snp_id] = temp_cluster_id
					if temp_cluster_id not in cluster_to_snp_arr:
						cluster_to_snp_arr[temp_cluster_id] = []
					cluster_to_snp_arr[temp_cluster_id].append(temp_snp_id)
				cluster_to_snp_arr.pop(cluster_id)
			elif cluster_id == temp_cluster_id:
				continue
			else:
				print('assumption: should not be here')
				pdb.set_trace()

	f.close()
	# Second pass
	f = open(input_file)
	t = open(output_file,'w')
	head_count = 0
	snp_counter = 0
	cluster_counter = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\t' + 'snp_cluster\n')
			continue
		snp_ids = np.asarray(data[4].split(',')).astype(int)
		cluster_ids = snp_to_cluster[snp_ids]
		unique_cluster_ids = np.unique(cluster_ids)
		if len(unique_cluster_ids) != 1:
			print('assumption error')
			pdb.set_trace()
		if unique_cluster_ids[0] == 'null':
			print('assumption error')
		t.write(line + '\t' + unique_cluster_ids[0] + '\n')
	f.close()
	t.close()
	return cluster_to_snp_arr, snp_to_cluster, num_snps


def print_cluster_mapping(dicti, snp_to_cluster, cluster_mapping_output_file, num_snps):
	t = open(cluster_mapping_output_file,'w')
	t.write('cluster_id\tnumber_of_snps\tsnp_ids\n')
	used_snps = {}
	for cluster_id in sorted(dicti.keys()):
		cluster_snps = dicti[cluster_id]
		t.write(cluster_id + '\t' + str(len(dicti[cluster_id])) + '\t' + ','.join(np.asarray(cluster_snps).astype(str)) + '\n')
		for snp in cluster_snps:
			if snp in used_snps:
				print('assumption error')
				pdb.set_trace()
			used_snps[snp] = 1
			if snp_to_cluster[snp] != cluster_id:
				print('assumption error')
				pdb.set_trace()
	if len(used_snps) != num_snps:
		print('assumption error')
		pdb.set_trace()
	t.close()


input_file = sys.argv[1]
output_file = sys.argv[2]
cluster_mapping_output_file = sys.argv[3]


cluster_to_snp_arr, snp_to_cluster, num_snps = add_cluster_labels(input_file, output_file)

print_cluster_mapping(cluster_to_snp_arr, snp_to_cluster, cluster_mapping_output_file, num_snps)