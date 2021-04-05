import numpy as np 
import os
import sys
import pdb
import pandas as pd
import numpy as np
from mofapy2.run.entry_point import entry_point

def test_mofa():
	D = [1000,1000] # Number of features per view
	M = len(D)      # Number of views
	K = 5           # Number of factors
	N = [100,100]   # Number of samples per group
	G = len(N)      # Number of groups
	data_dt = pd.read_csv("http://ftp.ebi.ac.uk/pub/databases/mofa/getting_started/data.txt.gz", sep="\t")
	ent = entry_point()
	ent.set_data_options(
		scale_groups = False, 
		scale_views = False
	)
	ent.set_data_df(data_dt, likelihoods = ["gaussian","gaussian"])
	ent.set_model_options(
    	factors = 10, 
    	spikeslab_weights = True, 
    	spikeslab_factors = True,
   		ard_factors = True,
    	ard_weights = True
	)
	ent.set_train_options(
	    iter = 1000, 
	    convergence_mode = "fast", 
	    startELBO = 1, 
	    freqELBO = 1, 
	    dropR2 = 0.001, 
	    gpu_mode = True, 
	   	verbose = False, 
   	 	seed = 1
	)
	ent.build()

	ent.run()
	pdb.set_trace()

def run_mofa_plus(Y):
	K=10
	G = 1
	M = 1
	D = [Y.shape[1]]
	N = [Y.shape[0]]
	data_mat = [[None for g in range(G)] for m in range(M)]
	data_mat[0][0] = np.copy(Y)
	ent = entry_point()
	ent.set_data_options(
		scale_groups = False, 
		use_float32 = True,
		scale_views = False
	)
	ent.set_data_matrix(data_mat, likelihoods = ["gaussian"])

	ent.set_model_options(
		factors = K*2, 
		spikeslab_weights = True, 
		ard_factors = True,
    	ard_weights = True
	)

	ent.set_train_options(
		iter = 300, 
		convergence_mode = "slow", 
		startELBO = 1, 
		freqELBO = 1, 
		dropR2 = 0.001, 
		nostop = True,
		gpu_mode = True, 
		startSparsity=1,
		verbose = False, 
		seed = 1
	)

	ent.build()
	ent.run()
	pdb.set_trace()
	learned_U = (ent.model.nodes['Z'].getExpectations())['EN']
	learned_S_U = (ent.model.nodes['Z'].getExpectations())['EB']
	learned_V = (ent.model.nodes['W'].getExpectations())[0]['EN']
	learned_S_V = (ent.model.nodes['W'].getExpectations())[0]['EB']
	theta_w = ent.model.nodes['ThetaW'].getExpectations()[0]['E']
	theta_z = ent.model.nodes['ThetaZ'].getExpectations()['E']
	fit = {'V': learned_V, 'V_S': learned_S_V, 'V_E': learned_V*learned_S_V, 'U': learned_U, 'U_S': learned_S_U, 'U_E': learned_U*learned_S_U, 'theta_V': theta_w, 'theta_U': theta_z}
	return fit


simulation_output_root = sys.argv[1]


# Load in simulated data
#Y_sim = np.load(simulation_output_root + 'sim_data_Y.npy')
#S_U_sim = np.load(simulation_output_root + 'sim_data_S_U.npy')
#S_V_sim = np.load(simulation_output_root + 'sim_data_S_V.npy')
#U_sim = np.load(simulation_output_root + 'sim_data_U.npy')
#V_sim = np.load(simulation_output_root + 'sim_data_V.npy')

num_snps=1000
num_studies=1000
simulation_k = 5
U_sim = np.random.randn(num_snps, simulation_k)*10.0
S_U_sim = np.random.randint(2, size=(num_snps, simulation_k))
U_E_sim = U_sim*S_U_sim
#S_U_sim = np.ones((num_snps, simulation_k))
V_sim = np.random.randn(simulation_k, num_studies)*10.0
S_V_sim = np.random.randint(2, size=(simulation_k, num_studies))
V_E_sim = V_sim*S_V_sim

tau_sim = 100.0
factor_predicted_mean_sim = np.dot(U_sim, V_sim*S_V_sim)
chi_squared_sim = np.zeros((num_snps, num_studies))

for study_num in range(num_studies):
	study_means = factor_predicted_mean_sim[:, study_num]
	chi_squared_sim[:, study_num] = np.random.normal(study_means, scale=np.sqrt(1.0/tau_sim))


#test_mofa()
mofa_fit = run_mofa_plus(chi_squared_sim)

pdb.set_trace()