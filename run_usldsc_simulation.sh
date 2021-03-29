#!/bin/bash -l

#SBATCH
#SBATCH --time=24:00:00
#SBATCH --mem=25GB
#SBATCH --partition=shared
#SBATCH --nodes=1


training_data_study_file="$1"
training_data_pairwise_ld_file="$2"
training_data_cluster_info_file="$3"
k="$4"
model_version="$5"
simulation_output_root="$6"
b_v="$7"


module load python/2.7-anaconda


python run_usldsc_simulation.py $training_data_study_file $training_data_pairwise_ld_file $training_data_cluster_info_file $k $model_version $simulation_output_root $b_v