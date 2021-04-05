#!/bin/bash -l

#SBATCH
#SBATCH --time=24:00:00
#SBATCH --mem=25GB
#SBATCH --partition=shared
#SBATCH --nodes=1

training_data_study_file="$1"
training_data_pairwise_ld_file="$2"
training_data_cluster_info_file="$3"
simulation_output_root="$4"


module load python/2.7-anaconda

python run_usldsc_simulation.py $training_data_study_file $training_data_pairwise_ld_file $training_data_cluster_info_file $simulation_output_root


source ~/.bash_profile

if false; then
python run_mofa_plus.py $simulation_output_root
fi