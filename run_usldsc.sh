#!/bin/bash -l

#SBATCH
#SBATCH --time=60:00:00
#SBATCH --mem=30GB
#SBATCH --partition=lrgmem
#SBATCH --nodes=1




training_data_study_file="$1"
training_data_pairwise_ld_file="$2"
trainnig_data_cluster_info_file="$3"
k="$4"
model_version="$5"
output_root="$6"
b_v="$7"

module load python/2.7-anaconda


python run_usldsc.py $training_data_study_file $training_data_pairwise_ld_file $trainnig_data_cluster_info_file $k $model_version $output_root $b_v