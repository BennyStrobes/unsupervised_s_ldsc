#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --mem=4GB
#SBATCH --partition=lrgmem
#SBATCH --nodes=1


model_dir="$1"
model_name="$2"
pairwise_ld_file="$3"
usldsc_visualize_results_dir="$4"
processed_ukbb_dir="$5"

module load python/2.7-anaconda
python organize_usldsc_results_for_visualization.py $model_dir $model_name $pairwise_ld_file $usldsc_visualize_results_dir



module load R/3.5.1
Rscript visualize_usldsc_results.R $model_name $usldsc_visualize_results_dir $processed_ukbb_dir
