#!/bin/bash -l

#SBATCH
#SBATCH --time=24:00:00
#SBATCH --mem=25GB
#SBATCH --partition=shared
#SBATCH --nodes=1



chromosomes_for_training="$1"
processed_ukbb_dir="$2"
processed_ld_score_dir="$3"
organized_training_data_dir="$4"


module load python/2.7-anaconda


python organize_usldsc_training_data.py $chromosomes_for_training $processed_ukbb_dir $processed_ld_score_dir $organized_training_data_dir