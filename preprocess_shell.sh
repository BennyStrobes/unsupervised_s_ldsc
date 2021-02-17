#!/bin/bash -l

#SBATCH
#SBATCH --time=5:00:00
#SBATCH --mem=4GB
#SBATCH --partition=shared
#SBATCH --nodes=1



one_k_genomes_dir="$1"
one_k_genomes_sample_annotation_file="$2"
ukbb_studies_file="$3"
centimorgan_map_dir="$4"
processed_1k_genomes_genotype_dir="$5"
processed_ukbb_dir="$6"
processed_ld_score_dir="$7"


# Generate file containing only european 1K_genomes sample names
eur_1k_genomes_samples_file=$processed_1k_genomes_genotype_dir"eur_1k_genomes_samples.txt"
python generate_european_1k_genomes_file.py $one_k_genomes_sample_annotation_file $eur_1k_genomes_samples_file


# ukbb_studies_file to contain column for sample size
updated_ukbb_studies_file=$processed_ukbb_dir"updated_ukbb_studies.txt"
python add_sample_size_column_to_ukbb_studies_file.py $ukbb_studies_file $updated_ukbb_studies_file


