# Ben Strober


#######################################
# Input data
#######################################

# Directory containing 1K genomes vcf files
# one file per chromosome of form ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
# I believe this was downloaded from: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
one_k_genomes_dir="/work-zfs/abattle4/lab_data/1k_genomes/"

# File containing which ancestry each 1K genomes sample belongs to
# Downloaded from: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
one_k_genomes_sample_annotation_file="/work-zfs/abattle4/lab_data/1k_genomes/integrated_call_samples_v3.20130502.ALL.panel"

# File with line for each UKBB study to be used containing location of all UKBB studies to be used
ukbb_studies_file="/work-zfs/abattle4/bstrober/unsupervised_s_ldsc/input_data/ukbb_studies.tsv"

# Directory containing centimorgan map file for each chromosome
# Files of format: genetic_map_chr19_combined_b37.txt
# Downloaded from https://mathgen.stats.ox.ac.uk/impute/1000GP%20Phase%203%20haplotypes%206%20October%202014.html
centimorgan_map_dir="/work-zfs/abattle4/lab_data/1k_genomes/cm_map/1000GP_Phase3/"

#######################################
# Output directories
#######################################
# root directory containing preprocessed data for unsupervised-s-ldsc
preprocess_root="/work-zfs/abattle4/bstrober/unsupervised_s_ldsc/preprocessed_data/"

# Directory containing preprocess 1K genomes genotype data
processed_1k_genomes_genotype_dir=$preprocess_root"processed_1k_genomes_genotype/"

# Directory containing preprocess UKBB dataa
processed_ukbb_dir=$preprocess_root"processed_ukbb/"

# Direcotry containing preprocess LD Score data
processed_ld_score_dir=$preprocess_root"processed_ld_scores/"
sh preprocess_shell.sh $one_k_genomes_dir $one_k_genomes_sample_annotation_file $ukbb_studies_file $centimorgan_map_dir $processed_1k_genomes_genotype_dir $processed_ukbb_dir $processed_ld_score_dir




if false; then
temp_genotype_dir="/work-zfs/abattle4/bstrober/tools/ldsc/1kg_eur/"

module load python/2.7-anaconda

python ldsc.py --bfile ${temp_genotype_dir}22 --l2 --l2-pairwise --chunk-size 1 --ld-wind-cm 1 --out 22
fi


