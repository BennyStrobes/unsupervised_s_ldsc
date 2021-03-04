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
ukbb_studies_file="/work-zfs/abattle4/bstrober/unsupervised_s_ldsc/input_data/independent_seed_1_heritable_ukbb_studies.tsv"

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

# Directory containing visualizations of processed data
visualize_processed_data_dir=$preprocess_root"visualize_processed_data/"


########################################
# Preprocess data
########################################
if false; then
sh preprocess_shell.sh $one_k_genomes_dir $one_k_genomes_sample_annotation_file $ukbb_studies_file $centimorgan_map_dir $processed_1k_genomes_genotype_dir $processed_ukbb_dir $processed_ld_score_dir
fi

if false; then
for chrom_num in {1..22}; do 
	echo $chrom_num
	sh preprocess_shell_parallel_per_chromosome.sh $one_k_genomes_dir $one_k_genomes_sample_annotation_file $ukbb_studies_file $centimorgan_map_dir $processed_1k_genomes_genotype_dir $processed_ukbb_dir $processed_ld_score_dir $chrom_num
done
fi
if false; then
chrom_num="4"
sh preprocess_shell_parallel_per_chromosome.sh $one_k_genomes_dir $one_k_genomes_sample_annotation_file $ukbb_studies_file $centimorgan_map_dir $processed_1k_genomes_genotype_dir $processed_ukbb_dir $processed_ld_score_dir $chrom_num
fi
if false; then
module load R/3.5.1
Rscript visualize_processed_data.R $processed_ld_score_dir $visualize_processed_data_dir
fi


########################################
# Run unsupervised S-LDSC
########################################


# Input data files
##################
chromosomes_for_training="/work-zfs/abattle4/bstrober/unsupervised_s_ldsc/input_data/chromosomes_for_unsupervised_s_ldsc.txt"

# Output directories
#####################
# root directory for analysis related to unsupervised-s-ldsc
usldsc_root="/work-zfs/abattle4/bstrober/unsupervised_s_ldsc/usldsc_training/"

# Directory containing organized data for usldsc training
organized_training_data_dir=$usldsc_root"organized_training_data/"

# Directory containing results of model training
trained_usldsc_model_dir=$usldsc_root"usldsc_results/"

# Directory containing results of model training
usldsc_visualize_results_dir=$usldsc_root"usldsc_results_visualization/"
if false; then
sh organize_usldsc_training_data.sh $chromosomes_for_training $processed_ukbb_dir $processed_ld_score_dir $organized_training_data_dir
fi

training_data_study_file=$organized_training_data_dir"usldsc_training_studies.txt"
training_data_pairwise_ld_file=$organized_training_data_dir"usldsc_training_pairwise_ld_files.txt"
training_data_cluster_info_file=$organized_training_data_dir"usldsc_training_snp_cluster_files.txt"

model_version="als"
k="7"
output_root=$trained_usldsc_model_dir"trained_usldsc_"$model_version"_k_"$k"_"
if false; then
sbatch run_usldsc.sh $training_data_study_file $training_data_pairwise_ld_file $training_data_cluster_info_file $k $model_version $output_root
fi



training_data_study_file=$organized_training_data_dir"usldsc_training_studies.txt"
training_data_pairwise_ld_file=$organized_training_data_dir"usldsc_training_pairwise_ld_files.txt"
training_data_cluster_info_file=$organized_training_data_dir"usldsc_training_snp_cluster_files.txt"

model_version="vi"
k="10"
output_root=$trained_usldsc_model_dir"trained_usldsc_"$model_version"_k_"$k"_"
echo "OPTIMIZE"
sh run_usldsc.sh $training_data_study_file $training_data_pairwise_ld_file $training_data_cluster_info_file $k $model_version $output_root



module load R/3.5.1
if false; then
Rscript visualize_usldsc_results.R $trained_usldsc_model_dir $processed_ukbb_dir $usldsc_results_dir
fi







#########################
# S-LDSC experimentation
##########################



if false; then
temp_genotype_dir="/work-zfs/abattle4/bstrober/tools/ldsc/1kg_eur/"

module load python/2.7-anaconda

python ldsc.py --bfile ${temp_genotype_dir}22 --l2 --l2-pairwise --chunk-size 1 --ld-wind-cm 1 --out 22
fi



module load python/2.7-anaconda

if false; then
python munge_sumstats.py --sumstats /work-zfs/abattle4/bstrober/tools/ldsc/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.txt --merge-alleles /work-zfs/abattle4/bstrober/tools/ldsc/w_hm3.snplist --out BMI --a1-inc
fi



if false; then
python ldsc.py --h2 BMI.sumstats.gz --ref-ld-chr /work-zfs/abattle4/bstrober/tools/ldsc/baseline/baseline. --w-ld-chr /work-zfs/abattle4/bstrober/tools/ldsc/weights_hm3_no_hla/weights. --overlap-annot --frqfile-chr /work-zfs/abattle4/bstrober/tools/ldsc/1000G_frq/1000G.mac5eur. --out BMI_baseline

fi