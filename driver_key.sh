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
echo "OPTIMIZE Non-sim"
if false; then
b_v="0"
output_root=$trained_usldsc_model_dir"trained_usldsc_"$model_version"_k_"$k"_b_v_prior_"$b_v"_shared_tau_"
sbatch run_usldsc.sh $training_data_study_file $training_data_pairwise_ld_file $training_data_cluster_info_file $k $model_version $output_root $b_v


b_v="50"
output_root=$trained_usldsc_model_dir"trained_usldsc_"$model_version"_k_"$k"_b_v_prior_"$b_v"_shared_tau_"
sbatch run_usldsc.sh $training_data_study_file $training_data_pairwise_ld_file $training_data_cluster_info_file $k $model_version $output_root $b_v


b_v="100"
output_root=$trained_usldsc_model_dir"trained_usldsc_"$model_version"_k_"$k"_b_v_prior_"$b_v"_shared_tau_"
sbatch run_usldsc.sh $training_data_study_file $training_data_pairwise_ld_file $training_data_cluster_info_file $k $model_version $output_root $b_v


b_v="150"
output_root=$trained_usldsc_model_dir"trained_usldsc_"$model_version"_k_"$k"_b_v_prior_"$b_v"_shared_tau_"
sbatch run_usldsc.sh $training_data_study_file $training_data_pairwise_ld_file $training_data_cluster_info_file $k $model_version $output_root $b_v

b_v="200"
output_root=$trained_usldsc_model_dir"trained_usldsc_"$model_version"_k_"$k"_b_v_prior_"$b_v"_shared_tau_"
sbatch run_usldsc.sh $training_data_study_file $training_data_pairwise_ld_file $training_data_cluster_info_file $k $model_version $output_root $b_v
fi

if false; then

b_v="20"
output_root=$trained_usldsc_model_dir"trained_usldsc_"$model_version"_k_"$k"_b_v_prior_"$b_v"_"
sbatch run_usldsc.sh $training_data_study_file $training_data_pairwise_ld_file $training_data_cluster_info_file $k $model_version $output_root $b_v

b_v="40"
output_root=$trained_usldsc_model_dir"trained_usldsc_"$model_version"_k_"$k"_b_v_prior_"$b_v"_"
sbatch run_usldsc.sh $training_data_study_file $training_data_pairwise_ld_file $training_data_cluster_info_file $k $model_version $output_root $b_v

b_v="60"
output_root=$trained_usldsc_model_dir"trained_usldsc_"$model_version"_k_"$k"_b_v_prior_"$b_v"_"
sbatch run_usldsc.sh $training_data_study_file $training_data_pairwise_ld_file $training_data_cluster_info_file $k $model_version $output_root $b_v

b_v="80"
output_root=$trained_usldsc_model_dir"trained_usldsc_"$model_version"_k_"$k"_b_v_prior_"$b_v"_"
sbatch run_usldsc.sh $training_data_study_file $training_data_pairwise_ld_file $training_data_cluster_info_file $k $model_version $output_root $b_v

b_v="100"
output_root=$trained_usldsc_model_dir"trained_usldsc_"$model_version"_k_"$k"_b_v_prior_"$b_v"_"
sbatch run_usldsc.sh $training_data_study_file $training_data_pairwise_ld_file $training_data_cluster_info_file $k $model_version $output_root $b_v

fi
b_v="0"
model_name="trained_usldsc_"$model_version"_k_"$k"_b_v_prior_"$b_v"_shared_tau_"
sh visualize_usldsc_results.sh $trained_usldsc_model_dir $model_name $training_data_pairwise_ld_file $usldsc_visualize_results_dir $processed_ukbb_dir















########################################
# Enrichment analysis for results of unsupervised S-LDSC
########################################

# Input data files
#####################

# from ashton
# baseline_sldsc_annotation_dir="/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/ldsc_annotations/baselineLF_v2.2.UKB/annots/"

# Directory containing baseline sldsc annotations
# Downloaded from https://alkesgroup.broadinstitute.org/LDSCORE/1000G_Phase3_EAS_baselineLD_v2.2_ldscores.tgz on 3/11/21
baseline_sldsc_annotation_dir="/work-zfs/abattle4/lab_data/ldsc_baseline_anno/"

# Directory containing sldsc cell type annotations
# Downloaded from https://alkesgroup.broadinstitute.org on 3/11/21
sldsc_cell_type_annotation_dir="/work-zfs/abattle4/lab_data/s_ldsc_cell_type_groups/"


# Directory containing CADD annotations for UKBB
cadd_annotation_dir="/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/cadd_annotations/imputed_cadd_annots_ukbb.tsv"




# Output directories
#####################
# root directory for analysis related to unsupervised-s-ldsc
enrichment_root="/work-zfs/abattle4/bstrober/unsupervised_s_ldsc/enrichment/"

# Directory containing genomic annotations assigned to variants used in U-LDSC
genomic_annotation_dir=$enrichment_root"genomic_annotations/"


model_name="trained_usldsc_vi_k_10_temp_gamma_"

chrom_num="4"
if false; then
sh extract_genomic_annotations_per_chromosome.sh $chrom_num $genomic_annotation_dir $processed_1k_genomes_genotype_dir $baseline_sldsc_annotation_dir $sldsc_cell_type_annotation_dir $trained_usldsc_model_dir $model_name
fi



