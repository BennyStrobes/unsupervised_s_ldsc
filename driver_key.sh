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
ukbb_studies_file="/work-zfs/abattle4/bstrober/unsupervised_s_ldsc/input_data/independent_.9_seed_2_observed_heritable_.1_ukbb_studies.tsv"

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
chrom_num="4"
if false; then
sbatch preprocess_shell_parallel_per_chromosome.sh $one_k_genomes_dir $one_k_genomes_sample_annotation_file $ukbb_studies_file $centimorgan_map_dir $processed_1k_genomes_genotype_dir $processed_ukbb_dir $processed_ld_score_dir $chrom_num
fi

module load R/3.5.1
if false; then
Rscript visualize_processed_data.R $processed_ld_score_dir $processed_ukbb_dir $visualize_processed_data_dir
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

# Directory containing results of Model simulation
usldsc_simulation_results_dir=$usldsc_root"usldsc_simulation_results/"
if false; then
sh organize_usldsc_training_data.sh $chromosomes_for_training $processed_ukbb_dir $processed_ld_score_dir $organized_training_data_dir
fi



training_data_study_file=$organized_training_data_dir"usldsc_training_studies.txt"
training_data_pairwise_ld_file=$organized_training_data_dir"usldsc_training_pairwise_ld_files.txt"
training_data_cluster_info_file=$organized_training_data_dir"usldsc_training_snp_cluster_files.txt"




model_version="vi"
simulation_output_root=$usldsc_simulation_results_dir"usldsc_simulation_"
if false; then
sh run_usldsc_simulation.sh $training_data_study_file $training_data_pairwise_ld_file $training_data_cluster_info_file $simulation_output_root
fi


model_version="vi"
k="10"
echo "OPTIMIZE Non-sim"
b_v="1"
output_root=$trained_usldsc_model_dir"trained_usldsc_"$model_version"_k_"$k"_b_v_prior_"$b_v"_"
if false; then
sbatch run_usldsc.sh $training_data_study_file $training_data_pairwise_ld_file $training_data_cluster_info_file $k $model_version $output_root $b_v
fi

k="10"
model_version="vi"
model_name="trained_usldsc_"$model_version"_k_"$k"_b_v_prior_"$b_v"_v_orig_"
if false; then
sh visualize_usldsc_results.sh $trained_usldsc_model_dir $model_name $training_data_pairwise_ld_file $usldsc_visualize_results_dir $processed_ukbb_dir
fi





















training_data_study_file=$organized_training_data_dir"usldsc_training_chromosome_standardized_studies.txt"
training_data_pairwise_ld_file=$organized_training_data_dir"usldsc_training_pairwise_ld_files.txt"
training_data_cluster_info_file=$organized_training_data_dir"usldsc_training_standardized_snp_cluster_files.txt"

model_version="weighted_component_gamma_vi"
k="10"
echo "OPTIMIZE Non-sim"
b_v="1000000000000000"
output_root=$trained_usldsc_model_dir"trained_standardized_usldsc_"$model_version"_k_"$k"_b_v_prior_"$b_v"_"
if false; then
sh run_usldsc.sh $training_data_study_file $training_data_pairwise_ld_file $training_data_cluster_info_file $k $model_version $output_root $b_v
fi


model_version="weighted_dirichlet_constrained_V_vi"
k="30"
echo "OPTIMIZE Non-sim"
b_v="1"
output_root=$trained_usldsc_model_dir"trained_standardized_usldsc_"$model_version"_k_"$k"_b_v_prior_"$b_v"_other_init2_"
if false; then
sbatch run_usldsc.sh $training_data_study_file $training_data_pairwise_ld_file $training_data_cluster_info_file $k $model_version $output_root $b_v
fi

k="30"
model_version="weighted_dirichlet_constrained_V_vi"
model_name="trained_standardized_usldsc_"$model_version"_k_"$k"_b_v_prior_"$b_v"_other_init2_"
if false; then
sh visualize_usldsc_results.sh $trained_usldsc_model_dir $model_name $training_data_pairwise_ld_file $usldsc_visualize_results_dir $processed_ukbb_dir
fi




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

# Directory containing visualizations of enrichment of geneomic annotations in U-LDSC variants
visualize_genomic_annotation_enrichment_dir=$enrichment_root"visualize_genomic_annotation_enrichment/"


k="10"
model_version="vi"
model_name="trained_usldsc_"$model_version"_k_"$k"_b_v_prior_"$b_v"_v_orig_"

chrom_num="4"
if false; then
sh extract_genomic_annotations_per_chromosome.sh $chrom_num $genomic_annotation_dir $processed_1k_genomes_genotype_dir $baseline_sldsc_annotation_dir $sldsc_cell_type_annotation_dir $trained_usldsc_model_dir $model_name $visualize_genomic_annotation_enrichment_dir
fi

























################### 
# Debug
###################
output_dir="/work-zfs/abattle4/bstrober/unsupervised_s_ldsc/debug_ldsc/"


study_name="KNEE_ARTHROSIS"
study_file="/work-zfs/abattle4/bstrober/unsupervised_s_ldsc/debug_ldsc/KNEE_ARTHROSIS.ldsc.imputed_v3.both_sexes.tsv.bgz"

output_root=$output_dir$study_name
module load python/2.7-anaconda
if false; then
python munge_sumstats.py --sumstats $study_file --out $output_root --merge-alleles "/work-zfs/abattle4/bstrober/unsupervised_s_ldsc/debug_ldsc/w_hm3.snplist"
fi

if false; then
python ldsc.py --h2 $study_file --ref-ld-chr /work-zfs/abattle4/bstrober/tools/ldsc/eur_w_ld_chr/ --w-ld-chr /work-zfs/abattle4/bstrober/tools/ldsc/eur_w_ld_chr/ --out $output_root"_h2_"
fi

module load python/2.7-anaconda
if false; then
python ldsc.py --h2 BMI.sumstats.gz --ref-ld-chr /work-zfs/abattle4/bstrober/tools/ldsc/baseline/baseline. --w-ld-chr /work-zfs/abattle4/bstrober/tools/ldsc/weights_hm3_no_hla/weights. --overlap-annot --frqfile-chr /work-zfs/abattle4/bstrober/tools/ldsc/1000G_frq/1000G.mac5eur. --out BMI_baseline
fi
if false; then
python ldsc.py --h2 BMI.sumstats.gz --ref-ld-chr /work-zfs/abattle4/bstrober/tools/ldsc/eur_w_ld_chr/ --w-ld-chr /work-zfs/abattle4/bstrober/tools/ldsc/eur_w_ld_chr/ --out BMI_h2
fi


