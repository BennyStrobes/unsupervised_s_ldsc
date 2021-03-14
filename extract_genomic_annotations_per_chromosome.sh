#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --mem=20GB
#SBATCH --partition=lrgmem
#SBATCH --nodes=1



chrom_num="$1"
genomic_annotation_dir="$2"
processed_1k_genomes_genotype_dir="$3"
baseline_sldsc_annotation_dir="$4"
sldsc_cell_type_annotation_dir="$5"
trained_usldsc_model_dir="$6"
model_name="$7"

# File containing list of variants used in u-ldsc analysis
uldsc_variant_file=$processed_1k_genomes_genotype_dir"chr_"$chrom_num"_1k_genomes_european_only_maf_thresh_05_intersect_ukbb.frq"
module load python/2.7-anaconda


allele_frequency_output_file=$genomic_annotation_dir"chr_"$chrom_num"_allele_frequencies.txt"
if false; then
python generate_allele_frequency_groupings.py $uldsc_variant_file $allele_frequency_output_file
fi

# Extract set of genomic annotations describing variants in $uldsc_variant_file
baseline_ldsc_input_file=$baseline_sldsc_annotation_dir"baselineLD."$chrom_num".annot"
baseline_sldsc_output_file=$genomic_annotation_dir"chr_"$chrom_num"_baseline_sldsc_genomic_annotations.txt"
if false; then
python extract_baseline_sldsc_genomic_annotations_per_chromosome.py $baseline_ldsc_input_file $baseline_sldsc_output_file $uldsc_variant_file
fi


# Extract sldsc cell type annotations describing variants in $uldsc_variant_file
cell_type_sldsc_output_file=$genomic_annotation_dir"chr_"$chrom_num"_sldsc_cell_type_genomic_annotations.txt"
if false; then
python extract_sldsc_cell_type_annotations_per_chromosome.py $sldsc_cell_type_annotation_dir $chrom_num $cell_type_sldsc_output_file $uldsc_variant_file
fi


U_S_npy_file=$trained_usldsc_model_dir$model_name"U_S.npy"
U_npy_file=$trained_usldsc_model_dir$model_name"U.npy"
enrichment_results_file_stem=$genomic_annotation_dir"chr_"$chrom_num"_baseline_sldsc_genomic_annotations_enrichment_within_"$model_name
if false; then
python perform_genomic_annotation_enrichment_analysis.py $allele_frequency_output_file $baseline_sldsc_output_file $U_S_npy_file $U_npy_file $enrichment_results_file_stem
fi


U_S_npy_file=$trained_usldsc_model_dir$model_name"U_S.npy"
U_npy_file=$trained_usldsc_model_dir$model_name"U.npy"
enrichment_results_file_stem=$genomic_annotation_dir"chr_"$chrom_num"_sldsc_cell_type_genomic_annotations_enrichment_within_"$model_name
python perform_genomic_annotation_enrichment_analysis.py $allele_frequency_output_file $cell_type_sldsc_output_file $U_S_npy_file $U_npy_file $enrichment_results_file_stem









