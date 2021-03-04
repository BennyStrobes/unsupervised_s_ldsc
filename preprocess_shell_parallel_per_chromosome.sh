#!/bin/bash -l

#SBATCH
#SBATCH --time=5:00:00
#SBATCH --mem=5GB
#SBATCH --partition=skylake
#SBATCH --nodes=1



one_k_genomes_dir="$1"
one_k_genomes_sample_annotation_file="$2"
ukbb_studies_file="$3"
centimorgan_map_dir="$4"
processed_1k_genomes_genotype_dir="$5"
processed_ukbb_dir="$6"
processed_ld_score_dir="$7"
chrom_num="$8"




# Generate file containing only european 1K_genomes sample names
eur_1k_genomes_samples_file=$processed_1k_genomes_genotype_dir"eur_1k_genomes_samples.txt"


# ukbb_studies_file to contain column for sample size
updated_ukbb_studies_file=$processed_ukbb_dir"updated_ukbb_studies.txt"


module load python/2.7-anaconda



echo $chrom_num
if false; then
# Use vcftools to filter to european ancestry samples
vcftools --gzvcf $one_k_genomes_dir"ALL.chr"$chrom_num".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" --keep $eur_1k_genomes_samples_file --out $processed_1k_genomes_genotype_dir"chr_"$chrom_num"_1k_genomes_european_only" --recode --recode-INFO-all
# Use vcftools to filter to sites that have MAF > .05 in european samples and those sites that dont pass filters
# further remove indels, as well as non bi-allelic snps
vcftools --vcf $processed_1k_genomes_genotype_dir"chr_"$chrom_num"_1k_genomes_european_only.recode.vcf" --remove-filtered-all --min-alleles 2 --max-alleles 2 --remove-indels --maf ".05" --out $processed_1k_genomes_genotype_dir"chr_"$chrom_num"_1k_genomes_european_only_maf_thresh_05" --recode --recode-INFO-all
# Generate file containing AF of each the above filtered, vcf file
vcftools --vcf $processed_1k_genomes_genotype_dir"chr_"$chrom_num"_1k_genomes_european_only_maf_thresh_05.recode.vcf" --freq --out $processed_1k_genomes_genotype_dir"chr_"$chrom_num"_1k_genomes_european_only_maf_thresh_05" 
fi

variant_file=$processed_ukbb_dir"chr_"$chrom_num"_variants_intersect_1k_genomes_ukbb.txt"
variant_position_file=$processed_ukbb_dir"chr_"$chrom_num"_variant_positions_intersect_1k_genomes_ukbb.txt"
if false; then
python filter_variants_to_those_in_both_1k_genomes_and_ukbb.py $chrom_num $processed_1k_genomes_genotype_dir $updated_ukbb_studies_file $variant_file $variant_position_file
fi
if false; then
vcftools --vcf $processed_1k_genomes_genotype_dir"chr_"$chrom_num"_1k_genomes_european_only_maf_thresh_05.recode.vcf" --positions $variant_position_file --out $processed_1k_genomes_genotype_dir"chr_"$chrom_num"_1k_genomes_european_only_maf_thresh_05_intersect_ukbb"  --recode --recode-INFO-all


vcftools --vcf $processed_1k_genomes_genotype_dir"chr_"$chrom_num"_1k_genomes_european_only_maf_thresh_05_intersect_ukbb.recode.vcf" --freq --out $processed_1k_genomes_genotype_dir"chr_"$chrom_num"_1k_genomes_european_only_maf_thresh_05_intersect_ukbb"


module load plink/1.90b6.4
plink --vcf $processed_1k_genomes_genotype_dir"chr_"$chrom_num"_1k_genomes_european_only_maf_thresh_05_intersect_ukbb.recode.vcf" --out $processed_1k_genomes_genotype_dir"chr_"$chrom_num"_1k_genomes_european_only_maf_thresh_05_intersect_ukbb"


plink --bfile $processed_1k_genomes_genotype_dir"chr_"$chrom_num"_1k_genomes_european_only_maf_thresh_05_intersect_ukbb" --cm-map $centimorgan_map_dir"genetic_map_chr"$chrom_num"_combined_b37.txt" $chrom_num --make-bed --out $processed_1k_genomes_genotype_dir"chr_"$chrom_num"_1k_genomes_european_only_maf_thresh_05_intersect_ukbb_with_cm_map"



python reorder_gwas_files_to_match_reference_genotype.py $chrom_num $processed_1k_genomes_genotype_dir"chr_"$chrom_num"_1k_genomes_european_only_maf_thresh_05_intersect_ukbb.frq" $updated_ukbb_studies_file $processed_ukbb_dir

fi

if false; then
python ldsc.py --bfile $processed_1k_genomes_genotype_dir"chr_"$chrom_num"_1k_genomes_european_only_maf_thresh_05_intersect_ukbb_with_cm_map" --l2 --l2-pairwise --ld-wind-cm 1 --chunk-size 1 --out $processed_ld_score_dir"chr_"$chrom_num"_ld_scores"
fi

if false; then
# Filter snp neighbors to exist if r_squared > thresh
r_squared_threshold=".1"
python filter_snp_neighbors.py $r_squared_threshold $processed_ld_score_dir"chr_"$chrom_num"_ld_scores_pairwise_l2.txt" $processed_ld_score_dir"chr_"$chrom_num"_ld_scores_pairwise_l2_neighbors_filtered"$r_squared_threshold".txt"


python cluster_snps.py $processed_ld_score_dir"chr_"$chrom_num"_ld_scores_pairwise_l2_neighbors_filtered"$r_squared_threshold".txt" $processed_ld_score_dir"chr_"$chrom_num"_ld_scores_pairwise_l2_neighbors_filtered"$r_squared_threshold"_clustered.txt" $processed_ld_score_dir"chr_"$chrom_num"_ld_scores_pairwise_l2_neighbors_filtered"$r_squared_threshold"_cluster_mapping.txt"
fi


# Filter snp neighbors to exist if r_squared > thresh
echo "FILTER"
r_squared_threshold=".2"
python filter_snp_neighbors.py $r_squared_threshold $processed_ld_score_dir"chr_"$chrom_num"_ld_scores_pairwise_l2.txt" $processed_ld_score_dir"chr_"$chrom_num"_ld_scores_pairwise_l2_neighbors_filtered"$r_squared_threshold".txt"

python cluster_snps.py $processed_ld_score_dir"chr_"$chrom_num"_ld_scores_pairwise_l2_neighbors_filtered"$r_squared_threshold".txt" $processed_ld_score_dir"chr_"$chrom_num"_ld_scores_pairwise_l2_neighbors_filtered"$r_squared_threshold"_clustered.txt" $processed_ld_score_dir"chr_"$chrom_num"_ld_scores_pairwise_l2_neighbors_filtered"$r_squared_threshold"_cluster_mapping.txt"


