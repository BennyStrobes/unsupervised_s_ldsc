args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(grid)
library(cowplot)
library(umap)
library(ggplot2)
library(RColorBrewer)
options(bitmapType = 'cairo', device = 'pdf')




gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}


generate_vector_of_number_of_neighbors_per_snp <- function(processed_ld_score_dir) {
	arr <- c()
	chrom_num <- 1

	file_name <- paste0(processed_ld_score_dir, "chr_", chrom_num, "_ld_scores_pairwise_l2_neighbors_filtered.2.txt")

	data = read.table(file_name, header=TRUE, sep="\t")
	vec = data$pairwise_snp_indices
	num_snps = length(vec)
	for (snp_num in 1:num_snps) {
		num_neighbors <- length(strsplit(as.character(vec[snp_num]), ",")[[1]])
		arr <- c(arr, num_neighbors)
	}
	print(summary(arr))
	return(arr)
}

make_num_neighbors_per_snp_histogram <- function(num_neighbors_per_snp) {
	df <- data.frame(num_neighbors_per_snp=num_neighbors_per_snp)
	p <- ggplot(df, aes(x=num_neighbors_per_snp))+
  		geom_histogram(color="darkblue", fill="lightblue") + 
  		labs(x="Number of neighbors / SNP", y = "Frequency") + 
  		gtex_v8_figure_theme()
  	return(p)
}

make_chi_squared_distribution_histogram <- function(df) {
	p <- ggplot(df, aes(x=chi_sq)) +
  			geom_histogram()+
  			facet_grid(study~.) +
  			labs(x="Chi sqared / sample size", y = "") + 
  			gtex_v8_figure_theme() +
  			 theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  	return(p)
}

 make_chi_squared_variance_distribution_histogram <- function(df) {
 	p <- ggplot(df, aes(x=variance)) +
  			geom_histogram()+
  			labs(x="Var(Chi sqared / sample size)", y = "") + 
  			gtex_v8_figure_theme() 
  	return(p)
 }

#############################
# Command line args
#############################
processed_ld_score_dir = args[1]
processed_ukbb_dir = args[2]
output_dir = args[3]


if (FALSE) {
# Generate vector of number of neighbors per snp
num_neighbors_per_snp = generate_vector_of_number_of_neighbors_per_snp(processed_ld_score_dir)

# Make histogram of number of neighbors per snp
output_file <- paste0(output_dir, "num_neighbors_per_snp_histogram.pdf")
num_neighbors_hist <- make_num_neighbors_per_snp_histogram(num_neighbors_per_snp)
ggsave(num_neighbors_hist, file=output_file, width=7.0, height=5.0, units="in")
}


study_file <- paste0(processed_ukbb_dir, "updated_ukbb_studies.txt")
studies_df <- read.table(study_file, header=TRUE, sep="\t")
nrow = length(studies_df$study)
chi_sq_arr <- c()
study_name_arr <- c()
avgs <- c()
varz <- c()
for (row_num in 1:nrow) {
	file_name <- paste0(processed_ukbb_dir, "chr4_", as.character(studies_df$study[row_num]), "_chi_squared.txt")
	sample_size <- studies_df$sample_size[row_num]
	print(sample_size)
	chi_sq_df <- read.table(file_name, header=TRUE)
	chi_sq_arr <- c(chi_sq_arr, (chi_sq_df$chi_squared)/sample_size)
	study_name_arr <- c(study_name_arr, rep(as.character(studies_df$study[row_num]), length(chi_sq_df$chi_squared)))
	avgs <- c(avgs, mean(chi_sq_df$chi_squared/sample_size))
	varz <- c(varz, var(chi_sq_df$chi_squared/sample_size))
}
print(studies_df$study[order(avgs)])
chi_df_final <- data.frame(chi_sq=chi_sq_arr, study=factor(study_name_arr, levels=studies_df$study[order(avgs)]))

output_file <- paste0(output_dir, "study_chi_squared_distributions.pdf")
chi_sq_hists <- make_chi_squared_distribution_histogram(chi_df_final)
ggsave(chi_sq_hists, file=output_file, width=7.0, height=20.0, units="in")

chi_df_var <- data.frame(variance=varz)
output_file <- paste0(output_dir, "study_chi_squared_variance_distribution.pdf")
chi_sq_variance_hists <- make_chi_squared_variance_distribution_histogram(chi_df_var)
ggsave(chi_sq_variance_hists, file=output_file, width=7.0, height=3.0, units="in")








