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

#############################
# Command line args
#############################
processed_ld_score_dir = args[1]
output_dir = args[2]



# Generate vector of number of neighbors per snp
num_neighbors_per_snp = generate_vector_of_number_of_neighbors_per_snp(processed_ld_score_dir)

# Make histogram of number of neighbors per snp
output_file <- paste0(output_dir, "num_neighbors_per_snp_histogram.pdf")
num_neighbors_hist <- make_num_neighbors_per_snp_histogram(num_neighbors_per_snp)
ggsave(num_neighbors_hist, file=output_file, width=7.0, height=5.0, units="in")