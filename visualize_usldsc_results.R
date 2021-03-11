args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(grid)
library(cowplot)
library(umap)
library(data.table)
library(ggplot2)
library(RColorBrewer)
options(bitmapType = 'cairo', device = 'pdf')





gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}

make_study_clustered_factor_heatmap <- function(factors, studies) {
	factors <- as.matrix(t(factors))

	# Initialize PVE heatmap
    factor_colnames <- paste0("usldsc_factor", 1:(dim(factors)[2]))
    factor_rownames <- studies
    colnames(factors) <- factor_colnames
    rownames(factors) <- factor_rownames

	ord <- hclust( dist(scale(factors), method = "euclidean"), method = "ward.D" )$order

    melted_mat <- melt(factors)
    colnames(melted_mat) <- c("Covariate", "Loading","factor_value")

    melted_mat$Covariate = factor(melted_mat$Covariate, levels=rownames(factors)[ord])
    melted_mat$Loading = factor(melted_mat$Loading, levels=factor_colnames)
	 #  Use factors to represent covariate and pc name
   	# melted_mat$Covariate 
    # melted_mat$Covariate <- factor(melted_mat$Covariate, levels = rownames(pve_map)[ord])
    #melted_mat$PC <- substr(as.character(melted_mat$PC),3,5)
    #melted_mat$PC <- factor(melted_mat$PC, levels=paste0("", 1:(length(unique(melted_mat$PC)))))

    
    #levels(melted_mat$PC) = paste0("PC", 1:(length(levels(melted_mat$PC))))
    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(y=Covariate, x=Loading)) + geom_tile(aes(fill=factor_value)) + scale_fill_gradient2(midpoint=0, guide="colorbar")
    heatmap <- heatmap + labs(y="",x="",fill="")
    heatmap <- heatmap + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
    # Save File
    return(heatmap)


}

factor_umap_colored_by_categorical <- function(factors, covariates, covariate_name) {
	umap_loadings = umap(t(factors))$layout
	df <- data.frame(loading_1=umap_loadings[,1], loading_2=umap_loadings[,2], covariate=factor(covariates))
	plotter <- ggplot(df) + 
	           geom_point(aes(x=loading_1, y=loading_2, color=covariate), size=.1) +
	           gtex_v8_figure_theme() + 
	           guides(colour = guide_legend(override.aes = list(size=2))) +
	           labs(x="UMAP 1", y = "UMAP 2", color=covariate_name) + 
	           guides(colour=guide_legend(nrow=3,byrow=TRUE, override.aes = list(size=2))) +
	           theme(legend.position="none") + 
	           theme(legend.text = element_text(size=8), legend.title = element_text(size=8))
	return(plotter)
}

make_loading_distribution_histogram_for_one_component <- function(loadings_vec, component_num, data_name) {
    df <- data.frame(loading=loadings_vec)
    p <- ggplot(df, aes(x=loading))+
            geom_histogram(color="darkblue", fill="lightblue") +
            gtex_v8_figure_theme() +
            labs(x=paste0(data_name, " ", component_num))
    return(p)

}

make_fraction_neighbor_loaded_histogram <- function(file_stem) {
    data_name <- "Fraction of neighbor variants loaded"
    k = 0
    file_name <- paste0(file_stem, "fraction_loaded_neighbors_component_", k, ".txt")
    data <- read.table(file_name, header=FALSE)$V1
    p0 <- make_loading_distribution_histogram_for_one_component(data, (k+1), data_name)

    k = 1
    file_name <- paste0(file_stem, "fraction_loaded_neighbors_component_", k, ".txt")
    data <- read.table(file_name, header=FALSE)$V1
    p1 <- make_loading_distribution_histogram_for_one_component(data, (k+1), data_name)

    k = 2
    file_name <- paste0(file_stem, "fraction_loaded_neighbors_component_", k, ".txt")
    data <- read.table(file_name, header=FALSE)$V1
    p2 <- make_loading_distribution_histogram_for_one_component(data, (k+1), data_name)

    k = 3
    file_name <- paste0(file_stem, "fraction_loaded_neighbors_component_", k, ".txt")
    data <- read.table(file_name, header=FALSE)$V1
    p3 <- make_loading_distribution_histogram_for_one_component(data, (k+1), data_name)

    k = 4
    file_name <- paste0(file_stem, "fraction_loaded_neighbors_component_", k, ".txt")
    data <- read.table(file_name, header=FALSE)$V1
    p4 <- make_loading_distribution_histogram_for_one_component(data, (k+1), data_name)

    k = 5
    file_name <- paste0(file_stem, "fraction_loaded_neighbors_component_", k, ".txt")
    data <- read.table(file_name, header=FALSE)$V1
    p5 <- make_loading_distribution_histogram_for_one_component(data, (k+1), data_name)

    k = 6
    file_name <- paste0(file_stem, "fraction_loaded_neighbors_component_", k, ".txt")
    data <- read.table(file_name, header=FALSE)$V1
    p6 <- make_loading_distribution_histogram_for_one_component(data, (k+1), data_name)

    k = 7
    file_name <- paste0(file_stem, "fraction_loaded_neighbors_component_", k, ".txt")
    data <- read.table(file_name, header=FALSE)$V1
    p7 <- make_loading_distribution_histogram_for_one_component(data, (k+1), data_name)


    k = 8
    file_name <- paste0(file_stem, "fraction_loaded_neighbors_component_", k, ".txt")
    data <- read.table(file_name, header=FALSE)$V1
    p8 <- make_loading_distribution_histogram_for_one_component(data, (k+1), data_name)

    k = 9
    file_name <- paste0(file_stem, "fraction_loaded_neighbors_component_", k, ".txt")
    data <- read.table(file_name, header=FALSE)$V1
    p9 <- make_loading_distribution_histogram_for_one_component(data, (k+1), data_name)

    p_merge = plot_grid(p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol=2)

    return(p_merge)

}

make_loadings_distributions_histograms <- function(loadings, data_name) {
    k = 1
    p1 <- make_loading_distribution_histogram_for_one_component(loadings[, k], k, data_name)
    k = 2
    p2 <- make_loading_distribution_histogram_for_one_component(loadings[, k], k, data_name)
    k = 3
    p3 <- make_loading_distribution_histogram_for_one_component(loadings[, k], k, data_name)
    k = 4
    p4 <- make_loading_distribution_histogram_for_one_component(loadings[, k], k, data_name)
    k = 5
    p5 <- make_loading_distribution_histogram_for_one_component(loadings[, k], k, data_name)
    k = 6
    p6 <- make_loading_distribution_histogram_for_one_component(loadings[, k], k, data_name)
    k = 7
    p7 <- make_loading_distribution_histogram_for_one_component(loadings[, k], k, data_name)
    k = 8
    p8 <- make_loading_distribution_histogram_for_one_component(loadings[, k], k, data_name)
    k = 9
    p9 <- make_loading_distribution_histogram_for_one_component(loadings[, k], k, data_name)
    k = 10
    p10 <- make_loading_distribution_histogram_for_one_component(loadings[, k], k, data_name)


    p_merge = plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, ncol=2)

    return(p_merge)
}


make_ld_score_vs_max_loading_scatter <- function(ld_scores, max_loading) {
    df <- data.frame(ld_score=ld_scores, max_loading=max_loading)

    print(cor.test(df$ld_score, df$max_loading, method="pearson"))

    p <- ggplot(df, aes(x=ld_score, y=max_loading)) + geom_point(size=.1,alpha=.07) +
        gtex_v8_figure_theme()
    return(p)
}

make_ld_score_vs_loading_boxplot <- function(ld_scores, loadings) {
    ld_score_bins = cut(ld_scores, breaks=c(0, 5, 20, 100, 250, Inf), include.lowest=TRUE, labels=c("<5", "5-20", "20-100", "100-250", ">250"))
    df <- data.frame(ld_score_bin=factor(ld_score_bins), loadings=loadings[,1])
    #df[, ld_score_bin:=cut(ld_score, breaks=c(0, 5, 20, 100, 250, Inf), include.lowest=TRUE, labels=c("<5", "5-20", "20-100", "100-250", ">250"))]

    p <- ggplot(df, aes(x=ld_score_bin, y=loadings)) + 
        geom_boxplot() + 
        gtex_v8_figure_theme()
    return(p)
}

#############################
# Command line args
#############################
model_name <- args[1]
visualization_dir <- args[2]
processed_ukbb_dir <- args[3]

############################
# Load in data
############################

factor_file = paste0(visualization_dir, model_name, "V_S.txt")
factors = read.table(factor_file, header=FALSE)

loadings_file = paste0(visualization_dir, model_name, "U_S.txt")
loadings = read.table(loadings_file, header=FALSE)

loadings_bernoulli_file = paste0(visualization_dir, model_name, "S_U.txt")
loadings_bernoulli = read.table(loadings_bernoulli_file, header=FALSE)

ld_scores_file = paste0(visualization_dir, model_name, "ld_scores.txt")
ld_scores = read.table(ld_scores_file, header=FALSE)$V1

max_loading = apply(abs(loadings), 1, FUN=max)
max_bernoulli = apply(loadings_bernoulli, 1, FUN=max)

studies_file <- paste0(processed_ukbb_dir, "updated_ukbb_studies.txt")
study_data = read.table(studies_file, header=TRUE)
studies = as.character(study_data$study_descriptor)


############################
# Visualize
############################

# Make histogram showing distribution of fraction of neighbors of loaded variant that are also loaded
output_file <- paste0(visualization_dir, model_name, "fraction_neighbor_loaded_histograms.pdf")
fraction_neighbor_loaded_hist <- make_fraction_neighbor_loaded_histogram(paste0(visualization_dir, model_name))
ggsave(fraction_neighbor_loaded_hist, file=output_file, width=7.0, height=6.0, units="in")

# Make scatter plot comparing ld scores with max loadings
output_file <- paste0(visualization_dir, model_name, "ld_score_max_loading_scatter.pdf")
ld_score_max_loading_scatter <- make_ld_score_vs_max_loading_scatter(ld_scores, max_loading)
ggsave(ld_score_max_loading_scatter, file=output_file, width=7.0, height=6.0, units="in")

# Make box plot comparing ld scores with max loadings
output_file <- paste0(visualization_dir, model_name, "ld_score_abs_loading_boxplot.pdf")
ld_score_loading_boxplot <- make_ld_score_vs_loading_boxplot(ld_scores, abs(loadings))
ggsave(ld_score_loading_boxplot, file=output_file, width=7.0, height=6.0, units="in")


# Make loadings distribution histograms
output_file <- paste0(visualization_dir, model_name, "loadings_distributions_histograms.pdf")
loading_histogram <- make_loadings_distributions_histograms(loadings, "Loading")
ggsave(loading_histogram, file=output_file, width=7.0, height=6.0, units="in")

# Make loadings distribution histograms
output_file <- paste0(visualization_dir, model_name, "abs_loadings_distributions_histograms.pdf")
loading_histogram <- make_loadings_distributions_histograms(abs(loadings), "Absolute Loading")
ggsave(loading_histogram, file=output_file, width=7.0, height=6.0, units="in")

# Make loadings distribution histograms
output_file <- paste0(visualization_dir, model_name, "bernoulli_prob_distributions_histograms.pdf")
loading_histogram <- make_loadings_distributions_histograms(loadings_bernoulli, "Loading bernoulli prob")
ggsave(loading_histogram, file=output_file, width=7.0, height=6.0, units="in")

# Make heatmap of factor matix
output_file <- paste0(visualization_dir, model_name, "study_clustered_factor_heatmap.pdf")
factor_heatmap <- make_study_clustered_factor_heatmap(factors, studies)
ggsave(factor_heatmap, file=output_file, width=7.0, height=30.0, units="in")

