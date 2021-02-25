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


#############################
# Command line args
#############################
model_results_dir <- args[1]
processed_ukbb_dir <- args[2]
visualization_dir <- args[3]

############################
# Load in data
############################

factor_file = paste0(model_results_dir, "V.txt")
factors = read.table(factor_file, header=FALSE)

studies_file <- paste0(processed_ukbb_dir, "updated_ukbb_studies.txt")
study_data = read.table(studies_file, header=TRUE)


studies = as.character(study_data$study)

for (study_iter in 1:length(studies)) {
	studies[study_iter] = as.character(paste0(studies[study_iter], "_", study_iter))
}

study_category <- as.character(read.table("/work-zfs/abattle4/bstrober/unsupervised_s_ldsc/input_data/studies_with_category.tsv",header=TRUE,sep="\t",quote="")$category)


############################
# Visualize
############################


# Make UMAP Scatter of factor matix colored by study category
output_file <- paste0(visualization_dir, "factor_umap_colored_by_study_category.pdf")
factor_umap <- factor_umap_colored_by_categorical(factors, study_category, "study_category")
ggsave(factor_umap, file=output_file, width=7.0, height=5.0, units="in")

# Make heatmap of factor matix
output_file <- paste0(visualization_dir, "study_clustered_factor_heatmap.pdf")
factor_heatmap <- make_study_clustered_factor_heatmap(factors, studies)
ggsave(factor_heatmap, file=output_file, width=7.0, height=30.0, units="in")

