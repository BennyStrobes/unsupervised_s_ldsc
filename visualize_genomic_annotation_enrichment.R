args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(grid)
library(cowplot)
library(umap)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
options(bitmapType = 'cairo', device = 'pdf')





gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}





make_genomic_anno_error_bar_plot <- function(df) {
	p <- ggplot(data=df, aes(x=annotation, y=beta, group=latent_factor)) +
		geom_errorbar(aes(ymin = beta_lb, ymax = beta_ub), width=0,position=position_dodge(width=0.5)) +
		geom_point(size=1.8, aes(color=latent_factor),position=position_dodge(width=0.5)) +
		gtex_v8_figure_theme() +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
		labs(y="log odds ratio",x="", color="latent factor") +
		geom_hline(yintercept=0.0) + 
		theme(legend.position="top")
	return(p)
}

make_genomic_anno_error_bar_plot_v2_for_one_factor <- function(df, factor_num, y_axis_label_bool) {
	print(dim(df))
	p <- ggplot(data=df, aes(y=annotation, x=beta)) +
		geom_errorbar(aes(xmin = beta_lb, xmax = beta_ub), width=0,position=position_dodge(width=0.5)) +
		geom_point(size=.8,position=position_dodge(width=0.5)) +
		gtex_v8_figure_theme() +
		geom_vline(xintercept=0.0) + 
		labs(y="", x=paste0("Factor ", factor_num))
	if (y_axis_label_bool == FALSE) {
		p <- p + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
	}
	return(p)

}

make_genomic_anno_error_bar_plot_v2 <- function(df) {
	factor <- 0
	p0 <- make_genomic_anno_error_bar_plot_v2_for_one_factor(df[df$latent_factor == factor, ], factor, TRUE)

	factor <- 1
	p1 <- make_genomic_anno_error_bar_plot_v2_for_one_factor(df[df$latent_factor == factor, ], factor, FALSE)


	factor <- 2
	p2 <- make_genomic_anno_error_bar_plot_v2_for_one_factor(df[df$latent_factor == factor, ], factor, FALSE)


	factor <- 3
	p3 <- make_genomic_anno_error_bar_plot_v2_for_one_factor(df[df$latent_factor == factor, ], factor, FALSE)


	factor <- 4
	p4 <- make_genomic_anno_error_bar_plot_v2_for_one_factor(df[df$latent_factor == factor, ], factor, FALSE)


	factor <- 5
	p5 <- make_genomic_anno_error_bar_plot_v2_for_one_factor(df[df$latent_factor == factor, ], factor, FALSE)


	factor <- 6
	p6 <- make_genomic_anno_error_bar_plot_v2_for_one_factor(df[df$latent_factor == factor, ], factor, FALSE)


	factor <- 7
	p7 <- make_genomic_anno_error_bar_plot_v2_for_one_factor(df[df$latent_factor == factor, ], factor, FALSE)


	factor <- 8
	p8 <- make_genomic_anno_error_bar_plot_v2_for_one_factor(df[df$latent_factor == factor, ], factor, FALSE)


	factor <- 9
	p9 <- make_genomic_anno_error_bar_plot_v2_for_one_factor(df[df$latent_factor == factor, ], factor, FALSE)


    p_merge = plot_grid(p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow=1, rel_widths=c(2.7, 1, 1, 1, 1, 1, 1, 1, 1, 1))

	return(p_merge)
}


make_annotation_beta_factor_heatmap <- function(df) {

	df2 <- data.frame(annotation=df$annotation, latent_factor=df$latent_factor, beta=df$beta)
	mat_t <- dcast(df2, annotation ~ latent_factor)
	anno_names <- mat_t[,1]
	mat <- as.matrix(mat_t[, 2:dim(mat_t)[2]])
	ord <- hclust( dist(scale(mat), method = "euclidean"), method = "ward.D" )$order
	#factors <- as.matrix(t(factors))

	# Initialize PVE heatmap
    #factor_colnames <- paste0("usldsc_factor", 1:(dim(factors)[2]))
    #factor_rownames <- studies
    #colnames(factors) <- factor_colnames
    #rownames(factors) <- factor_rownames

	#ord <- hclust( dist(scale(factors), method = "euclidean"), method = "ward.D" )$order

    #melted_mat <- melt(factors)
    #colnames(melted_mat) <- c("Covariate", "Loading","factor_value")

    #melted_mat$Covariate = factor(melted_mat$Covariate, levels=rownames(factors)[ord])
    #melted_mat$Loading = factor(melted_mat$Loading, levels=factor_colnames)
	 #  Use factors to represent covariate and pc name
   	# melted_mat$Covariate 
    # melted_mat$Covariate <- factor(melted_mat$Covariate, levels = rownames(pve_map)[ord])
    #melted_mat$PC <- substr(as.character(melted_mat$PC),3,5)
    #melted_mat$PC <- factor(melted_mat$PC, levels=paste0("", 1:(length(unique(melted_mat$PC)))))

    
    #levels(melted_mat$PC) = paste0("PC", 1:(length(levels(melted_mat$PC))))
    #  PLOT!
    df$annotation = factor(df$annotation, levels=anno_names[ord])

    heatmap <- ggplot(data=df, aes(y=annotation, x=latent_factor)) + geom_tile(aes(fill=beta)) + scale_fill_gradient2(midpoint=0, guide="colorbar")
    heatmap <- heatmap + labs(y="",x="",fill="")
    heatmap <- heatmap + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
    # Save File
    return(heatmap)


}


#############################
# Command line args
#############################
genomic_annotation_dir <- args[1]
chrom_num <- args[2]
model_name <- args[3]
output_dir <- args[4]


# File containing cell type genomic annotation enrichment results
ct_genomic_anno_enrichment_file <- paste0(genomic_annotation_dir, "chr_", chrom_num, "_sldsc_cell_type_genomic_annotations_enrichment_within_", model_name, "logstic_regression_results.txt")
ct_genomic_anno_enrichment_df <- read.table(ct_genomic_anno_enrichment_file, header=TRUE)
ct_genomic_anno_enrichment_df$latent_factor <- factor(ct_genomic_anno_enrichment_df$latent_factor)

sldsc_baseline_genomic_anno_enrichment_file <- paste0(genomic_annotation_dir, "chr_", chrom_num, "_baseline_sldsc_genomic_annotations_enrichment_within_", model_name, "logstic_regression_results.txt")
sldsc_baseline_genomic_anno_enrichment_df <- read.table(sldsc_baseline_genomic_anno_enrichment_file, header=TRUE)
sldsc_baseline_genomic_anno_enrichment_df$latent_factor <- factor(sldsc_baseline_genomic_anno_enrichment_df$latent_factor)
sldsc_baseline_genomic_anno_enrichment_df <- sldsc_baseline_genomic_anno_enrichment_df[complete.cases(sldsc_baseline_genomic_anno_enrichment_df), ]

############################
# Visualize
############################

# Make error-bar plot showing snp-factor enrichment within genomic annotations
output_file <- paste0(output_dir, model_name, "cell_type_genomic_anno_enrichment_error_bar.pdf")
#error_bar_plot <- make_genomic_anno_error_bar_plot(ct_genomic_anno_enrichment_df)
#ggsave(error_bar_plot, file=output_file, width=12.0, height=6.0, units="in")


# Maker error bar plot showing snp-factor enrichment with genomic annotations
output_file <- paste0(output_dir, model_name, "baseline_ldsc_genomic_anno_enrichment_error_bar_v2.pdf")
#error_bar_plot <- make_genomic_anno_error_bar_plot_v2(sldsc_baseline_genomic_anno_enrichment_df)
#ggsave(error_bar_plot, file=output_file, width=25, height=18.0, units="in")


# Maker error bar plot showing snp-factor enrichment with genomic annotations
output_file <- paste0(output_dir, model_name, "cell_type_genomic_anno_enrichment_beta_heatmap.pdf")
beta_heatmap <- make_annotation_beta_factor_heatmap(ct_genomic_anno_enrichment_df)
ggsave(beta_heatmap, file=output_file, width=7, height=9.0, units="in")

# Maker error bar plot showing snp-factor enrichment with genomic annotations
output_file <- paste0(output_dir, model_name, "baseline_ldsc_genomic_anno_enrichment_beta_heatmap.pdf")
beta_heatmap <- make_annotation_beta_factor_heatmap(sldsc_baseline_genomic_anno_enrichment_df)
ggsave(beta_heatmap, file=output_file, width=7, height=9.0, units="in")

