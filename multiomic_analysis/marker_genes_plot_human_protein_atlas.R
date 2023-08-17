rm(list = ls())


library(monocle3)
library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(biovizBase)
library(ggforce)
library(ggridges)
library(scatterpie)
library(scales)
library(ggtext)
library(ggnewscale)
library(stringr)

#
#
#
# # to test with highly variable genes follow bellow and skip the load data and info section
# load("~/scratch/hepatocyte_like_cells/multiomic_analysis/all_but_the_plot.RData")
# all_rna <- FindVariableFeatures(all_rna, selection.method = "vst", nfeatures = 400)
# master_marker_genes$ensembl_id = head(VariableFeatures(all_rna), 5)
# master_marker_genes$gene_short_name= gene_translations[match(head(VariableFeatures(all_rna), 5), gene_translations$id), "gene_short_name"]
#
#
# differentially_expressed_genes = read.table("~/scratch/hepatocyte_like_cells/multiomic_analysis/differentially_expressed/safety/all_markers_rna.txt", sep = "\t", header = T)
# for (i in 1:length(order_cell_types)){
#   cell_type_i_name = order_cell_types[i]
#   differentially_expressed_genes_filtered = differentially_expressed_genes[which(differentially_expressed_genes$cluster == cell_type_i_name),]
#   differentially_expressed_genes_filtered = differentially_expressed_genes_filtered[order(differentially_expressed_genes_filtered$p_val_adj, decreasing = F),]
#   master_marker_genes[i, "ensembl_id"] = differentially_expressed_genes_filtered[1, "gene"]
#   master_marker_genes[i, "gene_short_name"] = paste0(gene_translations[match(differentially_expressed_genes_filtered[1, "gene"], gene_translations$id), "gene_short_name"], " (", cell_type_i_name, ")")
# }
#


# =========================
# LOAD DATA AND INFO
# =========================
load("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/seurat_formatted_data.RData")


# save(master_marker_genes, file = "/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/master_marker_genes_test_delete_me_trash.RData")
#
#
# #dataframe of all unique cell type comparisons
# differential_expression_comparisons = as.data.frame(t(combn(order_cell_types, 2)))
#
# marker_genes_colnames = c("ensembl_id", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "gene_short_name", "cell_type_1", "cell_type_2", "ranking")
# master_marker_genes = data.frame(matrix(ncol = length(marker_genes_colnames), nrow = 0))
# colnames(master_marker_genes) = marker_genes_colnames
#
# #rna
# for (k in 1:nrow(differential_expression_comparisons)){
#   i = differential_expression_comparisons[k, 1]
#   j = differential_expression_comparisons[k, 2]
#
#   markers_rna_ij = read.table(paste0("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/differentially_expressed/rna/rna_differentially_expressed_", i, "_", j, ".txt"), sep = "\t", header = T)
#   markers_rna_ij = markers_rna_ij[order(markers_rna_ij$p_val_adj, decreasing = F), ]
#   markers_rna_ij$cell_type_1 = i
#   markers_rna_ij$cell_type_2 = j
#   markers_rna_ij$ranking = 1:nrow(markers_rna_ij)
#   colnames(markers_rna_ij) = marker_genes_colnames
#
#   master_marker_genes = rbind(master_marker_genes, markers_rna_ij)
# }

master_marker_genes = read.table("/home/jyu/scratch/hepatocyte_like_cells/human_protein_atlas_marker_genes/ce_enriched_liver_Hepatocyte_1_Hepatocyte_2_Very.tsv", header = T, sep = "\t")
# master_marker_genes = master_marker_genes[which(master_marker_genes$ranking %in% 1:2), ]


#
# master_marker_genes = master_marker_genes[master_marker_genes$p_val_adj == ave(master_marker_genes$p_val_adj, master_marker_genes$cluster, FUN=min),]
# master_marker_genes = master_marker_genes[abs(master_marker_genes$avg_log2FC) == ave(master_marker_genes$avg_log2FC, master_marker_genes$cluster, FUN=max),]




#
print("filter genes: prefiltered by only very high differrential score in human protein atlas")
print(paste0("number of genes: ", nrow(master_marker_genes)))
#
print("filter genes: genes that are contained in our seurat object's rownames")
master_marker_genes = master_marker_genes[which(master_marker_genes$Ensembl %in% rownames(all_rna)),]
print(paste0("number of genes: ", nrow(master_marker_genes)))

print("filter genes: genes with tissue expression clusters (i.e., not blank/empty)")
master_marker_genes = master_marker_genes[which(master_marker_genes$Tissue.expression.cluster != ""),]
print(paste0("number of genes: ", nrow(master_marker_genes)))

#
# print("filter genes: only genes with Enhanced reliability scores (both IF and IH)")
# master_marker_genes = master_marker_genes[which(master_marker_genes$Reliability..IH. == "Enhanced" & master_marker_genes$Reliability..IF. == "Enhanced"),]
# print(paste0("number of genes: ", nrow(master_marker_genes)))



colnames(master_marker_genes)[which(colnames(master_marker_genes) == "Gene")] = "gene_short_name"
colnames(master_marker_genes)[which(colnames(master_marker_genes) == "Ensembl")] = "gene"


tissue_clusters_table = as.data.frame(table(master_marker_genes$Tissue.expression.cluster))
tissue_clusters_table = tissue_clusters_table[order(tissue_clusters_table$Freq, decreasing = T),]
colnames(tissue_clusters_table) = c("tissue_cluster", "frequency")


print(tissue_clusters_table)






master_marker_genes_backup = master_marker_genes


for (tissue_cluster_index in 1:nrow(tissue_clusters_table)){

  tissue_cluster_i_name = tissue_clusters_table[tissue_cluster_index, "tissue_cluster"]
  tissue_cluster_i_filename = gsub(" ", "_", tissue_cluster_i_name)
  tissue_cluster_i_frequency = tissue_clusters_table[tissue_cluster_index, "frequency"]

  print(paste0("starting tissue cluster ", tissue_cluster_index, " of ", nrow(tissue_clusters_table), " (", tissue_cluster_i_frequency, " genes)"))

  master_marker_genes = master_marker_genes_backup[which(master_marker_genes_backup$Tissue.expression.cluster == tissue_cluster_i_name), ]

  # =========================
  # PLOTTIING
  # =========================

  genes_of_interest_colnames = c("counts_cell_type", "counts", "ensembl_id", "proportion_zero")

  #get the melted dataframe of all count information for genes of interest only
  master_all_genes_of_interest_counts = data.frame(matrix(ncol = length(genes_of_interest_colnames), nrow = 0))
  colnames(master_all_genes_of_interest_counts) = genes_of_interest_colnames

  for (i in 1:nrow(master_marker_genes)){
    gene_of_interest = master_marker_genes[i, "gene"]

    master_cell_type_j_gene_counts = data.frame(matrix(ncol = length(genes_of_interest_colnames), nrow = 0))
    colnames(master_cell_type_j_gene_counts) = genes_of_interest_colnames

    for (cell_type_j in order_cell_types){
      all_counts_vector = as.vector(all_rna[which(rownames(all_rna) == gene_of_interest), which(all_rna@meta.data[, "cell_type"] == cell_type_j)]@assays[["MAGIC_rna"]]@data)
      proportion_zero = length(which(all_counts_vector == 0))/length(all_counts_vector)
      cell_type_j_gene_counts = data.frame(counts_cell_type = cell_type_j,
                                           counts = all_counts_vector,
                                           ensembl_id = gene_of_interest,
                                           proportion_zero = proportion_zero)
      master_cell_type_j_gene_counts = rbind(master_cell_type_j_gene_counts, cell_type_j_gene_counts[which(cell_type_j_gene_counts$counts != 0),])
    }

    master_all_genes_of_interest_counts = rbind(master_all_genes_of_interest_counts, master_cell_type_j_gene_counts)


    print(paste0("done gene_of_interest: ", i, " of ", nrow(master_marker_genes)))
  }

  #merge for full info
  master_marker_genes_full = merge(master_marker_genes, master_all_genes_of_interest_counts, by.x = "gene", by.y = "ensembl_id")




  #proportion zero pie charts to plot
  #get the coordinates of the pie chart locations
  #determine where the center circle representing "all other cell types" should lie
  pie_chart_center_df<-data.frame(ID = 1,
                                  longitude = ceiling(max(all_rna@assays$MAGIC_rna@data, na.rm=T)),
                                  latitude = 0.5,
                                  cell_type = "all_others")

  make_circles <- function(centers, radius, nPoints = length(order_cell_types)+1){
    # centers: the data frame of centers with ID
    # radius: radius measured in kilometer
    #
    meanLat <- mean(centers$latitude)
    # length per longitude changes with lattitude, so need correction
    radiusLon <- radius /111 / cos(meanLat/57.3)
    radiusLat <- radius / 111
    circleDF <- data.frame(ID = rep(centers$ID, each = nPoints))
    angle <- seq(0,2*pi,length.out = nPoints)

    circleDF$longitude <- unlist(lapply(centers$longitude, function(x) x + radiusLon * cos(angle)))
    circleDF$latitude <- unlist(lapply(centers$latitude, function(x) x + radiusLat * sin(angle)))
    return(circleDF[1:nPoints-1, ])
  }


  #generate coordinates and append cell type name
  pie_chart_coordinates = make_circles(pie_chart_center_df, 25)
  pie_chart_coordinates$cell_type = order_cell_types
  pie_chart_coordinates = rbind(pie_chart_coordinates, pie_chart_center_df)

  #get information for pie charts such as proportion of 0s in the unimputed count matrix
  proportion_zero_plotting_genes_colnames = c("xloc", "yloc", "radius", "proportion_zero", "proportion_positive", "gene_short_name", "ensembl_id", "cell_type")
  proportion_zero_plotting_genes = data.frame(matrix(ncol = length(proportion_zero_plotting_genes_colnames), nrow = 0))


  master_marker_genes$gene_short_name = gene_translations[match(master_marker_genes$gene, gene_translations$id), "gene_short_name"]
  master_marker_genes_full$gene_short_name = gene_translations[match(master_marker_genes_full$gene, gene_translations$id), "gene_short_name"]


  for (marker_gene_i in 1:nrow(master_marker_genes)){
    print(paste0("starting gene ", marker_gene_i, " of ", nrow(master_marker_genes)))
    plotting_gene_ensembl_i = master_marker_genes[marker_gene_i, "gene"]
    plotting_gene_short_name_i = master_marker_genes[marker_gene_i, "gene_short_name"]

    for (cell_type_j in order_cell_types){

      number_zeros_ij = length(which(all_rna[plotting_gene_ensembl_i, which(Idents(all_rna) == cell_type_j)]@assays[["MAGIC_rna"]]@data == 0))
      number_positives_ij = length(which(all_rna[plotting_gene_ensembl_i, which(Idents(all_rna) == cell_type_j)]@assays[["MAGIC_rna"]]@data != 0))
      total_counts = sum(number_zeros_ij, number_positives_ij)

      plotting_gene_i_info = data.frame(xloc = pie_chart_coordinates[which(pie_chart_coordinates$cell_type == cell_type_j), "longitude"],
                                        yloc = pie_chart_coordinates[which(pie_chart_coordinates$cell_type == cell_type_j), "latitude"],
                                        radius = 0.1,
                                        proportion_zero = number_zeros_ij/total_counts,
                                        proportion_positive = number_positives_ij/total_counts,
                                        gene_short_name = plotting_gene_short_name_i,
                                        ensembl_id = plotting_gene_ensembl_i,
                                        cell_type = cell_type_j)

      proportion_zero_plotting_genes = rbind(proportion_zero_plotting_genes, plotting_gene_i_info)

      print(paste0("finished cell_type: ", cell_type_j))

    }

    print(paste0("finished binding gene ", marker_gene_i, " of ", nrow(master_marker_genes)))

  }

  master_marker_genes_full_substituted_zeros = master_marker_genes_full
  master_marker_genes_full_substituted_zeros[which(master_marker_genes_full_substituted_zeros$counts > 10), "counts"] = 1

  master_marker_genes_full_substituted_zeros_labels <- data.frame(gene_short_name=master_marker_genes[, "gene_short_name"], label=master_marker_genes[,"Gene.description"])




  #
  # marker_genes_plot = ggplot(master_marker_genes_full_substituted_zeros, aes(x = log10(counts), y = after_stat(density), color = counts_cell_type, fill = after_scale(color))) +
  #   geom_histogram(alpha = 0.35,
  #                  binwidth=1,
  #                  position = position_dodge(width = 0.4, preserve = "single")) +
  #   facet_wrap(~gene_short_name, ncol=1, strip.position = "left") +
  #   theme_bw() +
  #   theme(axis.line = element_line(colour = "black"),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         panel.border = element_blank(),
  #         panel.background = element_blank(),
  #         axis.text.y = element_text(colour = 'black'),
  #         strip.background = element_blank(),
  #         strip.text = element_text(colour = 'black'),
  #         strip.placement = "outside") +
  #   scale_x_continuous(breaks = c(0:ceiling(log10(max(all_rna@assays$MAGIC_rna@data, na.rm=T)))),
  #                      # limits = c(-0.5, ceiling(log10(max(all_rna@assays$MAGIC_rna@data, na.rm=T)))),
  #                      labels=c("< 1", 1:ceiling(log10(max(all_rna@assays$MAGIC_rna@data, na.rm=T)))),
  # 		     oob = scales::squish) +
  #   scale_y_continuous(breaks = seq(0,1,0.5), limits = c(0,1)) +
  #   xlab("log10(Non-Zero RNA Counts)") +
  #   ylab("Proportion") +
  #   labs(color = "Cell Type") +
  #   geom_label(data = master_marker_genes_full_substituted_zeros_labels, aes(label=paste0("top marker of: ", label)),
  #              x = 2.5, y = 0.7, hjust=0, vjust=0, size = 3,
  #              inherit.aes = FALSE)

  colour_scheme_vector = c("#ED5564","#FFCE54", "#A0D568", "#4FC1E8", "#AC92EB")


  marker_genes_plot = ggplot(master_marker_genes_full_substituted_zeros, aes(x = counts, y = after_stat(density), color = factor(counts_cell_type, order_cell_types), fill = after_scale(color))) +
    geom_histogram(alpha = 0.35,
                   binwidth=1,
                   position = position_dodge(width = 0.4, preserve = "single")) +
    facet_wrap(~factor(gene_short_name, levels = master_marker_genes[, "gene_short_name"]), ncol=1, strip.position = "left") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_text(colour = 'black'),
          strip.background = element_blank(),
          strip.text = element_text(colour = 'black', size = 8),
          strip.placement = "outside") +
    scale_x_continuous(breaks = c(0:ceiling(max(all_rna@assays$MAGIC_rna@data, na.rm=T))),
                       limits = c(-0.5, ceiling(max(all_rna@assays$MAGIC_rna@data, na.rm=T))+0.5),
                       labels=c("< 1", 1:ceiling(max(all_rna@assays$MAGIC_rna@data, na.rm=T))),
                       oob = scales::squish) +
    scale_y_continuous(breaks = seq(0,1,0.5), limits = c(0,1)) +
    xlab("Non-Zero RNA Counts") +
    ylab("Proportion") +
    labs(color = "Cell Type") +
    geom_richtext(data = master_marker_genes_full_substituted_zeros_labels, aes(label=str_wrap(label, 15)), # "top marker of: "
                  x = -0.5, y = 0, hjust=0, vjust=0, size = 1.5,
                  inherit.aes = FALSE,
                  label.r = unit(0, "pt"),
                  angle = 90,
                  label.padding = unit(0, "lines")) +# label.size = NA
    coord_equal() +
    scale_colour_manual(values = colour_scheme_vector) +
    geom_richtext(data = data.frame(x = ceiling(max(all_rna@assays$MAGIC_rna@data, na.rm=T)), y = 0.95, gene_short_name = master_marker_genes[, "gene_short_name"][1], label = "Non-Zero Fraction"),
                  aes(x = x, y = y, label = label),
                  fill = NA,
                  label.color = NA,
                  size = 3,
                  inherit.aes = FALSE,
                  label.r = unit(0, "pt")) +
    ggtitle(tissue_cluster_i_name)




  for (cell_type_index in 1:length(order_cell_types)){

    # colfunc <- c(colour_scheme_vector[cell_type_index], "#000000")
    # pie_chart_colors = c(colour_scheme_vector[cell_type_index], colfunc(1))

    marker_genes_plot = marker_genes_plot +
      new_scale_fill() +
      geom_scatterpie(aes(x=xloc, y=yloc, r = radius), data=proportion_zero_plotting_genes[which(proportion_zero_plotting_genes$cell_type == order_cell_types[cell_type_index]),],
                      cols=c("proportion_zero", "proportion_positive"),
                      color = "black",
                      alpha = 0.75,
                      show.legend = F) +
      scale_fill_manual(values = c("white", colour_scheme_vector[cell_type_index]))

    print(paste0("finished adding pie chart for cell type #", cell_type_index))
  }



  # +
  # scale_fill_manual(values = c("#ED5564","#FFCE54", "#A0D568", "#4FC1E8", "#AC92EB"))




  ggsave(paste0("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/plots/human_protein_atlas_hepatocyte_markers/rna_human_protein_atlas_hepatocyte_very_high_marker_gene_per_cell_type_histogram_tissue_cluster_", tissue_cluster_i_filename, ".pdf"), plot = marker_genes_plot, width = 8, height = 40, limitsize = FALSE)

  print(paste0("FINISHED ", tissue_cluster_i_name))

}



print("FINISHED ALL PLOTTING (THE END)")
