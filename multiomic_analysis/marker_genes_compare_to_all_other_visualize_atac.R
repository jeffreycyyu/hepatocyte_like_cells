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
library(fgsea)
library(GenomicRanges)
#



load("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/seurat_formatted_data.RData")


master_marker_genes = read.table("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/differentially_expressed/rna/all_markers_rna.txt", header = T, sep = "\t")
# master_marker_genes = master_marker_genes[which(master_marker_genes$ranking %in% 1:2), ]



master_marker_genes = master_marker_genes[master_marker_genes$p_val_adj == ave(master_marker_genes$p_val_adj, master_marker_genes$cluster, FUN=min),]
master_marker_genes = master_marker_genes[abs(master_marker_genes$avg_log2FC) == ave(master_marker_genes$avg_log2FC, master_marker_genes$cluster, FUN=max),]

master_marker_genes$gene_short_name = gene_translations[match(master_marker_genes$gene, gene_translations$id), "gene_short_name"]











#preprocsses atac seratassay















atac_peaks <- GetAssayData(all_atac, slot = "counts")


ccds$sequence_format_name = paste0("chr", ccds$chromosome, "-", ccds$cds_from, "-", ccds$cds_to)

master_marker_genes = merge(master_marker_genes, ccds, by.x = "gene_short_name", by.y = "gene", all.x = T, all.y = F)

all_atac_chromatin_assay = CreateChromatinAssay(atac_peaks,
                                       project = "hepatocyte_like_cells",
                                       assay = "peaks",
                                       sep = c("-", "-"),
                                       meta.data = as.data.frame(all_atac@meta.data))





for (i in 1:nrow(master_marker_genes)){


  cov_plot <- CoveragePlot(
    object = all_atac_chromatin_assay,
    region = master_marker_genes[i,"sequence_format_name"],
    annotation = F,
    peaks = F
  )

}
