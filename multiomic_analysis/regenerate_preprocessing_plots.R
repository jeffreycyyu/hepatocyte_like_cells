rm(list = ls())


library(monocle3)
library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(biovizBase)
library(ggforce)
library(foreach)
library(doParallel)
library(ggforce)



load("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/seurat_formatted_data.RData")
# load("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/seurat_formatted_data_intermediate_before_preprocessing_normalize_first_trash_me.RData")


#
#
#
# cell_type_translation = data.frame(short_name = order_cell_types,
#                                    full_name = c("induced pluripotent stem cells",
#                                                  "definitive endoderms",
#                                                  "hepatic endoderms",
#                                                  "immature hepatocytes",
#                                                  "mature hepatocytes"))
#
#
#
# DefaultAssay(all_rna) <- "MAGIC_rna"
#
# all_rna <- FindVariableFeatures(all_rna, selection.method = "vst", nfeatures = 400)
#
# #rna: scale data based on all genes (NOT just top x most variable)
# all_rna_genes <- rownames(all_rna)
# all_rna <- ScaleData(all_rna, features = all_rna_genes)
# all_rna <- RunPCA(all_rna, features = VariableFeatures(all_rna))
#
#
# #rna: add full name to metadata
# all_rna@meta.data$cell_type_full_name = cell_type_translation[match(all_rna@meta.data$cell_type, cell_type_translation$short_name), "full_name"]
#
#
# all_rna <- RunUMAP(all_rna, dims = 1:30, reduction.name = "umap.rna", reduction.key = "rnaUMAP_")
#





#
#
# #atac: remove those with less than x total regional counts
# all_atac = subset(all_atac, subset = nCount_peaks > 5000)
# #atac: normalize and log transform
# all_atac <- NormalizeData(all_atac, normalization.method = "LogNormalize", scale.factor = 10000)
# #atac: scale data based on all regions
# all_atac <- ScaleData(all_atac, features = rownames(all_atac))
#
# #atac: get all peaks that are greater than a certain numebr across all cells (NOT actually a "variable" feature)
# VariableFeatures(all_atac) <- names(which(Matrix::rowSums(all_atac) > 100))
# #atac: use TF-IDF to learn the structure of atac seq (Stuart & Butler et al. 2019)
# all_atac <- RunTFIDF(all_atac, method = 1, scale.factor = 10000)
# #atac: find top regions
# all_atac <- FindTopFeatures(all_atac, min.cutoff = "q0")
# #atac: svd
# all_atac <- RunSVD(all_atac)
# #atac: umap with lsi reduced components (do NOT include the first component as this is associated with sequencing depth for atac)
# all_atac <- RunUMAP(all_atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
#
#




#rna: add full name to metadata
all_atac@meta.data$cell_type_full_name = cell_type_translation[match(all_atac@meta.data$cell_type, cell_type_translation$short_name), "full_name"]




#=================

#rna: scatter plot of rna number of genes vs total gene counts
pdf("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/plots/feature_scatter_plot_rna_ncount_nfeature.pdf")
FeatureScatter(all_rna, feature1 = "nCount_rna", feature2 = "nFeature_rna") +
  ggtitle("") +
  xlab("Total RNA Counts") +
  ylab("Total Measurable Genes") +
  theme(legend.title=element_blank())
dev.off()



#rna: plot scatter plot of most variable genes based on standardized variance
pdf("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/plots/variable_features_plot_rna.pdf")
variable_features_plot <- VariableFeaturePlot(all_rna, selection.method = "vst")
top_variable_features <- head(VariableFeatures(all_rna), 10)
LabelPoints(plot = variable_features_plot,
            points = top_variable_features,
            labels = gene_translations[match(top_variable_features, gene_translations$id), "gene_short_name"],
            repel = TRUE)
dev.off()






#rna: pca plot
pdf("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/plots/pca_plot_rna.pdf")
buffer_space = 5
DimPlot(all_rna, reduction = "pca") +
  NoLegend() +
  xlab("PC 1") +
  ylab("PC 2") +
  geom_mark_ellipse(x = all_rna@reductions[["pca"]]@cell.embeddings[,"PC_1"],
                    y = all_rna@reductions[["pca"]]@cell.embeddings[,"PC_2"],
                    aes(color = all_rna@meta.data$cell_type,
                        label = all_rna@meta.data$cell_type,
                        description = all_rna@meta.data$cell_type_full_name),
                    label.fontsize = 5,
                    con.size = 0.5,
                    show.legend = F,
                    con.cap = 0) +
  xlim(c(min(all_rna@reductions[["pca"]]@cell.embeddings[,"PC_1"], na.rm=T)-buffer_space, max(all_rna@reductions[["pca"]]@cell.embeddings[,"PC_1"], na.rm=T)+buffer_space)) +
  ylim(c(min(all_rna@reductions[["pca"]]@cell.embeddings[,"PC_2"], na.rm=T)-buffer_space, max(all_rna@reductions[["pca"]]@cell.embeddings[,"PC_2"], na.rm=T)+buffer_space))
dev.off()




#rna: top x pcs heatmap for a subset of top genes and cells
pdf("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/plots/pc_heatmaps_plot_rna.pdf")
pcs_to_show = 9
DimHeatmap(all_rna, dims = 1:pcs_to_show, cells = 50, balanced = TRUE)
dev.off()



#atac: umap using lsi as reduction
pdf("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/plots/umap_plot_rna.pdf")
buffer_space=5
DimPlot(all_rna) +
  NoLegend() +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  geom_mark_ellipse(x = all_rna@reductions[["umap.rna"]]@cell.embeddings[,"rnaUMAP_1"],
                    y = all_rna@reductions[["umap.rna"]]@cell.embeddings[,"rnaUMAP_2"],
                    aes(color = all_rna@meta.data$cell_type,
                        label = all_rna@meta.data$cell_type,
                        description = all_rna@meta.data$cell_type_full_name),
                    label.fontsize = 5,
                    expand = unit(0, "mm"),
                    con.size = 0.5,
                    show.legend = F,
                    con.cap = 0) +
  xlim(c(min(all_rna@reductions[["umap.rna"]]@cell.embeddings[,"rnaUMAP_1"], na.rm=T)-buffer_space, max(all_rna@reductions[["umap.rna"]]@cell.embeddings[,"rnaUMAP_1"], na.rm=T)+buffer_space)) +
  ylim(c(min(all_rna@reductions[["umap.rna"]]@cell.embeddings[,"rnaUMAP_2"], na.rm=T)-buffer_space, max(all_rna@reductions[["umap.rna"]]@cell.embeddings[,"rnaUMAP_2"], na.rm=T)+buffer_space))
dev.off()


#
# #atac: umap using lsi as reduction
# pdf("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/plots/umap_lsi_plot_atac.pdf")
# buffer_space=5
# DimPlot(all_atac, reduction = "lsi") +
#   NoLegend() +
#   xlab("UMAP 1") +
#   ylab("UMAP 2") +
#   geom_mark_ellipse(x = all_atac@reductions[["lsi"]]@cell.embeddings[,"LSI_1"],
#                     y = all_atac@reductions[["lsi"]]@cell.embeddings[,"LSI_2"],
#                     aes(color = all_atac@meta.data$cell_type,
#                         label = all_atac@meta.data$cell_type,
#                         description = all_atac@meta.data$cell_type_full_name),
#                     label.fontsize = 5,
#                     expand = unit(0, "mm"),
#                     con.size = 0.5,
#                     show.legend = F,
#                     con.cap = 0) +
#   xlim(c(min(all_atac@reductions[["lsi"]]@cell.embeddings[,"LSI_1"], na.rm=T)-buffer_space, max(all_atac@reductions[["lsi"]]@cell.embeddings[,"LSI_1"], na.rm=T)+buffer_space)) +
#   ylim(c(min(all_atac@reductions[["lsi"]]@cell.embeddings[,"LSI_2"], na.rm=T)-buffer_space, max(all_atac@reductions[["lsi"]]@cell.embeddings[,"LSI_2"], na.rm=T)+buffer_space))
# dev.off()


#=====================

print("FINISHED EVERYTHING!!!")
