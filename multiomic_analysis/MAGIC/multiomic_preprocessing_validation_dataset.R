rm(list = ls())


# library(monocle3)
library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
# library(EnsDb.Hsapiens.v86)
# library(biovizBase)
library(ggforce)
library(foreach)
library(doParallel)
library(Rmagic)
library(ggforce)



# =========================
# LOAD DATA AND INFO
# =========================
#
# #RNA
# load("/home/jyu/scratch/hepatocyte_like_cells/rna_analysis/sample_cell_data_sets_to_transfer.RData")
# #ATAC
# load("~/scratch/hepatocyte_like_cells/multiomic_analysis/atac_data_to_transfer_0924_04july2023.RData")


#CCDS genes
ccds = read.table("/home/jyu/scratch/hepatocyte_like_cells/rna_analysis/ccds_genes/CCDS.current.txt", sep = "\t")
colnames(ccds) = c("chromosome", "nc_accession", "gene", "gene_id", "ccds_id", "ccds_status", "cds_strand", "cds_from", "cds_to", "cds_locations", "match_type")

#
# #referenced lists
# order_cell_types = c("ipsc", "dfendo", "hpendo", "immhp", "mathp")
# order_assays = c("rna", "atac")
# cell_type_translation = data.frame(short_name = order_cell_types,
#                                    full_name = c("induced pluripotent stem cells",
#                                                  "definitive endoderms",
#                                                  "hepatic endoderms",
#                                                  "immature hepatocytes",
#                                                  "mature hepatocytes"))





# =========================
# FILTERING, CONVERSION TO SEURAT OBJECT, AND DATA FORMATTING
# =========================








# The build is GRCh37, variant IDs are in the form chr_position_ref_alt. Not all variants are included for all donors - just heterozygous variants which had counts.
matrixdata_to_make = read.csv("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/validation_dataset/raw_counts.csv", row.names = 1)
coldata_to_make = read.table(file = '/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/validation_dataset/cell_metadata_cols.tsv', sep = '\t', header = TRUE, comment.char = "&", fill = T)


x = rownames(matrixdata_to_make)
s = lapply(strsplit(x, split = "_"), trimws)
positions = t(sapply(s, function(x) {length(x) <- max(lengths(s)); x}))
colnames(positions) = c("ensembl_id", "gene_short_name")

matrixdata_to_make = as.matrix(matrixdata_to_make)

row.names(coldata_to_make) = colnames(matrixdata_to_make)
row.names(positions) = rownames(matrixdata_to_make)

all_rna = CreateSeuratObject(matrixdata_to_make,
                                     project = "hepatocyte_like_cells",
                                     assay = "rna",
                                     names.delim = "",
                                     meta.data = coldata_to_make)



all_rna = all_rna[which(positions[,"gene_short_name"] %in% ccds$gene), ]




#stopping point for 24 july 2023

# #
# #
# #rna: get short names of all measured genes that are in ccds
# gene_translations = rbind(as.data.frame(rowData(ipsc_rna_matrix)),
#                           as.data.frame(rowData(dfendo_rna_matrix)),
#                           as.data.frame(rowData(hpendo_rna_matrix)),
#                           as.data.frame(rowData(immhp_rna_matrix)),
#                           as.data.frame(rowData(mathp_rna_matrix)))
#
# gene_translations = gene_translations[!duplicated(gene_translations), ]
#
# #
#
# #convert from mnocle to seurat object
# #rna
# for (i in 1:length(order_cell_types)){
#   cell_type_i_name = order_cell_types[i]
#
#   monocle_object_i = get(paste(order_cell_types[i], "_rna_matrix", sep = ""))
#
#   cell_metadata_i = colData(monocle_object_i)
#   cell_metadata_i$cell_type = cell_type_i_name
#
#   seurat_object_i = CreateSeuratObject(exprs(monocle_object_i),
#                                        project = "hepatocyte_like_cells",
#                                        assay = "rna",
#                                        names.delim = "",
#                                        meta.data = as.data.frame(cell_metadata_i))
#
#   Idents(seurat_object_i) = cell_type_i_name
#
#   assign(paste0(cell_type_i_name, "_rna"), seurat_object_i)
#
#
# }
# # atac
# for (i in 1:length(order_cell_types)){
#   cell_type_i_name = order_cell_types[i]
#
#   monocle_object_i = get(paste(order_cell_types[i], "_atac_input_cds", sep = ""))
#
#   cell_metadata_i = colData(monocle_object_i)
#   cell_metadata_i$cell_type = cell_type_i_name
#
#
#   signac_object_i = CreateChromatinAssay(exprs(monocle_object_i),
#                                          project = "hepatocyte_like_cells",
#                                          assay = "atac",
#                                          sep = c("_", "_"),
#                                          meta.data = as.data.frame(cell_metadata_i))
#
#
#   seurat_object_i = CreateSeuratObject(
#     counts = signac_object_i,
#     assay = "peaks",
#     meta.data = as.data.frame(cell_metadata_i))
#
#   Idents(seurat_object_i) = cell_type_i_name
#
#   assign(paste0(cell_type_i_name, "_atac"), seurat_object_i)
#
# }
#
#
# #rename cells in case duplicate names are present between cell types
# #rna
# ipsc_rna = RenameCells(ipsc_rna, add.cell.id = "ipsc")
# dfendo_rna = RenameCells(dfendo_rna, add.cell.id = "dfendo")
# hpendo_rna = RenameCells(hpendo_rna, add.cell.id = "hpendo")
# immhp_rna = RenameCells(immhp_rna, add.cell.id = "immhp")
# mathp_rna = RenameCells(mathp_rna, add.cell.id = "mathp")
# #atac
# ipsc_atac = RenameCells(ipsc_atac, add.cell.id = "ipsc")
# dfendo_atac = RenameCells(dfendo_atac, add.cell.id = "dfendo")
# hpendo_atac = RenameCells(hpendo_atac, add.cell.id = "hpendo")
# immhp_atac = RenameCells(immhp_atac, add.cell.id = "immhp")
# mathp_atac = RenameCells(mathp_atac, add.cell.id = "mathp")
#
#
# #merge all cell types into one seurat object
# #rna
# all_rna = merge(x = ipsc_rna,
#                 y = list(dfendo_rna,
#                          hpendo_rna,
#                          immhp_rna,
#                          mathp_rna))
#


#rna: remove those with less than x1 or more than x2 measured genes
all_rna = subset(all_rna, subset = (nFeature_rna > 200) & (nFeature_rna < 8000))
#rna: normalize and log transform
all_rna <- NormalizeData(all_rna, normalization.method = "LogNormalize", scale.factor = 10000)




#magic imputation
all_rna = magic(all_rna)



DefaultAssay(all_rna) <- "MAGIC_rna"

#
# #atac
# all_atac <- merge(x = ipsc_atac,
#                   y = list(dfendo_atac,
#                            hpendo_atac,
#                            immhp_atac,
#                            mathp_atac))
gene_translations = as.data.frame(positions)
colnames(gene_translations) = c("id", "gene_short_name")



all_rna@meta.data$plot_color = c("red", "yellow", "green", "blue")[match(all_rna@meta.data$day, unique(all_rna@meta.data$day))]


save(all_rna, ccds, gene_translations, file = "/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/validation_dataset/seurat_formatted_data_intermediate_before_preprocessing_trash_me_validation_dataset.RData")
print("FINISHED INTERMEDIATE SAVE!!!!!!!!!!!!!!!")


# =========================
# PREPROCESS RNA
# =========================

#rna: scatter plot of rna number of genes vs total gene counts
pdf("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/plots/validation_dataset/feature_scatter_plot_rna_ncount_nfeature.pdf")
FeatureScatter(all_rna, feature1 = "nCount_rna", feature2 = "nFeature_rna") +
  ggtitle("") +
  xlab("Total RNA Counts") +
  ylab("Total Measurable Genes") +
  theme(legend.title=element_blank())
dev.off()



#rna: find top x most variable genes
all_rna <- FindVariableFeatures(all_rna, selection.method = "vst", nfeatures = 400)


#rna: plot scatter plot of most variable genes based on standardized variance
pdf("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/plots/validation_dataset/variable_features_plot_rna.pdf")
variable_features_plot <- VariableFeaturePlot(all_rna, selection.method = "vst")
top_variable_features <- head(VariableFeatures(all_rna), 10)
x = top_variable_features
s = lapply(strsplit(x, split = "-"), trimws)
positions = t(sapply(s, function(x) {length(x) <- max(lengths(s)); x}))
colnames(positions) = c("ensembl_id", "gene_short_name")
LabelPoints(plot = variable_features_plot,
            points = top_variable_features,
            labels = positions[, "gene_short_name"],
            repel = TRUE)
dev.off()


#rna: scale data based on all genes (NOT just top x most variable)
all_rna_genes <- rownames(all_rna)
all_rna <- ScaleData(all_rna, features = all_rna_genes)
all_rna <- RunPCA(all_rna, features = VariableFeatures(all_rna))


# #rna: add full name to metadata
# all_rna@meta.data$cell_type_full_name = cell_type_translation[match(all_rna@meta.data$cell_type, cell_type_translation$short_name), "full_name"]


#rna: pca plot
pdf("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/plots/validation_dataset/pca_plot_rna.pdf")
buffer_space = 5
DimPlot(all_rna, reduction = "pca", group.by = "day") +
  NoLegend() +
  xlab("PC 1") +
  ylab("PC 2") +
  geom_mark_ellipse(x = all_rna@reductions[["pca"]]@cell.embeddings[,"PC_1"],
                    y = all_rna@reductions[["pca"]]@cell.embeddings[,"PC_2"],
                    aes(color = all_rna@meta.data$day,
                        label = all_rna@meta.data$day#, description = all_rna@meta.data$day
                        ),
                    label.fontsize = 5,
                    con.size = 0.5,
                    show.legend = F,
                    con.cap = 0) +
  xlim(c(min(all_rna@reductions[["pca"]]@cell.embeddings[,"PC_1"], na.rm=T)-buffer_space, max(all_rna@reductions[["pca"]]@cell.embeddings[,"PC_1"], na.rm=T)+buffer_space)) +
  ylim(c(min(all_rna@reductions[["pca"]]@cell.embeddings[,"PC_2"], na.rm=T)-buffer_space, max(all_rna@reductions[["pca"]]@cell.embeddings[,"PC_2"], na.rm=T)+buffer_space)) +
  ggtitle("")
dev.off()


#rna: top x pcs heatmap for a subset of top genes and cells
pdf("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/plots/validation_dataset/pc_heatmaps_plot_rna.pdf")
pcs_to_show = 9
DimHeatmap(all_rna, dims = 1:pcs_to_show, cells = 50, balanced = TRUE)
dev.off()


# pdf("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/plots/pc_jackstraw_plot_rna.pdf")
# all_rna <- JackStraw(all_rna, num.replicate = 100)
# all_rna <- ScoreJackStraw(all_rna, dims = 1:20)
# JackStrawPlot(all_rna, dims = 1:pcs_to_show)
# dev.off()

all_rna <- RunUMAP(all_rna, dims = 1:30, reduction.name = "umap.rna", reduction.key = "rnaUMAP_")


#atac: umap using lsi as reduction
pdf("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/plots/validation_dataset/umap_plot_rna.pdf")
buffer_space=5
DimPlot(all_rna, group.by = "day") +
  NoLegend() +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  geom_mark_ellipse(x = all_rna@reductions[["umap.rna"]]@cell.embeddings[,"rnaUMAP_1"],
                    y = all_rna@reductions[["umap.rna"]]@cell.embeddings[,"rnaUMAP_2"],
                    aes(color = all_rna@meta.data$day,
                        label = all_rna@meta.data$day#, description = all_rna@meta.data$day
                        ),
                    label.fontsize = 5,
                    expand = unit(0, "mm"),
                    con.size = 0.5,
                    show.legend = F,
                    con.cap = 0) +
  xlim(c(min(all_rna@reductions[["umap.rna"]]@cell.embeddings[,"rnaUMAP_1"], na.rm=T)-buffer_space, max(all_rna@reductions[["umap.rna"]]@cell.embeddings[,"rnaUMAP_1"], na.rm=T)+buffer_space)) +
  ylim(c(min(all_rna@reductions[["umap.rna"]]@cell.embeddings[,"rnaUMAP_2"], na.rm=T)-buffer_space, max(all_rna@reductions[["umap.rna"]]@cell.embeddings[,"rnaUMAP_2"], na.rm=T)+buffer_space)) +
  ggtitle("")
dev.off()



# =========================
# PREPROCESS ATAC
# =========================

# #atac: add gene annotation information generated from the file <atac_annotations.R
# load("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/atac_annotations.RData")
# Annotation(all_atac) = annotations


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
#
# #rna: add full name to metadata
# all_atac@meta.data$cell_type_full_name = cell_type_translation[match(all_atac@meta.data$cell_type, cell_type_translation$short_name), "full_name"]


#
#
# #atac: umap using lsi as reduction
# pdf("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/plots/validation_umap_lsi_plot_atac.pdf")
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
#



# =========================
# SAVE
# =========================

#save all objects
print("starting saving (normalization -> imputation -> variable features -> scaling -> dim reduction")
save(all_rna, gene_translations, ccds, file = "/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/validation_dataset/seurat_formatted_data_validation_dataset.RData")
print("done saving")



print("MASTER: done everything")





