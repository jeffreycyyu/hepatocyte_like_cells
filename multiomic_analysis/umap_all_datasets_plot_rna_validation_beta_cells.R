rm(list = ls())


# PACKAGES
# library(monocle3)
library(Seurat)
library(Signac)
library(Rmagic)
library(dplyr)
library(ggplot2)
library(ggforce)
library(foreach)
library(doParallel)


load("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/seurat_formatted_data.RData")


#CCDS genes
ccds = read.table("/home/jyu/scratch/hepatocyte_like_cells/rna_analysis/ccds_genes/CCDS.current.txt", sep = "\t")
colnames(ccds) = c("chromosome", "nc_accession", "gene", "gene_id", "ccds_id", "ccds_status", "cds_strand", "cds_from", "cds_to", "cds_locations", "match_type")


# LOAD TABULATED DATA
load("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/seurat_formatted_data.RData")
my_dataset_rna = all_rna
my_dataset_rna = my_dataset_rna[which(rownames(my_dataset_rna) != ""), ]
#===
load("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/validation_dataset/seurat_formatted_data_validation_dataset.RData")
validation_dataset_rna = all_rna
rownames(validation_dataset_rna@assays[["MAGIC_rna"]]@data) = substr(rownames(validation_dataset_rna), 1, regexpr("\\-", rownames(validation_dataset_rna))-1)
validation_dataset_rna = validation_dataset_rna[which(rownames(validation_dataset_rna) != ""), ]
Idents(validation_dataset_rna) = validation_dataset_rna@meta.data[["day"]]
#====
load("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/beta_cell_dataset/seurat_formatted_data_beta_cell_dataset.RData")
beta_cell_dataset_rna = all_rna
rownames(beta_cell_dataset_rna@assays[["MAGIC_rna"]]@data) = gene_translations[match(rownames(beta_cell_dataset_rna), gene_translations$gene_short_name), "id"]
beta_cell_dataset_rna = beta_cell_dataset_rna[which(rownames(beta_cell_dataset_rna) != ""), ]


# order_cell_types = c("ipsc", "dfendo", "hpendo", "immhp", "mathp")


rna_combined = merge(x = my_dataset_rna, y = c(validation_dataset_rna, beta_cell_dataset_rna), add.cell.ids = c("my_dataset_rna", "validation_dataset_rna", "beta_cell_dataset_rna"), project = "hepatocyte_like_cells")

save(rna_combined, gene_translations, file = "/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/compare_all_datasets_merged_non_preprocessed.RData")



rna_combined = rna_combined[which(rownames(rna_combined) %in% gene_translations[match(ccds$gene, gene_translations$gene_short_name), "id"]), ]




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
rna_combined = subset(rna_combined, subset = (nFeature_rna > 200) & (nFeature_rna < 8000))
#rna: normalize and log transform
rna_combined <- NormalizeData(rna_combined, normalization.method = "LogNormalize", scale.factor = 10000)




#magic imputation
rna_combined = magic(rna_combined)



DefaultAssay(rna_combined) <- "MAGIC_rna"







rna_combined <- FindVariableFeatures(rna_combined, selection.method = "vst", nfeatures = 400)
all_rna_genes <- rownames(rna_combined)
rna_combined <- ScaleData(rna_combined, features = all_rna_genes)
rna_combined <- RunPCA(rna_combined, features = VariableFeatures(rna_combined))

print("done pca")

rna_combined <- RunUMAP(rna_combined, dims = 1:30, reduction.name = "umap.rna", reduction.key = "rnaUMAP_")


pdf("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/plots/compare_datasets/umap_plot_compare_all_datasets_rna.pdf")
buffer_space=5
DimPlot(rna_combined) +#, group.by = "day") +
  NoLegend() +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  # geom_mark_ellipse(x = all_rna@reductions[["umap.rna"]]@cell.embeddings[,"rnaUMAP_1"],
  #                   y = all_rna@reductions[["umap.rna"]]@cell.embeddings[,"rnaUMAP_2"],
  #                   aes(color = all_rna@meta.data$day,
  #                       label = all_rna@meta.data$day#, description = all_rna@meta.data$day
  #                       ),
  #                   label.fontsize = 5,
  #                   expand = unit(0, "mm"),
  #                   con.size = 0.5,
  #                   show.legend = F,
  #                   con.cap = 0) +
  xlim(c(min(rna_combined@reductions[["umap.rna"]]@cell.embeddings[,"rnaUMAP_1"], na.rm=T)-buffer_space, max(rna_combined@reductions[["umap.rna"]]@cell.embeddings[,"rnaUMAP_1"], na.rm=T)+buffer_space)) +
  ylim(c(min(rna_combined@reductions[["umap.rna"]]@cell.embeddings[,"rnaUMAP_2"], na.rm=T)-buffer_space, max(rna_combined@reductions[["umap.rna"]]@cell.embeddings[,"rnaUMAP_2"], na.rm=T)+buffer_space)) +
  ggtitle("")
dev.off()

# save(rna_combined, file = "/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/compare_all_datasets_merged.RData")
# beta markers: CD20, CD19, CD45RO, CD80, CD86, CD27, CD24, CD79a

save(rna_combined, file = "/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/compare_all_datasets_merged.RData")


print("FINISHED EVERYTHING MASTER END!!!")
