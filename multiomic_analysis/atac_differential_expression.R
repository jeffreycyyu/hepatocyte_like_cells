
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






#differential expression analysis (total)
#rna: find all group markers compared to all other cells
all_markers_atac = Seurat::FindAllMarkers(all_atac, test.use = "bimod")
write.table(all_markers_atac, "/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/differentially_expressed/atac/all_markers_atac.txt", sep = "\t", row.names = T, col.names = T, quote = F)
print("atac: finished total differential expression analysis")
# #=====
# #atac: find all group markers compared to all other cells
# all_markers_atac = Seurat::FindAllMarkers(all_atac, test.use = "bimod")
# write.table(all_markers_atac, "/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/differentially_expressed/atac/all_markers_atac.txt", sep = "\t", row.names = T, col.names = T, quote = F)
# print("atac: finished total differential expression analysis")


print("FINISHED: differential expression analysis (total)")


#
#
# #==================================
#
#
#
# #differential expression analysis (pairwise)
# #rna
# print("rna: starting differential expression analysis")
# foreach(k = 1:nrow(differential_expression_comparisons)) %dopar% {
#   i = differential_expression_comparisons[k, 1]
#   j = differential_expression_comparisons[k, 2]
#
#   print(paste0("starting rna: ", i, " vs. ", j))
#
#   differentially_expressed_ij <- Seurat::FindMarkers(all_rna, ident.1 = i, ident.2 = j)
#   differentially_expressed_ij$gene_short_name = gene_translations[match(rownames(differentially_expressed_ij), gene_translations$id), "gene_short_name"]
#   write.table(differentially_expressed_ij, paste0("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/differentially_expressed/rna/rna_differentially_expressed_", i, "_", j, ".txt"), sep = "\t", row.names = T, col.names = T, quote = F)
#
#   print(paste0("finished rna: ", i, " vs. ", j))
# }
# print("rna: finished pairwise differential expression analysis")
# #======
# #atac
# print("atac: starting differential expression analysis")
# foreach(k = 1:nrow(differential_expression_comparisons)) %dopar% {
#
#   i = differential_expression_comparisons[k, 1]
#   j = differential_expression_comparisons[k, 2]
#
#   print(paste0("starting atac: ", i, " vs. ", j))
#
#   differentially_expressed_ij <- Seurat::FindMarkers(all_atac, ident.1 = i, ident.2 = j)
#   write.table(differentially_expressed_ij, paste0("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/differentially_expressed/atac/atac_differentially_expressed_", i, "_", j, ".txt"), sep = "\t", row.names = T, col.names = T, quote = F)
#
#   print(paste0("finished atac: ", i, " vs. ", j))
# }
# print("atac: finished pairwise differential expression analysis")
#
# print("FINISHED: rna and atac differential expression analysis (pairwise)")
#

#=============

print("MASTER: done everything")

