

rm(list = ls())


# PACKAGES
# library(monocle3)
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(ggforce)
library(foreach)
library(doParallel)
library(stringr)
library(annotatr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)



atac_differentially_expressed = read.table("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/differentially_expressed/atac/all_markers_atac.txt", header = T, sep = "\t")
colnames(atac_differentially_expressed) = c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "region")
atac_differentially_expressed[, c("seqnames", "start", "end")] = str_split_fixed(atac_differentially_expressed$region, "-", 3)
atac_differentially_expressed$stand = NA
# atac_differentially_expressed$chr = gsub("chr", "", atac_differentially_expressed$chr)



annotation_shortcuts <- c("hg19_genes_1to5kb",
                          "hg19_genes_promoters",
                          "hg19_genes_cds",
                          "hg19_genes_5UTRs",
                          "hg19_genes_exons",
                          "hg19_genes_firstexons",
                          "hg19_genes_introns",
                          "hg19_genes_3UTRs",
                          "hg19_genes_intergenic") # head(builtin_annotations()[grepl("hg19", builtin_annotations())], n=11)[c(1:7,10:11)]


annotations = build_annotations(genome = 'hg19', annotations = annotation_shortcuts)
region_sites <- as.data.frame(annotations)


chromosomes_to_keep = c(1:22, "X", "Y")







