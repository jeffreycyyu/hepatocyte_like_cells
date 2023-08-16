

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



master_gene_check = as.data.frame(matrix(NA, 0, 23))

colnames(master_gene_check) = c(paste0("atac_", colnames(atac_differentially_expressed)), paste0("annotation_", colnames(region_sites)))

for (chr in chromosomes_to_keep){

  print(paste0("starting chromsome: ", chr))

  atac_differentially_expressed_chr_i <- makeGRangesFromDataFrame(atac_differentially_expressed[which(atac_differentially_expressed[, "seqnames"] == paste0("chr", chr)),], keep.extra.columns = T)
  region_sites_chr_i <- makeGRangesFromDataFrame(region_sites[region_sites[,"seqnames"] == paste0('chr',chr),], keep.extra.columns = T)


  hits <- findOverlaps(query = atac_differentially_expressed_chr_i, subject = region_sites_chr_i, ignore.strand = TRUE)
  #
  #
  # sequences_differential_atac_regions <- Map(function(num1, num2) seq(as.numeric(num1), as.numeric(num2)),
  #                              atac_differentially_expressed_chr_i[,which(colnames(atac_differentially_expressed_chr_i) == "start")],
  #                              atac_differentially_expressed_chr_i[,which(colnames(atac_differentially_expressed_chr_i) == "stop")])
  #
  #
  # sequences_gene_regions <- Map(function(num1, num2) seq(as.numeric(num1), as.numeric(num2)),
  #                               region_sites_chr_i[,which(colnames(region_sites_chr_i) == "start")],
  #                               region_sites_chr_i[,which(colnames(region_sites_chr_i) == "end")])
  #
  #
  #
  # intersected_positions <- do.call("cbind", lapply(sequences_gene_regions, function(x) lapply(sequences_differential_atac_regions, function(y) length(intersect(x,y)))))
  #
  #
  # intersection_matrix_indices <- which(intersected_positions != 0, arr.ind = TRUE)
  #
  # colnames(intersection_matrix_indices) <- c("DMR index (use these as row indices for atac_differentially_expressed_chr_i)", "gene region index (use these as row indices for region_sites_chr_i)")
  #
  #
  # gene_check <- cbind(rep(chr, nrow(intersection_matrix_indices)), #chromosome
  #                     region_sites_chr_i[unlist(intersection_matrix_indices[,2]), which(colnames(region_sites_chr_i) == "start")], #gene region start site
  #                     region_sites_chr_i[unlist(intersection_matrix_indices[,2]), which(colnames(region_sites_chr_i) == "end")], #gene region stop site
  #                     region_sites_chr_i[unlist(intersection_matrix_indices[,2]), which(colnames(region_sites_chr_i) == "symbol")], #gene symbol
  #                     region_sites_chr_i[unlist(intersection_matrix_indices[,2]), which(colnames(region_sites_chr_i) == "type")], #gene region type
  #                     atac_differentially_expressed_chr_i[unlist(intersection_matrix_indices[,1]), which(colnames(atac_differentially_expressed_chr_i) == "start")], #DMR start site
  #                     atac_differentially_expressed_chr_i[unlist(intersection_matrix_indices[,1]), which(colnames(atac_differentially_expressed_chr_i) == "stop")], #DMR stop site
  #                     atac_differentially_expressed_chr_i[unlist(intersection_matrix_indices[,1]), which(colnames(atac_differentially_expressed_chr_i) == "p_val_adj")], #DMR disease p value
  #                     data.frame(n_bp_intersected = unlist(intersected_positions[cbind(intersection_matrix_indices[,1], intersection_matrix_indices[,2])]))) #number of intersected positions (i.e., overlap size, NOTE not representative of number of CpG regions)

  # colnames(gene_check) <- c("chr", "searched_gene_start", "searched_gene_end", "searched_gene_symbol", "searched_gene_region_type", "atac_start", "atac_end", "atac_p_value_adj", "n_bp_intersected")

  atac_differentially_expressed_chr_i = as.data.frame(atac_differentially_expressed_chr_i)
  region_sites_chr_i = as.data.frame(region_sites_chr_i)

  gene_check = cbind(atac_differentially_expressed_chr_i[queryHits(hits),], region_sites_chr_i[subjectHits(hits),])


  master_gene_check = rbind(master_gene_check, gene_check)

  print(paste0("finished chromsome: ", chr))

}



print(paste0("number of distinct genic regions overlapping differnetially accessed chromatin regions: ", nrow(master_gene_check)))

print("starting save")

save.image(file = "/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/differentially_expressed/atac/all_markers_atac_annotated_hg19.RData")

print("ending save")

print("FINISHED ALL MASTER !!!")
