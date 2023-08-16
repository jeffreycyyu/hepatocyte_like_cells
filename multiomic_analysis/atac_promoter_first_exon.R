#######
#
#RUN ON RSTUDIO, NO .SH SCRIPT
#
#######


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
library(fgsea)
library(ggrepel)
library(data.table)
library(gridExtra)

load("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/differentially_expressed/atac/all_markers_atac_annotated_hg19.RData")

order_cell_types = c("ipsc", "dfendo", "hpendo", "immhp", "mathp")

promoter_first_exons_genes = master_gene_check[which(master_gene_check$type %in% c("hg19_genes_promoters", "hg19_genes_firstexons")), ]
promoter_first_exons_genes = promoter_first_exons_genes[-is.na(promoter_first_exons_genes$symbol),]
promoter_first_exons_genes = na.omit(promoter_first_exons_genes[, -which(colnames(promoter_first_exons_genes) == "stand")])
promoter_first_exons_genes$cluster = factor(promoter_first_exons_genes$cluster, levels=order_cell_types)





#START GSEA========================
# promoter_first_exons_genes <- aggregate( . ~ symbol, promoter_first_exons_genes, mean)
promoter_first_exons_genes_formatted = promoter_first_exons_genes$avg_log2FC
names(promoter_first_exons_genes_formatted) = promoter_first_exons_genes$symbol

gene_set_info = gmtPathways("/home/jyu/scratch/hepatocyte_like_cells/msigdb_gene_sets/c2.all.v2023.1.Hs.symbols.gmt")

hpa_tissue_cell_type_markers_filtered_high <- data.frame(Tissue = rep(names(gene_set_info), sapply(gene_set_info, length)),
                                                         Gene.name = unlist(gene_set_info))
# hpa_tissue_cell_type_markers_filtered_high$Gene = gene_translations[match(hpa_tissue_cell_type_markers_filtered_high$Gene.name, gene_translations$gene_short_name), "id"]
# hpa_tissue_cell_type_markers_filtered_high$Tissue = gsub(pattern="HALLMARK_",
#                                                          replacement="",
#                                                          x=hpa_tissue_cell_type_markers_filtered_high$Tissue)

gsea_results = fgsea(pathways = gene_set_info,
                     stats    = promoter_first_exons_genes_formatted,
                                 eps      = 0.0,
                                 minSize  = 15,
                                 maxSize  = 500)

gsea_results = gsea_results[order(gsea_results$padj, decreasing = F),]
pdf("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/plots/atac_differentially_expressed_annotated_table_gsea.pdf", height=4, width=12)

grid.table(gsea_results[1:10, c("pathway", "padj", "log2err", "NES", "size")])

dev.off()


#END GSEA========================

# labels
promoter_first_exons_genes_labels = as.data.frame(matrix(NA, 0, ncol(promoter_first_exons_genes)))
colnames(promoter_first_exons_genes_labels) = colnames(promoter_first_exons_genes)
for (cluster_i in order_cell_types){
  promoter_first_exons_genes_labels_i = promoter_first_exons_genes[which(promoter_first_exons_genes$cluster == cluster_i), ]
  promoter_first_exons_genes_labels_i_min = promoter_first_exons_genes_labels_i[which.min(promoter_first_exons_genes_labels_i$avg_log2FC), ]
  promoter_first_exons_genes_labels_i_max = promoter_first_exons_genes_labels_i[which.max(promoter_first_exons_genes_labels_i$avg_log2FC), ]
  promoter_first_exons_genes_labels = rbind(promoter_first_exons_genes_labels, rbind(promoter_first_exons_genes_labels_i_min, promoter_first_exons_genes_labels_i_max))
}



volcano_plot <- ggplot(data=promoter_first_exons_genes,
            aes(x=avg_log2FC, y=-log10(p_val_adj),
                col=type)) +
  geom_point() +
  facet_wrap(~cluster) +
  theme_minimal() +
  geom_vline(xintercept=0, col="grey") +
  geom_hline(yintercept=-log10(0.05), col="grey") +
  ylab("Adjusted -log10(P_Value)") +
  xlab("Average log2(Fold-Change)") +
  xlim(-max(abs(promoter_first_exons_genes$avg_log2FC)+0.5), max(abs(promoter_first_exons_genes$avg_log2FC)+0.5)) +
  scale_color_discrete(name = "hg19 Genic Region", labels = c("first_exon", "promoter")) +
  geom_label_repel(data = promoter_first_exons_genes_labels, aes(label = symbol), min.segment.length = 0, label.size = 0.1, color = "black")



ggsave(plot = volcano_plot, file = "/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/plots/atac_differentially_expressed_annotated_volcano_plot.pdf", height = 6, width = 10)
