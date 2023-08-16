
#======
#======
#======
# RUN HERE, NO ".sh" FILE PRESENT
#======
#======
#======



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





#CCDS genes
ccds = read.table("/home/jyu/scratch/hepatocyte_like_cells/rna_analysis/ccds_genes/CCDS.current.txt", sep = "\t")
colnames(ccds) = c("chromosome", "nc_accession", "gene", "gene_id", "ccds_id", "ccds_status", "cds_strand", "cds_from", "cds_to", "cds_locations", "match_type")


# LOAD TABULATED DATA
load("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/seurat_formatted_data.RData")
my_dataset_rna = all_rna
load("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/validation_dataset/seurat_formatted_data_validation_dataset.RData")
validation_dataset_rna = all_rna
rownames(validation_dataset_rna@assays[["MAGIC_rna"]]@data) = substr(rownames(validation_dataset_rna), 1, regexpr("\\-", rownames(validation_dataset_rna))-1)
Idents(validation_dataset_rna) = validation_dataset_rna@meta.data[["day"]]
load("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/beta_cell_dataset/seurat_formatted_data_beta_cell_dataset.RData")
beta_cell_dataset_rna = all_rna
rownames(beta_cell_dataset_rna@assays[["MAGIC_rna"]]@data) = gene_translations[match(rownames(beta_cell_dataset_rna), gene_translations$gene_short_name), "id"]

# order_cell_types = c("ipsc", "dfendo", "hpendo", "immhp", "mathp")



# UREA CYCLE GENES
target_genes_translation <- rbind(data.frame(ensembl_id = "ENSMUSG00000025991",
                                             gene_name = "carbamoyl-phosphate synthetase 1 (CPS1)",
                                             gene_abbreviation = "CPS1"),
                                  data.frame(ensembl_id = "ENSG00000036473",
                                             gene_name = "ornithine transcarbamylase (OTC)",
                                             gene_abbreviation = "OTC"),
                                  data.frame(ensembl_id = "ENSG00000130707",
                                             gene_name = "argininosuccinate synthase 1 (ASSN1)",
                                             gene_abbreviation = "ASSN1"),
                                  data.frame(ensembl_id = "ENSG00000126522",
                                             gene_name = "argininosuccinate lyase (ASL)",
                                             gene_abbreviation = "ASL"),
                                  data.frame(ensembl_id = "ENSG00000118520",
                                             gene_name = "arginase 1 (AGR1)",
                                             gene_abbreviation = "AGR1"),
                                  data.frame(ensembl_id = "ENSG00000161653",
                                             gene_name = "N-acetylglutamate synthase (NAGS)",
                                             gene_abbreviation = "NAGS"),
                                  data.frame(ensembl_id = "ENSMUSG00000031482",
                                             gene_name = "mitochondrial carrier ornithine transporter (SLC25A15)",
                                             gene_abbreviation = "SLC25A15"),
                                  data.frame(ensembl_id = "ENSG00000004864",
                                             gene_name = "citrin (SLC25A13)",
                                             gene_abbreviation = "SLC25A13"))





all_datasets_names = c("this_study", "validation", "beta_cell")
all_datasets_data = list(my_dataset_rna, validation_dataset_rna, beta_cell_dataset_rna)

for (i in 1:length(all_datasets_data)){

  dataset_i_name = all_datasets_names[i]
  rna_dataset_to_loop = all_datasets_data[[i]]

  # gene_short_name_id_interest_i =  target_genes_translation[which(target_genes_translation$ensembl_id %in% rownames(rna_dataset_to_loop)), "gene_abbreviation"]
  gene_full_name_id_interest_i =  target_genes_translation[which(target_genes_translation$ensembl_id %in% rownames(rna_dataset_to_loop)), "gene_name"]
  ensembl_id_interest_i =  target_genes_translation[which(target_genes_translation$ensembl_id %in% rownames(rna_dataset_to_loop)), "ensembl_id"]
  plot_list = FeaturePlot(object = rna_dataset_to_loop, reduction = "umap.rna", features = ensembl_id_interest_i, label = T, repel = T)
  for (j in 1:length(ensembl_id_interest_i)){
    plot_list[[j]] = plot_list[[j]] +
      labs(title = gene_full_name_id_interest_i[j]) +
      xlab("UMAP 1") +
      ylab("UMAP 2")
  }
  ggsave(paste0("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/plots/compare_datasets/amonia_genes_feature_plots_", dataset_i_name, ".pdf"), plot = plot_list, width = 12, height = 18)
  print(paste0("done dataset: ", dataset_i_name))
}





print("DONE ALL PLOTS !!!!!")
