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
library(reshape2)



load("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/seurat_formatted_data.RData")


# atac_regions_locations = colsplit(rownames(all_atac), "\\-", names=c("chromosome", "start", "end"))
# atac_regions_locations$chromosome <- gsub("^chr", "", atac_regions_locations$chromosome)





hpa_tissue_cell_type_markers = read.table("/home/jyu/scratch/hepatocyte_like_cells/human_protein_atlas_marker_genes/normal_tissue.tsv", quote = "", header = T, sep = "\t")


#filter: Reliability = Approved
hpa_tissue_cell_type_markers_filtered = hpa_tissue_cell_type_markers[which(hpa_tissue_cell_type_markers$Reliability == "Approved"), ]
#filter: Level = High OR Low
hpa_tissue_cell_type_markers_filtered = hpa_tissue_cell_type_markers_filtered[c(which(hpa_tissue_cell_type_markers_filtered$Level == "High"),
                                                                                which(hpa_tissue_cell_type_markers_filtered$Level == "Low")), ]
#seperate high and lows
hpa_tissue_cell_type_markers_filtered_high = hpa_tissue_cell_type_markers_filtered[which(hpa_tissue_cell_type_markers_filtered$Level == "High"), ]
hpa_tissue_cell_type_markers_filtered_low = hpa_tissue_cell_type_markers_filtered[which(hpa_tissue_cell_type_markers_filtered$Level == "Low"), ]


# all_rna = all_rna[c(1:100, 5000:5100), c(1:100, 5000:5100, 10500:10600)]

ungrouped_rna_matrix = as.data.frame(t(as.matrix(GetAssayData(object = all_rna, assay = "MAGIC_rna", slot = "data"))))
ungrouped_rna_matrix$cell_type = all_rna@meta.data$cell_type
grouped_rna_matrix = ungrouped_rna_matrix %>%
  group_by(cell_type) %>%
  summarize_all(sum) %>%
  as.data.frame()


#construct plotting dataframe
rna_counts_all_markers_all_cells_high = data.frame(Gene = rep(hpa_tissue_cell_type_markers_filtered_high$Gene, nrow(cell_type_translation)),
                                                   Gene.name = rep(hpa_tissue_cell_type_markers_filtered_high$Gene.name, nrow(cell_type_translation)),
                                                   marker_Tissue = rep(hpa_tissue_cell_type_markers_filtered_high$Tissue, nrow(cell_type_translation)),
                                                   marker_Cell.type = rep(hpa_tissue_cell_type_markers_filtered_high$Cell.type, nrow(cell_type_translation)),
                                                   measurement_cell_type = rep(cell_type_translation$short_name, each=nrow(hpa_tissue_cell_type_markers_filtered_high)),
                                                   measurement_expression = NA)
high_row_indices = match(rna_counts_all_markers_all_cells_high$measurement_cell_type, grouped_rna_matrix$cell_type)
high_col_indices = match(rna_counts_all_markers_all_cells_high$Gene, colnames(grouped_rna_matrix))
rna_counts_all_markers_all_cells_high$measurement_expression = grouped_rna_matrix[cbind(high_row_indices, high_col_indices)]
rna_counts_all_markers_all_cells_high_tissue = na.omit(rna_counts_all_markers_all_cells_high[, c("measurement_cell_type", "marker_Tissue", "measurement_expression")]) %>%
  group_by(marker_Tissue, measurement_cell_type) %>%
  summarise(summed_measurement_expression = sum(as.numeric(measurement_expression), na.rm = T)) %>%
  as.data.frame()
rna_counts_all_markers_all_cells_high_tissue$measurement_cell_type <- factor(rna_counts_all_markers_all_cells_high_tissue$measurement_cell_type, levels = cell_type_translation$short_name)

#plotting
proportion_rna_counts_all_tissue_markers_all_cells_high_plot = ggplot(rna_counts_all_markers_all_cells_high_tissue,
                                                    aes(fill = marker_Tissue,
                                                        colour = marker_Tissue,
                                                        y = summed_measurement_expression,
                                                        x = measurement_cell_type,
                                                        alpha = 0.9)) +
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank() #, legend.title = element_text(size = 3), legend.text = element_text(size = 3)
        ) +
  xlab("Cell Type") +
  ylab("Proportion of Positive Marker Gene RNA Expressed") +
  guides(fill=guide_legend(title="Tissue", override.aes = list(size = 0.1)),
         color = "none",
         alpha = "none")
#save plot
ggsave("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/plots/proportion_rna_counts_all_tissue_markers_all_cells_high_plot.pdf", plot = proportion_rna_counts_all_tissue_markers_all_cells_high_plot, width = 10, height = 10)
print("FINISHED positive marker genes !!!!!")



#construct plotting dataframe
rna_counts_all_markers_all_cells_low = data.frame(Gene = rep(hpa_tissue_cell_type_markers_filtered_low$Gene, nrow(cell_type_translation)),
                                                  Gene.name = rep(hpa_tissue_cell_type_markers_filtered_low$Gene.name, nrow(cell_type_translation)),
                                                  marker_Tissue = rep(hpa_tissue_cell_type_markers_filtered_low$Tissue, nrow(cell_type_translation)),
                                                  marker_Cell.type = rep(hpa_tissue_cell_type_markers_filtered_low$Cell.type, nrow(cell_type_translation)),
                                                  measurement_cell_type = rep(cell_type_translation$short_name, each=nrow(hpa_tissue_cell_type_markers_filtered_low)),
                                                  measurement_expression = NA)
low_row_indices = match(rna_counts_all_markers_all_cells_low$measurement_cell_type, grouped_rna_matrix$cell_type)
low_col_indices = match(rna_counts_all_markers_all_cells_low$Gene, colnames(grouped_rna_matrix))
rna_counts_all_markers_all_cells_low$measurement_expression = grouped_rna_matrix[cbind(low_row_indices, low_col_indices)]
rna_counts_all_markers_all_cells_low_tissue = na.omit(rna_counts_all_markers_all_cells_low[, c("measurement_cell_type", "marker_Tissue", "measurement_expression")]) %>%
  group_by(marker_Tissue, measurement_cell_type) %>%
  summarise(summed_measurement_expression = sum(as.numeric(measurement_expression), na.rm = T)) %>%
  as.data.frame()
rna_counts_all_markers_all_cells_low_tissue$measurement_cell_type <- factor(rna_counts_all_markers_all_cells_low_tissue$measurement_cell_type, levels = cell_type_translation$short_name)

#plotting
proportion_rna_counts_all_tissue_markers_all_cells_low_plot = ggplot(rna_counts_all_markers_all_cells_low_tissue,
                                                          aes(fill = marker_Tissue,
                                                              colour = marker_Tissue,
                                                              y = summed_measurement_expression,
                                                              x = measurement_cell_type,
                                                              alpha = 0.9)) +
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank() #, legend.title = element_text(size = 3), legend.text = element_text(size = 3)
  ) +
  xlab("Cell Type") +
  ylab("Proportion of Negative Marker Gene RNA Expressed") +
  guides(fill=guide_legend(title="Tissue", override.aes = list(size = 0.1)),
         color = "none",
         alpha = "none")
#save plot
ggsave("/home/jyu/scratch/hepatocyte_like_cells/multiomic_analysis/plots/proportion_rna_counts_all_tissue_markers_all_cells_low_plot.pdf", plot = proportion_rna_counts_all_tissue_markers_all_cells_low_plot, width = 10, height = 10)
print("FINISHED negative marker genes !!!!!")



print("FINISHED EVERYTHING !!!!!")
