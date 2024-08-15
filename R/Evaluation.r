# This script is used to evaluate the cell type proportions and marker genes of Myeloid/T subtypes in the spatial transcriptomics data using scRNA-seq data as a reference.




library(Seurat)
subtype <- temp@cell_metadata$cell$rna$cell_type_isML_updated.sub
# T_cell <- c("Progenitor_like_T","CD8_TRM","CD8_Tc","CD8_Tn","Plasma","Tcm","Tem","Treg","ILC1/NK","Th2/ILC2")
# select_id <- which(subtype %in% T_cell)
count <- as.matrix(temp@expression$cell$rna$normalized@exprMat)



T_data = read.csv("subtype_T_marker_expression_filtered_matrix.csv")
T_cell = T_data$X
T_cell

Myeloid_data = read.csv("subtype_Myeloid_marker_expression_filtered_matrix.csv")
Myeloid_cell = Myeloid_data$X
Myeloid_cell



metadata <- data.frame(cell_id = colnames(count), cell_types = subtype)
rownames(metadata) <- metadata$`cell_id`
seurat <- CreateSeuratObject(counts = count, meta.data = metadata)

seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
seurat@meta.data$cell_types <- as.factor(seurat@meta.data$cell_types)
Idents(seurat) <- seurat@meta.data$cell_types

markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Extracting top 5 marker genes for each cell type, sorted by avg_log2FC in descending order
top5_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# Creating a new column to rank genes within each cell type, sorting from high to low
top5_markers <- top5_markers %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  mutate(gene_rank = row_number())

# Formatting the data for the bubble plot
# Extracting log fold change and proportion of cells expressing the marker genes
bubble_data <- top5_markers %>% 
  select(gene, cell_type = cluster, logFC = avg_log2FC, proportion = pct.1, gene_rank)

# Sorting genes for left-diagonal layout
bubble_data <- bubble_data %>%
  arrange(desc(cell_type), desc(gene_rank))

# Converting gene names to factors to maintain the order in the plot
bubble_data$gene <- factor(bubble_data$gene, levels = rev(unique(bubble_data$gene)))

# Creating the bubble plot
bubble_plot <- ggplot(bubble_data, aes(x = cell_type, y = gene, size = proportion, color = logFC)) +
  geom_point() +
  scale_color_gradient2(low = "white", mid = "grey", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Marker Genes of Myeloid Subtypes",
       x = "Cell Type",
       y = "Gene",
       size = "Proportion of Cells Expressing",
       color = "Log Fold Change")

# Display the bubble plot
#print(bubble_plot)

# Save the bubble plot
ggsave("Myeloid_bubble_plot_filter.png", plot = bubble_plot, width = 10, height = 8)




