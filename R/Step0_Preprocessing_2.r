# This script uses the Seurat package to identify marker genes for each cell type in the spatial transcriptomics dataset.


# parse input parameters from linux scripts
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Please supply the results_folder, path of scRNA-seq RDS file and giotto_object_path for tiling as command line arguments!", call.=FALSE)
} 

results_folder <- args[1]
scRNA_path <- args[2]
giotto_object_path <- args[3]

## Load the scRNA-seq reference data
spatial = readRDS(scRNA_path)

## Load the Giotto object
if(!"Giotto" %in% installed.packages()) {
  devtools::install_github("drieslab/Giotto@suite")
}
library(Giotto)
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table,RColorBrewer,ggplot2,dplyr,hrbrthemes,viridis,ggpubr,sjmisc) 
genv_exists = checkGiottoEnvironment()
if(!genv_exists){
  installGiottoEnvironment()
}
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
python_path = NULL
xenium_gobj = loadGiotto(giotto_object_path)
setwd(results_folder)
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  return_plot = FALSE)
xenium_gobj = replaceGiottoInstructions(xenium_gobj, instructions = instrs)



# Obtain the marker_gene * main_celltype matrix
## select the overlapped genes
xenium_gene <- rownames(xenium_gobj@expression$cell$rna$raw@exprMat)
scrna_gene <- rownames(spatial@assays$RNA@data)
intersect_genes <- intersect(xenium_gene, scrna_gene)
scrna_subset <- subset(spatial, features = intersect_genes)
library(Seurat)
all_markers_global <- FindAllMarkers(scrna_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, group.by = "celltype.global")

## Extract the unique marker genes
marker_genes <- unique(all_markers_global$gene)
## Get the expression data (assay: RNA, slot: data)
expression_matrix <- scrna_subset@assays$RNA@data
## Subset the expression matrix to include only marker genes
expression_matrix_markers <- as.matrix(expression_matrix)
expression_matrix_markers <- expression_matrix_markers[marker_genes,]
## Extract cell type information
cell_types <- scrna_subset@meta.data$celltype.global
## Create a matrix to store the bulk-level expression
bulk_expression_matrix <- matrix(0, nrow = length(unique(cell_types)), ncol = length(marker_genes))
rownames(bulk_expression_matrix) <- unique(cell_types)
colnames(bulk_expression_matrix) <- marker_genes
## Aggregate expression data by cell type
for (cell_type in unique(cell_types)) {
  cells_in_type <- which(cell_types == cell_type)
  bulk_expression_matrix[cell_type, ] <- rowMeans(expression_matrix_markers[, cells_in_type])
}
## Convert to data frame for easier saving
bulk_expression_df <- as.data.frame(bulk_expression_matrix)
bulk_expression_df$celltype <- rownames(bulk_expression_matrix)
## Reorder columns to have cell types as the first column
bulk_expression_df <- bulk_expression_df %>%
  select(celltype, everything())
## Save the matrix as a CSV file
write.csv(bulk_expression_matrix_df, "global_ct_marker_expression_matrix.csv", row.names = TRUE)





cell_types_with_subtypes <- scrna_subset@meta.data %>%
  group_by(`celltype.global`) %>%
  summarise(subtype_count = n_distinct(`celltype`)) %>%
  filter(subtype_count > 1)
cell_types_with_subtypes

# only T cells, Myeloid cells, and keratinocytes cells have subtypes
keratinocytes_obj <- subset(scrna_subset, subset = `celltype.global` == "Keratinocytes")
Keratinocytes_markers_subtype <- FindAllMarkers(keratinocytes_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, group.by = "celltype")
# Extract the unique marker genes
marker_genes <- unique(Keratinocytes_markers_subtype$gene)
# Get the expression data (assay: RNA, slot: data)
expression_matrix <- keratinocytes_obj@assays$RNA@data
# Subset the expression matrix to include only marker genes
expression_matrix_markers <- as.matrix(expression_matrix)
expression_matrix_markers <- expression_matrix_markers[marker_genes,]
# Extract cell type information
cell_types <- keratinocytes_obj@meta.data$celltype
# Create a matrix to store the bulk-level expression
bulk_expression_matrix <- matrix(0, nrow = length(unique(cell_types)), ncol = length(marker_genes))
rownames(bulk_expression_matrix) <- unique(cell_types)
colnames(bulk_expression_matrix) <- marker_genes
rownames(bulk_expression_matrix)[5] <- "Empty"
# Aggregate expression data by cell type
for (cell_type in unique(cell_types)) {
  if(cell_type == ""){
    cells_in_type <- which(cell_types == cell_type)
    bulk_expression_matrix["Empty", ] <- rowMeans(expression_matrix_markers[, cells_in_type])
  }
  else{
    cells_in_type <- which(cell_types == cell_type)
    bulk_expression_matrix[cell_type, ] <- rowMeans(expression_matrix_markers[, cells_in_type])
  }
}
# Convert to data frame for easier saving
bulk_expression_df <- as.data.frame(bulk_expression_matrix)
bulk_expression_df$celltype <- rownames(bulk_expression_matrix)
# Reorder columns to have cell types as the first column
bulk_expression_df <- bulk_expression_df %>%
  select(celltype, everything())
write.csv(bulk_expression_df, "subtype_keratinocyte_marker_expression_matrix.csv", row.names = TRUE)



Myeloid_obj <- subset(scrna_subset, subset = `celltype.global` ==  "Myeloid Cells")
Myeloid_obj <- subset(Myeloid_obj, cells = which(Myeloid_obj@meta.data$celltype != ""))
Myeloid_obj@meta.data$celltype <- gsub(" ", "_", Myeloid_obj@meta.data$celltype)
Myeloid_obj@meta.data$celltype <- gsub("-", "_", Myeloid_obj@meta.data$celltype)

# Set the cell identities to the standardized celltype column
Idents(Myeloid_obj) <- Myeloid_obj@meta.data$celltype

# Verify the standardized levels (unique identities) in the object
print(levels(Myeloid_obj))
Myeloid_markers_subtype <- FindAllMarkers(Myeloid_obj, 
                                          only.pos = TRUE, 
                                          min.pct = 0.25, 
                                          logfc.threshold = 0.25, 
                                          group.by = "celltype")


# Extract the unique marker genes
marker_genes <- unique(Myeloid_markers_subtype$gene)

# Get the expression data (assay: RNA, slot: data)
expression_matrix <- Myeloid_obj@assays$RNA@data

# Subset the expression matrix to include only marker genes
expression_matrix_markers <- as.matrix(expression_matrix)
expression_matrix_markers <- expression_matrix_markers[marker_genes,]

# Extract cell type information
cell_types <- Myeloid_obj@meta.data$celltype

# Create a matrix to store the bulk-level expression
bulk_expression_matrix <- matrix(0, nrow = length(unique(cell_types)), ncol = length(marker_genes))
rownames(bulk_expression_matrix) <- unique(cell_types)
colnames(bulk_expression_matrix) <- marker_genes
# rownames(bulk_expression_matrix)[5] <- "Empty"

# Aggregate expression data by cell type
for (cell_type in unique(cell_types)) {
  if(cell_type == ""){
    cells_in_type <- which(cell_types == cell_type)
    bulk_expression_matrix["Empty", ] <- rowMeans(expression_matrix_markers[, cells_in_type])
  }
  else{
    cells_in_type <- which(cell_types == cell_type)
    bulk_expression_matrix[cell_type, ] <- rowMeans(expression_matrix_markers[, cells_in_type])
  }
}

# Convert to data frame for easier saving
bulk_expression_df <- as.data.frame(bulk_expression_matrix)
bulk_expression_df$celltype <- rownames(bulk_expression_matrix)

# Reorder columns to have cell types as the first column
bulk_expression_df <- bulk_expression_df %>%
  select(celltype, everything())

write.csv(bulk_expression_df, "subtype_Myeloid_marker_expression_matrix.csv", row.names = TRUE)



# Count the occurrences of each gene
gene_counts <- table(Myeloid_markers_subtype$gene)
# Filter genes that occur more than once
filtered_genes <- names(gene_counts[gene_counts == 1])

bulk_expression_filter_df = bulk_expression_df[,filtered_genes]

write.csv(bulk_expression_filter_df, "subtype_Myeloid_marker_expression_filtered_matrix.csv", row.names = TRUE)





T_obj <- subset(scrna_subset, subset = `celltype.global` == "T Cells")
T_obj <- subset(T_obj, cells = which(T_obj@meta.data$celltype != ""))
T_obj@meta.data$celltype <- gsub(" ", "_", T_obj@meta.data$celltype)
T_obj@meta.data$celltype <- gsub("-", "_", T_obj@meta.data$celltype)
Idents(T_obj) <- T_obj@meta.data$celltype
T_markers_subtype <- FindAllMarkers(T_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, group.by = "celltype")

# Extract the unique marker genes
marker_genes <- unique(T_markers_subtype$gene)

# Get the expression data (assay: RNA, slot: data)
expression_matrix <- T_obj@assays$RNA@data

# Subset the expression matrix to include only marker genes
expression_matrix_markers <- as.matrix(expression_matrix)
expression_matrix_markers <- expression_matrix_markers[marker_genes,]

# Extract cell type information
cell_types <- T_obj@meta.data$celltype

# Create a matrix to store the bulk-level expression
bulk_expression_matrix <- matrix(0, nrow = length(unique(cell_types)), ncol = length(marker_genes))
rownames(bulk_expression_matrix) <- unique(cell_types)
colnames(bulk_expression_matrix) <- marker_genes
# rownames(bulk_expression_matrix)[5] <- "Empty"

# Aggregate expression data by cell type
for (cell_type in unique(cell_types)) {
  if(cell_type == ""){
    cells_in_type <- which(cell_types == cell_type)
    bulk_expression_matrix["Empty", ] <- rowMeans(expression_matrix_markers[, cells_in_type])
  }
  else{
    cells_in_type <- which(cell_types == cell_type)
    bulk_expression_matrix[cell_type, ] <- rowMeans(expression_matrix_markers[, cells_in_type])
  }
}

# Convert to data frame for easier saving
bulk_expression_df <- as.data.frame(bulk_expression_matrix)
bulk_expression_df$celltype <- rownames(bulk_expression_matrix)

# Reorder columns to have cell types as the first column
bulk_expression_df <- bulk_expression_df %>%
  select(celltype, everything())

write.csv(bulk_expression_df, "subtype_T_marker_expression_matrix.csv", row.names = TRUE)
# Count the occurrences of each gene
gene_counts <- table(T_markers_subtype$gene)
# Filter genes that occur more than once
filtered_genes <- names(gene_counts[gene_counts == 1])

bulk_expression_filter_df = bulk_expression_df[,filtered_genes]

write.csv(bulk_expression_filter_df, "subtype_T_marker_expression_filtered_matrix.csv", row.names = TRUE)