# This script uses the Seurat package to identify marker genes for each cell type in the spatial transcriptomics dataset.


# parse input parameters from linux scripts
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=5) {
  stop("Please supply the results_folder, path of scRNA-seq RDS file, giotto_object_path, main_celltype name in scRNA-seq RDS file, and sub_celltype name in scRNA-seq RDS file as command line arguments!", call.=FALSE)
} 

results_folder = args[1]
scRNA_path = args[2]
giotto_object_path = args[3]
main_celltype = args[4]
sub_celltype = args[5]

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






generate_marker_gene_matrix <- function(scrna, target_main_celltype){
    T_obj <- subset(scrna_subset, subset = eval(parse(text = paste(main_celltype, "==", target_main_celltype))))

    T_obj <- subset(T_obj, cells = which(T_obj@meta.data[[sub_celltype]] != ""))
    T_obj@meta.data[[sub_celltype]] <- gsub(" ", "_", T_obj@meta.data[[sub_celltype]])
    T_obj@meta.data[[sub_celltype]] <- gsub("-", "_", T_obj@meta.data[[sub_celltype]])
    Idents(T_obj) <- T_obj@meta.data[[sub_celltype]]
    T_markers_subtype <- FindAllMarkers(T_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, group.by = sub_celltype)

    # Extract the unique marker genes
    marker_genes <- unique(T_markers_subtype$gene)

    # Get the expression data (assay: RNA, slot: data)
    expression_matrix <- T_obj@assays$RNA@data

    # Subset the expression matrix to include only marker genes
    expression_matrix_markers <- as.matrix(expression_matrix)
    expression_matrix_markers <- expression_matrix_markers[marker_genes,]

    # Extract cell type information
    cell_types <- T_obj@meta.data[[sub_celltype]]

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

    exp_file_name1 = paste0("subtype_", target_main_celltype, "_marker_expression_matrix.csv")
    write.csv(bulk_expression_df, exp_file_name1, row.names = TRUE)
    # Count the occurrences of each gene
    gene_counts <- table(T_markers_subtype$gene)
    # Filter genes that occur more than once
    filtered_genes <- names(gene_counts[gene_counts == 1])


    bulk_expression_filter_df = bulk_expression_df[,filtered_genes]

    exp_file_name2 = paste0("subtype_", target_main_celltype, "_marker_expression_filtered_matrix.csv")
    write.csv(bulk_expression_filter_df, exp_file_name2, row.names = TRUE)
}



library(dplyr)
# main_celltype = "celltype.global"
# sub_celltype = "celltype"
# unique(spatial@meta.data[[sub_celltype]])

vec_main_celltype = spatial@meta.data[[main_celltype]]
vec_sub_celltype = spatial@meta.data[[sub_celltype]]

data = data.frame(main_celltype = vec_main_celltype, sub_celltype = vec_sub_celltype)
data[data == ""] <- NA
data = na.omit(data)
data = data[!duplicated(data),]

result <- data %>%
  group_by(main_celltype) %>%
  summarise(sub_celltypes = paste(unique(sub_celltype), collapse = ", "))



for (i in 1:nrow(result)){
  vec_sub_celltype = unlist(strsplit(result$sub_celltypes[i], ", "))
  if (length(vec_sub_celltype) > 1){
    #print(i)
    print(result$main_celltype[i])
    #print(vec_sub_celltype)
    generate_marker_gene_matrix(spatial, result$main_celltype[i])
  }
}
