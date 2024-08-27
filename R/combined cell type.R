# This script is used to consolidate the prediction of cell type annotations derived from two different methods that we developed.

results_folder = "/home/zhanh/Janssen_AD_xenium/results_aug12"
setwd(results_folder)

# Install and load additional necessary packages
library(Giotto)
library(ggplot2)
library(ggthemes)
library(ggrastr)
library(RColorBrewer)
library(patchwork)
library(ggpubr)
install.packages("Polychrome")
library(Polychrome)
install.packages("symphony")
library(symphony)
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
library(Seurat)

combine_celltype_data <- function(symphony_path, giotto_paths) {
  # Load Symphony data
  data_symphony <- readRDS(symphony_path)
  data_symphony$cell_ID <- rownames(data_symphony)
  data_symph <- data_symphony[, c("cell_ID", "celltype.pred.combined")]
  rownames(data_symph) <- NULL
  data_symph$cell_ID <- sub("_", "-", data_symph$cell_ID)
  
  # Load Giotto objects
  giotto_objects <- lapply(giotto_paths, loadGiotto)
  
  # Extract and format metadata
  data_list <- lapply(giotto_objects, function(gobj) {
    metaDT <- gobj@cell_metadata$cell$rna@metaDT
    if ("cell_type_isML_subcelltyping" %in% colnames(metaDT)) {
      data_insitu <- metaDT[, c("cell_ID", "cell_type_isML_subcelltyping")]
    } else {
      data_insitu <- metaDT[, c("cell_ID", "sub_cell_type_isML_updated")]
    }
    data_insitu$cell_ID <- sub("^obj", "", data_insitu$cell_ID)
    return(data_insitu)
  })
  
  # Combine all data into a single data frame
  df <- c(data_list, list(data_symph))
  combined_results <- Reduce(function(x, y) merge(x, y, by = "cell_ID"), df)
  
  # Rename columns
  colnames(combined_results) <- c("cell_ID", "celltype_3steps_3rd_bass", "celltype_3steps", "celltype_1step", "celltype_symphony")
  
  # Save to CSV
  write.csv(combined_results, "celltype_combined1.csv", row.names = FALSE)
  
  # Return the combined results
  return(combined_results)
}

# Example usage of the function
symphony_path <- "/home/zhanh/Janssen_AD_xenium/results_aug12/symphony_celltype_results.rds"
giotto_paths <- c(
  "/home/zhanh/Janssen_AD_xenium/results_aug12/gobj_updated_final",
  "/home/zhanh/Janssen_AD_xenium/results_aug12/gobj_3_step_change_name",
  "/home/zhanh/Janssen_AD_xenium/results_aug12/gobj_1_step_change_name"
)

# Call the function
combined_results <- combine_celltype_data(symphony_path, giotto_paths)
print(head(combined_results))
