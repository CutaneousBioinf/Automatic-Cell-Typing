library(dplyr)
library(tidyr)
library(data.table)

#### function of majority vote

majority_vote <- function(gobj_path) {
  # Load Giotto object
  xenium_gobj <- loadGiotto(gobj_path, reconnect_giottoImage = FALSE)
  
  # Extract cell metadata
  cell_metadata <- pDataDT(xenium_gobj)
  cell_types <- cell_metadata$cell_type_isML
  clusters <- cell_metadata$leiden_clus
  
  # Clusters by cell type
  cluster_table <- as.data.table(cell_metadata[, .N, by = .(leiden_clus, cell_type_isML)])
  
  # Clusters by cell type percentage
  cluster_table[, percentage := N / sum(N) * 1, by = leiden_clus]
  percentage_table <- dcast(cluster_table, leiden_clus ~ cell_type_isML, value.var = "percentage", fill = 0)
  results_folder <- dirname(gobj_path)
  write.csv(percentage_table, file.path(results_folder, "8.1 cluster_by_cell_type_percentage.csv"), row.names = FALSE)
  
  # Overall percentage from spatial insitutype
  cell_type_totals <- colSums(percentage_table[,-1])
  total_sum <- sum(cell_type_totals)
  cell_type_percentages <- cell_type_totals / total_sum * 1
  overall_percentage_isML <- data.frame(
    cell_type = names(cell_type_percentages),
    overall_percentage_isML = cell_type_percentages
  )
  write.csv(overall_percentage_isML, file.path(results_folder, "8.2 overall_percentage_isML.csv"), row.names = FALSE)
  
  ##### Function to assign weighted main cell types and flag clusters needing subclustering #####
  cluster_data <- percentage_table
  overall_data <- overall_percentage_isML
  
  # Convert to Long Format and Merge
  cluster_data_long <- pivot_longer(cluster_data, cols = -leiden_clus, names_to = "cell_type", values_to = "percentage") %>%
    left_join(overall_data, by = "cell_type")
  
  # Calculate Weighted Proportion
  cluster_data_long <- cluster_data_long %>%
    mutate(weighted_proportion = percentage / overall_percentage_isML)
  
  # Normalize the Weighted Proportion
  cluster_data_long <- cluster_data_long %>%
    group_by(leiden_clus) %>%
    mutate(total_weighted = sum(weighted_proportion),
           normalized_weighted_proportion = weighted_proportion / total_weighted) %>%
    ungroup()
  
  # Calculate Cumulative Percentage
  cluster_data_long <- cluster_data_long %>%
    arrange(leiden_clus, desc(normalized_weighted_proportion)) %>%
    group_by(leiden_clus) %>%
    mutate(cumulative_percentage = cumsum(normalized_weighted_proportion)) %>%
    ungroup()
  write.csv(cluster_data_long, file.path(results_folder, "8.3 cluster_data_long.csv"), row.names = FALSE)
  
  # Assign main cell types and flag clusters needing sub-clustering
  assigned_clusters <- cluster_data_long %>%
    group_by(leiden_clus) %>%
    summarize(Assigned_Cell_Type = {
      top_cell_types <- if (cumulative_percentage[1] >= 0.95) {
        cell_type[1]
      } else {
        cell_type[1:which(cumulative_percentage >= 0.95)[1]]
      }
      # Apply the additional filter for weighted_proportion >= 1
      top_cell_types <- top_cell_types[weighted_proportion[cell_type %in% top_cell_types] >= 1]
      paste(top_cell_types, collapse = ", ")
    }) %>%
    ungroup() %>%
    mutate(Assigned_Cell_Type = sapply(Assigned_Cell_Type, function(x) paste(x, collapse = ", ")),
           Subclustering_Status = ifelse(grepl(",", Assigned_Cell_Type), "Need Subclustering", Assigned_Cell_Type))
  write.csv(assigned_clusters, file.path(results_folder, "8.4 assigned_clusters.csv"), row.names = FALSE)
  
  # Clusters need to sub
  needs_subclustering_clusters <- assigned_clusters %>%
    filter(Subclustering_Status == "Need Subclustering") %>%
    pull(leiden_clus)
  print(needs_subclustering_clusters)
  
  # Function to update the original Giotto object with new cell type assignments
  update_giotto_cell_types <- function(gobject, merged_results, spat_unit = "cell") {
    # Get original cell metadata
    cell_metadata <- pDataDT(gobject)
    
    # Remove any duplicate 'cell_ID' column if it exists
    if ("cell_ID" %in% colnames(cell_metadata)) {
      cell_metadata <- cell_metadata[, !("cell_ID"), with = FALSE]
    }
    
    # Update cell types based on the merged results
    for (i in 1:nrow(merged_results)) {
      cluster_id <- merged_results$leiden_clus[i]
      cell_type <- merged_results$Assigned_Cell_Type[i]
      cell_metadata[cell_metadata$leiden_clus == cluster_id, cell_type_majorvote := cell_type]
    }
    
    # Update cell type annotations in the Giotto object
    gobject <- addCellMetadata(gobject, new_metadata = cell_metadata, spat_unit = spat_unit)
    
    return(gobject)
  }
  
  # Update the original Giotto object with the new cell type assignments
  xenium_gobj_updated <- update_giotto_cell_types(xenium_gobj, assigned_clusters, spat_unit = "cell")
  
  # save the updated cell metadata
  saveGiotto(xenium_gobj_updated, file.path(results_folder, "gobj_reassign"), overwrite = TRUE)
}

# example
majority_vote("/home/zhanh/Janssen_AD_xenium/results_jul29/gobj")


# Plot UMAP using the updated Giotto object
plotUMAP(
  gobject = xenium_gobj_updated,
  cell_color = 'cell_type_majorvote',
  cell_color_code = col_vector,
  show_NN_network = FALSE,
  show_center_label = FALSE,
  nn_network_to_use = 'sNN',
  network_name = "sNN.umap",
  save_param = list(save_name = '1.4_UMAP_cell.types_majorvote', base_width = 15, base_height = 12)
)

# Plot UMAP using the cluster number
plotUMAP(
  gobject = xenium_gobj,
  cell_color = 'leiden_clus',
  cell_color_code = col_vector,
  show_NN_network = FALSE,
  show_center_label = FALSE,
  nn_network_to_use = 'sNN',
  network_name = "sNN.umap",
  save_param = list(save_name = '1.5_UMAP_clusters_isML_org', base_width = 15, base_height = 12)
)

