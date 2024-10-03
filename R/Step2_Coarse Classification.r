# parse input parameters from linux scripts
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Please supply the results_folder, path of scRNA-seq reference csv file generated before, and giotto_object_path as command line arguments!", call.=FALSE)
} 

results_folder = args[1]
scRNA_path = args[2]
giotto_object_path = args[3]

setwd(results_folder)

# load giotto object
xenium_gobj = loadGiotto(giotto_object_path, reconnect_giottoImage=F) 
instrs <- createGiottoInstructions(
        save_dir = results_folder,
        save_plot = TRUE,
        show_plot = FALSE,
        return_plot = FALSE,
        python_path = python_path)

xenium_gobj <- replaceGiottoInstructions(xenium_gobj, instructions = instrs)

    

sub_celltype_insitutypeML <- function(cluster, sub_celltype){
  #########subclustering for mixed main celltypes
  c = cluster
  temp <- subsetGiotto(xenium_gobj_updated, cell_ids=xenium_gobj_updated@cell_ID$cell[pDataDT(xenium_gobj_updated)$cell_type_majorvote==c])
  
  temp = runPCA(gobject = temp,spat_unit = 'cell',expression_values = 'normalized',feats_to_use = NULL,scale_unit = F,center = F)
  
  # run Harmony
  temp = runGiottoHarmony(temp, vars_use = "fov",dim_reduction_to_use = "pca",dimensions_to_use = 1:25)
  
  # Generate UMAP from Harmony
  temp = runUMAP(temp,dim_reduction_to_use = "harmony",dim_reduction_name = "harmony",dimensions_to_use = 1:25,n_threads = 4)
  
  
  # sNN and Leiden clustering
  temp = createNearestNetwork(temp,
                              type = "sNN",
                              dim_reduction_to_use = "umap",
                              dim_reduction_name = "umap",
                              dimensions_to_use = 1:2, name="sNN.umap")
  
  # the resolution here defaultly to be 0.1, but can be customized by the user
  # the higher the resolution is, the more clusters there will be
  temp = doLeidenCluster(temp,
                         network_name = "sNN.umap",
                         resolution = 0.05,
                         n_iterations = 1000)
  
  
  ### in situ type for T-cell subclustering
  temp.cell_info = pDataDT(xenium_gobj_updated)
  prb_ind <- rownames(xenium_gobj_updated@expression$cell$neg_probe$raw@exprMat)
  
  # input 1: A spatial matrix of counts data, cells x genes
  # input 2: A vector giving each cell’s mean negative control value
  temp.exp_mat <- list()
  temp.negmean <- list()
  
  ### giotto normalized matrix
  temp.exp_mat[["gio_norm"]] <- as.matrix(temp@expression$cell$rna$normalized@exprMat) + 1e-7 ### gene x cell matrix  
  temp.negmean[["gio_norm"]]  <- colMeans(as.matrix(temp@expression$cell$neg_probe$normalized@exprMat))### 
  #temp.negmean[["gio_norm"]]  <- weighted_sum
  
  ### giotto normalized matrix and batch corrected
  a <- adjustGiottoMatrix(temp, expression_values="normalized", feat_type="rna", batch_columns="list_ID")
  a <- a@expression$cell$rna$custom@exprMat
  temp.exp_mat[["gio_norm_batch"]] <- a
  
  a <- adjustGiottoMatrix(temp, expression_values="normalized", feat_type="neg_probe", batch_columns="list_ID")
  a <- a@expression$cell$neg_probe$custom@exprMat
  temp.negmean[["gio_norm_batch"]] <- Matrix::colMeans(a)
  
  ### remove cells with zero counts for all genes
  for (i in names(temp.exp_mat)){
    temp.exp_mat[[i]] <- temp.exp_mat[[i]][,colSums(temp.exp_mat[[i]])!=0]
    mode(temp.exp_mat[[i]]) <- "integer"
  }
  
  for (i in names(temp.exp_mat)){
    temp.negmean[[i]] <- temp.negmean[[i]][colnames(temp.exp_mat[[i]])]
  }
  
  # input 3: A "reference matrix" giving the expected expression profile of each cell type,
  # with genes in rows and cell types in columns. The reference matrix must be in linear-scale, not log-scale.
  temp.ref_mat <- list()
  temp.sc <- as.matrix(read.table(scRNA_path,row.names=1,header=T,sep=",",quote=NULL,comment.char="",check.names=F))
  row.names(temp.sc)<- gsub("\"","",row.names(temp.sc) )
  colnames(temp.sc)<- gsub("\"","",colnames(temp.sc) )
  temp.sc = t(temp.sc)

  #### !!!! need to select the sub-celltypes
  positions <- match(sub_celltype, colnames(temp.sc))
  if (length(positions) == 0){
    stop("sub_celltype not found in the reference matrix")
  }
  temp.ref_mat[["ref.scRNA"]] <- temp.sc[, positions]
  temp.ref_mat[["ref.scRNA"]] <- temp.sc
  
  
  # input 4: 
  ## Incorporating additional data types: use immunofluorescence and spatial coordinates data for cohorting
  ## immunofluorescence
  
  spatial_location <- get_spatial_locations(temp,spat_unit = NULL,spat_loc_name = NULL,output = c("spatLocsObj", "data.table"),copy_obj = TRUE,verbose = TRUE,set_defaults = TRUE) 
  
  xy_coordinates <- spatial_location@coordinates 
  
  temp.umap <- temp@dimension_reduction$cells$cell$rna$umap$umap@coordinates
  temp.meta_info <- cbind(xy_coordinates[,-1], temp.umap)
  rownames(temp.meta_info) <- xy_coordinates$cell_ID
  
  temp.cohort <- list()
  temp.cohort.value <- fastCohorting(temp.meta_info, gaussian_transform = T) 
  
  ### if directly using the expression cluster
  ### temp.cohort.value <- pDataDT(temp)$leiden_clus
  
  ##########################################################
  ### do the in situtype supervised clustering one fov at a time
  
  fs <- as.data.frame(unique(temp.cell_info[,3]))[,1]
  temp.is.ML <- list()
  
  for (i in names(temp.exp_mat)[1]){
    print(i)
    temp.is.ML[[i]] <- list()
    for (f in fs){
      print(f)
      templogical <- grep(colnames(temp.exp_mat[[i]]),pattern=f)
      
      temp.is.ML[[i]][[f]] <- insitutypeML (x=t(temp.exp_mat[[i]][,templogical]), neg=temp.negmean[[i]][templogical], cohort=temp.cohort.value[templogical],reference_profiles=temp.ref_mat[["ref.scRNA"]])
    }
  }
  
  ########################
  # project to original space
  
  tempcell = xenium_gobj_updated@cell_metadata$cell$rna$cell_type_isML_updated
  names(tempcell) = xenium_gobj_updated@cell_ID$cell
  
  i=1
  for (f in names(temp.is.ML[[i]])){
    a <- temp.is.ML[[i]][[f]]$clust
    tempcell[names(a)] <- a
  }
  
  xenium_gobj_updated = annotateGiotto(gobject = xenium_gobj_updated,
                                       annotation_vector = tempcell,
                                       cluster_column = 'cell_ID',
                                       name = 'cell_type_isML_updated')
  
  return(xenium_gobj_updated)
  
}




xenium_gobj_updated = xenium_gobj

main_celltype_list = unique(xenium_gobj_updated@cell_metadata$cell$rna$cell_type_majorvote) 




for (main_celltype in main_celltype_list) {
  contains_comma = grepl(",", main_celltype)
  if (contains_comma){
    sub_celltype_list = as.vector(strsplit(main_celltype, ", ")[[1]])
    print(sub_celltype_list)
    xenium_gobj_updated = sub_celltype_insitutypeML(main_celltype, sub_celltype_list)
  }
}

saveGiotto(xenium_gobj_updated, save_dir = results_folder, save_name = "gobj_updated")


plotUMAP(gobject = xenium_gobj_updated,
         cell_color = 'cell_type_isML_updated',
         cell_color_code = col_vector,
         show_NN_network = F, 
         show_center_label=F,
         nn_network_to_use='sNN', 
         network_name="sNN.umap",
         point_size = 1,
         dim_reduction_name="umap",
         save_param = list(save_name = '8.1_UMAP_cell.types.updated_isML',base_width=15,base_height=12)
)
