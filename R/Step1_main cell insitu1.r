
### in situ type analysis

###### Supervised cell typing using Insitutype ######
library(InSituType)
library(data.table)

# parse input parameters from linux scripts
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Please supply the results_folder, path of scRNA-seq reference csv file generated before, and giotto_object_path as command line arguments!", call.=FALSE)
} 

results_folder = args[1]
scRNA_path = args[2]
giotto_object_path = args[3]
    
# set working directory
# results_folder = "/home/zhanh/Janssen_AD_xenium/results_jul29"
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


xenium_gobj = calculateOverlapRaster(xenium_gobj, feat_info = 'neg_probe')
### can check through get_polygon_info(fov_join, polygon_name="cell",polygon_overlap="rna")
# Convert the overlap information into a cell by feature expression matrix which
# is then stored in the Giotto object's expression slot
#xenium_gobj = overlapToMatrix(xenium_gobj, feat_info = 'rna')
xenium_gobj = overlapToMatrix(xenium_gobj, feat_info = 'neg_probe')
# standard method of normalization (log normalization based)
xenium_gobj <- normalizeGiotto(gobject = xenium_gobj,
                                feat_type = 'neg_probe',
                                norm_methods = 'standard',
                                library_size_norm = FALSE,
                                verbose = TRUE)


###########################
### sanity check !!!! !!!!
showGiottoExpression(xenium_gobj)
showGiottoSpatLocs(xenium_gobj)
Giotto:::list_expression(xenium_gobj)

temp.cell_info = pDataDT(xenium_gobj)
prb_ind <- rownames(xenium_gobj@expression$cell$neg_probe$raw@exprMat)


# input 1: A spatial matrix of counts data, cells x genes
# input 2: A vector giving each cell’s mean negative control value
temp.exp_mat <- list()
temp.negmean <- list()

### giotto normalized matrix
temp.exp_mat[["gio_norm"]] <- as.matrix(xenium_gobj@expression$cell$rna$normalized@exprMat) ### gene x cell matrix  
temp.negmean[["gio_norm"]]  <- colMeans(as.matrix(xenium_gobj@expression$cell$neg_probe$normalized@exprMat))### 
#temp.negmean[["gio_norm"]]  <- weighted_sum

### giotto normalized matrix and batch corrected
a <- adjustGiottoMatrix(xenium_gobj, expression_values="normalized", feat_type="rna", batch_columns="list_ID")
a <- a@expression$cell$rna$custom@exprMat
temp.exp_mat[["gio_norm_batch"]] <- a

a <- adjustGiottoMatrix(xenium_gobj, expression_values="normalized", feat_type="neg_probe", batch_columns="list_ID")
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
##########################################################
### temp.ref_mat[["ref.scRNA"]] <- as.matrix(read.table("~/db/Normal.Control.Skin.scRNA/control.skin_G.by.CellType.txt",row.names=1,header=T,sep="\t",quote=NULL,comment.char=""))
# temp.ref_mat[["ref.scRNA"]] <- as.matrix(read.table("/home/zhanh/Xenium_psoriasis/test3/pso.skin_G.by.MainCellType.2.txt",row.names=1,header=T,sep="\t",quote=NULL,comment.char=""))
temp.ref_mat[["ref.scRNA"]] <- as.matrix(read.csv(scRNA_path, row.names=1, header=TRUE, sep=",", quote=NULL, comment.char=""))
row.names(temp.ref_mat[[1]]) <- gsub("\"", "", row.names(temp.ref_mat[[1]]))
colnames(temp.ref_mat[[1]]) <- gsub("\"", "", colnames(temp.ref_mat[[1]]))

row.names(temp.ref_mat[[1]]) <- gsub("^X\\.", "", row.names(temp.ref_mat[[1]]))
colnames(temp.ref_mat[[1]]) <- gsub("^X\\.", "", colnames(temp.ref_mat[[1]]))
row.names(temp.ref_mat[[1]]) <- gsub("\\.$", "", row.names(temp.ref_mat[[1]]))
colnames(temp.ref_mat[[1]]) <- gsub("\\.$", "", colnames(temp.ref_mat[[1]]))
temp.ref_mat[[1]] <- temp.ref_mat[[1]][, colnames(temp.ref_mat[[1]]) != "V15"]
temp.ref_mat[[1]] <- temp.ref_mat[[1]][, colnames(temp.ref_mat[[1]]) != "V1"]
### df <- read.table("/home/alextsoi/db/psoriasis.scRNA/pso.skin_G.by.CellType.txt",row.names=1,header=T,sep="\t",quote=NULL,comment.char="")


# input 4: 
## Incorporating additional data types: use immunofluorescence and spatial coordinates data for cohorting
## immunofluorescence

spatial_location <- get_spatial_locations(xenium_gobj,
    spat_unit = NULL,
    spat_loc_name = NULL,
    output = c("spatLocsObj", "data.table"),
    copy_obj = TRUE,
    verbose = TRUE,
    set_defaults = TRUE
    ) 

xy_coordinates <- spatial_location@coordinates 


### if using umap_harmony to store harmony corrected umap
### temp.umap <- xenium_gobj@dimension_reduction$cells$cell$rna$umap$umap_harmony@coordinates
temp.umap <- xenium_gobj@dimension_reduction$cells$cell$rna$umap$umap@coordinates
temp.meta_info <- cbind(xy_coordinates[,-1], temp.umap)
rownames(temp.meta_info) <- xy_coordinates$cell_ID

temp.cohort <- list()
temp.cohort.value <- fastCohorting(temp.meta_info, gaussian_transform = T) 


### !!!or directly using the clusters from UMAP
### temp.cohort.value <- xenium_gobj@cell_metadata$cell$rna@metaDT$leiden_clus

##########################################################
### do the in situtype supervised clustering one fov at a time

fs <- as.data.frame(unique(temp.cell_info[,3]))[,1]
is.ML <- list()

for (i in names(temp.exp_mat)[1:2]){
    print(i)
    is.ML[[i]] <- list()
    for (f in fs){
        print(f)
        templogical <- grep(colnames(temp.exp_mat[[i]]),pattern=f)
        
        is.ML[[i]][[f]] <- insitutypeML (x=t(temp.exp_mat[[i]][,templogical]), neg=temp.negmean[[i]][templogical], cohort=temp.cohort.value[templogical],reference_profiles=temp.ref_mat[["ref.scRNA"]])
    }
}

### add the in situ cluster label to giotto
tempcell <- character(length(xenium_gobj@cell_ID$cell))
names(tempcell) <- xenium_gobj@cell_ID$cell


### looks like the giotto normalized matrix provides reasonable results
i=1
for (f in names(is.ML[[i]])){
    a <- is.ML[[i]][[f]]$clust
    tempcell[names(a)] <- a
}

xenium_gobj = annotateGiotto(gobject = xenium_gobj,
                                annotation_vector = tempcell,
                                cluster_column = 'cell_ID',
                                name = 'cell_type_isML')
cell_types <- pDataDT(xenium_gobj)$cell_type_isML

saveGiotto(xenium_gobj,"gobj",overwrite=T)
# xenium_gobj = loadGiotto("/home/zhanh/Janssen_AD_xenium/results/gobj", reconnect_giottoImage=F)
### !!!!!!!!!!!!!!!!!!!!!
### output the pData including cell type label in case we lose the giotto object
write.table(pDataDT(xenium_gobj), file="xenium_gobj.pData",row.names=T,col.names=T,sep="\t",quote=F)

#########################################################
plotUMAP(gobject = xenium_gobj,
            cell_color = 'cell_type_isML',
            cell_color_code = col_vector,
            show_NN_network = F, 
            show_center_label=F,
            nn_network_to_use='sNN', 
            network_name="sNN.umap",
            point_size = 1,
            dim_reduction_name="umap",
            save_param = list(save_name = '7.1_UMAP_cell.types_isML_ori',base_width=15,base_height=12)
)

    
spatInSituPlotPoints(xenium_gobj,
                        show_image = FALSE,
                        feats = NULL,
                        point_size = 0.05,
                        show_polygon = TRUE,
                        polygon_feat_type = 'cell',
                        polygon_alpha = 1,
                        polygon_color = 'black',
                        polygon_line_size = 0.01,
                        polygon_fill = 'cell_type_isML',
                        polygon_fill_as_factor = TRUE,
                        polygon_fill_code = col_vector,
                        coord_fix_ratio = TRUE,
                        save_para = list(dpi = 1000,
                                        save_name = '7.1_insitu_cell.types_isML'))


    

