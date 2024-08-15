library(BASS)

#xenium_gobj=loadGiotto("/home/alextsoi/Researches/Xenium_skin_UMcustom.panel/analysis/gobj")
#obj_ids = paste0("obj",1:length(xenium_folders))
# psoriasis
xenium_gobj = loadGiotto("/home/zhaixt/Xenium_Janssen_0715/gobj_updated_0719")

# set working directory
#results_folder = '/path/to/save/directory/'
#results_folder = getwd()
results_folder = "/home/zhaixt/Xenium_Janssen_0715/"
setwd(results_folder)

instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  return_plot = FALSE)

xenium_gobj = replaceGiottoInstructions(xenium_gobj, instructions = instrs)

c1="Myeloid.Cells"
# c2="T.Cells"
# temp <- subsetGiotto(xenium_gobj, cell_ids=xenium_gobj@cell_ID$cell[(pDataDT(xenium_gobj)$cell_type_isML_updated==c1) | (pDataDT(xenium_gobj)$cell_type_isML_updated==c2)])
temp <- subsetGiotto(xenium_gobj, cell_ids=xenium_gobj@cell_ID$cell[(pDataDT(xenium_gobj)$cell_type_isML_updated==c1)])


# get cell info
cell_info = pDataDT(temp)[,c(1,3,6,7,10,11)]
cell_info = as.data.frame(cell_info)
head(cell_info)

# position
spatial_location = get_spatial_locations(temp,spat_unit = NULL,spat_loc_name = NULL,output = c("spatLocsObj", "data.table"),copy_obj = TRUE,verbose = TRUE,set_defaults = TRUE) 
pos_file <- spatial_location@coordinates 
pos_file = as.data.frame(pos_file)
rownames(pos_file) <- cell_info$cell_ID
pos_file <- pos_file[,c(2,3)]
colnames(pos_file) <- c("x","y")

# raw count matrix
raw_cnts = temp@expression$cell$rna$raw@exprMat + 1e-8
raw_cnts = as.matrix(raw_cnts)

sum(rownames(pos_file) != colnames(raw_cnts))

# create lists of counts and xy coordinates
fov_names = unique(cell_info$list_ID)
cnts = list()
xy = list()
clust = list()
clust_ref = 1:length(unique(cell_info$leiden_clus)) - 1
names(clust_ref) = unique(cell_info$leiden_clus)
for (i in 1:length(fov_names)){
  cell_ids = cell_info[cell_info$list_ID == fov_names[i], "cell_ID"]
  cnts[[i]] = as.matrix(raw_cnts[,cell_ids])
  xy[[i]] = pos_file[cell_ids,]
  clust[[i]] = cell_info[cell_info$list_ID == fov_names[i], "leiden_clus"] - 1
  clust[[i]] = clust_ref[cell_info[cell_info$list_ID == fov_names[i], "leiden_clus"]]
}

names(cnts) = fov_names
names(xy) = fov_names
names(clust) = fov_names

# create a BASS object
BASS <- createBASSObject(cnts, xy, C = 7, R = 5, beta_method = "SW")
listAllHyper(BASS)
BASS <- BASS.preprocess(BASS)
# run BASS algorithm using self-provided clustering labels
#BASS <- BASS.run(BASS, clust)
BASS <- BASS.run(BASS)
# Post-process posterior samples
BASS <- BASS.postprocess(BASS)
clabels <- BASS@results$c # cell type clusters
zlabels <- BASS@results$z # spatial domain labels
pi_est <- BASS@results$pi # cell type composition matrix

giotto_c = c()
giotto_z = c()
for (i in 1:length(clabels)){
  giotto_c = c(giotto_c,clabels[[i]])
  giotto_z = c(giotto_z,zlabels[[i]])
}
names(giotto_c) = pDataDT(temp)$cell_ID
names(giotto_z) = pDataDT(temp)$cell_ID

temp = annotateGiotto(gobject = temp,
                          annotation_vector = giotto_c,
                          cluster_column = 'cell_ID',
                          name = 'bass_clus_M')

temp = annotateGiotto(gobject = temp,
                          annotation_vector = giotto_z,
                          cluster_column = 'cell_ID',
                          name = 'bass_domain_M')

plotUMAP(temp,
         cell_color = 'bass_clus_M',
         cell_color_code = col_vector,
         point_size = 1,
         save_param = list(save_name = 'UMAP_bass_clust_M'))

plotUMAP(temp,
         cell_color = 'bass_domain_M',
         cell_color_code = col_vector,
         point_size = 1,
         save_param = list(save_name = 'UMAP_bass_domain_M'))

spatInSituPlotPoints(xenium_gobj,
                     show_polygon = TRUE,
                     polygon_color = 'white',
                     polygon_line_size = 0.01,
                     polygon_fill = 'bass_clus',
                     polygon_fill_as_factor = TRUE,
                     polygon_fill_code = col_vector,
                     polygon_alpha = 1, # transparency of polygon colors
                     save_param = list(base_height = 10, dpi = 1000,
                                       save_name = 'inSitu_bass_cluster'))
spatInSituPlotPoints(xenium_gobj,
                     show_polygon = TRUE,
                     polygon_color = 'white',
                     polygon_line_size = 0.01,
                     polygon_fill = 'bass_domain',
                     polygon_fill_as_factor = TRUE,
                     polygon_fill_code = col_vector,
                     polygon_alpha = 1, # transparency of polygon colors
                     save_param = list(base_height = 10, dpi = 1000,
                                       save_name = 'inSitu_bass_domain'))


########## sub-celltyping

### add the in situ cluster label to giotto
temp.cell_info = pDataDT(xenium_gobj)
prb_ind <- rownames(xenium_gobj@expression$cell$neg_probe$raw@exprMat)

tempcell = xenium_gobj@cell_metadata$cell$rna$cell_type_isML_updated.sub
names(tempcell) = xenium_gobj@cell_ID$cell

temp.exp_mat <- list()
temp.negmean <- list()

### giotto normalized matrix
temp.exp_mat[["gio_norm"]] <- as.matrix(temp@expression$cell$rna$raw@exprMat) ### gene x cell matrix  
temp.negmean[["gio_norm"]]  <- colMeans(as.matrix(temp@expression$cell$neg_probe$raw@exprMat))### 
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
temp.sc <- as.matrix(read.table("/home/zhaixt/Myeloid Cells_subset.csv",row.names=1,header=T,sep=",",quote=NULL,comment.char="",check.names=F))
row.names(temp.sc) <- gsub("\"","",row.names(temp.sc) )
colnames(temp.sc) <- gsub("\"","",colnames(temp.sc) )
  
temp.ref_mat[["ref.scRNA"]] <- temp.sc
temp.cohort.value <- giotto_c

fs <- as.data.frame(unique(temp.cell_info[,3]))[,1]
temp.is.ML <- list()

for (i in names(temp.exp_mat)[1]){
  print(i)
  temp.is.ML[[i]] <- list()
  for (f in fs){
    print(f)
    templogical <- grep(colnames(temp.exp_mat[[i]]),pattern=f)
    temp.is.ML[[i]][[f]] <- insitutypeML(x=t(temp.exp_mat[[i]][,templogical]), neg=temp.negmean[[i]][templogical], cohort=temp.cohort.value[templogical],reference_profiles=temp.ref_mat[["ref.scRNA"]])
  }
}




### looks like the giotto normalized matrix provides reasonable results
i=1
for (f in names(temp.is.ML[[i]])){
  a <- temp.is.ML[[i]][[f]]$clust
  tempcell[names(a)] <- a
}

sup = temp.is.ML[[i]]
heatmap(sweep(sup$profiles, 1, pmax(apply(sup$profiles, 1, max), .2), "/"), scale = "none",
        main = "Mean cell type expression profiles")

sup = temp.is.ML[[i]]
cols <-
  c(
    '#8DD3C7',
    '#BEBADA',
    '#FB8072',
    '#80B1D3',
    '#FDB462',
    '#B3DE69',
    '#FCCDE5',
    '#D9D9D9',
    '#BC80BD',
    '#CCEBC5',
    '#FFED6F',
    '#E41A1C',
    '#377EB8',
    '#4DAF4A',
    '#984EA3',
    '#FF7F00',
    '#FFFF33',
    '#A65628',
    '#F781BF',
    '#999999'
  )

cols <- cols[seq_along(unique(sup$clust))]
names(cols) <- unique(sup$clust)

# make the flightpath plot
png("output.png", width = 500, height = 700, res = 300)
fp <- flightpath_plot(flightpath_result = NULL, insitutype_result = sup, col = cols[sup$clust])
dev.off()

print(fp)




#xenium_gobj@cell_metadata$cell$rna$sub_cell_type_isML_updated = xenium_gobj@cell_metadata$cell$rna$cell_type_isML_updated.sub

xenium_gobj = annotateGiotto(gobject = xenium_gobj,
                                       annotation_vector = tempcell,
                                       cluster_column = 'cell_ID',
                                       name = 'sub_cell_type_isML_updated')

plotUMAP(gobject = xenium_gobj,
         cell_color = 'sub_cell_type_isML_updated',
         cell_color_code = col_vector,
         show_NN_network = F, 
         show_center_label=F,
         nn_network_to_use='sNN', 
         network_name="sNN.umap",
         point_size = 1,
         dim_reduction_name="umap",
         save_param = list(save_name = '8.1_UMAP_M_T.sub.cell.types.updated_isML',base_width=15,base_height=12)
)

saveGiotto(xenium_gobj, save_dir = results_folder, save_name = "gobj_updated_0723")
