library('Seurat')
library('symphony')
library('tibble')
library('dplyr')
library('irlba')


########## Set Parameters ##########
set.seed(0)

## Build reference
skip_build_ref = TRUE
ref_exp_path = '/home/yulicai/symphony/exprs_norm.rds'	# janssen ad
ref_metadata_path = '/home/yulicai/symphony/metadata.rds'
ref_exp = readRDS(ref_exp_path)	
ref_metadata = readRDS(ref_metadata_path)
maintype_col_name = 'celltype.global'  # column in metadata representing main celltype
subtype_col_name = 'celltype'   # column in metadata representing sub celltype

#seurat_obj = readRDS('/home/alextsoi/db/psoriasis.scRNA/seurat_with_rachael_subtype.RDS')
#ref_exp = seurat_obj$RNA@data    # gene x cell, should be log(CP10K + 1) normalized 
#ref_metadata = seurat_obj@meta.data
#maintype_col_name = 'celltype'  # column in metadata representing main celltype
#subtype_col_name = 'subtypes'   # column in metadata representing sub celltype

downsample = TRUE    # You may need downsampling when you have too many cells.
downsample_to = 100000  # The number of cells you want to keep. Final number may be slightly different. Ignore it when downsample=FALSE.
vars_use = c('orig.ident') # Harmony parameter. If meta_data is dataframe, this defined which variable(s) to remove (character vector).
save_main_ref_dir = '/home/yulicai/symphony/test_pipeline/reference/main'
save_main_uwot_dir = '/home/yulicai/symphony/test_pipeline/reference/main'
save_sub_ref_dir = '/home/yulicai/symphony/test_pipeline/reference/sub'
save_sub_uwot_dir = '/home/yulicai/symphony/test_pipeline/reference/sub'


## Query
output_dir <- '/home/yulicai/symphony/test_pipeline/outputs'
query_paths <- c('/hits/AGC/Spatial/Xenium/hSkin_480g/10265-JF_and_10628-JF/20240517__214754__10265_10628-JF/20240517__214754__10265_10628-JF/output-XETG00077__0022387__10628-JF-1_ROI_A1__20240517__215149/cell_feature_matrix.h5',
                 '/hits/AGC/Spatial/Xenium/hSkin_480g/10265-JF_and_10628-JF/20240517__214754__10265_10628-JF/20240517__214754__10265_10628-JF/output-XETG00077__0022387__10628-JF-1_ROI_A2__20240517__215149/cell_feature_matrix.h5',
                 '/hits/AGC/Spatial/Xenium/hSkin_480g/10265-JF_and_10628-JF/20240517__214754__10265_10628-JF/20240517__214754__10265_10628-JF/output-XETG00077__0022387__10628-JF-1_ROI_B__20240517__215149/cell_feature_matrix.h5',
                 '/hits/AGC/Spatial/Xenium/hSkin_480g/10265-JF_and_10628-JF/20240517__214754__10265_10628-JF/20240517__214754__10265_10628-JF/output-XETG00077__0022387__10628-JF-1_ROI_C1__20240517__215149/cell_feature_matrix.h5'
                )
k = 5   # the k of KNN in mapping.



########## Build the Reference for Main Celltypes ##########
if (skip_build_ref) {print('Skip building reference.')} else {
print('Start building reference for main celltypes...')
for (p in c(save_main_ref_dir, save_main_uwot_dir, save_sub_ref_dir, save_sub_uwot_dir)) {
    if (!dir.exists(p)){
        dir.create(p, recursive = TRUE)
    }
}
# Preprocess
print('Preprocessing...')
# remove missing celltypes and downsample
if (downsample) {
    ref_metadata <- ref_metadata %>% 
                    filter(get(subtype_col_name) != '') %>% 
                    rownames_to_column() %>% 
                    group_by(get(subtype_col_name)) %>% 
                    sample_frac(downsample_to/dim(ref_exp)[2]) %>% 
                    column_to_rownames()
} else {
    ref_metadata <- ref_metadata %>% 
                    filter(get(subtype_col_name) != '') %>% 
                    rownames_to_column() %>% 
                    column_to_rownames()
}
ref_exp <- ref_exp[,rownames(ref_metadata)]
# remove unexpressed genes, otherwise they may cause errors
ref_exp <- ref_exp[rowSums(ref_exp)>0,]

# Calculate and save the mean and standard deviations for each gene
print('Calculating mean and sd...')
genes_means_sds = tibble(symbol = rownames(ref_exp), mean = Matrix::rowMeans(ref_exp))
genes_means_sds$stddev = singlecellmethods::rowSDs(ref_exp, genes_means_sds$mean)

# Scale data using calculated gene means and standard deviations
print('Scaling data...')
ref_exp_scaled = singlecellmethods::scaleDataWithStats(ref_exp, genes_means_sds$mean, genes_means_sds$stddev, 1)

# Run SVD, save gene loadings (s$u)
print('Running SVD...')
s = irlba(ref_exp_scaled, nv = 20, fastpath=FALSE)
Z_pca_ref = diag(s$d) %*% t(s$v) # [pcs by cells]
loadings = s$u

# Run Harmony integration
ref_harmObj = harmony::HarmonyMatrix(
        data_mat = t(Z_pca_ref),  ## PCA embedding matrix of cells
        meta_data = ref_metadata, ## dataframe with cell labels
        theta = c(2),             ## cluster diversity enforcement
        vars_use = vars_use,    ## variable to integrate out
        nclust = 100,             ## number of clusters in Harmony model
        max.iter.harmony = 20,
        return_object = TRUE,     ## return the full Harmony model object
        do_pca = FALSE            ## don't recompute PCs
)

# Compress a Harmony object into a Symphony reference
reference = symphony::buildReferenceFromHarmonyObj(
                           ref_harmObj,            # output object from HarmonyMatrix()
                           ref_metadata,           # reference cell metadata
                           genes_means_sds,     # gene names, means, and std devs for scaling
                           loadings,               # genes x PCs matrix
                           verbose = TRUE,         # verbose output
                           do_umap = TRUE,         # Set to TRUE only when UMAP model was saved for reference
                           save_uwot_path = paste(save_main_uwot_dir,'/uwot_main', sep=''))

# specify which normalization method was used to build the reference 
reference$normalization_method = 'log(CP10k+1)'

# Save Symphony reference
saveRDS(reference, paste(save_main_ref_dir, '/ref_main.rds', sep=''))



########## Build the Reference for Subtypes ##########
print('Start building reference for sub celltypes...')
print('Counts of each main cell type:')
print(table(ref_metadata[,maintype_col_name]))

for (main_type in names(table(ref_metadata[,maintype_col_name]))){
	print(paste('Start building reference for', main_type, '...'))
	ref_metadata_sub <- ref_metadata %>% filter(get(maintype_col_name) == main_type) %>% filter(get(subtype_col_name) != '') %>% rownames_to_column() %>% column_to_rownames()
	if (length(table(ref_metadata_sub[,subtype_col_name])) <= 1){
		print('No subtype found. Skip.')
		next
	}
	ref_exp_sub <- ref_exp[,rownames(ref_metadata_sub)]
	ref_exp_sub <- ref_exp_sub[rowSums(ref_exp_sub)>0,]

	# Calculate and save the mean and standard deviations for each gene
	print('Calculating mean and sd...')
	genes_means_sds = tibble(symbol = rownames(ref_exp_sub), mean = Matrix::rowMeans(ref_exp_sub))
	genes_means_sds$stddev = singlecellmethods::rowSDs(ref_exp_sub, genes_means_sds$mean)

	# Scale data using calculated gene means and standard deviations
	print('Scaling data...')
	ref_exp_scaled = singlecellmethods::scaleDataWithStats(ref_exp_sub, genes_means_sds$mean, genes_means_sds$stddev, 1)
	print(str(ref_exp_scaled))

	# Run SVD, save gene loadings (s$u)
	print('Running SVD...')
	s = irlba(ref_exp_scaled, nv = 20, fastpath=FALSE)
	Z_pca_ref = diag(s$d) %*% t(s$v) # [pcs by cells]
	loadings = s$u

	# Run Harmony integration
	ref_harmObj = harmony::HarmonyMatrix(
	        data_mat = t(Z_pca_ref),  ## PCA embedding matrix of cells
	        meta_data = ref_metadata_sub, ## dataframe with cell labels
	        theta = c(2),             ## cluster diversity enforcement
	        vars_use = vars_use,    ## variable to integrate out
	        nclust = 100,             ## number of clusters in Harmony model
	        max.iter.harmony = 20,
	        return_object = TRUE,     ## return the full Harmony model object
	        do_pca = FALSE            ## don't recompute PCs
	)

	# Compress a Harmony object into a Symphony reference
	reference = symphony::buildReferenceFromHarmonyObj(
	                           ref_harmObj,            # output object from HarmonyMatrix()
	                           ref_metadata_sub,           # reference cell metadata
	                           genes_means_sds,     # gene names, means, and std devs for scaling
	                           loadings,               # genes x PCs matrix
	                           verbose = TRUE,         # verbose output
	                           do_umap = TRUE,         # Set to TRUE only when UMAP model was saved for reference
	                           save_uwot_path = paste(save_sub_uwot_dir,'/uwot_sub_',gsub(' ','-',main_type),sep='')
	                           )

	# specify which normalization method
	reference$normalization_method = 'log(CP10k+1)'

	# Save Symphony reference
	saveRDS(reference, paste(save_sub_ref_dir,'/ref_sub_',gsub(' ','-',main_type),'.rds',sep=''))
}
}


########## Set up Functions for plot ##########
library(ggplot2)
library(ggthemes)
library(ggrastr)
library(RColorBrewer)
library(patchwork)
library(ggpubr)
library(Polychrome)

plotBasic <- function(umap_labels,                # metadata, with UMAP labels in UMAP1 and UMAP2 slots
                      title = 'Query',         # Plot title
                      color.by = 'cell_type',  # metadata column name for coloring
                      facet.by = NULL,         # (optional) metadata column name for faceting
                      color.mapping = NULL,    # custom color mapping
                      legend.position = 'right', # Show cell type legend
					  save_path = './plot.png',
                      size = c(8,10)) {  
    
    p = umap_labels %>%
            dplyr::sample_frac(1L) %>% # permute rows randomly
            ggplot(aes(x = UMAP1, y = UMAP2)) + 
            geom_point_rast(aes(col = get(color.by)), size = 0.3, stroke = 0.2, shape = 16)
        if (!is.null(color.mapping)) { p = p + scale_color_manual(values = color.mapping) }
    
    # Default formatting
    p = p + theme_bw() +
            labs(title = title, color = color.by) + 
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(legend.position=legend.position) +
            theme(legend.text = element_text(size=8), legend.title=element_text(size=12)) + 
            guides(colour = guide_legend(override.aes = list(size = 4))) + guides(alpha = 'none')

    if(!is.null(facet.by)) {
        p = p + facet_wrap(~get(facet.by)) +
                theme(strip.text.x = element_text(size = 12)) }    
    
    ggsave(save_path, plot=p, width=size[2], height=size[1])
}




########## Predict Main Celltype ##########
for (p in c(output_dir)) {
    if (!dir.exists(p)){
        dir.create(p, recursive = TRUE)
    }
}
print('Start working on main celltypes.')
## load reference for main celltype
reference_main <- readRDS(paste(save_main_ref_dir, '/ref_main.rds', sep=''))
# If it cannot automatically find uwot model, manually set the path.
#reference_main$save_uwot_path <- ''

# plot distribution of celltypes in reference space
plotBasic(cbind(reference_main$meta_data, reference_main$umap$embedding), 
          title = 'Reference Cells in Reference UMAP Space of Main celltypes', 
          color.by = maintype_col_name,
          save_path = paste(output_dir, '/1-main_celltype_refer.png', sep=''),
          size = c(8,10))
#plotBasic(cbind(reference_main$meta_data, reference_main$umap$embedding), 
#          title = 'Reference Cells in Reference UMAP Space of Main celltypes', 
#          color.by = subtype_col_name,
#          save_path = paste(output_dir, '/1-02-sub_celltype_refer.png', sep=''),
#          size = c(8,10.5))

## load query
print('Loading queries...')
seurat_objs <- c()
for (path in query_paths){
    h5f <- Read10X_h5(path)
    seurat_obj <- CreateSeuratObject(h5f$`Gene Expression`)
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_objs <- c(seurat_objs, seurat_obj)
}
if (length(seurat_objs)>1) {
seurat_merged <- merge(seurat_objs[[1]], y=seurat_objs[2:length(seurat_objs)], add.cell.ids=as.character(c(1:length(seurat_objs))))
} else {
	seurat_merged <- seurat_objs[1]
}

print(str(seurat_merged))

# Map main celltype
print('Mapping main celltypes...')
query_exp = seurat_merged$RNA@data
query_metadata= seurat_merged@meta.data
query = mapQuery(query_exp,             # query gene expression (genes x cells)
                 query_metadata,        # query metadata (cells x attributes)
                 reference_main,             # Symphony reference object
                 do_normalize = FALSE,  # perform log(CP10k+1) normalization on query
                 do_umap = TRUE)        # project query cells into reference UMAP
query = knnPredict(query, 
                   reference_main, 
                   reference_main$meta_data[,maintype_col_name], 
                   k = k,
                   save_as = paste(maintype_col_name,'.pred',sep=''))
plotBasic(cbind(query$meta_data, query$umap), 
          title = 'Query Cells in Reference UMAP Space', 
          color.by = paste(maintype_col_name,'.pred',sep=''),
          save_path = paste(output_dir, '/2-main_celltype_pred.png', sep=''),
          size = c(8,10))


########## Predict Sub Celltype ##########
print('Start working on sub celltypes.')
celltypes_with_sub = list()
for (main_type in names(table(reference_main$meta_data[,maintype_col_name]))){
	ref_metadata_sub <- reference_main$meta_data %>% filter(get(maintype_col_name) == main_type) %>% filter(get(subtype_col_name) != '') %>% rownames_to_column() %>% column_to_rownames()
	if (length(table(ref_metadata_sub[,subtype_col_name])) <= 1){next}
    celltypes_with_sub[main_type] <- as.data.frame(names(table(ref_metadata_sub[,subtype_col_name])))
}
print("Main celltypes with sub celltypes:")
print(celltypes_with_sub)

# plot distribution of celltypes in reference space
for (main_type in names(celltypes_with_sub)){
    ## load reference for main celltype
    reference_sub <- readRDS(paste(save_sub_ref_dir,'/ref_sub_',gsub(' ','-',main_type),'.rds', sep=''))
    # If it cannot automatically find uwot model, manually set the path.
    #reference_sub$save_uwot_path <- ''

    # plot distribution of celltypes in reference space
    plotBasic(cbind(reference_sub$meta_data, reference_sub$umap$embedding), 
              title = paste('Reference Cells in Reference UMAP Space of', main_type), 
              color.by = subtype_col_name,
              save_path = paste(output_dir, '/3-sub_celltype_refer_',gsub(' ','-',main_type),'.png', sep=''),
              size = c(8,10))
}

# Map sub celltype
print('Mapping main celltypes...')
i = 1
for (main_type in names(celltypes_with_sub)){
    reference_sub <- readRDS(paste(save_sub_ref_dir,'/ref_sub_',gsub(' ','-',main_type),'.rds', sep=''))
	query_metadata_sub = query$meta_data[query$meta_data[paste(maintype_col_name,'.pred',sep='')] == main_type,]
	query_exp_sub = as.array(query$exp)[,rownames(query_metadata_sub)]
	query_sub = mapQuery(query_exp_sub,             # query gene expression (genes x cells)
	                     query_metadata_sub,        # query metadata (cells x attributes)
	                     reference_sub,             # Symphony reference object
	                     do_normalize = FALSE,  # perform log(CP10k+1) normalization on query
	                     do_umap = TRUE)        # project query cells into reference UMAP
	query_sub = knnPredict(query_sub, 
	                       reference_sub, 
	                       reference_sub$meta_data[,subtype_col_name], 
	                       k = k,
	                       save_as = paste(main_type,'.pred',sep=''))
	umap_combined_labels = cbind(query_sub$meta_data, query_sub$umap)
    plotBasic(umap_combined_labels, 
              title = paste('Query Cells in Reference UMAP Space of', main_type), 
              color.by = paste(main_type,'.pred',sep=''),
              save_path = paste(output_dir, '/4-sub_celltype_pred_',gsub(' ','-',main_type),'.png', sep=''),
              size = c(8,10))
    query_sub$meta_data$celltype.pred.final <- query_sub$meta_data[,paste(main_type,'.pred',sep='')]
    query_sub$meta_data[paste(main_type,'.UMAP1',sep='')] <- query_sub$umap['UMAP1']
    query_sub$meta_data[paste(main_type,'.UMAP2',sep='')] <- query_sub$umap['UMAP2']
    if (i==1){
        sub_results = query_sub$meta_data
        i <- i+1
    } else {
        sub_results = bind_rows(sub_results, query_sub$meta_data)
    }
}

query$meta_data <- cbind(query$meta_data, query$umap)
query$meta_data$celltype.pred.final <- query$meta_data[,paste(maintype_col_name,'.pred',sep='')]
label_main <- query$meta_data %>% filter(!(celltype.pred.final %in% names(celltypes_with_sub)))
print(str(label_main))
label_final <- bind_rows(label_main, sub_results)
label_final$celltype.pred.final <- as.factor(as.character(label_final$celltype.pred.final))

saveRDS(label_final, file=paste(output_dir, '/symphony_celltype_results.rds', sep=''))

print(str(label_final))
print(table(label_final$celltype.pred.final))