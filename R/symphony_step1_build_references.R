library('Seurat')
library('symphony')
library('tibble')
library('dplyr')
library('irlba')


########## Set Parameters ##########
set.seed(0)

# Build reference
ref_exp = readRDS('/home/yulicai/symphony/exprs_norm.rds')	# janssen ad
ref_metadata = readRDS('/home/yulicai/symphony/metadata.rds')
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


for (p in c(save_main_ref_dir, save_main_uwot_dir, save_sub_ref_dir, save_sub_uwot_dir)) {
    if (!dir.exists(p)){
        dir.create(p, recursive = TRUE)
    }
}



########## Build the Reference for Main Celltypes ##########
print('Start building reference for main celltypes...')
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















