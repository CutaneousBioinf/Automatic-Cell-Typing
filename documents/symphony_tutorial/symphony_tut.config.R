library(Matrix, lib.loc="/home/alextsoi/R/R-4.4/lib/")
library('Seurat',lib.loc="/home/alextsoi/R/R-4.4/lib/")
library('Giotto')
set.seed(0)

########## Parameters ##########
skip_build_ref_main = FALSE
skip_build_ref_sub = FALSE

maintype_col_name = 'main_cell_type'  # column in metadata representing main celltype
subtype_col_name = 'cell_type'   # column in metadata representing sub celltype

downsample = FALSE    # You may need downsampling when you have too many cells. If input is too large, it will cause an error in library Matrix.
downsample_to = 10000  # The number of cells you want to keep. Final number may be slightly different. Ignore it when downsample=FALSE.
vars_use = c('sample_id') # Covariants to remove

save_main_ref_dir = './references' 
save_main_uwot_dir = save_main_ref_dir
save_sub_ref_dir = './references'
save_sub_uwot_dir = save_sub_ref_dir

k = 5   # the k of KNN in mapping.

output_dir <- './outputs'   # directory for saving analyzing results



########## Read Inputs For References ##########

## if provide a giotto object directly
giotto_object = readRDS('./example_data/example_refer_giotto.rds')

giotto_object = normalizeGiotto(giotto_object, norm_methods = "standard", logbase=exp(1), scalefactor = 10000, scale_genes = FALSE, scale_cells = FALSE)
ref_exp = giotto_object@norm_expr
ref_metadata = as.data.frame(giotto_object@cell_metadata)
rownames(ref_metadata) = ref_metadata[,1]  # use first column as cell id


## if provide expression matrix and metadata seperately
#ref_exp_path = ''
#ref_metadata_path = ''
#
#ref_exp = readRDS(ref_exp_path)	
#ref_metadata = readRDS(ref_metadata_path)
## If the expression matric is raw counts, do normalization.
#seurat_obj = CreateSeuratObject(ref_exp)
#seurat_obj = NormalizeData(seurat_obj)     # log(CP10K + 1) normalization


## if provide a seurat object directly
#seurat_obj = readRDS('')
#
#seurat_obj = NormalizeData(seurat_obj)     # log(CP10K + 1) normalization
#ref_exp = seurat_obj$RNA@data    # gene x cell, should be log(CP10K + 1) normalized 
#ref_metadata = seurat_obj@meta.data



########## Read Inputs For Queries ##########

### If provide full paths of all cell_features_matrix.h5 file of the xenium data
#query_paths <- as.matrix(read.table("xxxxx"))[,1]
##query_paths <- c('h5_file1', 'h5_file2')
#for (path in query_paths){
#    h5f <- Read10X_h5(path)
#    seurat_obj <- CreateSeuratObject(h5f$`Gene Expression`) # Not sure if it is applied to other h5 files
#    seurat_obj <- NormalizeData(seurat_obj)
#    seurat_objs <- c(seurat_objs, seurat_obj)
#}


## If provide a giotto object
library('Giotto')
giotto_object <- readRDS('./example_data/example_refer_giotto.rds')
seurat_obj <- CreateSeuratObject(giotto_object@raw_exprs)
seurat_obj <- NormalizeData(seurat_obj)
seurat_objs <- c()
seurat_objs <- c(seurat_objs, seurat_obj)




if (length(seurat_objs)>1) {
    seurat_merged <- merge(seurat_objs[[1]], y=seurat_objs[2:length(seurat_objs)], add.cell.ids=paste('obj',as.character(c(1:length(seurat_objs))),sep=''))
    colnames(seurat_merged) <- sub("_", "-", colnames(seurat_merged))
} else {
	seurat_merged <- seurat_objs[[1]]
}



