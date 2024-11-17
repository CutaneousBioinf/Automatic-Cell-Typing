library(Matrix, lib.loc="/home/alextsoi/R/R-4.4/lib/")
library('Seurat',lib.loc="/home/alextsoi/R/R-4.4/lib/")
library('Giotto',lib.loc="/home/alextsoi/R/R-4.4/lib/")
set.seed(0)

########## Parameters ##########
skip_build_ref_main = FALSE
skip_build_ref_sub = FALSE

maintype_col_name = 'celltype'  # column in metadata representing main celltype
subtype_col_name = 'celltype'   # column in metadata representing sub celltype

remove_maintype = c()   # You may want to exclude some celltypes, like missing celltypes. If not, leave it empty, c(). 
remove_subtype = c('')    # It happens at the very beginning, before downsampling and calculating cell proportions. 

downsample = FALSE    # You may need downsampling when you have too many cells. If input is too large, it will cause an error in library Matrix.
downsample_to = 10000  # The number of cells you want to keep. Final number may be slightly different. Ignore it when downsample=FALSE.
vars_use = c('orig.ident') # Covariants to remove
num_PC = 20     # The number of principle components used to build references. Default=20
theta = c(2)   # Diversity clustering penalty parameter. Specify for each variable in vars_use Default theta=2. theta=0 does not encourage any diversity. Larger theta results in more diverse clusters.

save_main_ref_dir = '/home/yulicai/symphony/mucosa/references'
save_main_uwot_dir = save_main_ref_dir
save_sub_ref_dir = '/home/yulicai/symphony/mucosa/references'
save_sub_uwot_dir = save_sub_ref_dir

k = 5   # the k of KNN in mapping.

output_dir <- '/home/yulicai/symphony/mucosa/outputs'   # directory for saving analyzing results



########## Read Inputs For References ##########

## if provide a giotto object directly

## This is for older Giotto version
#giotto_object = readRDS('./example_data/example_refer_giotto.rds')
#
#giotto_object = normalizeGiotto(giotto_object, norm_methods = "standard", logbase=exp(1), scalefactor = 10000, scale_genes = FALSE, scale_cells = FALSE)
#ref_exp = giotto_object@norm_expr
#ref_metadata = as.data.frame(giotto_object@cell_metadata)
#rownames(ref_metadata) = ref_metadata[,1]  # use first column as cell id

## This is for newer Giotto version
#giotto_object = readRDS('./example_data/example_refer_giotto_v4.rds')
#
#giotto_object = normalizeGiotto(giotto_object, norm_methods = "standard", logbase=exp(1), scalefactor = 10000, scale_feats = FALSE, scale_cells = FALSE)
#ref_exp = giotto_object@expression$cell$rna$normalized@exprMat
#ref_metadata = as.data.frame(giotto_object@cell_metadata$cell$rna@metaDT)
#rownames(ref_metadata) = ref_metadata[,1]  # use first column as cell id


## if provide expression matrix and metadata seperately
## Because reading whole scRNA data is time-consuming, I extract the required parts when testing my code.
ref_exp_path = '/home/yulicai/symphony/mucosa/mucosa_expr_raw.rds'
ref_metadata_path = '/home/yulicai/symphony/mucosa/mucosa_metadata.rds'

ref_exp = readRDS(ref_exp_path)	
ref_metadata = readRDS(ref_metadata_path)
# If the expression matric is raw counts, do normalization.
seurat_obj = CreateSeuratObject(ref_exp)
seurat_obj = NormalizeData(seurat_obj)     # log(CP10K + 1) normalization
ref_exp = seurat_obj$RNA$data    # gene x cell, should be log(CP10K + 1) normalized 


## if provide a seurat object directly
#seurat_obj = readRDS('/hits/wasikowr/Novartis/mucosa/seurat.RDS')
#
#seurat_obj = NormalizeData(seurat_obj)     # log(CP10K + 1) normalization
#ref_exp = seurat_obj$RNA@data    # gene x cell, should be log(CP10K + 1) normalized 
#ref_metadata = seurat_obj@meta.data



########## Read Inputs For Queries ##########

### If provide full paths of all cell_features_matrix.h5 file of the xenium data
#query_paths <- as.matrix(read.table("xxxxx"))[,1]
##query_paths <- c('/hits/AGC/Spatial/Xenium/hSkin_481g/10265-JF_and_10628-JF/20240517__214754__10265_10628-JF/20240517__214754__10265_10628-JF/output-XETG00077__0022387__10628-JF-1_ROI_A1__20240517__215149/cell_feature_matrix.h5', '/hits/AGC/Spatial/Xenium/hSkin_480g/10265-JF_and_10628-JF/20240517__214754__10265_10628-JF/20240517__214754__10265_10628-JF/output-XETG00077__0022387__10628-JF-1_ROI_A2__20240517__215149/cell_feature_matrix.h5','/hits/AGC/Spatial/Xenium/hSkin_480g/10265-JF_and_10628-JF/20240517__214754__10265_10628-JF/20240517__214754__10265_10628-JF/output-XETG00077__0022387__10628-JF-1_ROI_B__20240517__215149/cell_feature_matrix.h5','/hits/AGC/Spatial/Xenium/hSkin_480g/10265-JF_and_10628-JF/20240517__214754__10265_10628-JF/20240517__214754__10265_10628-JF/output-XETG00077__0022387__10628-JF-1_ROI_C1__20240517__215149/cell_feature_matrix.h5')
#seurat_objs <- c()
#for (path in query_paths){
#    h5f <- Read10X_h5(path)
#    seurat_obj <- CreateSeuratObject(h5f$`Gene Expression`) # Not sure if it is applied to other h5 files
#    seurat_obj <- NormalizeData(seurat_obj)
#    seurat_objs <- c(seurat_objs, seurat_obj)
#}


## If provide a giotto object
library('Giotto')
giotto_object <- readRDS('/home/yulicai/symphony/mucosa/data/gobject.RDS')
#seurat_obj <- CreateSeuratObject(giotto_object@raw_exprs)   # This is for older Giotto version
seurat_obj <- CreateSeuratObject(as.data.frame(giotto_object@expression$cell$rna$raw@exprMat))    # This is for newer Giotto version (tested for v4.1.0)
#seurat_obj <- CreateSeuratObject(giotto_object@norm_expr)
seurat_obj <- NormalizeData(seurat_obj)
seurat_objs <- c()
seurat_objs <- c(seurat_objs, seurat_obj)




if (length(seurat_objs)>1) {
    seurat_merged <- merge(seurat_objs[[1]], y=seurat_objs[2:length(seurat_objs)], add.cell.ids=paste('obj',as.character(c(1:length(seurat_objs))),sep=''))
    colnames(seurat_merged) <- sub("_", "-", colnames(seurat_merged))
} else {
	seurat_merged <- seurat_objs[[1]]
}



