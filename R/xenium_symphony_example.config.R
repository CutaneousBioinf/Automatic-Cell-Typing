########## Set Parameters ##########
set.seed(0)

## Build reference
skip_build_ref_main = FALSE
skip_build_ref_sub = FALSE # nolint


# if provide a giotto object directly (the Insitutype pipeline will generate a QCed giotto object)
giotto_object = readRDS('/home/alextsoi/Researches/Novartis_Xenium_LP/analysis_mucosa/gobj')
maintype_col_name = 'celltype.global'  # column in metadata representing main celltype
subtype_col_name = 'celltype.sub'   # column in metadata representing sub celltype

giotto_object = normalizeGiotto(giotto_object, norm_methods = "standard", logbase=exp(1), scalefactor = 10000)
ref_exp = giotto_object@expression$cell$rna$normalized@exprMat
ref_metadata = giotto_object$cell$rna@metaDT


## if provide expression matrix and metadata seperately
#ref_exp_path = '/home/yulicai/symphony/exprs_norm.rds'	# gene x cell, should be log(CP10K + 1) normalized 
#ref_metadata_path = '/home/yulicai/symphony/metadata.rds'
#maintype_col_name = 'celltype.global'  # column in metadata representing main celltype
#subtype_col_name = 'celltype.sub'   # column in metadata representing sub celltype
#
#ref_exp = readRDS(ref_exp_path)	
#ref_metadata = readRDS(ref_metadata_path)


## if provide a seurat object directly
#seurat_obj = readRDS('/home/alextsoi/db/psoriasis.scRNA/seurat_with_rachael_subtype.RDS')
#maintype_col_name = 'celltype.global'  # column in metadata representing main celltype
#subtype_col_name = 'celltype.sub'   # column in metadata representing sub celltype
#
#seurat = NormalizeData(seurat_obj)
#ref_exp = seurat_obj$RNA@data    # gene x cell, should be log(CP10K + 1) normalized 
#ref_metadata = seurat_obj@meta.data


downsample = TRUE    # You may need downsampling when you have too many cells. If input is too large, it will cause an error in library Matrix.
downsample_to = 10000  # The number of cells you want to keep. Final number may be slightly different. Ignore it when downsample=FALSE.
vars_use = c('orig.ident') # Harmony parameter. If meta_data is dataframe, this defined which variable(s) to remove (character vector).

save_main_ref_dir = '/home/alextsoi/test/symphony'
save_main_uwot_dir = '/home/alextsoi/test/symphony'
save_sub_ref_dir = '/home/alextsoi/test/symphony'
save_sub_uwot_dir = '/home/alextsoi/test/symphony'


## Query
output_dir <- '/home/alextsoi/test/symphony'

### provide full paths of all cell_features_matrix.h5 file of the xenium data
query_paths <- as.matrix(read.table("xxxxx"))[,1]
###query_paths <- c('/hits/AGC/Spatial/Xenium/hSkin_481g/10265-JF_and_10628-JF/20240517__214754__10265_10628-JF/20240517__214754__10265_10628-JF/output-XETG00077__0022387__10628-JF-1_ROI_A1__20240517__215149/cell_feature_matrix.h5',              '/hits/AGC/Spatial/Xenium/hSkin_480g/10265-JF_and_10628-JF/20240517__214754__10265_10628-JF/20240517__214754__10265_10628-JF/output-XETG00077__0022387__10628-JF-1_ROI_A2__20240517__215149/cell_feature_matrix.h5','/hits/AGC/Spatial/Xenium/hSkin_480g/10265-JF_and_10628-JF/20240517__214754__10265_10628-JF/20240517__214754__10265_10628-JF/output-XETG00077__0022387__10628-JF-1_ROI_B__20240517__215149/cell_feature_matrix.h5','/hits/AGC/Spatial/Xenium/hSkin_480g/10265-JF_and_10628-JF/20240517__214754__10265_10628-JF/20240517__214754__10265_10628-JF/output-XETG00077__0022387__10628-JF-1_ROI_C1__20240517__215149/cell_feature_matrix.h5')


k = 5   # the k of KNN in mapping.

