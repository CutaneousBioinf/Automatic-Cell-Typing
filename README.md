# Automatic-Cell-Typing
An automatic cell type identification pipeline for spatial transcriptomics(CosMx and Xenium). 




## Step 0: Preprocessing

Preprocessing imports `Seurat` package to find the marker genes of spatial trnscriptomics overlapped with scRNA-seq reference, and use our `spatial-pipeline` to get the preprocessed spatial transcriptomics in the format of giotto object.


## Step 1: Majority Vote



## Step 2: Coarse Classification



## Step 3: Precise Classification


## Evaluation

The purpose of evaluation is to verify whether this pipeline reveals informative biological insights or not. Two types of evaluation will be implemented automatically: 

1) Cell type proportions comparison between scRNA-seq reference and spatial transcriptomics
2) Marker gene identification of each sub-celltype
