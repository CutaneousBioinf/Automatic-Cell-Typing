# Automatic-Cell-Typing
An automatic cell type identification pipeline for spatial transcriptomics(CosMx and Xenium). 

The motivation of developing this pipeline is using the rich gene expression information from scRNA-seq to address the issue of sub-celltype classification in spatial data.

This pipeline includes three main steps, one preprocessing and evaluation. All the things showed here can be done automatically.


## Step 0: Preprocessing

Preprocessing imports `Seurat` package to find the marker genes of spatial trnscriptomics overlapped with scRNA-seq reference, and use our `spatial-pipeline` to get the preprocessed spatial transcriptomics in the format of giotto object.

The input of this step is the file path of scRNA-seq reference and spatial transcriptomics, and the output contains two kinds of files, a giotto object and expression matrix. Specifically, the expression matrix for main celltypes and all sub-celltypes contained in each main celltype. 

## Step 1: Majority Vote

Majority Vote includes two parts, the first is using `Insitutype` package to perform supervised-celltype-classification for main celltypes. And the second is using majority voting to give every cluster a label inlcuding all the possible main celltypes. 

## Step 2: Coarse Classification

Coarse Classification aims at assigning each cell to one specific main celltype. 
For each leiden cluster that has score criteria larger than 95%, which is defined by more than one celltype in Step 1, a subset Giotto object of main celltypes included in the mixture cluster will be extracted.

Run PCA, UMAP, and Harmony to get the renewed umap info for each selected cell
InsitutypeML:
(1) Prior info: umap and global coordinates
(2) Reference matrix: Marker genes* main cell types, delete all the overlapped marker genes

## Step 3: Precise Classification



## Evaluation

The purpose of evaluation is to verify whether this pipeline reveals informative biological insights or not. Two types of evaluation will be implemented automatically: 

1) Cell type proportions comparison between scRNA-seq reference and spatial transcriptomics
2) Marker gene identification of each sub-celltype
