# Automatic-Cell-Typing
An automatic cell type identification pipeline for spatial transcriptomics(CosMx and Xenium). 

The motivation of developing this pipeline is using the rich gene expression information from scRNA-seq to address the issue of sub-celltype classification in spatial data.

This pipeline includes three main steps, one preprocessing and evaluation. All the things showed here can be done automatically. The workflow is demonstrated as below:

<div style="text-align:center;">
    <img src="pic/image.png" alt="Workflow" width="800px">
</div>

The input of this pipeline includes the path to the scRNA-seq reference, and the path to the spatial transcriptomics dataset. The output of this pipeline is a renewed Giotto object with two kinds of annotations: automatic-cell-typing using [Insitutype](https://github.com/Nanostring-Biostats/InSituType) and that using [Symphony](https://github.com/immunogenomics/symphony?tab=readme-ov-file). We also provided a spearman correlation matrix to make comparisons between these two methods.

The umap plot below demonstrates the final performance of our pipeline in multi-sample problem. 

<div style="text-align:center;">
    <img src="pic/8.18_8.1_UMAP_sequenced.cell.types.subcelltyping_isML.png" alt="Final result" width="480px">
</div>

## Step 0: Preprocessing

Preprocessing imports `Seurat` package to find the marker genes of spatial trnscriptomics overlapped with scRNA-seq reference, and use the [spatial-pipeline](https://github.com/CutaneousBioinf/spatial-pipeline/tree/main/Xenium) developed by our team to get the preprocessed spatial transcriptomics in the format of a Giotto object.

Users should first tell us whether the dataset they are going to analyze is CosMx or Xenium. Because the preprocessing in `spatial-pipeline` will be slightly different.

The input of this step includes the file path of scRNA-seq reference and spatial transcriptomics, and the output contains two kinds of files, a Giotto object and expression matrix. Specifically, the expression matrix for main celltypes and all sub-celltypes contained in each main celltype. 

## Step 1: Majority Vote

Majority Vote involves clustering cells by celltype, calculating and normalizing the percentages of each cell type within clusters, and identifying main cell types based on cumulative percentage thresholds. Clusters needing further sub-clustering are flagged. Majority Vote applies both the thoughts from the school of probability and the school of Bayesian. If users have any prior knowledge about the celltypes, they are approved to slightly adjust the results.

The input of both CosMx and Xenium includes the file path of scRNA-seq reference and generated Giotto object after preprocessing. While for CosMx datasets, users should also tell us the path of the immunofluorescence data.
The method then updates the Giotto object with new cell type assignments and saves it, while identifying clusters that require additional sub-clustering.

## Step 2: Coarse Classification

Coarse Classification aims at assigning each cell to one specific main celltype. 
For each leiden cluster that has score criteria larger than 95%, which is defined by more than one celltype in Step 1, a subset Giotto object of main celltypes included in the mixture cluster will be extracted. Then insitutypeML is applied to 

## Step 3: Precise Classification

Precise Classification further classifies the cells within each main cell type and identify each sub-celltype.
Clustering results from [BASS](https://github.com/zhengli09/BASS) are used as prior info. 

The output of this step contains a renewed giotto object, and a umap plot of whole dataset.



## Evaluation

The purpose of evaluation is to verify whether this pipeline reveals informative biological insights or not. Two types of evaluation will be implemented automatically: 

1) Cell type proportions comparison between scRNA-seq reference and spatial transcriptomics
2) Marker gene identification of each sub-celltype
