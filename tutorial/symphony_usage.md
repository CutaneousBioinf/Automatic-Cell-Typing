# 2-Step Symphony Pipeline

## Contents

1. [Basic Usage](#Basic-Usage)
2. [Dependencies](#Dependencies)
3. [Configuration](#Configuration)

    1. [Parameters](#Parameters)
    2. [Inputs For References](#inputs-for-references)
    3. [Inputs For Queries](#inputs-for-queries)


## Basic Usage

First, clone the repository. 

```bash
git clone https://github.com/CutaneousBioinf/Automatic-Cell-Typing.git
```

Enter the directory `R/`.

```bash
cd Automatic-Cell-Typing/R
ls
```

You can find two files in this folder:

`xenium_symphony_pipeline.R`: 2-Step Symphony Pipeline

`xenium_symphony_example.config.R`: Configuration of 2-Step Symphony Pipeline

To run the 2-Step Symphony Pipeline, modifiy the configuration file (see details [below](#Modify-Configuration)). Then use the command:

```bash
Rscript ./xenium_symphony_pipeline.R {configuration_file}
```

Here, it is:

```bash
Rscript ./xenium_symphony_pipeline.R ./xenium_symphony_example.config.R
```

## Dependencies <a name="dependencies"></a>

If some libraries are missing, you can use these commands to install libraries in R.

```R
install.packages("Seurat")
install.packages("symphony")
install.packages("tibble")
install.packages("dplyr")
install.packages("irlba")
install.packages("ggplot2")
install.packages("ggthemes")
install.packages("ggrastr")
install.packages("RColorBrewer")
install.packages("patchwork")
install.packages("ggpubr")
install.packages("cowplot")
devtools::install_github('immunogenomics/singlecellmethods')
devtools::install_github("drieslab/Giotto@suite")
```

## Configuration

Current version of the pipeline requires users to manually modify the configuration file. The configuration file contains 3 parts.

### Parameters

Here are explanations of these parameters.

`skip_build_ref_main`: Skip building references for main celltypes or not.

`skip_build_ref_sub`: Skip building references for sub-celltypes or not.

`maintype_col_name`: The name of the column representing main celltype in metadata.

`subtype_col_name`: The name of the column representing sub-celltype in metadata. If you only have main celltypes, set `subtype_col_name = maintype_col_name`.

`downsample`: Use a subset of cells to build references or not.

`downsample_to`: The number of cells in the subset. Ignore it when `downsample=FALSE`. 

NOTICE: The pipeline first calculates `downsample_to` divided by total number of cells in the data for building reference. Then use the ratio to do downsampling for each sub-celltype. Therefore, the proportions of main celltypes and sub-celltypes will not change a lot after the downsampling, and the number of cells in the subset may not be exactly `downsample_to`. When building references for sub-celltypes, if the number of cells in the main celltype is larger than `downsample_to` and `downsample=TRUE`, it will also downsample the cells to `downsample_to`.

`vars_use`: Column in meta_data that defines dataset for each cell. If meta_data is dataframe, this defined which variable(s) to remove (character vector).

`save_main_ref_dir`: The path of the directory saving references for main celltypes.

`save_main_uwot_dir`: The path of the directory saving uwot models for main celltypes.

`save_sub_ref_dir`: The path of the directory saving references for sub-celltypes.

`save_sub_uwot_dir`: The path of the directory saving uwot models for sub-celltypes.

NOTICE: Each reference has an absolute path of corresponding uwot modle. If the path of a uwot mode changes, the reference cannot find it.

`k`: The K of KNN. After projecting reference cells and queries, KNN is used to assign celltypes.

`output_dir`: The path of the directory saving analyzing results.

### Inputs For References

Here we provide some examples for reading inputs:

1. A Giotto object.
2. Expression matrix and metadata are provided seperately.
3. A Seurat object.

To use these examples, users only need to modify paths to inputs.

Users can use their own codes, make sure that it generates these two variables:

`ref_exp`: Expression matrix for references, gene x cell, log(CP10K + 1) normalized. 

log(CP10K + 1) normalzation: for each count, calculate

$ln(\frac{count}{total\ counts\ in\ the\ cell} \cdot 1000+1).$

It is the defult normalization method of Seurat, but it is not the case for Giotto.

`ref_metadata`: Metadata matrix for references, cell x feature, including main celltypes and sub-celltypes.

### Inputs For Queries

Here we provide some examples for reading inputs:

1. Multiple h5 files.
2. A Giotto object.

Users can use their own codes, make sure that it generates this variable:

`seurat_objs`: A list with one or multiple Seurat object(s).
