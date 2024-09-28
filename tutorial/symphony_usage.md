# 2-Step Symphony Pipeline



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

To run the 2-Step Symphony Pipeline, modifiy the configuration file (see details below). Then use the command:

```bash
Rscript ./xenium_symphony_pipeline.R {configuration_file}
```

Here, it is:

```bash
Rscript ./xenium_symphony_pipeline.R ./xenium_symphony_example.config.R
```

## Dependencies

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

## Modify Configuration File