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

## Modify Configuration File