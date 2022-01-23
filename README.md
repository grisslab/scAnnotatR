# scAnnotatR

The `scAnnotatR` package automatically classifies cells in scRNA-seq datasets. It is simple to use with a clear infrastructure to easily add additional cell type classification models. `scAnnotatR` support both `Seurat` and `SingleCellExperiment` objects as input.

## Installation

You can install the latest version directly from GitHub using the `devtools` package:

```r
# install devtools if needed
if (!require(devtools)) {
  install.packages("devtools")
}

if (!require(scAnnotatR)) {
  install_github("grisslab/scAnnotatR")
}
```

## Help

The complete usage is shown in the vignettes:

  * [Basic classification of cells](vignettes/classifying-cells.Rmd)
  * [Basic training of a new cell classification model](vignettes/training-basic-model.Rmd)
  * [Training of child-celltype models](vignettes/training-child-model.Rmd)

For more questions / feedback please simply post an [Issue](https://github.com/grisslab/scAnnotatR/issues/new).

## Citation

If you used scAnnotatR in your research, we would be grateful if you could cite the following manuscript:

Nguyen, V., Griss, J. scAnnotatR: framework to accurately classify cell types in single-cell RNA-sequencing data. BMC Bioinformatics 23, 44 (2022). https://doi.org/10.1186/s12859-022-04574-5
