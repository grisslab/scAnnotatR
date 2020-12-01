# scTypeR

The `scTypeR` package automatically classifies cells in scRNA-seq datasets. It is simple to use with a clear infrastructure to easily add additional cell type classification models. `scTypeR` support both `Seurat` and `SingleCellExperiment` objects as input.

## Installation

You can install the latest version directly from GitHub using the `devtools` package:

```r
# install devtools if needed
if (!require(devtools)) {
  install.packages("devtools")
}

# install the ReactomeGSA package
if (!require(scTypeR)) {
  install_github("grisslab/SingleCellClassR")
}
```

## Help

The complete usage is shown in the vignettes:

  * [Basic classification of cells](vignettes/classifying-cells.Rmd)
  * [Basic training of a new cell classification model](vignettes/training-basic-model.Rmd)
  * [Training of child-celltype models](vignettes/training-child-model.Rmd)

For more questions / feedback please simply post an [Issue](https://github.com/grisslab/SingleCellClassR/issues/new/choose).

