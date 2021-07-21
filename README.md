# `plotmics`  <img src="logo.png" align="right" alt="" width="350" />

Visualize omics and sequencing data in R.

## Information

The goal of this package is to provide simple functions to visualize several omics and sequencing data:

* `chromReads()`: draws a `ggplot2`-based horizontal barplot with the number and the percentage of reads mapped to each chromosome.
* `barDEGs()`: draws a `ggplot2`-based horizontal barplot with the upregulated and downregulated genes coming from DESeq2.
* `volcanoPlot()`: draws a `ggplot2`-based volcano plot with the DE data coming out from DESeq2.
* `DEcompare()`: draws a `ggplot2`-based scatter plot comparing the Log2FCs of two different contrasts.
* `barAnno()`: takes a list of modified outputs from `ChIPseeker::annotatePeak()`, computes the proportion of peaks in distal or promoter regions and draws a `ggplot2`-based barplot with the correseponding proportion for each sample.
* `getVennCounts()`: helper function that calls `ChIPpeakAnno::makeVennDiagram()` to intersect different sets of peaks and returns the Venn counts and a list of the peaks present in each set of peaks.
* `upsetPeaks()`: calls `getVennCounts()` and draws an UpSet plot using the Venn counts and the `UpSetR` package.
* `ggUpsetPeaks()`: calls `getVennCounts()` and draws a `ggplot2`-based UpSet plot using the Venn counts.
* `ggVennPeaks()`: calls `getVennCounts()` and draws a Venn diagram using the package `ggvenn`.
* `plotDendogram()`: helper function to draw a dendogram for the heatmaps in `expressionHeatmap()` and `expressionHeatmap2()`.
* `expressionHeatmap()`: function that takes a data frame of expression data (including *Geneid*) and plots a heatmap of the selected genes. It can be any expression data, but it is better for expression values, counts, normalized counts in which the frist column is the *Geneid* and the other columns are each of the samples that are to be plotted.
* `expressionHeatmap2()`: function that takes a list of data frames with the columns *Geneid* and *log2FoldChange* (it can be another type data, such as FPKMs, but its easier to maintain the name of *log2FoldChange* since it is how it output from `DESeq2`) and plots a heatmap with the selected genes. It is better for *log2FoldChange* heatmaps in which each element in the list is a data frame that outputs from `DESeq2`.

## Install `plotmics` 

To install `plotmics` you have to run the following command in R:

```
# if not installed, install the devtools package from CRAN 
if(!require(devtools)) { install.packages("devtools") }

# install plotmics from this Github repository 
devtools::install_github("amitjavilaventura/plotmics")
```

## Contributors

This package has been developed by [Adrià Mitjavila Ventura](https://amitjavilaventura.github.io), with some contributions from [dfernandezperez](https://github.com/dfernandezperez)

If you want to contribute to this package, make a post in the issues section in this repository or fork this repository adding your code and do a pull request.

## Cite

If you use this package, please cite [this repository](https://github.com/amitjavilaventura/plotmics) and give it a star.

## Updates

* `v1.0.0`: first version.
* `v1.1.0`: 
  + Rescale `ggVennPeaks()` output to remove blank space around the Venn diagram.
  + Add possibility to scale (`scale()`) data by rows or columns in `expressionHeatmap()` and `expressionHeatmap2()`.
