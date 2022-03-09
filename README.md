# `plotmics`  <img src="logo.png" align="right" alt="" width="350" />

_pl**O**t**MICS**_: Visualization of omics and sequencing data in R.

## Information

_pl**O**t**MICS**_ is an R package to visualize omics and sequencing data in R using simple but useful functions:

* `chromReads()`: draws a `ggplot2`-based horizontal barplot with the number and the percentage of reads mapped to each chromosome.
* `barDEGs()`: draws a `ggplot2`-based horizontal barplot with the upregulated and downregulated genes coming from DESeq2.
* `volcanoPlot()`: draws a `ggplot2`-based volcano plot with the DE data coming out from DESeq2.
* `DEcompare()`: draws a `ggplot2`-based scatter plot comparing the Log2FCs of two different contrasts.
* `barAnno()`: takes a list of modified outputs from `ChIPseeker::annotatePeak()`, computes the proportion of peaks in distal or promoter regions and draws a `ggplot2`-based barplot with the correseponding proportion for each sample.
* `getVennCounts()`: helper function that calls `ChIPpeakAnno::makeVennDiagram()` to intersect different sets of peaks and returns the Venn counts and a list of the peaks present in each set of peaks. Note that the intersection of peaks may have a smaller number than expected to prevent the intersection being bigger than one of the sets.
* `upsetPeaks()`: calls `getVennCounts()` and draws an UpSet plot using the Venn counts and the `UpSetR` package.
* `ggUpsetPeaks()`: calls `getVennCounts()` and draws a `ggplot2`-based UpSet plot using the Venn counts.
* `ggVennPeaks()`: calls `getVennCounts()` and draws a Venn diagram using the package `ggvenn`.
* `plotDendogram()`: helper function to draw a dendogram for the heatmaps in `expressionHeatmap()` and `expressionHeatmap2()`.
* `expressionHeatmap()`: function that takes a data frame of expression data (including *Geneid*) and plots a heatmap of the selected genes. It can be any expression data, but it is better for expression values, counts, normalized counts in which the frist column is the *Geneid* and the other columns are each of the samples that are to be plotted.
* `expressionHeatmap2()`: function that takes a list of data frames with the columns *Geneid* and *log2FoldChange* (it can be another type data, such as FPKMs, but its easier to maintain the name of *log2FoldChange* since it is how it output from `DESeq2`) and plots a heatmap with the selected genes. It is better for *log2FoldChange* heatmaps in which each element in the list is a data frame that outputs from `DESeq2`.
* `chromRegions()`: function that takes the chromosome size information and a list of regions (e.g. BED files) to draw the positions of that regions in the genome.
* `circleRegions()`: function that takes a list of chromosome size information and a list of regions (e.g. BED files) to draw the positions of that regions in the genome. Similar to `chromRegions()`, but in a circular plot, allowing different assemblies, pairing regions, etc. 
* `expressionCor()`: function that computes the correlation of expression values between samples.
* `compareCounts()`: function that draws a scatter plot of gene expression values between 2 samples.
* `ggVennBed()`: draws a Venn diagram with the intersection of two BED files, obtained with `bedtoolsr::bt.intersect()`.

## Install `plotmics` 

To install `plotmics` you have to run the following command in R:

```
# if not installed, install the devtools package from CRAN 
if(!require(devtools)) { install.packages("devtools") }

# install plotmics from this Github repository 
devtools::install_github("amitjavilaventura/plotmics")
```

## Contributors

This package has been developed by [Adrià Mitjavila Ventura](https://amitjavilaventura.github.io), with some contributions from [dfernandezperez](https://github.com/dfernandezperez).

If you want to contribute to this package, make a post in the issues section in this repository or fork this repository adding your code and do a pull request.

## Cite

If you use this package, please cite [this repository](https://github.com/amitjavilaventura/plotmics) and give it a star.

## Versions

`plotmics` versions have the structure of `1.2.3`. The first number (`1`, *major*) implies the addition of a function and/or major changes in the packages; the second number (`2`,*minor*) implies the addition of new features to a function and possible corrections; the third number (`3`, *micro*) implies the correction of minor bugs or addition of minor features.

`plotmics` version history is shown below:

* `v1.0.0`: 

  + First version. 
  
<br>

* `v1.1.0`: 
  + `ggVennPeaks()`: Rescale output to remove blank space around the Venn diagram.
  + `expressionHeatmap()` and `expressionHeatmap2()`: Add possibility to scale (`scale()`) data by rows or columns. 
  
<br>

* `v1.1.1`: 
  + `expressionHeatmap()` and `expressionHeatmap2()`: Add minor formatting options (remove the gene names, change sizes of texts and titles, change the color of the cell border, etc).
  + `ggVennPeaks()`: Add minor changes in order to make it easier to visualize more peaks sets.  
  
<br>
  
* `v1.1.2`:
  + `chromReads()`: Change chromosome filtering method.
  + `volcanoPlot()`: Allow dataframes without `DEG` column as input.
  + `volcanoPlot()`: Change `scale = FALSE` for `scale = "none"`. 
  
<br>
  
* `v1.1.3`:
  + `expressionHeatmap()` and `expressionHeatmap2()`: Fix error in labelling.
  
<br>
  
* `v1.1.4`:
  + `expressionHeatmap()` and `expressionHeatmap2()`: Add possibility to color the NA values.
  + `barDEGs()`: Change title format in. 
  
<br>
  
* `v2.0.0` *(2021-09-22)*:
  + Add new function: `chromRegions()`. 
  
<br>
  
* `v2.1.0` *(2021-09-26)*:
  + `chromRegions()`: Allow to take a list of regions as input.
  + `chromRegions()`: Allow to order the region sets.
  
<br>
  
* `v3.0.0` *(2021-09-28)*:
  + Add new function: `circleRegions()`.
  + `chromRegions()`: Allow to color by different parameters.
  + `chromRegions()`: Allow to add extra info.
  + `chromRegions()`: Allow to remove or change size of text in the Y axis.
  
<br>
  
* `v3.0.1` *(2021-09-30)*:
  + `circleRegions()`: Fix a minor bug about plotting the chromosome labels.
  + `ggVennPeaks()`: Fix a minor bug that caused an intersection with one region more than expected.
  + `getVennCounts()`: Remove `pkgcond` from required packages.
  
<br>

* `v3.1.0` *(2021-10-04)*:
  + `ggVennPeaks()`and `getVennCounts()`: Allow to consider strand information.
  
<br>
  
* `v3.1.1` *(2021-10-07)*:
  + `barDEGs()` and `volcanoPlot()`: Change minor features.
  + `circleRegions()`: Exclude the chromosomes before looking if the chromosome names in regions are also in the chrom.sizes files.
  
<br>

* `v3.1.2` *(2021-10-12)*:
  + `circleRegions()`: Write the chromosome labels after drawing all the lines and points, so the labels won't be under many layers.
  + `circleRegions()`: Fix error to allow different assemblies.
  + `barDEGs()`: Add the possibility to add title and subtitle, and also to count the total number of genes in each contrast.
  + `barAnno()`: Add the possibility to add percentage/counts.
  
<br>

* `v4.0.0` *(2021-10-27)*:
  + Add new function: `expressionCor()`.
  + `ggVennPeaks()`: Add the possibility to annotate the number of true overlaps for each set of peaks.
  + `barDEGs()`: Add the possibility to do a `prop.test()` and add the p-value to the plot.

<br>
  
* `v5.0.0` *(2021-11-29)*:
  + Add new function: `compareCounts()`.
  + `volcanoPlot()`: Add the possibility to draw grid lines.
  + `ggVennPeaks()`: Add the possibility to return the lists of overlapping peaks.
  
<br>

* `v5.1.0` *(2021-12-03)*:
  + `chromRegions()`: Add the possibility to color the cytogenetic bands.
  + `expressionHeatmap()`: Add the possibility to have more fields in the data frame so it will be possible to use facets and scale.
  + `volcanoPlot()`: Fix plotting of `degLabels`.
  
<br>

* `v6.0.0` *(2022-02-11)*:
  + Add new function: `ggVennBed()`
  
<br>

* `v6.0.1` *(2022-03-09)*:
  + `volcanoPlot()`: change how DEG labels are written and colored, from `geom_text()` to `annotate()`.
  
<br>
