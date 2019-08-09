---
layout: article
title: "GF-IGF"
tags: gf-icf
sidebar:
  nav: docs-en
aside:
  toc: true
---

<!-- Global site tag (gtag.js) - Google Analytics -->
<script async src="https://www.googletagmanager.com/gtag/js?id=UA-144257957-1"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-144257957-1');
</script>

## What is GF-IGF
An R implementation of the Gene Frequency - Inverse Cell Frequency method for single cell data
normalization [(Gambardella et al. 2019)](https://www.frontiersin.org/articles/10.3389/fgene.2019.00734/abstract){:target="_blank"}.
The package also includes [Phenograph](https://www.cell.com/cell/fulltext/S0092-8674(15)00637-6){:target="_blank"}
[Louvain method](https://sites.google.com/site/findcommunities/){:target="_blank"}
clustering using [RcppAnnoy](https://cran.r-project.org/package=RcppAnnoy){:target="_blank"} library
from [uwot](https://github.com/jlmelville/uwot){:target="_blank"}.
The package also include data reduction with either Principal Component Analisys (PCA) or
Latent Semantic Anlisys (LSA) before to apply t-SNE or UMAP for single cell data visualization.
  
## Installing
  
### From github
  
  `gficf` makes use of annoy library in `uwot`. So you may have to carry out
a few extra steps before being able to build this package like for `uwot` installation:
  
  **Windows**: install 
[Rtools](https://cran.r-project.org/bin/windows/Rtools/){:target="_blank"} and ensure 
`C:\Rtools\bin` is on your path.

**Mac OS X**: using a custom `~/.R/Makevars` 
[may cause linking errors](https://github.com/jlmelville/uwot/issues/1){:target="_blank"}.
This sort of thing is a potential problem on all platforms but seems to bite
Mac owners more.
[The R for Mac OS X FAQ](https://cran.r-project.org/bin/macosx/RMacOSX-FAQ.html#Installation-of-source-packages){:target="_blank"}
may be helpful here to work out what you can get away with. To be on the safe
side, I would advise building `uwot` without a custom `Makevars`.
                         
```R
if(!require(devtools)){
  install.packages("devtools") # If not already installed
}
devtools::install_github("dibbelab/gficf")

# Install also edgeR
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")
```
## Phenograph Implementation Details

In the package `gficf` the function `clustcells` implement the [Phenograph](https://www.cell.com/cell/fulltext/S0092-8674(15)00637-6) algorithm,
which is a clustering method designed for high-dimensional single-cell data analysis. It works by creating a graph ("network") representing phenotypic similarities between cells by calculating the Jaccard coefficient between nearest-neighbor sets, and then identifying communities using the well known [Louvain method](https://sites.google.com/site/findcommunities/) in this graph. 

In this particular implementation of Phenograph we use approximate nearest neighbors found using [RcppAnnoy](https://cran.r-project.org/package=RcppAnnoy)
libraries present in the `uwot` package. The supported distance metrics for KNN (set by the `dist.method` parameter) are:

* Euclidean (default)
* Cosine
* Manhattan
* Hamming

Please note that the Hamming support is a lot slower than the
other metrics. It is not recomadded to use it if you have more than a few hundred
features, and even then expect it to take several minutes during the index 
building phase in situations where the Euclidean metric would take only a few
seconds.

After computation of Jaccard distances among cells, the Louvain community detection is instead performed using `igraph` implementation.
All supported communities detection algorithm (set by the `community.algo` parameter) are:

* Louvain (default)
* Louvian with modularity optimization (c++ function imported from Seurat)
* Louvain algorithm with multilevel refinement (c++ function imported from Seurat)
* Walktrap
* Fastgreedy
