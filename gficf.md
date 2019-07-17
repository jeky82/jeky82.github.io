---
layout: article
title: "GF-IGF"
tags: gf-icf
sidebar:
  nav: docs-en
aside:
  toc: true
---

## What is GF-IGF
An R implementation of the Gene Frequency - Inverse Cell Frequency method for single cell data
normalization [(Gambardella et al. 2019)](https://www.frontiersin.org/articles/10.3389/fgene.2019.00734/abstract).
The package also includes [Phenograph](https://www.cell.com/cell/fulltext/S0092-8674(15)00637-6)
[Louvain method](https://sites.google.com/site/findcommunities/)
clustering using [RcppAnnoy](https://cran.r-project.org/package=RcppAnnoy) library
from [uwot](https://github.com/jlmelville/uwot).
The package also include data reduction with either Principal Component Analisys (PCA) or
Latent Semantic Anlisys (LSA) before to apply t-SNE or UMAP for single cell data visualization.
  
## Installing
  
### From github
  
  `gficf` makes use of annoy library in `uwot`. So you may have to carry out
a few extra steps before being able to build this package like for `uwot` installation:
  
  **Windows**: install 
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) and ensure 
`C:\Rtools\bin` is on your path.

**Mac OS X**: using a custom `~/.R/Makevars` 
[may cause linking errors](https://github.com/jlmelville/uwot/issues/1).
This sort of thing is a potential problem on all platforms but seems to bite
Mac owners more.
[The R for Mac OS X FAQ](https://cran.r-project.org/bin/macosx/RMacOSX-FAQ.html#Installation-of-source-packages)
may be helpful here to work out what you can get away with. To be on the safe
side, I would advise building `uwot` without a custom `Makevars`.
                         
```R
if(!require(devtools)){
  install.packages("devtools") # If not already installed
}
devtools::install_github("dibbelab/gficf")
```
