---
layout: article
title: "How to use GF-IGF"
tags: gf-icf
sidebar:
  nav: docs-en
aside:
  toc: true
---

## Tabula Muris
Download Tabula Muris dataset from [HERE](https://drive.google.com/open?id=1yX8IQ7DiWG8PCmYieFFS7vj53Hf1OfT2) and
annotation from [HERE](https://drive.google.com/open?id=10ixOOsqZqf6GgwQP1okwoe_TMP_ZTzn5).

```R
library(gficf)
library(ggplot2)

# See function man page for help
?gficf

# Common pipeline to use that goes from normalization to clustering

# Step 1: Nomrmalize data with gficf
data = gficf::gficf(M = readRDS("path/to/TabulaMuris.10x.mouse.RAW.rds"),cell_proportion_max = 1,cell_proportion_min = .05,storeRaw = F,normalize = T)

# Step 2: Reduce data with Latent Semantic Anlysis before to apply t-SNE or UMAP
data = gficf::runLSA(data = data,dim = 50)

# Alternative Step 2: Reduce data with Principal Component Analysis before to apply t-SNE or UMAP
# data = gficf::runPCA(data = data,dim = 50)

# Step 3: Applay t-SNE on reduced data and plot cells
data = gficf::runReduction(data = data,reduction = "tsne",seed = 0,nt=4)
gficf::plotCells(data = data)

# Alternative Step 3: Applay UMAP on reduced data
# data = gficf::runReduction(data = data,reduction = "umap",seed = 0,nt=4)
# gficf::plotCells(data = data)

# Step 4: Cell clustering using Phenograph algorithm
data = gficf::clustcells(data = data,from.embedded = F,dist.method = "manhattan",nt = 4,k = 50,community.algo = "louvian",seed = 0)

# Step 5: Visualize cells by identified clusters
gficf::plotCells(data = data,colorBy="cluster",pointSize = .05) + xlab("t-SNE1") + ylab("t-SNE2") + ggtitle("Cells colored by Clusters") 

```
![tabula_clusters.png](https://github.com/dibbelab/gficf/blob/master/img/tabula_clusters.png?raw=true)

```R
# Additional steps: add annotation to cells and plot it.
info = readRDS("/path/to/TabulaMuris.10x.mouse.annotation.rds")
data$embedded$tissue = info$tissue
data$embedded$subtissue = info$subtissue
data$embedded$cell_ontology_class = info$cell_ontology_class
gficf::plotCells(data = data,colorBy="cell_ontology_class",pointSize = .05) + xlab("t-SNE1") + ylab("t-SNE2") + ggtitle("Cells colored by Clusters") 

```
![tabula_annotated.png](https://github.com/dibbelab/gficf/blob/master/img/tabula_annotated.png?raw=true)

## PBMCs from 10X
