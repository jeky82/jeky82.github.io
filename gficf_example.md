---
layout: article
title: "How to use GF-IGF"
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

## Normalization and Clustering (Tabula Muris)
Download Tabula Muris dataset from [HERE](https://drive.google.com/open?id=1yX8IQ7DiWG8PCmYieFFS7vj53Hf1OfT2){:target="_blank"} and
annotation from [HERE](https://drive.google.com/open?id=10ixOOsqZqf6GgwQP1okwoe_TMP_ZTzn5){:target="_blank"}.

```R
library(gficf)
library(ggplot2)

# See function man page for help
?gficf

# Common pipeline to use that goes from normalization to clustering

# Step 1: Nomrmalize data with gficf
data = gficf::gficf(M = readRDS("path/to/TabulaMuris.10x.mouse.RAW.rds"),cell_proportion_max = 1,cell_proportion_min = .05,storeRaw = T,normalize = T)

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

```R
# Plot the relative expression of selected genes
p = plotGenes(data = data,genes = c("Cd34","Cd8"))
p[[1]] + xlab("t-SNE1") + ylab("t-SNE2") + ggtitle("Cd34")
p[[2]] + xlab("t-SNE1") + ylab("t-SNE2") + ggtitle("Cd8a")

```

|                                   |                                 |
|-----------------------------------|---------------------------------|
![Cd8a_expression.png](https://github.com/dibbelab/gficf/blob/master/img/Cd8a_expression.png?raw=true) | ![Cd34_expression.png](https://github.com/dibbelab/gficf/blob/master/img/Cd34_expression.png?raw=true)

## How to embedd new cells in an existing space

Download PBMCs dataset from [HERE](https://drive.google.com/open?id=13cuTP7cjV62Ma4aV9jkzFpR4VoyBBmKj){:target="_blank"}.

```R
library(gficf)
library(ggplot2)

# Common pipeline to use that goes from normalization to clustering

# Step 1: load data and split in training and test set
M = readRDS("/path/to/Purified.PBMC.RAW.rds")

set.seed(0)
ix = sample(x = 1:ncol(M),size = 1000)
M.test = M[,ix] 
M.train = M[,-ix]
rm(M,ix);gc()


# Step 2: Nomrmalize the training set data with gficf
data = gficf::gficf(M = M.train,cell_proportion_max = 1,cell_proportion_min = .05,storeRaw = F,normalize = T)

# Step 3: Reduce data with Latent Semantic Anlysis before to apply t-SNE or UMAP
data = gficf::runPCA(data = data,dim = 50)

# Step 4: Applay UMAP on reduced data and plot cells
# Note: You can pass more parameters directly to umap. a and b are specific parameters controlling the embendding.
# use ?umap for details and additional parameters to use.
data = gficf::runReduction(data = data,reduction = "umap",seed = 0,nt=4,a=2,b=2)

# cell type is contained in the name of each cell
data$embedded$cell.type = unlist(sapply(strsplit(x = rownames(data$embedded),split = ".",fixed = T),function(x) x[1]))
gficf::plotCells(data = data,colorBy = "cell.type")
```

![pbmc.umap.png](https://github.com/dibbelab/gficf/blob/master/img/pbmc.umap.png?raw=true)


```R
# Step 5: We can now embed the new cell in the already existing space and predct thei type 
data = gficf::embedNewCells(data = data,x = M.test,nt = 6,seed = 0)

# Let's add the know cell type for each predicted cell
data$embedded$cell.type = unlist(sapply(strsplit(x = rownames(data$embedded),split = ".",fixed = T),function(x) x[1])) 

# Step 6: Plot results. Embededd cell are shown as triangle and colored according to their original cell type.
ggplot(data = data$embedded,aes(x=X,y=Y,color=cell.type)) + geom_point(aes(shape=predicted,size=predicted)) + theme_bw() + scale_shape_manual(values = c(20,17)) + scale_size_manual(values = c(.1,3))

```

![pbmc.predicted.png](https://github.com/dibbelab/gficf/blob/master/img/pbmc.predicted.png?raw=true)
