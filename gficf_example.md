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
Download Tabula Muris dataset from [HERE](https://drive.google.com/file/d/1rBClWTzRtxLGJ8JUT3MM_p9w7Ec80yCf/view?usp=sharing){:target="_blank"} and
annotation from [HERE](https://drive.google.com/file/d/1wDaF6ONd59cDKbIBijpE2Cc2ViE7iGG9/view?usp=sharing){:target="_blank"}.

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
data$embedded$tissue = info$tissue[match(rownames(data$embedded),info$id)]
data$embedded$subtissue = info$subtissue[match(rownames(data$embedded),info$id)]
data$embedded$cell_ontology_class = info$cell_ontology_class[match(rownames(data$embedded),info$id)]
gficf::plotCells(data = data,colorBy="cell_ontology_class",pointSize = .05) + xlab("t-SNE1") + ylab("t-SNE2") + ggtitle("Cells colored by Types") 

```
![tabula_annotated.png](https://github.com/dibbelab/gficf/blob/master/img/tabula_annotated.png?raw=true)

```R
# Plot the relative expression of selected genes
p = gficf::plotGenes(data = data,genes = c("Cd34","Cd8a"))
p[[1]] + xlab("t-SNE1") + ylab("t-SNE2") + ggtitle("Cd34")
p[[2]] + xlab("t-SNE1") + ylab("t-SNE2") + ggtitle("Cd8a")

```

|                                   |                                 |
|-----------------------------------|---------------------------------|
![Cd8a_expression.png](https://github.com/dibbelab/gficf/blob/master/img/Cd8a_expression.png?raw=true) | ![Cd34_expression.png](https://github.com/dibbelab/gficf/blob/master/img/Cd34_expression.png?raw=true)

## How to embedd new cells in an existing space

Download PBMCs dataset from [HERE](https://drive.google.com/file/d/15pW1JNFz7TjBXuT5Z90h2yq-SO9xrj77/view?usp=sharing){:target="_blank"}.

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

## How to perform GSEA to identify active pathways in each cluster

Download PBMCs dataset from [HERE](https://drive.google.com/file/d/15pW1JNFz7TjBXuT5Z90h2yq-SO9xrj77/view?usp=sharing){:target="_blank"}.   
Download gmt file containing the gene sets from [HERE](https://drive.google.com/file/d/1_N8-LCXJrPgGyuZwdEfxQ7bdLWdzZdi3/view?usp=sharing){:target="_blank"}.

```R
library(gficf)
library(ggplot2)

# Example on how use GSEA for the identification
# of active pathways for each identified cluster

# Step 1: load PBMC dstsdet
pbmc = readRDS("/path/to/Purified.PBMC.RAW.rds")

# Step 2: Nomrmalize the training set data with gficf
data = gficf::gficf(M = pbmc,cell_proportion_max = 1,cell_proportion_min = .05,storeRaw = F,normalize = T)
rm(pbmc);gc()

# Step 3: Reduce data with PCA before to apply t-SNE or UMAP
data = gficf::runPCA(data = data,dim = 50)

# Step 4: Applay umap on reduced data and plot cells
# see ?umap for what a,b,n_neighbors and metric parameters are.
data = gficf::runReduction(data = data,reduction = "umap",seed = 0,nt = 2,a=2,b=2,n_neighbors=30,metric="manhattan")

# cell type is contained in the name of each cell
data$embedded$cell.type = unlist(sapply(strsplit(x = rownames(data$embedded),split = ".",fixed = T),function(x) x[1]))
p1 = gficf::plotCells(data = data,colorBy = "cell.type")
print(p1)

# Step 5: Cluster cells with phenograph method and visualize the results
# see ?clustcells for more details
data = gficf::clustcells(data = data,dist.method = "manhattan",nt = 2,k = 50,community.algo = "louvian 2",seed = 0,resolution = .75,n.start = 25,n.iter = 50)
p2 = gficf::plotCells(data = data,colorBy = "cluster")
print(p2)
```

|        Plot p1 (by cell type)     |       Plot p2 (by clusters)     |
|-----------------------------------|---------------------------------|
|![PBMC_umap.png](https://github.com/jeky82/jeky82.github.io/blob/master/img/PBMC_umap.png?raw=true)|![PBMC_clusters.png](https://github.com/jeky82/jeky82.github.io/blob/master/img/PBMC_clusters.png?raw=true)|


Now that we have identified clusters we can use GSEA to identify pathway activity
across cells of the same cluster. Briefly Gene Ranks are first summedd across cells
of the same cluster and then GSEA is performed.   
In this example we use the 50 hallmarks pathways from mSigDB. However You can find additional gene sets files (gmt) on
mSigDB [HERE](http://software.broadinstitute.org/gsea/msigdb/collections.jsp){:target="_blank"}. **Only the one with official gene symbols are
actually supported.**   
   
**Tips:** Use different combination of convertToEns and convertHu2Mm in the function runGSEA to convert gene set files in the one you need before to perform GSEA. Choose the right combination to match identifier you are using to rapresent genes.   

**Examples:**   
convertToEns = **T** and convertHu2Mm = **F** **-->** human symbols are converted to human ensamble id   
convertToEns = **T** and convertHu2Mm = **T** **-->** human symbols are converted to mouse ensamble id   
convertToEns = **F** and convertHu2Mm = **T** **-->** human symbols are converted to mouse symbols   
convertToEns = **F** and convertHu2Mm = **F** **-->** no conversions, original symbols in the gmt file are used   

```R
# Step 6: Run GSEA to identify active pathways in each group of cells
# ?runGSEA for details

gmt.file.path = "/path/to/h.all.v6.2.symbols.gmt"
data = gficf::runGSEA(data = data,gmt.file = gmt.file.path,nsim = 10000,convertToEns = T,convertHu2Mm = F,minSize = 15,nt = 4)

# Step 7: Plot GSEA results to show significant pathways in eac cluster
# Note: To be comparable across pathways the Normalized Enrichement Scores (NES) are plotted.
p3 = gficf::plotGSEA(data = data,fdr = .1)
print(p3)

# Note: Info about pathways and number of used genes are in data$gsea$stat dataframe
print(head(data$gsea$stat))

# Step 8: Plot pathway activity across cells in the embedded space
p4 = gficf::plotPathway(data = data, pathwayName = "HALLMARK_PI3K_AKT_MTOR_SIGNALING",fdr = .1)
print(p4)
```

|        Plot p3 (GSEA results)     |Plot p4 (PI3K AKT MTOR pathway activity)|
|-----------------------------------|---------------------------------|
|![PBMC_gsea.png](https://github.com/jeky82/jeky82.github.io/blob/master/img/PBMC_gsea.png?raw=true)|![PBMC_mtor.png](https://github.com/jeky82/jeky82.github.io/blob/master/img/PBMC_mtor.png?raw=true)|

