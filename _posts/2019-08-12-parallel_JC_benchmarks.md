---
title: "Parallel Jaccard with RcppParallel"
tags: jaccard-coefficient RcppParallel gf-icf
sidebar:
  nav: docs-en
---

<!-- Global site tag (gtag.js) - Google Analytics -->
<script async src="https://www.googletagmanager.com/gtag/js?id=UA-144257957-1"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-144257957-1');
</script>
  
## Parallel Jaccard Benchmarks
We have implemented a naive but fast parallel version of Jaccard Coefficient estimation for the Phenograph method, thus obtaing a speed boot of about 20X when comparade with previous serial implementation.

|     test       |  elapsed | relative | cell neighbors number|
|----------------|---------------------|-------------|
Parallel Jaccard |  54.822  |  1.000   |   250       |
**Serial Jaccard**   | 1293.185 |  **23.589**  |   250       |
Parallel Jaccard |  7.192   |  1.000   |   100       |
**Serial Jaccard**   | 161.533  |  **22.460**  |   100       |
Parallel Jaccard | 1.614    |  1.000   |   50        |
**Serial Jaccard**   | 38.642   |  **23.942**  |   50        |
Parallel Jaccard | 0.372    |  1.000   |   25        |
**Serial Jaccard**   | 9.546    |  **25.661**  |   25        |
Parallel Jaccard | 0.114    |  1.000   |   15        |
**Serial Jaccard**   | 3.366    |  **29.526**  |   15        |


Below the code used to generate the above table using PBMCs dataset from [HERE](https://drive.google.com/file/d/15pW1JNFz7TjBXuT5Z90h2yq-SO9xrj77/view?usp=sharing){:target="_blank"}.

```R
library(gficf)
library(uwot)
library(RcppParallel)
library(rbenchmark)

# Step 1: load PBMC dstsdet
pbmc = readRDS("/path/to/Purified.PBMC.RAW.rds")
data = gficf::gficf(M = pbmc,cell_proportion_max = 1,cell_proportion_min = .05,storeRaw = F,normalize = F)
rm(pbmc);gc()

# Step 2: Reduce data with PCA
data = gficf::runPCA(data = data,dim = 10)

# Benchmarks
RcppParallel::setThreadOptions(numThreads = 6)
res = NULL
for (neigh.number in c(15,25,50,100,250))
{
  message(paste("Testing neigh = ",neigh.number))
  neigh.idx = uwot:::find_nn(data$pca$cells,k=neigh.number+1,include_self = T,n_threads = 6,verbose = F,method = "annoy")$idx
  neigh.idx = neigh.idx[,-1]

  # compare performance of serial and parallel
  tmp <- benchmark(gficf:::jaccard_coeff(neigh.idx,T),
                   gficf:::rcpp_parallel_jaccard_coef(neigh.idx,T),
                   replications=2,
                   order="elapsed")
  tmp$neigh_number = neigh.number
  res = rbind(tmp,res)
}

print(res[,c("test","elapsed","relative","neigh_number")])
```