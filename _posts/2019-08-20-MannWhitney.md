---
title: "RcppParallel Mann–Whitney U test"
tags: gf-icf RcppParallel Rcpp Mann–Whitney
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
  
## Embed new cells in an already existing space
I have implemented an RcppParalle and Rcpp version of Mann–Whitney U test with continuity correction, that using 6 threads was about **100X faster then native R function** and **5X faster then corresponding serial C++ function**. Test were performed using i7 with 6 cores.   

The implemented C++ Mann–Whitney U test correspond to the R command wilcox.test(x,y,alternative = "two.sided", paired = F,exact = F,correct = T). 

Below the code used to generate the above table.

```R
set.seed(0)
ngenes = 3000;
ncells = 5000;
expr = matrix( rnorm(ngenes*ncells,mean=0,sd=1), ngenes, ncells)


# Benchmarks
RcppParallel::setThreadOptions(numThreads = 6)

ix1 = 1:round(ncells/2)
ix2 = (round(ncells/2)+1):ncells

res <- benchmark(apply(X = expr, MARGIN = 1, FUN = function(x,i1=ix1,i2=ix2) wilcox.test(x = x[i1],y = x[i2],alternative = "two.sided", paired = F,exact = F,correct = T)),
                 gficf:::rcpp_WMU_test(M = expr,idx1 = ix1,idx2 = ix2),
                 gficf:::rcpp_parallel_WMU_test(matX = expr[,ix1],matY = expr[,ix2],printOutput = F),
                 replications=3,
                 order="elapsed"
                 )

print(res[,c("test","elapsed","relative")])
```

