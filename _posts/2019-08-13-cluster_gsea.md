---
title: "New Functionality: Predict active pathways"
tags: gf-icf
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
We can now use GSEA (Gene Set Enrichement Analysis) to predict pathway activity across cells of the same cluster. Briefly gf-icf gene ranks are first summedd across cells of the same cluster and then GSEA is performed. This analysis can be now performed via `runGSEA` function. [See the example here](https://jeky82.github.io/gficf_example.html#how-to-perform-gsea-to-identify-active-pathways-in-each-cluster) about this new functionality. 
