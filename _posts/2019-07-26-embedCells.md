---
title: "New Functionality: Emmbed new cells"
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
You can now use an existing embedding to add new cells via `embedNewCells`. Briefly new cells are first normalized with GF-ICF method but using as ICF weigth estimated on the training set and than projected in the existing PCA/LSA space before to be embedded in the already existing UMAP space via `uwot::umap_transform`. [See the example here](https://jeky82.github.io/gficf_example.html#how-to-embedd-new-cells-in-an-existing-space) about this new functionality. 
