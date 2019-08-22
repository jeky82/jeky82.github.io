---
title: "New Functionality: Find Marker Genes"
tags: gf-icf marker-genes DE-genes
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
GFICF package try to identify marker genes across clusters performing [Mann-Whitney U test with continuity correction](https://jeky82.github.io/2019/08/20/MannWhitney.html){:target="_blank"}. Briefly differentially expressed (DE) genes of each cluster are identified comparing the expression of each gene in each cluster versus the all the others. [See the example here](https://jeky82.github.io/gficf_example.html#find-marker-genes) about this new functionality. 
