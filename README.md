# GREP2 : GEO RNA-seq Experiments Processing Pipeline

The Gene Expression Omnibus ([GEO](https://www.ncbi.nlm.nih.gov/geo/)) is a public repository of gene expression data 
that hosts more than 6,000 RNA-Seq datasets and this number is increasing. Most of these datasets are deposited in raw sequencing 
format which needs to be downloaded and processed. With an aim to transform all these datasets in an analysis-ready format, 
we have developed a comprehensive pipeline to simultaneously download and process RNA-Seq data sets from GEO. 
This R-based automated pipeline can process the available RNA-Seq data of human, mouse, and rat from GEO.

---
## Installation

You can install `GREP2` using `devtools`:

```
library(devtools)
install_github("uc-bd2k/GREP2")
```
