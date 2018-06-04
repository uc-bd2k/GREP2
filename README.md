# GREP2 : GEO RNA-seq Experiments Processing Pipeline

The Gene Expression Omnibus ([GEO](https://www.ncbi.nlm.nih.gov/geo/)) is a public repository of gene expression data 
that hosts more than 6,000 RNA-seq datasets and this number is increasing. Most of these datasets are deposited in raw sequencing 
format which needs to be downloaded and processed. With an aim to transform all these datasets in an analysis-ready format, 
we have developed a comprehensive pipeline to simultaneously download and process RNA-seq data sets from GEO. 
This R-based automated pipeline can process the available RNA-seq data of human, mouse, and rat from GEO. This package is 
recommended to use in the unix environment as many of the features are not available in windows.

---
## Installation

Before installing `GREP2`, you need to install the following software packages first:

1. [SRA toolkit](http://www.sthda.com/english/wiki/install-sra-toolkit)
2. [Aspera-connect](http://download.asperasoft.com/download/docs/connect/2.3/aspera-connect-linux.html#installation)
3. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/seqmonk/INSTALL.txt)
4. [Salmon](http://salmon.readthedocs.io/en/latest/building.html)
5. [MultiQC](https://github.com/ewels/MultiQC/blob/master/docs/installation.md)

You will also need to install the following R and Bioconductor packages:
```
install.packages(c("devtools", "XML", "parallel", "utils", "rentrez", "RCurl")
source("https://bioconductor.org/biocLite.R")
biocLite(c("GEOquery", "Biobase", "tximport", "EnsDb.Hsapiens.v86", "EnsDb.Rnorvegicus.v79", "EnsDb.Mmusculus.v79",
    "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db"))
``` 

Once you install the above packages, you can now install `GREP2` using `devtools`:

```
library(devtools)
install_github("uc-bd2k/GREP2")
```

