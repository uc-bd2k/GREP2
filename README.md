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

Once you install the above packages, you can now install `GREP2` using `CRAN`:
```
install.packages("GREP2")
```

You can also install `GREP2` using `devtools`:
```
library(devtools)
install_github("uc-bd2k/GREP2")
```



## GREP2 pipeline workflow

To consistently process GEO RNA-seq datasets through a robust and uniform 
system, we have built GEO RNA-seq evenly processing pipeline (GREP2). 
To demonstrate the usage of the package, we demonstrate the processing steps
with a small dataset from GEO:
[GSE102170](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102170).
The whole processing workflow can be summarized in the following steps:

1. The pipeline starts with a valid GEO series accession ID. Currently 
the pipeline works for human, mouse, and rat species. We then retrieve 
metadata for the GEO series accession using Bioconductor package 
GEOquery. We also download metadata file from 
the sequence read archive (SRA) to get corresponding run information.
```
library(GREP2)
metadata <- get_metadata(geo_series_acc="GSE102170",destdir=tempdir(),
geo_only=FALSE,download_method="auto")
```

2. Download corresponding experiment run files from the SRA using `ascp` 
utility of [Aspera Connect](http://download.asperasoft.com/download/docs/connect/2.3/aspera-connect-linux.html#installation) or 
regular download. All the downloaded files are stored in the local repository until processed. You can skip this
step by downloading fastq files directly.
```
srr_id <- metadata$metadata_sra$Run
for(i in 1:length(srr_id)){
	get_srr(srr_id=srr_id[i], destdir=tempdir(), ascp=FALSE,
	prefetch_workspace=NULL,ascp_path=NULL)
}
```

3. Convert SRA files to fastq format using [NCBI SRA toolkit](http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software) 
or download fastq files directly. 
```
library_layout <- metadata$metadata_sra$LibraryLayout
for(i in 1:length(srr_id)){
	get_fastq(srr_id=srr_id[i],library_layout=library_layout[i],
	use_sra_file=FALSE,sra_files_dir=NULL,n_thread=2,
	destdir=tempdir())
}
```

4. Run [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
on each fastq file to generate quality control (QC) reports.
```
run_fastqc(destdir=tempdir(),fastq_dir=tempdir(),
n_thread=2)
```

5. Remove adapter sequences if necessary using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).
```
for(i in 1:length(srr_id)){
	trim_fastq(srr_id=srr_id[i],fastq_dir=tempdir(),
	instrument="MiSeq",library_layout=library_layout[i],destdir=tempdir(),n_thread=2)
}
```

6. Quantify transcript abundances using [Salmon](http://salmon.readthedocs.io/en/latest/building.html). 
Transcript level estimates are then summarized to gene level using 
Bioconductor package [tximport](http://bioconductor.org/packages/release/bioc/html/tximport.html).
We obtained gene annotation for Homo sapiens (GRCh38), 
Mus musculus (GRCm38), and Rattus norvegicus (Rnor_6.0)
from Ensemble.
```
# Before running Salmon, you will have to build index first.
build_index(species="human",kmer=31,ens_release=92,
destdir=tempdir())
# Run Salmon
for(i in 1:length(srr_id)){
	run_salmon(srr_id=srr_id[i],library_layout=library_layout[i],
	index_dir=tempdir(),destdir=tempdir(),
	fastq_dir=tempdir(),use_trimmed_fastq=FALSE,
	other_opts=NULL,n_thread=2)
}
# Run tximport
counts_list <- run_tximport(srr_id=srr_id, species="human",
salmon_dir=paste0(tempdir(),"/salmon"),countsFromAbundance="lengthScaledTPM")
```

7. Compile FastQC reports and Salmon log files into a single
interactive HTML report using [MultiQC](http://multiqc.info/).  
```
run_multiqc(fastqc_dir=tempdir(),salmon_dir=tempdir(),
destdir=tempdir())
```

---
You can run the above individual functions for each step or run the 
whole pipeline using the following `process_geo_rnaseq` function. 
All of the above steps are combined into the following single function. 
We would recommend using this function for processing GEO RNA-seq data.
```
process_geo_rnaseq (geo_series_acc=geo_series_acc,destdir=tempdir(),
download_method="auto",
ascp=FALSE,prefetch_workspace=NULL,
ascp_path=NULL,use_sra_file=FALSE,trim_fastq=FALSE,
index_dir=tempdir(),species="human",
countsFromAbundance="lengthScaledTPM",n_thread=1)
```



