#' Trim fastq files using Trimmomatic
#'
#' \code{trim_fastq} trim fastq files based on the illumina instruments
#' using Trimmomatic.
#'
#' The following parameters are used as default in the trimmoatic function:
#' \enumerate{
#' \item Remove leading low quality or N bases (below quality 3) 
#' (LEADING:3)
#' \item Remove trailing low quality or N bases (below quality 3) 
#' (TRAILING:3)
#' \item Scan the read with a 4-base wide sliding window, cutting when
#' the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
#' \item Drop reads below the 36 bases long (MINLEN:36)
#' }
#'
#' @param srr_id SRA run accession ID.
#' @param fastq_dir directory of the fastq files.
#' @param instrument name of the illumina sequencing platform.
#' For example, \code{'HiSeq'}. 
#' @param trimmomatic_path path to the Trimmomatic software.
#' @param library_layout layout of the library used. Either \code{'SINGLE'}
#' or \code{'PAIRED'}.
#' @param destdir directory where the trimmed fastq files will be saved.
#' @param n_thread number of cores.
#'
#' @return trimmed fastq files. 
#' 
#' @references 
#' 
#' Anthony M. Bolger, Marc Lohse, and Bjoern Usadel (2014):
#' Trimmomatic: a flexible trimmer for Illumina sequence data.
#' Bioinformatics, 30(15), 2114-2120.
#' \url{https://doi.org/10.1093/bioinformatics/btu170}
#' 
#' @examples
#'
#' \donttest{
#' fastq_dir=system.file("extdata","", package="GREP2")
#' trimmomatic_path=system.file("java","trimmomatic-0.36.jar", package="GREP2")
#' trim_fastq(srr_id="SRR5890521",fastq_dir=fastq_dir,
#' instrument="MiSeq",trimmomatic_path=trimmomatic_path,
#' library_layout="SINGLE",destdir=tempdir(),n_thread=2)
#' }
#'
#' @export
trim_fastq <- function(srr_id,fastq_dir,instrument,trimmomatic_path,
library_layout=c("SINGLE","PAIRED"),destdir,n_thread){

    adapters<- function(instrument) {
        if (grepl("HiSeq|MiSeq",instrument)) {	
            if (library_layout=="SINGLE") {
                return(paste0(system.file("extdata","TruSeq3-SE.fa", package="GREP2")))
            } else{
                return(paste0(system.file("extdata","TruSeq3-PE.fa", package="GREP2")))
            }
        }
        else if (grepl("GA|Genome Analyzer",instrument)) {
            if (library_layout=="SINGLE") {
                return(paste0(system.file("extdata","TruSeq2-SE.fa", package="GREP2")))
            } else{
                return(paste0(system.file("extdata","TruSeq2-PE.fa", package="GREP2")))
            }
        }
        else {
            if (library_layout=="SINGLE") {
                return(paste0(system.file("extdata","TruSeq3-SE.fa", package="GREP2")))
            } else{
                return(paste0(system.file("extdata","TruSeq3-PE.fa", package="GREP2")))
            }
        }
    }
    library_layout <- match.arg(library_layout,c("SINGLE","PAIRED"))
    if (library_layout=="SINGLE") {
        fq=paste(fastq_dir,"/",srr_id,"*.fastq",sep="")
        fq_se=paste(destdir,"/",srr_id,"_trimmed.fastq",sep="")
        system(paste0("java -jar ",trimmomatic_path,
        " SE -phred33 -threads ",n_thread," ",
        fq," ",fq_se," ","ILLUMINACLIP:",
        adapters(instrument),
        ":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 ",
        "MINLEN:36"))
    } else {
        fq1=paste(fastq_dir,"/",srr_id,"*_1.fastq",sep="")
        fq2=paste(fastq_dir,"/",srr_id,"*_2.fastq",sep="")
        fq1_paired=paste(destdir,"/",srr_id,"_1_trimmed.fastq",sep="")
        fq2_paired=paste(destdir,"/",srr_id,"_2_trimmed.fastq",sep="")
        fq1_unpaired=paste(destdir,"/",srr_id,"_1_unpaired.fastq",sep="")
        fq2_unpaired=paste(destdir,"/",srr_id,"_2_unpaired.fastq",sep="")

        system(paste0("java -jar ",trimmomatic_path,
        " PE -phred33 -threads ",n_thread," ",
        fq1," ",fq2," ",fq1_paired," ",fq1_unpaired," ",
        fq2_paired," ",fq2_unpaired," ","ILLUMINACLIP:",
        adapters(instrument),
        ":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 ",
        "MINLEN:36"))
    }
}
