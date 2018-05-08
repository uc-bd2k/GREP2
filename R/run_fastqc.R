#' QC report for each fastq files using FastQC
#'
#' \code{run_fastqc} HTML report of each fastq files using FastQC. You need to install FastQC from \url{https://www.bioinformatics.babraham.ac.uk/projects/fastqc/} 
#'
#' @param destdir directory where all the results will be saved.
#' @param fastq_dir directory of the fastq files.
#' @param n_thread number of cores to use.
#'
#' @return HTML report of the fastq files under fastqc directory.
#' 
#' @examples
#'
#' run_fastqc(destdir="/home", fastq_dir="/home/SRR6324192/", n_thread=2)
#'
#' @export
run_fastqc <- function(destdir, fastq_dir, n_thread ) {
	
	cat(paste("Running FastQC... ",Sys.time(),"\n",sep=""))
	setwd(destdir)
	fastq_files = list.files(fastq_dir, pattern=".fastq$", full=TRUE)
	
	if(!dir.exists("fastqc")){
		system(paste0("mkdir ",destdir,"/fastqc"))
	}
	setwd(paste0(destdir,"/fastqc/"))
	system(paste0("fastqc -o ",getwd()," --threads ",n_thread," ", fastq_files))
}