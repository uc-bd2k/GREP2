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
#' destdir="/mnt/raid/test"
#' fastq_dir="/mnt/raid/test/GSE107363/SRR6324192"
#' \dontrun{
#' run_fastqc(destdir=destdir, fastq_dir=fastq_dir, n_thread=2)
#' }
#'
#' @export
run_fastqc <- function(destdir, fastq_dir, n_thread ) {
	
	cat(paste("Running FastQC... ",Sys.time(),"\n",sep=""))
	#setwd(destdir)
	fastq_files = list.files(fastq_dir, pattern=".fastq$", full.names=TRUE)
	
	if(!dir.exists("fastqc")){
		system(paste0("mkdir ",destdir,"/fastqc"))
	}
	#setwd(paste0(destdir,"/fastqc/"))
	system(paste0("fastqc -o ",destdir,"/fastqc/ --threads ",n_thread," ", fastq_files))
}