#' Download fastq files
#'
#' \code{get_fastq} downloads fastq files using SRA toolkit. 
#' We recommend using Aspera for fast downloading. You need to 
#' install Aspera for using \code{ascp} option. 
#'
#' @param srr_id SRA run accession ID.
#' @param library_layout layout of the library used. Either 
#' \code{'SINGLE'} or \code{'PAIRED'}.
#' @param use_sra_file logical, whether to use downloaded SRA
#' files to get fastq files or directly download
#' fastq files.
#' @param sra_files_dir directory where SRA files are saved.
#' If you use \code{use_sra_file=FALSE} then \code{sra_files_dir=NULL}.
#' @param n_thread number of cores to use.
#' @param destdir directory where all the results will be saved.
#'
#' @return A single fastq file will be generated for SINGLE end
#' reads and two files for PAIRED end reads.
#' 
#' @examples
#'
#' \donttest{
#' get_fastq(srr_id="SRR5890521",library_layout="SINGLE",
#' use_sra_file=FALSE,sra_files_dir=NULL,n_thread=2,
#' destdir=tempdir())
#' }
#'
#' @export 
get_fastq <- function(srr_id,library_layout=c("SINGLE","PAIRED"),
use_sra_file=FALSE,sra_files_dir=NULL,n_thread,destdir) {

    library_layout <- match.arg(library_layout,c("SINGLE","PAIRED"))
    if (library_layout=="SINGLE") {
        if(length(list.files(destdir,pattern=".fastq"))
            ==1) {
            warning("Fastq file exist. Processing next sample...")
        } else {
            cat("Converting sra to fastq...")
            if(use_sra_file){
                system (paste0("fastq-dump --outdir ",destdir,
                " --skip-technical  --readids --read-filter pass",
                " --dumpbase --split-spot --clip ",
                sra_files_dir,"/",srr_id,".sra"))
            } else {
                system (paste0("fastq-dump --outdir ",destdir,
				" --skip-technical  --readids --read-filter pass",
				" --dumpbase --split-spot --clip ",srr_id))
			}
        }
    } else {
        if(length(list.files(paste0(destdir,"/",srr_id),
            pattern=".fastq"))==2) {
            warning("Fastq file exist. Processing next sample...")
        } else {
            cat("Converting sra to fastq...")
            if(use_sra_file){
                system (paste0("fastq-dump --outdir ",destdir,
                " --skip-technical  --readids --read-filter pass",
                " --dumpbase --split-files --clip ",sra_files_dir,"/",
                srr_id,".sra"))
            } else {
                system (paste0("fastq-dump --outdir ",destdir,
                " --skip-technical  --readids --read-filter pass",
                " --dumpbase --split-files --clip ",srr_id))
            }
        }
    }
    
    n_fastq <- if(library_layout=="PAIRED") {2} else {1}
    fastq_dumped <- length(list.files(destdir,
        pattern="\\.fastq$",recursive=FALSE,full.names=FALSE))
    if(n_fastq!=fastq_dumped){
        warning("Incomplete fastq download...")
    } else {
        cat(paste("All fastq files are generated successfully. ",
        Sys.time(),"\n",sep=""))
    }
}
