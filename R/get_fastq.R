#' Download fastq files
#'
#' \code{get_fastq} downloads fastq files using SRA toolkit. 
#' We recommend using Aspera for fast downloading. You need to 
#' install Aspera(\url{http://www.asperasoft.com/}) for using 
#' \code{ascp} option. 
#'
#' @param srr_id SRA run accession ID.
#' @param library_layout layout of the library used. Either 
#' \code{'SINGLE'} or \code{'PAIRED'}.
#' @param get_sra_file logical, whether to download SRA file
#' first and get fastq files afterwards or directly download
#' fastq files.
#' @param sra_files_dir directory where SRA files are saved.
#' If you use \code{get_sra_file=FALSE} then \code{sra_files_dir=NULL}.
#' @param n_thread number of cores to use.
#' @param destdir directory where all the results will be saved.
#'
#' @return A single fastq file will be generated for SINGLE end
#' reads and two files for PAIRED end reads.
#' 
#' @examples
#' srr_id="SRR6324192"
#' \dontrun{
#' get_fastq(srr_id=srr_id, library_layout="SINGLE",
#' get_sra_file=FALSE,
#' sra_files_dir=NULL, n_thread=2,
#' destdir="/mnt/raid/test")
#' }
#'
#' @export 
get_fastq <- function(srr_id, library_layout=c("SINGLE","PAIRED"),
get_sra_file=FALSE,sra_files_dir=NULL,n_thread,destdir) {

    if(!dir.exists(paste0(destdir,"/",srr_id))){
        system(paste0("mkdir ",destdir,"/",srr_id))
    }
    library_layout <- match.arg(library_layout, c("SINGLE","PAIRED"))

    if (library_layout=="SINGLE") {
        if(length(list.files(paste0(destdir,"/",srr_id), pattern=".fastq"))
            ==1) {
            warning("Fastq file exist. Processing next sample...")
        } else {
            cat("Converting sra to fastq...")
            if(get_sra_file){
                system (paste0("fastq-dump --outdir ",destdir,"/",srr_id,
                " --skip-technical  --readids --read-filter pass 
                --dumpbase --split-spot --clip ",
                sra_files_dir,"/",srr_id,".sra"))
            } else {
                system (paste0("fastq-dump --outdir ",destdir,"/",srr_id,
                " --skip-technical  --readids --read-filter pass 
                --dumpbase --split-spot --clip ",srr_id))
            }
        }
    } else {
        if(length(list.files(paste0(destdir,"/",srr_id),
            pattern=".fastq"))==2) {
            warning("Fastq file exist. Processing next sample...")
        } else {
            cat("Converting sra to fastq...")
            if(get_sra_file){
                system (paste0("fastq-dump --outdir ",destdir,"/",srr_id,
                " --skip-technical  --readids --read-filter pass 
                --dumpbase --split-files --clip ", sra_files_dir,"/",
                srr_id,".sra"))
            } else {
                system (paste0("fastq-dump --outdir ",destdir,"/",srr_id,
                " --skip-technical  --readids --read-filter pass 
                --dumpbase --split-files --clip ",srr_id))
            }
        }
    }
    
    n_fastq <- if(library_layout=="PAIRED") {2} else {1}
    fastq_dumped <- length(list.files(paste0(destdir,"/",srr_id),
        pattern = "\\.fastq$",recursive=TRUE,full.names=FALSE))
    if(n_fastq!=fastq_dumped){
        warning("Incomplete fastq download...")
    } else {
        cat(paste("All fastq files are generated successfully. ",
        Sys.time(),"\n",sep=""))
    }
}
