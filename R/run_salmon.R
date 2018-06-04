#' Quantify transcript abundances using Salmon
#'
#' \code{run_salmon} is a wrapper function for mapping reads to quantify
#' transcript abundances using Salmon. You need to install Salmon and
#' build index to run this function. 
#' For index building see function \code{'build_index'}. 
#'
#' \code{run_salmon} We use default options of Salmon. This function
#' works for a single sample. You can use this function in a loop
#' for multiple samples. 
#' For other options from Salmon use \code{'other_opts'}.
#'
#' @param srr_id SRA run accession ID.
#' @param library_layout layout of the library used. Either \code{'SINGLE'}
#' or \code{'PAIRED'}.
#' @param index_dir directory of the indexing files needed for read 
#' mapping using Salmon. See function \code{'build_index'}.
#' @param destdir directory where all the results will be saved.
#' @param fastq_dir directory of the fastq files.
#' @param use_trimmed_fastq logical, whether to use trimmed fastq files. 
#' @param other_opts other options to use. See Salmon documentation for
#' the available options.
#' @param n_thread number of cores to use.
#'
#' @return The following items will be returned and saved in the salmon
#' directory:
#' \enumerate{
#' \item quant_new.sf: plain-text, tab-separated quantification file
#' that contains 5 column: Name,Length,EffectiveLength,TPM, and NumReads.
#' \item cmd_info.json: A JSON format file that records the main command
#' line parameters with which Salmon was invoked for the run that produced
#' the output in this directory.
#' \item aux_info: This directory will have a number of files (and
#' subfolders) depending on how salmon was invoked.
#' \item meta_info.json: A JSON file that contains meta information about
#' the run, including stats such as the number of observed and mapped
#' fragments, details of the bias modeling etc. 
#' \item ambig_info.tsv: This file contains information about the number
#' of uniquely-mapping reads as well as the total number of
#' ambiguously-mapping reads for each transcript. 
#' \item lib_format_counts.json: This JSON file reports the number of
#' fragments that had at least one mapping compatible with the designated
#' library format, as well as the number that didn't. 
#' \item libParams: The auxiliary directory will contain a text file
#' called flenDist.txt. This file contains an approximation of the
#' observed fragment length distribution.
#' }
#' 
#' @references 
#' 
#' Rob Patro, Geet Duggal, Michael I. Love, Rafael A. Irizarry, and
#' Carl Kingsford (2017): Salmon provides fast and bias-aware
#' quantification of transcript expression. Nature methods, 14(4), 417.
#' \url{https://www.nature.com/articles/nmeth.4197}
#'
#' @examples
#'
#' #You will have to build index first to run this function
#' fastq_dir=system.file("extdata","", package="GREP2")
#' \donttest{
#' build_index(species="human",kmer=31,ens_release=92,
#' destdir=tempdir())
#' run_salmon(srr_id="SRR5890521",library_layout="SINGLE",
#' index_dir=tempdir(),destdir=tempdir(),
#' fastq_dir=fastq_dir,use_trimmed_fastq=FALSE,
#' other_opts=NULL,n_thread=2)
#' }
#'
#' @export 
run_salmon <- function(srr_id, library_layout=c("SINGLE","PAIRED"),
index_dir, destdir, fastq_dir, use_trimmed_fastq=FALSE,
other_opts=NULL, n_thread ) {

    if(!dir.exists(paste0(destdir,"/salmon"))){
        system(paste0("mkdir ",destdir,"/salmon"))
    }
    library_layout <- match.arg(library_layout, c("SINGLE","PAIRED"))

    if (library_layout=="SINGLE") {
        if(use_trimmed_fastq){
            system(paste0("salmon quant -i ",index_dir," -p ",n_thread,
            " ",other_opts, " -l A -r ",fastq_dir,"/",srr_id, 
            "_trimmed.fastq -o ",destdir,"/salmon/",srr_id,
            "_transcripts_quant"))
        } else {
            system(paste0("salmon quant -i ",index_dir," -p ",n_thread,
            " ",other_opts, " -l A -r ",fastq_dir,"/",srr_id, 
            "_pass.fastq -o ",destdir,"/salmon/",srr_id,
            "_transcripts_quant"))
        }
    } else {
        if(use_trimmed_fastq){
            system(paste0("salmon quant -i ",index_dir," -p ",n_thread,
            " ",other_opts, " -l A -1 ",fastq_dir,"/",srr_id,
            "_trimmed_1.fastq ", "-2 ",fastq_dir,"/",srr_id,
            "_trimmed_2.fastq -o ",destdir,"/salmon/",srr_id,
            "_transcripts_quant"))
        } else {
            system(paste0("salmon quant -i ",index_dir," -p ",n_thread,
            " ",other_opts, " -l A -1 ",fastq_dir,"/",srr_id,
            "_pass_1.fastq ", "-2 ",fastq_dir,"/",srr_id,
            "_pass_2.fastq -o ",destdir,"/salmon/",srr_id,
            "_transcripts_quant"))
        }
    }
    if (file.exists(paste0(destdir,"/salmon/",srr_id,
        "_transcripts_quant/quant.sf"))) {
        system(paste("cat ",destdir,"/salmon/",srr_id,
        "_transcripts_quant/quant.sf","| sed -E 's/\\.[0-9]+//' > ",
        destdir,"/salmon/",srr_id,"_transcripts_quant","/",srr_id,
        "_quant_new.sf", sep=""))
    } else {
        cat("quant.sf doesn't exist. Processing next sample.")
    }
}
