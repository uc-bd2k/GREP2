#' Generate combined QC report for Salmon and FastQC
#'
#' \code{run_fastqc} generates a single HTML report from the fastQC
#' reports and salmon read mapping results using MultiQC.
#'
#' @param fastqc_dir directory where all the FastQC files are saved.
#' @param salmon_dir directory of the salmon files.
#' @param destdir directory where you want to save the combined QC report.
#'
#' @return HTML report.
#' 
#' @references 
#' 
#' Philip Ewels, Mans Magnusson, Sverker Lundin, and Max Kaller (2016):
#' MultiQC: summarize analysis results for multiple tools and samples 
#' in a single report. Bioinformatics, 32(19), 3047-3048.
#' \url{https://doi.org/10.1093/bioinformatics/btw354}
#'
#' @examples
#'
#' \dontrun{
#' run_multiqc(fastqc_dir=tempdir(),salmon_dir=tempdir(),
#' destdir=tempdir())
#' }
#'
#' @export 
run_multiqc <- function(fastqc_dir, salmon_dir, destdir) {
    cat(paste("Creating MultiQC report.\n",sep=""))
    system(paste0("multiqc ",fastqc_dir," ",salmon_dir," -o ",destdir))
    system(paste0("sed -i 's/A modular tool to aggregate results 
        from bioinformatics analyses across many samples into a single 
        report./''/g' ", destdir,"/","multiqc_report.html"))
}
