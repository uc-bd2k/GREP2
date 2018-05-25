#' Download SRA run files
#'
#' \code{get_srr} downloads SRA files using Aspera or FTP. 
#' We recommend using Aspera for fast downloading. 
#' You need to install Aspera for using \code{ascp} option. 
#'
#' @param srr_id SRA run accession ID.
#' @param destdir directory where all the results will be saved.
#' @param ascp logical, whether to use Aspera for downloading SRA files.
#' @param prefetch_workspace directory where SRA run files will be downloaded.
#' This parameter is needed if \code{ascp=TRUE}. 
#' The location of this directory can be found by going to the aspera 
#' directory (/.aspera/connect/bin/) and typing \code{'vdb-config -i'}.
#' A new window will pop-up and under the \code{'Workspace Name'}, 
#' you will find the location. Usually the default is 
#' \code{'/home/username/ncbi/public'}.
#' @param ascp_path path to the Aspera software. 
#'
#' @return SRA run accession file with extension ".sra". If you use 
#' \code{ascp=TRUE}, then downloaded files will be saved under 
#' \code{'/prefetch_workspace/sra'} directory. 
#' If \code{ascp=FALSE}, then files will be saved in the 
#' \code{'destdir'}
#' 
#' @examples
#'
#' \donttest{
#' get_srr(srr_id="SRR5890521", destdir=tempdir(),ascp=FALSE,
#' prefetch_workspace=NULL, ascp_path=NULL)
#' }
#'
#' @export
get_srr <- function(srr_id,destdir,ascp,prefetch_workspace,ascp_path){

    if(ascp){
        system(paste0("prefetch -X 500G --ascp-path '",
            ascp_path,"/connect/bin/ascp|",ascp_path,
            "/connect/etc/asperaweb_id_dsa.openssh' ",
            "--ascp-options '-k 1 -QT -l 400m' ", srr_id))
    } else {
        #system(paste0("mkdir ",destdir,"/",srr_id))
        system(paste0("wget -O ", destdir,"/",
        srr_id,".sra ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/",
        "ByRun/sra/", substr(srr_id, 1,3),"/", substr(srr_id, 1,6),"/",
        srr_id, "/",srr_id,".sra"))
    }

    # Check if the run file is not downloaded for any technical reason:
    if(ascp){
        if(file.exists(paste0(prefetch_workspace,"/sra","/",srr_id,".sra"))){
            cat(paste("All SRA files are downloaded successfully. ",
                Sys.time(),"\n",sep=""))
        } else {
            repeat {
                system(paste0("prefetch -X 500G --ascp-path '",
                    ascp_path,"/connect/bin/ascp|",ascp_path,
                    "/connect/etc/asperaweb_id_dsa.openssh' --ascp-options",
                    " '-k 1 -QT -l 400m' ", srr_id))
                if (file.exists(paste0(prefetch_workspace,"/sra","/",
                    srr_id,".sra"))) break
                cat(paste("All SRA files are downloaded successfully."
                    ,Sys.time(),"\n",sep=""))
            }
        }
    } else {
        if(file.exists(paste0(destdir,"/",srr_id,".sra"))){
            cat(paste("All SRA files are downloaded successfully. ",
                Sys.time(),"\n",sep=""))
        } else {
            repeat {
                system(paste0("wget -O ", destdir,"/",srr_id,
                    ".sra ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/",
                    "reads/ByRun/sra/", substr(srr_id, 1,3),"/",
                    substr(srr_id, 1,6),"/",srr_id, "/",srr_id,".sra"))
                if (file.exists(paste0(destdir,"/",srr_id,
                    ".sra"))) break
                cat(paste("All SRA files are downloaded successfully. ",
                Sys.time(),"\n",sep=""))
            }
        }
    }
}
