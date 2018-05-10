#' Download metadata from GEO and SRA
#'
#' @param geo_series_acc GEO series accession ID.
#'
#' @return a list of GEO and SRA metadata. 
#' 
#' @examples
#' geo_series_acc="GSE107363"
#' \dontrun{
#' get_metadata(geo_series_acc=geo_series_acc)
#' }
#'
#' @importFrom rentrez entrez_search
#' @importFrom XML xmlRoot xmlValue xmlTreeParse getNodeSet
#' @importFrom RCurl getURL
#' @importFrom GEOquery getGEO
#' @importFrom Biobase phenoData pData 
#' @importFrom utils read.csv
#'
#' @export 
get_metadata <- function(geo_series_acc) {

    options(warn=-1)
    geo_id <- rentrez::entrez_search(db="gds",term=geo_series_acc)$ids[[1]]
    geo_summary <- XML::xmlRoot(XML::xmlTreeParse(RCurl::getURL(
        paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
        esummary.fcgi?db=gds&id=", geo_id))))
    sra_study_acc <- XML::xmlValue(XML::getNodeSet(geo_summary[[1]],
        "//DocSum//Item//Item//Item[@Name='TargetObject']")[1][[1]])

    # GEO metadata
    x <- GEOquery::getGEO(geo_series_acc, AnnotGPL=FALSE,GSEMatrix=FALSE,
        destdir = getwd() ,getGPL=FALSE)
    nm <- names(x@gpls)
    stub=gsub("\\d{1,3}$", "nnn", geo_series_acc, perl = TRUE)
    gseurl <- "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s"

    if(length(nm)>1){
        for (i in 1: length(nm)){
            myurl <- sprintf(gseurl, stub, geo_series_acc, 
                paste0(geo_series_acc,"-",nm[i],"_series_matrix.txt.gz"))
            system(paste0("wget ", myurl))
        }
        geo_list_prim <- lapply(lapply(nm, function(x) 
            GEOquery::getGEO(geo_series_acc, filename = 
            paste0("./",geo_series_acc,"-",x ,"_series_matrix.txt.gz"), 
            GSEMatrix=TRUE,destdir = getwd() ,getGPL=FALSE)
            ),function(x) Biobase::pData(phenoData(x)))
        geo_list_sec <- lapply(geo_list_prim, function(x) 
            x[,Reduce(intersect, lapply(geo_list_prim, function(x) 
            colnames(x)))])
        geo_df <- do.call(rbind, geo_list_sec)
        closeAllConnections()
        metadata_geo <- data.frame(lapply(geo_df, as.character), 
            stringsAsFactors=FALSE)
    } else {
        myurl <- sprintf(gseurl, stub, geo_series_acc, 
            paste0(geo_series_acc,"_series_matrix.txt.gz"))
        system(paste0("wget ", myurl))
        geo_df <- Biobase::pData(phenoData(getGEO(geo_series_acc, 
            filename = paste0("./",geo_series_acc, "_series_matrix.txt.gz"), 
            GSEMatrix=TRUE,destdir = getwd() ,getGPL=FALSE)))
        metadata_geo <- data.frame(lapply(geo_df, as.character), 
            stringsAsFactors=FALSE)
    }

    # SRA metadata
    system(paste("wget -O ",geo_series_acc,
    "_metadata.csv 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=
    efetch&db=sra&rettype=runinfo&term=",sra_study_acc,"'", sep=""))
    metadata_sra <- data.frame(lapply(utils::read.csv(file=
        paste(geo_series_acc,"_metadata.csv", sep=""), header=TRUE), 
        as.character), stringsAsFactors=FALSE)

    #system(paste("rm -rf *matrix.txt.gz *_metadata.csv *.gz"))
    metadata <- list(metadata_geo=metadata_geo, metadata_sra=metadata_sra)
    return(metadata)
    options(warn=0)
}
