#' Build index for mapping using Salmon
#'
#' \code{build_index} for mapping reads using Salmon.  
#'
#' @param species name of the species. Only \code{'human'}, \code{'mouse'}, 
#' and \code{'rat'} are allowed to use.
#' @param kmer k-mer size for indexing. default is 31. See \code{'Salmon'} 
#' for details.
#' @param ens_release version of Ensembl release.
#'
#' @return directory of index files
#' 
#' @references 
#' 
#' Rob Patro, Geet Duggal, Michael I. Love, Rafael A. Irizarry, and 
#' Carl Kingsford (2017): Salmon provides fast and bias-aware 
# quantification of transcript expression. Nature methods, 14(4),417.
#' \url{https://www.nature.com/articles/nmeth.4197}
#'
#' @examples
#' build_index(species="rat", kmer=31,
#' ens_release=92)
#'
#' @export
build_index <- function(species=c("human","mouse","rat"),kmer=31,
    ens_release=92){
    
    species <- match.arg(species, c("human","mouse","rat"))

    if(species=="human"){
        system(paste0("wget ftp://ftp.ensembl.org/pub/release-",ens_release,
            "/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"))
        system(paste0("wget ftp://ftp.ensembl.org/pub/release-",ens_release,
            "/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz"))
        system(paste0("gunzip -c Homo_sapiens.GRCh38.cdna.all.fa.gz 
            Homo_sapiens.GRCh38.ncrna.fa.gz > Homo_sapiens.GRCh38.release",
            ens_release,".cdna.ncrna.fa"))
        system(paste0("salmon index -t Homo_sapiens.GRCh38.release",
            ens_release,".cdna.ncrna.fa -i human_transcripts_release",
            ens_release,"_index"))
    }
    if(species=="mouse"){
        system(paste0("wget ftp://ftp.ensembl.org/pub/release-",ens_release,
            "/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz"))
        system(paste0("wget ftp://ftp.ensembl.org/pub/release-",ens_release,
            "/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz"))
        system(paste0("gunzip -c Mus_musculus.GRCm38.cdna.all.fa.gz 
            Mus_musculus.GRCm38.ncrna.fa.gz > Mus_musculus.GRCm38.release",
            ens_release,".cdna.ncrna.fa"))
        system(paste0("salmon index -t Mus_musculus.GRCm38.release",
            ens_release,".cdna.ncrna.fa -i mouse_transcripts_release",
            ens_release,"_index"))
    }
    if(species=="rat"){
        system(paste0("wget ftp://ftp.ensembl.org/pub/release-",ens_release,
            "/fasta/rattus_norvegicus/cdna/
            Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz"))
        system(paste0("wget ftp://ftp.ensembl.org/pub/release-",ens_release,
            "/fasta/rattus_norvegicus/ncrna/
            Rattus_norvegicus.Rnor_6.0.ncrna.fa.gz"))
        system(paste0("gunzip -c Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz 
            Rattus_norvegicus.Rnor_6.0.ncrna.fa.gz > 
            Rattus_norvegicus.Rnor_6.0.release",ens_release,".cdna.ncrna.fa"))
        system(paste0("salmon index -t Rattus_norvegicus.Rnor_6.0.release",
            ens_release,".cdna.ncrna.fa -i rat_transcripts_release",
            ens_release,"_index"))
    }
}
