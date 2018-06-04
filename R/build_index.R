#' Build index for mapping using Salmon
#'
#' \code{build_index} for mapping reads using Salmon.  
#'
#' @param species name of the species. Only \code{'human'}, \code{'mouse'}, 
#' and \code{'rat'} are allowed to use.
#' @param kmer k-mer size for indexing. default is 31. See \code{'Salmon'} 
#' for details.
#' @param ens_release version of Ensembl release.
#' @param destdir directory where all the files will be saved.
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
#'
#' #Running this function will take some time.
#' \donttest{
#' build_index(species="human",kmer=31,
#' ens_release=92,destdir=tempdir())
#' }
#'
#' @export
#'
build_index<-function(species=c("human","mouse","rat"),kmer=31,
ens_release=92,destdir){

    species <- match.arg(species, c("human","mouse","rat"))

    if(species=="human"){
        system(paste0("wget -O ",destdir,"/Homo_sapiens.GRCh38.cdna.",
			"all.fa.gz ftp://ftp.ensembl.org/pub/release-",ens_release,
            "/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"))
        system(paste0("wget -O ",destdir,"/Homo_sapiens.GRCh38.ncrna.fa.gz",
			" ftp://ftp.ensembl.org/pub/release-",ens_release,
            "/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz"))
        system(paste0("gunzip -c ",destdir,"/Homo_sapiens.GRCh38.cdna.all.",
			"fa.gz ",destdir,"/Homo_sapiens.GRCh38.ncrna.fa.gz > ",destdir,
			"/Homo_sapiens.GRCh38.release",ens_release,".cdna.ncrna.fa"))
        system(paste0("salmon index -t ",destdir,"/Homo_sapiens.GRCh38.",
			"release",ens_release,".cdna.ncrna.fa -i ",destdir,
			"/human_transcripts_release",ens_release,"_index"))
    } else if(species=="mouse"){
        system(paste0("wget -O ",destdir,"/Mus_musculus.GRCm38.cdna.",
			"all.fa.gz ftp://ftp.ensembl.org/pub/release-",ens_release,
            "/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz"))
        system(paste0("wget -O ",destdir,"/Mus_musculus.GRCm38.ncrna.fa.gz",
			" ftp://ftp.ensembl.org/pub/release-",ens_release,
            "/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz"))
        system(paste0("gunzip -c ",destdir,"/Mus_musculus.GRCm38.cdna.all.",
			"fa.gz ",destdir,"/Mus_musculus.GRCm38.ncrna.fa.gz > ",destdir,
			"/Mus_musculus.GRCm38.","release",ens_release,".cdna.ncrna.fa"))
        system(paste0("salmon index -t ",destdir,"/Mus_musculus.GRCm38.",
			"release",ens_release,".cdna.ncrna.fa -i ",destdir,"/mouse_",
			"transcripts_","release",ens_release,"_index"))
    } else if(species=="rat"){
        system(paste0("wget -O ",destdir,"/Rattus_norvegicus.Rnor_6.0.",
			"cdna.all.fa.gz ftp://ftp.ensembl.org/pub/release-",ens_release,
            "/fasta/rattus_norvegicus/cdna/",
            "Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz"))
        system(paste0("wget -O ",destdir,"/Rattus_norvegicus.Rnor_6.0.ncrna.",
			"fa.gz ftp://ftp.ensembl.org/pub/release-",ens_release,
            "/fasta/rattus_norvegicus/ncrna/",
            "Rattus_norvegicus.Rnor_6.0.ncrna.fa.gz"))
        system(paste0("gunzip -c ",destdir,"/Rattus_norvegicus.Rnor_6.0.cdna.",
			"all.fa.gz ",destdir,"/Rattus_norvegicus.Rnor_6.0.ncrna.fa.gz > ",
            destdir,"/Rattus_norvegicus.Rnor_6.0.release",ens_release,
			".cdna.ncrna.fa"))
        system(paste0("salmon index -t ",destdir,"/Rattus_norvegicus.Rnor_6.0.",
			"release",ens_release,".cdna.ncrna.fa -i ",destdir,"/rat_",
			"transcripts_release",ens_release,"_index"))
    } else {
        warning("A valid species name is missing")
    }
}
