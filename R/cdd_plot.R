#' Statistical sequence length
#' 
#' @title fastaleng
#' @param fasta_file The path of protein fasta file.
#' @export
#' @author Shiqi Zhao
#' @return data.frame
#' @examples
#' fasta_path <- system.file("extdata", "idpep.fa", package = "BioVizSeq") 
#' fastaleng(fasta_path)
#' 
#' @importFrom seqinr read.fasta

fastaleng <- function(fasta_file) {
  seqs <- seqinr::read.fasta(fasta_file, seqtype = "AA", as.string = TRUE, whole.header = TRUE)
  
  ids <- names(seqs)
  lengths <- numeric()
  for(i in seq_along(seqs)){
    lengths[i] <- nchar(seqs[[i]][1])
  }
  
  result_df <- data.frame(ID = ids, length = lengths)
  
  return(result_df)
}

#' Extract the location information of domain from cdd file 
#' 
#' @title cdd_to_loc
#' @param cdd_file CDD file.
#' @export
#' @author Shiqi Zhao
#' @return data.frame
#' @examples
#' hitdata_path <- system.file("extdata", "hitdata.txt", package = "BioVizSeq")
#' cdd_file <- readLines(hitdata_path)
#' domain_loc <- cdd_to_loc(cdd_file)
cdd_to_loc <- function(cdd_file){
  data <- cdd_file[!grepl("^#", cdd_file)]
  data <- data[data != ""]
  data <- data[!grepl("^Query", data)]
  data <- as.data.frame(data)
  data <- separate(data, 1, into = c("ID", "Hit_type", "PSSM", "start","end","evalue","bitscore","accession","Superfamily","Incomplete","i"), sep='\t')
  domain_loc <- data[c(1,9,4,5)]
  domain_loc[,1] <- gsub("Q#.* - >", "", domain_loc[,1])
  domain_loc[,2] <- gsub(" superfamily", "", domain_loc[,2])
  domain_loc[,3] <- as.numeric(domain_loc[,3])
  domain_loc[,4] <- as.numeric(domain_loc[,4])
  return(domain_loc)
}

#' Visualization of domain in CDD file
#' 
#' @title cdd_plot
#' @param cdd_file The path of cdd file.
#' @param fasta_file The path of fasta file.
#' @param the_order The path of order file. A List of Gene ID , One ID Per Line.
#' @param domain_select The domain ID which you want to align with.
#' @param shape RoundRect or Rect.
#' @param r The radius of rounded corners.
#' @param legend_size The size of legend.
#' @param domain_color The color set of domain.
#' @export
#' @author Shiqi Zhao
#' @return p
#' @examples
#' hitdata_path <- system.file("extdata", "hitdata.txt", package = "BioVizSeq")
#' fa_path <- system.file("extdata", "idpep.fa", package = "BioVizSeq")
#' cdd_plot(hitdata_path, fa_path)
#' 
#' order_path <- system.file("extdata", "order.csv", package = "BioVizSeq")
#' cdd_plot(hitdata_path, fa_path, the_order = order_path)

cdd_plot <- function(cdd_file, fasta_file, the_order = NULL, 
                     domain_select = NULL, shape = "RoundRect", 
                     r = 0.3, legend_size= 15, domain_color = NULL
                     ){
  cdd_data <- readLines(cdd_file)
  domain_loc <- cdd_to_loc(cdd_data)
  
  gene_length <- fastaleng(fasta_file)
  
  if (is.null(the_order)) {
    the_order=NULL
  }else{
    the_order <- readLines(the_order)
  }
  
  if (is.null(domain_select)) {
    motif_select=NULL
  }else{
    motif_select = domain_select
  }
  
  if (is.null(domain_color)) {
    domain_color = NULL
  }

  motif_plot(domain_loc, gene_length, the_order = the_order, 
             motif_select = motif_select, shape = shape, r=r, 
             legend_size= legend_size, motif_color =domain_color)
}
