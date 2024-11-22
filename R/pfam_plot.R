#' Extract the location information of domain from pfam result
#' 
#' @title pfam_to_loc
#' @param pfam_data The result file (.tsv) of pfam (via InterPro).
#' @export
#' @author Shiqi Zhao
#' @return list
#' @examples
#' pfam_path <- system.file("extdata", "iprscan.tsv", package = "BioVizSeq")
#' pfam_file <- read.table(pfam_path, sep='\t', header = FALSE)
#' motif_loc <- pfam_to_loc(pfam_file)

pfam_to_loc <- function(pfam_data){
  gene_length <- pfam_data[,c(1,3)]
  colnames(gene_length) <- c("ID", "length")
  domain_loc <- pfam_data[,c(1,6,7,8)]
  colnames(domain_loc) <- c("ID", "Domain", "start", "end")
  data_list <- list(gene_length = gene_length, table_loc = domain_loc)
  return(data_list)
}

#' Visualization of domain in pfam result file
#' 
#' @title pfam_plot
#' @param pfam_file The path of meme file or mast file.
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
#' pfam_path <- system.file("extdata", "iprscan.tsv", package = "BioVizSeq")
#' order_path <- system.file("extdata", "order.csv", package = "BioVizSeq")
#' pfam_plot(pfam_path)
#' pfam_plot(pfam_path, the_order=order_path)
#' 
#' @importFrom utils read.table

pfam_plot <- function(pfam_file, the_order = NULL, domain_select = NULL, 
                      shape = "RoundRect", r = 0.3, legend_size= 15,
                      domain_color = NULL){
  pfam_data <- read.table(pfam_file, sep='\t', header = F)
  domain_loc <- pfam_to_loc(pfam_data)
  
  if (is.null(the_order)) {
    the_order=NULL
  }else{
    the_order <- readLines(the_order)
  }
  
  if (is.null(domain_select)) {
    motif_select=NULL
  }else{
    motif_select=domain_select
  }
  
  if (is.null(domain_color)) {
    domain_color = NULL
  }

  motif_plot(domain_loc$table_loc, domain_loc$gene_length, 
             the_order = the_order, motif_select = motif_select, 
             shape = shape, r=r, legend_size=legend_size,
             motif_color=domain_color)
}

