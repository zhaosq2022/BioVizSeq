
#' Extract the location information of motif from mast or meme file
#' 
#' @title meme_to_loc
#' @param motif_file The motif data of mast or meme file.
#' @export
#' @author Shiqi Zhao
#' @return list
#' @examples
#' 
#' mast_path <- system.file("extdata", "mast.xml", package = "BioVizSeq") 
#' mast_file <- readLines(mast_path)
#' motif_loc <- meme_to_loc(mast_file)
#' 
#' @importFrom tidyr separate
#' @importFrom stringr str_detect
#' 

meme_to_loc <- function(motif_file){
  get_mast_location <- function(mast_file){
    
    process_mast <- function(mast_file){
      data <- gsub('"', "", mast_file)
      data <- grep('alt=MEME|sequence db=|idx=', data, value = TRUE)
      data <- gsub(".*alt=", "", data)
      data <- gsub(" nsite.*", "", data)
      data <- gsub(".*name=", "name=", data)
      data <- gsub(".*pos=", "pos=", data)
      data <- gsub(" pvalue.*", "", data)
      data <- gsub(">", "", data)
      data <- gsub("comment= ", "", data)
      return(data)
    }
    processed_mast <- process_mast(mast_file)
    
    seqlen_mast <- function(processed_data){
      data <- grep("name", processed_data, value = TRUE)
      data <- gsub(".*name=", "", data )
      data <- gsub(" length=", "\t", data)
      data <- as.data.frame(data)
      
      gene_length <- tidyr::separate(data, 1, into = c("ID", "length"), sep='\t')
      gene_length$length <- as.numeric(gene_length$length)
      return(gene_length)
    }
    gene_length <- seqlen_mast(processed_mast)
    
    data_raw <- processed_mast
    data <- grep("MEME", data_raw, value = TRUE, invert = TRUE)
    data <- as.data.frame(data)
    data <- tidyr::separate(data, 1, into = c("ID", "length"), sep=' ')
    
    #get the start position of motif in gene
    newfunction <- function(x){
      data_x <- x
      #table_null <- data.frame(ID=as.character(), motif=as.character(), start=as.numeric())
      table_loc <- data.frame(ID=as.character(), motif=as.numeric(), start=as.numeric())
      for(i in 1:nrow(data_x)){
        if(stringr::str_detect(data_x[i,1],"name") == TRUE){ #str_detect in stringr package
          id<-data_x[i,1]
        }else{
          tab_loc_i<-data.frame(ID=id, motif=data_x[i,2],start=data_x[i,1])
          table_loc<- rbind(table_loc, tab_loc_i)
        }
      }
      return(table_loc)
    }
    table_loc <- newfunction(data)
    
    table_loc$ID <- gsub("name=", "", table_loc$ID)
    table_loc$start <- gsub("pos=", "", table_loc$start)
    table_loc$motif <- gsub("idx=", "", table_loc$motif)
    table_loc$motif <- as.numeric(table_loc$motif) + 1
    table_loc$motif <- as.character(table_loc$motif)
    table_loc$start <- as.numeric(table_loc$start)
    
    # get the length of MEME motif
    get_meme_length <- function(x){
      data <- grep("MEME", x, value = TRUE)
      data <- gsub("MEME-", "", data )
      data <- gsub(" length=", "\t", data)
      data <- as.data.frame(data)
      
      ##separate in tidyr package
      meme_length <- tidyr::separate(data, 1, into = c("name", "length"), sep='\t')
      meme_length$length <- as.numeric(meme_length$length)
      return(meme_length)
    }
    meme_length <- get_meme_length(data_raw)
    
    #get the end position of motif in gene
    for(i in 1: nrow(table_loc)){
      table_loc[i,4] <- table_loc[i,3] + meme_length[as.numeric(table_loc[i,2]),2] - 1
    }
    colnames(table_loc)[4] <- "end"
    table_loc[,2] <- as.numeric(table_loc[,2])
    colnames(table_loc)[2] <- "Motif"
    data_list <- list(gene_length = gene_length, table_loc = table_loc)
    return(data_list)
  }
  get_meme_location <- function(meme_file){
    
    process_meme <- function(meme_file){
      data <- grep('<sequence id=|<motif id=|<scanned_sites sequence_id|<scanned_site motif_id', meme_file, value = TRUE)
      data <- gsub('"', "", data)
      data <- gsub("<sequence | weight=.*|.*alt=| sites=.*", "", data)
      data <- gsub("><scanned_site ", "length=\n", data)
      data <- unlist(strsplit(data, "\n"))
      data <- gsub("pvalue=.*length=", "length=", data)
      data <- gsub("<scanned_sites | pvalue=.*|strand=none |<scanned_site |motif_", "", data)
      return(data)
    }
    processed_meme <- process_meme(meme_file)
    
    seqlen_mast <- function(processed_data){
      data <- grep("name", processed_data, value = TRUE)
      data <- gsub(".*name=", "", data )
      data <- gsub(" length=", "\t", data)
      data <- as.data.frame(data)
      
      gene_length <- tidyr::separate(data, 1, into = c("ID", "length"), sep='\t')
      gene_length$length <- as.numeric(gene_length$length)
      return(gene_length)
    }
    gene_length <- seqlen_mast(processed_meme)
    
    data_raw <- processed_meme
    data <- grep("MEME|name=", data_raw, value = TRUE, invert = TRUE)
    data<-as.data.frame(data)
    data <- tidyr::separate(data, 1, into = c("ID", "length"), sep=' ')
    
    ID_convert <- grep("^id=sequence_", data_raw, value = TRUE)
    ID_convert <- gsub(" length=.*|id=|name=", "", ID_convert)
    ID_convert <- as.data.frame(ID_convert)
    ID_convert <- tidyr::separate(ID_convert, 1, into = c("ID", "id"), sep=' ')
    
    #get the start position of motif in gene
    newfunction <- function(x){
      data_x <- x
      #table_null <- data.frame(ID=as.character(), motif=as.character(), start=as.numeric())
      table_loc <- data.frame(ID=as.character(), motif=as.numeric(), start=as.numeric())
      for(i in 1:nrow(data_x)){
        if(str_detect(data_x[i,1],"sequence") == TRUE){ #str_detect in stringr package
          id<-data_x[i,1]
        }else{
          tab_loc_i<-data.frame(ID=id, motif=data_x[i,1],start=data_x[i,2])
          table_loc<- rbind(table_loc, tab_loc_i)
        }
      }
      return(table_loc)
    }
    table_loc <- newfunction(data)
    
    table_loc$ID <- gsub("sequence_id=", "", table_loc$ID)
    table_loc$start <- gsub("position=", "", table_loc$start)
    table_loc$motif <- gsub("id=", "", table_loc$motif)
    table_loc$start <- as.numeric(table_loc$start)+1
    
    # get the length of MEME motif
    get_meme_length <- function(x){
      data <- grep("MEME", x, value = TRUE)
      data <- gsub("MEME-", "", data )
      data <- gsub(" width=", "\t", data)
      data <- as.data.frame(data)
      
      ##separate in tidyr package
      meme_length <- tidyr::separate(data, 1, into = c("name", "length"), sep='\t')
      meme_length$name <- as.numeric(meme_length$name)
      meme_length$name <- as.character(meme_length$name)
      meme_length$length <- as.numeric(meme_length$length)
      return(meme_length)
    }
    meme_length <- get_meme_length(data_raw)
    
    #get the end position of motif in gene
    for(i in 1: nrow(table_loc)){
      table_loc[i,4] <- table_loc[i,3] + meme_length[as.numeric(table_loc[i,2]),2] - 1
    }
    colnames(table_loc)[4] <- "end"
    table_loc[,2] <- as.numeric(table_loc[,2])
    table_motif_loc <- merge(ID_convert, table_loc[, c("ID", "motif", "start", "end")], by = "ID", all.x = TRUE)
    table_motif_loc <- table_motif_loc[,-1]
    colnames(table_motif_loc)[1] <- "ID"
    colnames(table_motif_loc)[2] <- "Motif"
    data_list <- list(gene_length = gene_length, table_loc = table_motif_loc)
    return(data_list)
  }
  
  if(grepl("mast", motif_file[2])){
    motif_loc <- get_mast_location(motif_file)
  }else{
    motif_loc <- get_meme_location(motif_file)
  }
  return(motif_loc)
}


#' Visualization of motif in meme file or mast file
#' 
#' @title meme_plot
#' @param meme_file The path of meme file or mast file.
#' @param the_order The path of order file. A List of Gene ID , One ID Per Line.
#' @param motif_select The motif ID which you want to align with.
#' @param shape RoundRect or Rect.
#' @param show_motif_id Display the name of the motif.
#' @param r The radius of rounded corners.
#' @param legend_size The size of legend.
#' @param motif_color The color set of motif.
#' @export
#' @author Shiqi Zhao
#' @return p
#' @examples
#' 
#' mast_path <- system.file("extdata", "mast.xml", package = "BioVizSeq")
#' meme_plot(mast_path)
#' 
#' meme_plot(mast_path, motif_select="1", show_motif_id = TRUE)
#' 
#' order_path <- system.file("extdata", "order.csv", package = "BioVizSeq")
#' meme_plot(mast_path, the_order=order_path, motif_select="1")

meme_plot <- function(meme_file, the_order = NULL, motif_select = NULL, 
                      shape = "RoundRect", show_motif_id = FALSE, r = 0.3, 
                      legend_size= 15, motif_color = NULL){
  raw_file <- readLines(meme_file)
  meme_motif_loc <- meme_to_loc(raw_file)

  if (is.null(the_order)) {
    the_order=NULL
  }else{
    the_order <- readLines(the_order)
  }
  
  if (is.null(motif_select)) {
    motif_select=NULL
  }
  
  if (is.null(show_motif_id)|show_motif_id == FALSE|show_motif_id == F) {
    show_motif_id = FALSE
  }else if (show_motif_id == TRUE|show_motif_id == T) {
    show_motif_id = TRUE
  }
  
  if (is.null(motif_color)) {
    motif_color = NULL
  }
  
  motif_plot(meme_motif_loc$table_loc, meme_motif_loc$gene_length,
             the_order=the_order, motif_select = motif_select, 
             shape = shape, show_motif_id = show_motif_id, r=r, 
             legend_size= legend_size, motif_color =motif_color)
}

#' Get motif sequence from meme file or mast file
#' 
#' @title meme_seq
#' @param meme_file The path of meme file or mast file.
#' @export
#' @author Shiqi Zhao
#' @return data.frame
#' @examples
#' 
#' mast_path <- system.file("extdata", "mast.xml", package = "BioVizSeq")
#' mast_file <- readLines(mast_path)
#' motifseq<- meme_seq(mast_file)
#' 
meme_seq <- function(meme_file){
  raw_file <- meme_file
  if(grepl("mast", raw_file[2])){
    data <- gsub('"', "", raw_file)
    data <- grep('motif db=', data, value = T)
    data <- gsub(".*id=|>| bad=.*", "", data)
    data <- gsub("alt=MEME-|length=|nsites=|evalue=", "", data)
    data <- as.data.frame(data)
    data <- tidyr::separate(data, 1, into = c("Sequence", "Motif", "Length", "Nsites", "Evalue"), sep=' ')
    data[,3] <- as.numeric(data[,3])
    data[,4] <- as.numeric(data[,4])
  }else{
    data <- gsub('"', "", raw_file)
    data <- grep('motif id=', data, value = T)
    data <- gsub("ic=.*p_value", "p_value", data)
    data <- gsub(".*name=| bayes_threshold.*", "", data)
    data <- gsub("alt=MEME-|width=|sites=|p_value=|e_value=", "", data)
    data <- as.data.frame(data)
    data <- tidyr::separate(data, 1, into = c("Sequence", "Motif", "Width", "sites", "p_value", "e_value"), sep=' ')
    data[,3] <- as.numeric(data[,3])
    data[,4] <- as.numeric(data[,4])
  }
  return(data)
}
