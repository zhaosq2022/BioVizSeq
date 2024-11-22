#' Extract the location information of element from gff or gtf file 
#' 
#' @title gff_to_loc
#' @param gff_data gff file.
#' @param mRNA_ID The mRNA you selected. If NULL, it means selecting all mRNAs.
#' @export
#' @author Shiqi Zhao
#' @return list
#' @examples
#' 
#' gff_path <- system.file("extdata", "test.gff", package = "BioVizSeq")
#' gff_data <- read.table(gff_path, header = FALSE, sep = '\t')
#' gff_loc <- gff_to_loc(gff_data)
#' 
#' ID_path <- system.file("extdata", "ID_select.csv", package = "BioVizSeq")
#' mRNA_ID <- readLines(ID_path) 
#' gff_loc <- gff_to_loc(gff_data, mRNA_ID=mRNA_ID)

gff_to_loc <- function(gff_data, mRNA_ID = NULL){
  
  gff_mRNA <- gff_data[grepl("mRNA", gff_data[,3], ignore.case = T),]
  
  mRNA_length <- gff_mRNA[,c(9,4,5)]
  
  mRNA_length[,1] <- gsub(";.*", "", mRNA_length[,1])
  mRNA_length[,1] <- gsub(".*=", "", mRNA_length[,1])
  mRNA_length[,3] <- mRNA_length[,3] - mRNA_length[,2]
  mRNA_length <- mRNA_length[,-2]
  
  if(!is.null(mRNA_ID)){
    mRNA_length <- mRNA_length[mRNA_length[,1] %in% mRNA_ID,]
  }
  
  gff_data <- gff_data[grepl("exon|cds", gff_data[,3], ignore.case = T),]
  
  gff_data[,9] <- gsub(".*Parent=", "", gff_data[,9])
  
  gff_data[,9] <- gsub(";.*", "", gff_data[,9])
  
  if(!is.null(mRNA_ID)){
    gff_data <- gff_data[gff_data[,9] %in% mRNA_ID,]
  }
  
  table_loc <- data.frame(ID=as.character(), element=as.character(), start=as.numeric(), end=as.numeric())
  for(i in 1:nrow(mRNA_length)){
    data_mRNA <- gff_data[gff_data[,9] == mRNA_length[i,1],]
    
    if(data_mRNA[1,7] == "+"){
      data_mRNA <- data_mRNA[order(data_mRNA[, 4]), ]
      mRNA_start <- data_mRNA[1,4]
      data_mRNA[,4] <- data_mRNA[,4] - mRNA_start
      data_mRNA[,5] <- data_mRNA[,5] - mRNA_start
      data_mRNA_new <- data_mRNA[,c(9,3,4,5)]
      
      data_CDS <- data_mRNA_new[grepl("CDS", data_mRNA_new[,2], ignore.case = T),]
      data_CDS <- data_CDS[order(data_CDS[, 3]), ]
      CDS_UTR <- data_CDS
      data_exon <- data_mRNA_new[grepl("exon", data_mRNA_new[,2], ignore.case = T),]
      data_CDS <- data_CDS[order(data_CDS[, 3]), ]
      
      for(i in 1:nrow(data_exon)){
        if(data_exon[i,4] < data_CDS[1,4]){
          CDS_UTR <- rbind(CDS_UTR, data_exon[i,])
        }else if(data_exon[i,4] == data_CDS[1,4] && data_exon[i,3] < data_CDS[1,3]){
          data_exon[i,4] <- data_CDS[1,3]
          CDS_UTR <- rbind(CDS_UTR, data_exon[i,])
        }else if(data_exon[i,3] > data_CDS[nrow(data_CDS),3]){
          CDS_UTR <- rbind(CDS_UTR, data_exon[i,])
        }else if(data_exon[i,3] == data_CDS[nrow(data_CDS),3] && data_exon[i,4] > data_CDS[nrow(data_CDS),4]){
          data_exon[i,3] <- data_CDS[nrow(data_CDS),4]
          CDS_UTR <- rbind(CDS_UTR, data_exon[i,])
        }else if(data_exon[1,3] < data_CDS[1,3] && data_exon[1,4] > data_CDS[1,4]){
          utr5 <- data_exon[1,]
          utr5[1,4] <- data_CDS[1,3]
          utr3 <- data_exon[1,]
          utr3[1,3] <-data_CDS[1,4]
          CDS_UTR <- rbind(CDS_UTR, utr5, utr3)
        }
      }
      table_loc <- rbind(table_loc, CDS_UTR)
    }else if(data_mRNA[1,7] == "-"){
      data_mRNA <- data_mRNA[order(data_mRNA[, 5]), ]
      mRNA_start <- data_mRNA[nrow(data_mRNA),5]
      data_mRNA[,4] <- abs(data_mRNA[,4] - mRNA_start)
      data_mRNA[,5] <- abs(data_mRNA[,5] - mRNA_start)
      data_mRNA_new <- data_mRNA[,c(9,3,4,5)]
      data_mRNA_new[,3] <- data_mRNA_new[,4]
      data_mRNA_new[,4] <- data_mRNA[,4]
      
      data_CDS <- data_mRNA_new[grepl("CDS", data_mRNA_new[,2], ignore.case = T),]
      data_CDS <- data_CDS[order(data_CDS[, 3]), ]
      CDS_UTR <- data_CDS
      data_exon <- data_mRNA_new[grepl("exon", data_mRNA_new[,2], ignore.case = T),]
      data_CDS <- data_CDS[order(data_CDS[, 3]), ]
      
      for(i in 1:nrow(data_exon)){
        if(data_exon[i,4] < data_CDS[1,4]){
          CDS_UTR <- rbind(CDS_UTR, data_exon[i,])
        }else if(data_exon[i,4] == data_CDS[1,4] && data_exon[i,3] < data_CDS[1,3]){
          data_exon[i,4] <- data_CDS[1,3]
          CDS_UTR <- rbind(CDS_UTR, data_exon[i,])
        }else if(data_exon[i,3] > data_CDS[nrow(data_CDS),3]){
          CDS_UTR <- rbind(CDS_UTR, data_exon[i,])
        }else if(data_exon[i,3] == data_CDS[nrow(data_CDS),3] && data_exon[i,4] > data_CDS[nrow(data_CDS),4]){
          data_exon[i,3] <- data_CDS[nrow(data_CDS),4]
          CDS_UTR <- rbind(CDS_UTR, data_exon[i,])
        }else if(data_exon[1,3] < data_CDS[1,3] && data_exon[1,4] > data_CDS[1,4]){
          utr5 <- data_exon[1,]
          utr5[1,4] <- data_CDS[1,3]
          utr3 <- data_exon[1,]
          utr3[1,3] <-data_CDS[1,4]
          CDS_UTR <- rbind(CDS_UTR, utr5, utr3)
        }
        
      }
      table_loc <- rbind(table_loc, CDS_UTR)
    }
  }
  table_loc[,2] <- gsub("exon", "UTR", table_loc[,2], ignore.case = TRUE)
  colnames(table_loc) <- c("ID", "element", "start", "end")
  
  data_list <- list(mRNA_length = mRNA_length, table_loc = table_loc)
  return(data_list)
}


#' Visualization of element in gff or gtf file
#' 
#' @title gff_plot
#' @param gff_file The path of gff file.
#' @param the_order The path of order of mRNA. It is also the mRNA you want to showcase. A List of Gene ID , One ID Per Line.
#' @param shape RoundRect or Rect.
#' @param r The radius of rounded corners.
#' @param element_color The color set of element.
#' @param legend_size The size of legend.
#' @export
#' @author Shiqi Zhao
#' @return p
#' @examples
#' gff_path <- system.file("extdata", "test.gff", package = "BioVizSeq")
#' gff_plot(gff_path)

gff_plot <- function(gff_file, the_order = NULL, shape = "Rect", 
                     r = 0.3, legend_size= 15, 
                     element_color = NULL){
  gff_data <- read.table(gff_file, header = FALSE, sep = '\t')
  if(!is.null(the_order)){
    mRNA_ID <- readLines(the_order)
  }else{
    mRNA_ID = NULL
  }
  
  if (is.null(element_color)) {
    element_color = NULL
  }
  
  gff_to_loc_data <- gff_to_loc(gff_data, mRNA_ID = mRNA_ID)
  
  motif_plot(gff_to_loc_data$table_loc, gff_to_loc_data$mRNA_length, 
             the_order = mRNA_ID, shape = shape, r = r, 
             legend_size= legend_size, motif_color=element_color) +
    labs(x="", y="")
}
