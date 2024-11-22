#' Global Variables
#'
#' This block declares global variables to avoid R CMD check warnings.
#'
#' @name globalVariables
#' @title Global Variables Declaration
#' @keywords internal
utils::globalVariables(c("motif", "start", "ID", "length", "seq_order", "end", "name", "class", "class.x"))


#' Draws multiple rounded rectangle.
#' 
#' motif_plot() draws multiple rounded rectangle to represent the above elements of biosequences, 
#' but not limited to biosequences
#' 
#' @title motif_plot
#' @param motif_loc A data.frame contains for columuns: ID, motif, start, end.
#' @param gene_length A data.fram of the length of biosequences. Two columns: ID, length.
#' @param the_order A List of Gene ID , One ID Per Line.
#' @param motif_select The motif ID which you want to align with.
#' @param shape RoundRect or Rect.
#' @param show_motif_id Display the name of the motif.
#' @param r The radius of rounded corners.
#' @param legend_size The size of legend.
#' @param motif_color The color set of motif.
#' 
#' @export
#' @author Shiqi Zhao
#' @return p
#' @examples
#' df <- data.frame(
#'  ID = rep(c("geneA", "geneB", "geneC"), each = 3),
#'  motif = rep(c("1", "2", "3"), times = 3),
#'  start = c(1, 3, 6, 1, 6, 10, 10, 7, 17),
#'  end = c(3, 5, 11, 3, 8, 15, 12, 9, 22)
#'  )
#'  
#'  length_data <- data.frame(
#'  ID = c("geneA", "geneB", "geneC"),
#'  length = c(15, 27, 30)
#'  )
#'  
#'  order_data <- c("geneB", "geneA", "geneC")
#'  
#'  motif_plot(df, length_data)
#'  motif_plot(df, length_data, the_order = order_data)
#' 
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal


motif_plot <- function(motif_loc, gene_length, the_order = NULL, 
                       motif_select = NULL, shape = "RoundRect", 
                       show_motif_id = FALSE, r = 0.3, legend_size= 15, 
                       motif_color = NULL){
  if (!shape %in% c("RoundRect", "Rect")) {
    message("The parameter 'shape' must be either 'RoundRect' or 'Rect'")
  }else{
    legend_name <- colnames(motif_loc)[2]
    colnames(motif_loc) <- c("ID", "motif", "start", "end")
    colnames(gene_length) <- c("ID", "length")
    
    legend_number <- as.character(levels(as.factor(motif_loc[,2])))
    
    motif_loc[,2]<-as.factor(motif_loc[,2])
    
    if (is.null(the_order)) {
      the_order <- as.vector(unique(motif_loc$ID))
    }else{
      the_order <- as.data.frame(the_order)
      the_order <- as.vector(the_order[,1])
    }
    
    motif_loc$ID<-factor(motif_loc$ID,levels=rev(the_order))
    gene_length$ID<-factor(gene_length$ID,levels=rev(the_order))
    gene_length$start <- 0
    
    for(i in 1:nrow(motif_loc)){
      ID_levels<-levels(motif_loc$ID)
      motif_loc[i,5]<- which(ID_levels == motif_loc[i,1])
    }
    
    colnames(motif_loc)[5]  <- "seq_order"
    
    if(is.null(motif_select)){
      motif_loc <- motif_loc
      gene_length <- gene_length
    }else{
      select_motif <- filter(motif_loc, motif == motif_select)
      max_loc <- max(select_motif$start)
      
      for(i in 1:nrow(motif_loc)){
        for(j in 1:nrow(select_motif)){
          if(motif_loc$ID[i] == select_motif$ID[j]){
            motif_loc$start[i] <- motif_loc$start[i] - select_motif$start[j]
            motif_loc$end[i] <- motif_loc$end[i] - select_motif$start[j]
          }
        }
      }
      
      for(i in 1:nrow(gene_length)){
        for(j in 1:nrow(select_motif)){
          if(gene_length$ID[i] == select_motif$ID[j]){
            gene_length$start[i] <- gene_length$start[i] - select_motif$start[j]
            gene_length$length[i] <- gene_length$length[i] - select_motif$start[j]
          }
        }
      }
    }

    if(shape == "RoundRect"){
      if(is.null(motif_color)){
        p <- ggplot() +
          geom_segment(data = gene_length, aes(x=start,xend=length,y=ID,yend=ID), colour="black") + 
          geom_rrect(data = motif_loc, 
            aes(
              xmin = start,
              xmax = end,
              ymin = seq_order-0.3,
              ymax = seq_order+0.3,
              fill = motif
            ),
            r = r
          ) + 
          labs(x = "Protein length (5'-3')", y = 'Protein name') + 
          theme(panel.grid.major = element_blank(),
                panel.grid.minor=element_blank(),
                legend.key.size= unit(legend_size, "pt"),
                panel.background = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.x = element_line(linewidth = 1, colour = "black", linetype = 1)
          ) + 
          scale_fill_discrete(limits = legend_number) + 
          scale_x_continuous(expand = c(0, 0)) +
          labs(fill = legend_name)
      }else{
        p <- ggplot() +
          geom_segment(data = gene_length, aes(x=start,xend=length,y=ID,yend=ID),colour="black") + 
          geom_rrect(data = motif_loc,
            aes(
              xmin = start,
              xmax = end,
              ymin = seq_order-0.3,
              ymax = seq_order+0.3,
              fill = motif
            ),
            r = r
          ) + 
          scale_fill_manual(values = motif_color) + 
          labs(x = "Protein length (5'-3')", y = 'Protein name') + 
          theme(panel.grid.major = element_blank(),
                panel.grid.minor=element_blank(),
                legend.key.size= unit(legend_size, "pt"),
                panel.background = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.x = element_line(linewidth = 1, colour = "black", linetype = 1)
          ) +
          scale_x_continuous(expand = c(0, 0)) +
          labs(fill = legend_name)
      }
      
    }else{
      if(is.null(motif_color)){
        p <- ggplot() +
          geom_segment(data = gene_length, aes(x=start,xend=length,y=ID,yend=ID),colour="black") + 
          geom_rect(data = motif_loc,
            aes(
              xmin = start,
              xmax = end,
              ymin = seq_order-0.3,
              ymax = seq_order+0.3,
              fill = motif
            )
          ) + 
          labs(x = "Protein length (5'-3')", y = 'Protein name') + 
          theme(panel.grid.major = element_blank(),
                panel.grid.minor=element_blank(),
                legend.key.size= unit(legend_size, "pt"),
                panel.background = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.x = element_line(linewidth = 1, colour = "black", linetype = 1)
          ) + 
          scale_fill_discrete(limits = legend_number) +
          scale_x_continuous(expand = c(0, 0)) +
          labs(fill = legend_name)
      }else{
        p <- ggplot() +
          geom_segment(data = gene_length, aes(x=start,xend=length,y=ID,yend=ID),colour="black") + 
          geom_rect(data = motif_loc,
            aes(
              xmin = start,
              xmax = end,
              ymin = seq_order-0.3,
              ymax = seq_order+0.3,
              fill = motif
            )
          ) + 
          scale_fill_manual(values = motif_color) + 
          labs(x = "Protein length (5'-3')", y = 'Protein name') + 
          theme(panel.grid.major = element_blank(),
                panel.grid.minor=element_blank(),
                legend.key.size= unit(legend_size, "pt"),
                panel.background = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.x = element_line(linewidth = 1, colour = "black", linetype = 1)
          ) + 
          scale_x_continuous(expand = c(0, 0)) +
          labs(fill = legend_name)
      }
    }
    
    if(is.null(show_motif_id)|show_motif_id == FALSE|show_motif_id == F){
      p
    }else if (show_motif_id == TRUE|show_motif_id == T){
      p + geom_text(data = motif_loc,aes(x = (start + end) / 2, y = (seq_order + seq_order) / 2, label = motif), color = "black")
    }
    
  }
}
