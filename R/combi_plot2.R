#' Get ggplot2 files to facilitate free combination in patchwork
#' 
#' @title combi_p2
#' @param order_path The path of order file (.txt/.csv).
#' @param gff_path The path of .gff/gtf file.
#' @param meme_path The path of .meme/mast file.
#' @param pfam_path The path of pfam result file (.tsv).
#' @param cdd_path The path of cdd result file (.txt).
#' @param fa_path The path of protein file (.fa/fasta).
#' @param smart_path Do SMART or not. (TRUE or FALSE) 
#' @param plantcare_path The path of plantcare file (.tab).
#' @param promoter_length The length of promoter.
#' @param renamefile Rename file. Two cols: new_name and old_name.
#' @param shape RoundRect or Rect.
#' @param r The radius of rounded corners.
#' @param legend_size The size of legend.
#' @export
#' @author Shiqi Zhao
#' @return list
#' 
#' @import ggtree
#' @importFrom treeio read.newick 
#' @importFrom utils write.table
#' @importFrom stats setNames
#' @examples
#' order_path <- system.file("extdata", "order.csv", package = "BioVizSeq")
#' gff_path <- system.file("extdata", "idpro.gff3", package = "BioVizSeq")
#' plot_file <- combi_p2(order_path, gff_path = gff_path)

combi_p2 <- function(order_path, gff_path = NULL, 
                    meme_path = NULL, pfam_path = NULL, 
                    cdd_path = NULL, fa_path = NULL, smart_path = FALSE, 
                    plantcare_path = NULL, promoter_length = NULL, 
                    renamefile = NULL, shape = "RoundRect", 
                    r = 0.3, legend_size= 6
                       ){
  
  if(is.null(gff_path)){
    p_gff = NULL
  }else{
    p_gff <- gff_plot(gff_path, the_order = order_path,
                      shape = shape, r=r) +
      labs(y="") + 
      theme(axis.title.y = element_blank(), 
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), 
            plot.margin = margin(10, 10, 10, 0),
            legend.position = "right",
            legend.margin = margin(0, 0, 0, 0),
            legend.text = element_text(size = 10),
            legend.key.size = unit(legend_size, "pt")) 
  }
  
  #meme
  if(is.null(meme_path)){
    p_meme = NULL
  }else{
    p_meme <- meme_plot(meme_path, the_order = order_path,
                        shape = shape, r=r) +
      labs(y="") + 
      theme(axis.title.y = element_blank(), 
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), 
            plot.margin = margin(10, 10, 10, 0),
            legend.position = "right",
            legend.margin = margin(0, 0, 0, 0),
            legend.text = element_text(size = 10),
            legend.key.size = unit(legend_size, "pt")) 
  }
  
  #CDD
  if(!is.null(cdd_path) && !is.null(fa_path)){
    p_cdd <- cdd_plot(cdd_path, fa_path, the_order = order_path,
                      shape = shape, r=r) +
      labs(y="") + 
      theme(axis.title.y = element_blank(), 
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), 
            plot.margin = margin(10, 10, 10, 0),
            legend.position = "right",
            legend.margin = margin(0, 0, 0, 0),
            legend.text = element_text(size = 10),
            legend.key.size = unit(legend_size, "pt"))
  }else{
     p_cdd = NULL
  }
  
  #SMART
  if(smart_path == TRUE && !is.null(fa_path)){
    p_smart <- smart_plot(fa_path, the_order = order_path,
                          shape = shape, r=r) +
      labs(y="") + 
      theme(axis.title.y = element_blank(), 
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), 
            plot.margin = margin(10, 10, 10, 0),
            legend.position = "right",
            legend.margin = margin(0, 0, 0, 0),
            legend.text = element_text(size = 10),
            legend.key.size = unit(legend_size, "pt"))
  }else{
    p_smart = NULL
  }
  
  
  #pfam
  if(is.null(pfam_path)){
    p_pfam = NULL
  }else{
    p_pfam <- pfam_plot(pfam_path, the_order = order_path,
                        shape = shape, r=r) +
      labs(y="") + 
      theme(axis.title.y = element_blank(), 
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), 
            plot.margin = margin(10, 10, 10, 0),
            legend.position = "right",
            legend.margin = margin(0, 0, 0, 0),
            legend.text = element_text(size = 10),
            legend.key.size = unit(legend_size, "pt"))
  }
  
  #plantcare
  if(!is.null(plantcare_path) && !is.null(promoter_length)){
    p_plantcare <- plantcare_plot(plantcare_path, 
                                  promoter_length = promoter_length, 
                                  the_order = order_path,
                                  shape = shape, r=r) +
      labs(y="") + 
      theme(axis.title.y = element_blank(), 
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), 
            plot.margin = margin(10, 10, 10, 0),
            legend.position = "right",
            legend.margin = margin(0, 0, 0, 0),
            legend.text = element_text(size = 10),
            legend.key.size = unit(legend_size, "pt"))
    
    p_plantcare1 <- plantcare_plot1(plantcare_path, the_order = order_path) +
      labs(y="") + 
      theme(axis.title.y = element_blank(), 
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), 
            plot.margin = margin(10, 10, 10, 0),
            legend.position = "right",
            legend.justification = c(0.5, 1),
            legend.margin = margin(0, 0, 0, 0),
            legend.text = element_text(size = 8),
            legend.key.size = unit(legend_size, "pt"))
    
    p_plantcare2 <- plantcare_plot2(plantcare_path, the_order = order_path) +
      labs(y="") + 
      theme(axis.title.y = element_blank(), 
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), 
            plot.margin = margin(10, 10, 10, 0),
            legend.position = "right",
            legend.justification = c(0.5, 1),
            legend.margin = margin(0, 0, 0, 0),
            legend.text = element_text(size = 8),
            legend.key.size = unit(legend_size, "pt"))
    
    p_plantcare3 <- plantcare_plot3(plantcare_path, the_order = order_path) +
      labs(y="") + 
      theme(axis.title.y = element_blank(), 
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), 
            plot.margin = margin(10, 10, 10, 0),
            legend.position = "right",
            legend.justification = c(0.5, 1),
            legend.margin = margin(0, 0, 0, 0),
            legend.text = element_text(size = 8),
            legend.key.size = unit(legend_size, "pt"))
  }else{
    p_plantcare = NULL
    p_plantcare1 = NULL
    p_plantcare2 = NULL
    p_plantcare3 = NULL
  }
  
  plot_list <- list(p_gff = p_gff, p_pfam = p_pfam, 
                    p_meme = p_meme, p_smart = p_smart, p_cdd = p_cdd,
                    p_plantcare = p_plantcare, p_plantcare1 = p_plantcare1,
                    p_plantcare2 = p_plantcare2, p_plantcare3 = p_plantcare3
                    )
  
  return(plot_list)
  
}