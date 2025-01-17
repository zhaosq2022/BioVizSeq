#' Get ggplot2 files to facilitate free combination in patchwork
#' 
#' @title combi_p
#' @param tree_path The path of tree file (.newick).
#' @param gff_path The path of .gff/gtf file.
#' @param meme_path The path of .meme/mast file.
#' @param pfam_path The path of pfam result file (.tsv).
#' @param cdd_path The path of cdd result file (.txt).
#' @param fa_path The path of protein file (.fa/fasta).
#' @param smart_path Do SMART or not. (TRUE or FALSE) 
#' @param plantcare_path The path of plantcare file (.tab).
#' @param promoter_length The length of promoter.
#' @param renamefile Rename file. Two cols: new_name and old_name.
#' @param groupfile Group information. Two cols: label and Group.
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
#' tree_path <- system.file("extdata", "idpep.nwk", package = "BioVizSeq")
#' plot_file <- combi_p(tree_path)

combi_p <- function(tree_path, gff_path = NULL, 
                    meme_path = NULL, pfam_path = NULL, 
                    cdd_path = NULL, fa_path = NULL, smart_path = FALSE, 
                    plantcare_path = NULL, promoter_length = NULL, 
                    renamefile = NULL, groupfile = NULL, shape = "RoundRect", 
                    r = 0.3, legend_size= 6
                       ){
  tree_file <- read.newick(tree_path, node.label = "support")
  p_tree <- ggtree(tree_file, branch.length = 'none') + theme_tree() + geom_tiplab(size = 4) + 
    xlim(NA,12) + geom_nodelab(aes(label = support), hjust=1.2, vjust= -0.3, size=3)
  
  the_order <- get_taxa_name(p_tree)
  
  if(is.null(renamefile)){
    p_tree = p_tree
  }else{
    name_map <- setNames(renamefile$new_name, renamefile$old_name)
    p_tree$data$label <- name_map[p_tree$data$label]
  }
  
  if(is.null(groupfile)){
    p_tree = p_tree
  }else{
    group_loc <- merge(p_tree$data, groupfile, by= "label")
    group_loc$yend <- group_loc$y
    group_loc$yend <- ifelse(group_loc$yend == min(group_loc$y), group_loc$yend + 0.5, group_loc$yend)
    group_loc$y <- ifelse(group_loc$y == max(group_loc$yend), group_loc$y - 0.5, group_loc$y)
    
    p_tree <- p_tree + geom_rect(data = group_loc, aes(
      xmin = 0,
      xmax = x,
      ymin = y+0.5,
      ymax = yend-0.5,
      fill = Group
    ),alpha = 0.5) + 
      geom_tree() + 
      theme(legend.position = "right",
            legend.margin = margin(0, 0, 0, 0),
            legend.title = element_blank(),
            legend.text = element_text(size = 10),
            legend.key.size = unit(legend_size, "pt")
      ) 
  }

  
  file_dir <- file.path(tempdir(), "the_order")
  
  if (!dir.exists(file_dir)) {
    dir.create(file_dir)
  }
  
  order_path <- file.path(file_dir, "the_order.csv")
  write.table(the_order, order_path, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

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
            #panel.spacing = unit(0, "lines"),
            #legend.position = "top",
            legend.position = "right",
            #legend.justification = c(0.5, 1),
            legend.margin = margin(0, 0, 0, 0),
            legend.title = element_blank(),
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
            #panel.spacing = unit(0, "lines"),
            #legend.position = "top",
            legend.position = "right",
            #legend.justification = c(0.5, 1),
            legend.margin = margin(0, 0, 0, 0),
            legend.title = element_blank(),
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
            #legend.justification = c(0.5, 1),
            legend.margin = margin(0, 0, 0, 0),
            legend.title = element_blank(),
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
            #legend.position = "top",
            legend.position = "right",
            #legend.justification = c(0.5, 1),
            legend.margin = margin(0, 0, 0, 0),
            legend.title = element_blank(),
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
            #legend.justification = c(0.5, 1),
            legend.margin = margin(0, 0, 0, 0),
            legend.title = element_blank(),
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
            #legend.justification = c(0.5, 1),
            legend.margin = margin(0, 0, 0, 0),
            legend.title = element_blank(),
            legend.text = element_text(size = 10),
            legend.key.size = unit(legend_size, "pt"))
  }else{
    p_plantcare = NULL
  }

  unlink(file_dir, recursive = TRUE)
  
  plot_list <- list(p_tree = p_tree, p_gff = p_gff, p_pfam = p_pfam, 
                    p_meme = p_meme, p_smart = p_smart, p_cdd = p_cdd,
                    p_plantcare = p_plantcare
                    )
  
  return(plot_list)
  
}