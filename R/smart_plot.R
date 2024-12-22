#library(httr)
#library(Biostrings)

#' Extract the location information of domain from SMART result
#' 
#' @title smart_to_loc
#' @param input_file The path of potein fasta file.
#' @param do_pfam Include the pfam domain or not.
#' @export
#' @author Shiqi Zhao
#' @return list
#' 
#' @import httr
#' @importFrom dplyr filter 
#' @importFrom seqinr read.fasta

smart_to_loc <- function(input_file, do_pfam = TRUE){
  
  file_dir <- file.path(tempdir(), "SMART_results_task")
  dir.create(file_dir)

  submit_to_smart <- function(input_file, do_pfam = FALSE, do_signalp = FALSE, 
                              do_rep = FALSE, do_disembl = FALSE, do_schnipsel = FALSE) {

    submit_url <- "http://smart.embl.de/smart/show_motifs.pl"
    job_status_url <- "http://smart.embl.de/smart/job_status.pl"
    
    if (!file.exists(input_file)) {
      stop("Input file does not exist.")
    }
    
    output_directory <- file_dir
    
    seqs <- seqinr::read.fasta(input_file, seqtype = "AA", as.string = TRUE, whole.header = TRUE)
    
    for (i in seq_along(seqs)) {
      seq_id <- names(seqs)[i]
      output_file <- file.path(output_directory, paste0(seq_id, "_SMART_results.txt"))
      
      message("Submitting sequence ", seq_id, "...\n")
      
      post_content <- list(
        SEQUENCE = as.character(seqs[[i]][1]),
        TEXTONLY = 1
      )
      
      if (do_pfam) post_content$DO_PFAM <- 'DO_PFAM'
      if (do_signalp) post_content$INCLUDE_SIGNALP <- 'INCLUDE_SIGNALP'
      if (do_rep) post_content$DO_PROSPERO <- 'DO_PROSPERO'
      if (do_disembl) post_content$DO_DISEMBL <- 'DO_DISEMBL'
      if (do_schnipsel) post_content$INCLUDE_BLAST <- 'INCLUDE_BLAST'
      
      response <- POST(submit_url, body = post_content, encode = "form")
      
      if (http_status(response)$category == "Success") {
        res_content <- content(response, "text", encoding = "UTF-8")
        if (grepl("^--\\ SMART\\ RESULT", res_content)) {
          
          writeLines(res_content, output_file)
          #cat("Results saved to", output_file, "\n")
        } else {
          job_id <- sub(".* <a href='job_status\\.pl\\?jobid=(\\d+.+?)'>click here</a>.*", "\\1", res_content)
          if (nchar(job_id) == 0) {
            error_file <- file.path(output_directory, paste0(seq_id, "_SMART_error.html"))
            writeLines(res_content, error_file)
            #cat("SMART returned an error page, saved to", error_file, "\n")
            stop("Aborting further submissions.")
          } else {
            message("Job entered the queue with ID", job_id, ". Waiting for results.\n")
            repeat {
              Sys.sleep(5)
              job_status_response <- GET(paste0(job_status_url, "?jobid=", job_id))
              if (http_status(job_status_response)$category == "Success") {
                job_status_content <- content(job_status_response, "text", encoding = "UTF-8")
                if (grepl("-- SMART RESULT", job_status_content)) {
                  writeLines(job_status_content, output_file)
                  #cat("Results saved to", output_file, "\n")
                  break
                }
              } else {
                stop("SMART returned a web server error.")
              }
            }
          }
        }
      } else {
        stop("SMART returned a web server error.")
      }
      
      Sys.sleep(5) # Be nice to other users
    }
    #unlink(file_dir, recursive = TRUE)
  }
  
  smart_file_convert <- function(smart_txt){
    ID <- gsub(".*\\/(.*)_SMART_results\\.txt", "\\1", smart_txt)
    lines <- readLines(smart_txt)
    lines <- lines[lines != "" & !grepl("--|CRC_PASSED=|NUMBER_OF_FEATURES_FOUND=", lines)]
    cleaned_text <- paste(lines, collapse = "")
    cleaned_text <- gsub("DOMAIN=", "\n", cleaned_text)
    cleaned_text <- gsub("START=|END=|EVALUE=|TYPE=|STATUS=", "\t", cleaned_text)
    cleaned_text <- sub("\n", "", cleaned_text)
    text_connection <- textConnection(cleaned_text)
    cleaned_text <- read.table(text_connection, sep='\t')
    filtered_text <- filter(cleaned_text, .data$V1 != "low_complexity_region" & .data$V6 == "visible|OK")
    filtered_text2 <- data.frame(ID=ID, Domain =filtered_text[,1], start=filtered_text[,2], end=filtered_text[,3])
    filtered_text2[,1] <- sub(".*:", "", filtered_text2[,1])
    return(filtered_text2)
  }
  
  submit_to_smart(input_file, do_pfam = do_pfam)
   
   file_list <- list.files(path = file_dir, full.names = TRUE)
   
   data_list <- lapply(file_list, smart_file_convert)
   
   combined_data <- do.call(rbind, data_list)
   
   gene_length <- fastaleng(input_file)
   
   data_list <- list(gene_length = gene_length, table_loc = combined_data)
   unlink(file_dir, recursive = TRUE)
   return(data_list)
 }

#' Visualization of domain in SMART result file
#' 
#' @title smart_plot
#' @param fasta_file The path of protein fasta file.
#' @param the_order The path of order file. A List of Gene ID , One ID Per Line.
#' @param domain_select The domain ID which you want to align with.
#' @param shape RoundRect or Rect.
#' @param r The radius of rounded corners.
#' @param legend_size The size of legend.
#' @param domain_color The color set of domain.
#' @export
#' @author Shiqi Zhao
#' @return p


smart_plot <- function(fasta_file, the_order = NULL, domain_select = NULL, 
                       shape = "RoundRect", r = 0.3, legend_size= 15, 
                       domain_color = NULL){

  smart_domain_loc <- smart_to_loc(fasta_file)
  if (is.null(the_order)) {
    the_order=NULL
  }
  
  if (is.null(domain_select)) {
    motif_select=NULL
  }else{
    motif_select=domain_select
  }
  
  if (is.null(domain_color)) {
    domain_color = NULL
  }

  motif_plot(smart_domain_loc$table_loc, smart_domain_loc$gene_length, 
             the_order = the_order, motif_select = motif_select, 
             shape = shape, r=r, legend_size= legend_size, 
             motif_color=domain_color)
}

