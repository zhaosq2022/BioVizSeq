#' Extract the information of protein sequence
#' 
#' @title ProtParam_calc
#' @param input_file The path of potein fasta file.
#' @export
#' @author Shiqi Zhao
#' @return data.frame
#' 
#' @import httr
#' @importFrom seqinr read.fasta
#' @importFrom stats runif
#' 

ProtParam_calc <- function(input_file){
  
  file_dir <- file.path(tempdir(), "ProtParam_results")
  dir.create(file_dir)
  submit_to_protparam <- function(input_file){
    submit_url <- "https://web.expasy.org/cgi-bin/protparam/protparam"
    if (!file.exists(input_file)) {
      stop("Input file does not exist.")
    }
    
    output_directory <- file_dir
    
    seqs <- seqinr::read.fasta(input_file, seqtype = "AA", as.string = TRUE, whole.header = TRUE)
    
    for (i in seq_along(seqs)) {
      seq_id <- names(seqs)[i]
      output_file <- file.path(output_directory, paste0(seq_id, "_results.txt"))
      
      message("Submitting sequence ", seq_id, "...\n")
      
      post_content <- list(
        sequence = as.character(seqs[[i]][1]),
        TEXTONLY = 1
      )
      response <- httr::POST(submit_url, body = post_content, encode = "form")
      res_content <- httr::content(response, "text", encoding = "UTF-8")  
      writeLines(res_content, output_file)
      Sys.sleep(runif(1,1,3))
    }
  }
    
    #protparam_txt:the path of protparam file
  protparam_file_convert <- function(protparam_txt){
      ID <- gsub(".*\\/(.*)_results\\.txt", "\\1", protparam_txt)
      lines <- readLines(protparam_txt)
      lines <- lines[grepl("Number of amino acids|Molecular weight|Theoretical pI|The instability index|Aliphatic index|Grand average of hydropathicity", lines)]
      cleaned_text <- gsub("<.*?>", "", lines)
      cleaned_text <- gsub(": ", "\t", cleaned_text)
      cleaned_text <- gsub(" \\(GRAVY\\):", "\t", cleaned_text)
      cleaned_text <- gsub(" \\(II\\) is computed to be ", "\t", cleaned_text)
      id_text <- paste("ID", ID, sep = '\t')
      id_text <- paste(id_text, "Number of amino acids", sep = '\n')
      cleaned_text <- gsub("Number of amino acids", id_text, cleaned_text)
      
      text_connection <- textConnection(cleaned_text)
      cleaned_text <- read.table(text_connection, row.names = 1,  sep='\t')
      filtered_text <- as.data.frame(t(cleaned_text))
      row.names(filtered_text) <- NULL
      return(filtered_text)
  }
  
  submit_to_protparam(input_file)
  
  file_list <- list.files(path = file_dir, full.names = TRUE)
  
  data_list <- lapply(file_list, protparam_file_convert)
  
  combined_data <- do.call(rbind, data_list)
  
  unlink(file_dir, recursive = TRUE)
  
  return(combined_data)
    
}