#' Upload the promoter file to Plantcare database
#' 
#' 
#' Due to the file size limitation of plantcare on fasta, 
#' upload_fa_to_plantcare() first splits fasta file. 
#' Then uploads the splited fasta files to the plantcare database, 
#' and automatically returns the results to the email provided by the user.
#' 
#' upload_fa_to_plantcare("the path/test.fasta", "your e-mail address")
#' 
#' @title upload_fa_to_plantcare
#' @param fasta_file The path of promoter file.
#' @param email e-mail address.
#' @export
#' @author Shiqi Zhao
#' @return plantcare_result
#'
#' @importFrom seqinr read.fasta
#' @importFrom seqinr write.fasta

#library(httr)
#library(seqinr)
upload_fa_to_plantcare <- function(fasta_file, email) {
  # Read the FASTA file
  fasta_sequences <- read.fasta(fasta_file, seqtype = "DNA", as.string = TRUE, whole.header = TRUE)
  
  # Get the total number of sequences
  total_sequences <- length(fasta_sequences)
  
  # If the number of sequences exceeds 50, split into multiple files
  if (total_sequences > 50) {
    message("Number of sequences exceeds 50, splitting...\n")
    
    # Calculate the number of split files
    split_files <- split(fasta_sequences, ceiling(seq_along(fasta_sequences) / 50))
    file_names <- vector()
    
    # Save the split files
    for (i in seq_along(split_files)) {
      split_file_name <- paste0("split_fasta_", i, ".fa")
      write.fasta(split_files[[i]], names = names(split_files[[i]]), file.out = split_file_name)
      file_names <- c(file_names, split_file_name)
      message("Saved file:", split_file_name, "\n")
    }
  } else {
    file_names <- fasta_file  # If sequences are 50 or fewer, just upload the original file
  }
  
  # Upload the files
  for (file in file_names) {
    url <- "https://bioinformatics.psb.ugent.be/webtools/plantcare/cgi-bin/CallMat_onCluster.htpl"
    
    fasta_content <- upload_file(file)
    
    response <- POST(
      url,
      body = list(
        Field_UserEmail = email,
        Field_File = fasta_content,
        MatInspector = "Search"
      ),
      encode = "multipart"
    )
    
    # Check response result
    if (status_code(response) == 200) {
      message("File", file, "submitted successfully!\n")
    } else {
      message("File", file, "submission failed, status code:", status_code(response), "\n")
    }
    
    # Random wait time (between 10 and 60 seconds)
    wait_time <- sample(10:60, 1)
    message("Waiting for", wait_time, "seconds before uploading the next file...\n")
    Sys.sleep(wait_time)
  }
}

#' Classify the functions of cis element from Plantcare
#' 
#' @title plantcare_classify
#' @param plantcare_file The result file (.tab) of Plantcare.
#' @export
#' @author Shiqi Zhao
#' @return data.frame
#' @examples
#' plantcare_path <- system.file("extdata", "plantCARE_output.tab", package = "BioVizSeq") 
#' plantcare_file <- read.table(plantcare_path, header = FALSE, sep = '\t', quote="")
#' plantcare_data <- plantcare_classify(plantcare_file)
plantcare_classify <- function(plantcare_file){
  plantcare_data <- plantcare_file
  plantcare_data <- plantcare_data[plantcare_data[, 2] != "", ]
  plantcare_data <- plantcare_data[!(plantcare_data[, 8] %in% c("", "?")), ]
  plantcare_data$V9 <- NA
  plantcare_data$V10 <- NA
  plantcare_data[,2] <- gsub("G-Box", "G-box", plantcare_data[,2])
  
  plantcare_data$V9 <- ifelse(grepl("light respon", plantcare_data[, 8]), "light", plantcare_data[, 9])
  plantcare_data$V9 <- ifelse(grepl("heat", plantcare_data[, 8]), "heat", plantcare_data[, 9])
  plantcare_data$V9 <- ifelse(grepl("low-temp", plantcare_data[, 8]), "low-temperature", plantcare_data[, 9])
  plantcare_data$V9 <- ifelse(grepl("circadian", plantcare_data[, 8]), "circadian", plantcare_data[, 9])
  plantcare_data$V9 <- ifelse(grepl("defense and stress responsivenes", plantcare_data[, 8]), "mixed stress", plantcare_data[, 9])
  plantcare_data$V9 <- ifelse(grepl("wound", plantcare_data[, 8]), "wound", plantcare_data[, 9])
  plantcare_data$V9 <- ifelse(grepl("drought", plantcare_data[, 8]), "drought", plantcare_data[, 9])
  plantcare_data$V9 <- ifelse(grepl("anaerobic", plantcare_data[, 8]), "anaerobic", plantcare_data[, 9])
  plantcare_data$V9 <- ifelse(grepl("anoxic", plantcare_data[, 8]), "anoxic", plantcare_data[, 9])
  plantcare_data$V9 <- ifelse(grepl("anoxic", plantcare_data[, 8]), "anoxic", plantcare_data[, 9])
  plantcare_data$V9 <- ifelse(grepl("dehydration, low-temp, salt stresses", plantcare_data[, 8]), "mixed stress", plantcare_data[, 9])
  plantcare_data$V10 <- ifelse(grepl("light respon|low-temp|heat|circadian|defense and stress responsivenes|wound|drought|anaerobic|anoxic", plantcare_data[, 8]), "Environment", plantcare_data[, 10])
  
  plantcare_data$V9 <- ifelse(grepl("auxin", plantcare_data[, 8]), "auxin", plantcare_data[, 9])
  plantcare_data$V9 <- ifelse(grepl("gibberellin", plantcare_data[, 8]), "gibberellin", plantcare_data[, 9])
  plantcare_data$V9 <- ifelse(grepl("salicylic acid", plantcare_data[, 8]), "salicylic acid", plantcare_data[, 9])
  plantcare_data$V9 <- ifelse(grepl("abscisic acid", plantcare_data[, 8]), "abscisic acid", plantcare_data[, 9])
  plantcare_data$V9 <- ifelse(grepl("JA", plantcare_data[, 8]), "JA", plantcare_data[, 9])
  plantcare_data$V10 <- ifelse(grepl("auxin|gibberellin|salicylic acid|abscisic acid|JA", plantcare_data[, 8]), "Phytohormone", plantcare_data[, 10])
  
  plantcare_data$V9 <- ifelse(grepl("differentiation", plantcare_data[, 8]), "tissue", plantcare_data[, 9])
  plantcare_data$V9 <- ifelse(grepl("endosperm", plantcare_data[, 8]), "tissue", plantcare_data[, 9])
  plantcare_data$V9 <- ifelse(grepl("meristem", plantcare_data[, 8]), "tissue", plantcare_data[, 9])
  plantcare_data$V9 <- ifelse(grepl("metabolism", plantcare_data[, 8]), "metabolism", plantcare_data[, 9])
  plantcare_data$V9 <- ifelse(grepl("seed|root", plantcare_data[, 8]), "tissue", plantcare_data[, 9])
  plantcare_data$V9 <- ifelse(grepl("cell cycle", plantcare_data[, 8]), "cell cycle", plantcare_data[, 9])
  plantcare_data$V9 <- ifelse(grepl("biosynth", plantcare_data[, 8]), "biosynthesis", plantcare_data[, 9])
  plantcare_data$V10 <- ifelse(grepl("differentiation|endosperm|meristem|metabolism|seed|root|cell cycle|biosynth", plantcare_data[, 8]), "Growth and development", plantcare_data[, 10])
  
  plantcare_data <- plantcare_data[!is.na(plantcare_data[, 9]), ]
  
  colnames(plantcare_data)[1] <- "ID"
  colnames(plantcare_data)[2] <- "name"
  colnames(plantcare_data)[4] <- "start"
  
  colnames(plantcare_data)[9] <- "type"
  colnames(plantcare_data)[10] <- "class"
  return(plantcare_data)
}

#' Extract the location information of cis-element from Plantcare
#' 
#' @title plantcare_to_loc
#' @param plantcare_data The result of plantcare_classify().
#' @export
#' @author Shiqi Zhao
#' @return data.frame
#' @examples
#' plantcare_path <- system.file("extdata", "plantCARE_output.tab", package = "BioVizSeq") 
#' plantcare_file <- read.table(plantcare_path, header = FALSE, sep = '\t', quote="")
#' plantcare_data <- plantcare_classify(plantcare_file)
#' plantcare_loc <- plantcare_to_loc(plantcare_data)

plantcare_to_loc <- function(plantcare_data){
  element_loc <- plantcare_data[,c(1,9,4,5)]
  element_loc[,4] <- element_loc[,3] + element_loc[,4]-1
  colnames(element_loc) <- c("ID", "element", "start", "end")
  return(element_loc)
}

#' Count the number of cis element from Plantcare for heatmap
#' 
#' @title plantcare_statistic1
#' @param plantcare_data The result of plantcare_classify().
#' @export
#' @author Shiqi Zhao
#' @return data.frame
#' @examples
#' plantcare_path <- system.file("extdata", "plantCARE_output.tab", package = "BioVizSeq") 
#' plantcare_file <- read.table(plantcare_path, header = FALSE, sep = '\t', quote="")
#' plantcare_data <- plantcare_classify(plantcare_file)
#' statistic_data1 <- plantcare_statistic1(plantcare_data)
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr n
#' @importFrom dplyr left_join
#' @importFrom dplyr select
#' @importFrom tidyr complete 

plantcare_statistic1 <- function(plantcare_data){
  df <- plantcare_data
  df_combined <- df %>%
    dplyr::group_by(ID, name, class) %>% 
    dplyr::summarise(frequency = n(), .groups = 'drop') 
  
  description_df <- unique(df_combined[,2:3])
  
  df_complete <- df_combined %>%
    tidyr::complete(ID, name, fill = list(frequency = 0)) %>%
    dplyr::left_join(description_df, by = "name") %>%
    dplyr::select(-class.x)
  colnames(df_complete)[4] <- "class"
  return(df_complete)
}

#' Count the number of cis element from Plantcare for Bar chart
#' 
#' @title plantcare_statistic2
#' @param plantcare_data The result of plantcare_classify().
#' @export
#' @author Shiqi Zhao
#' @return data.frame
#' @examples
#' plantcare_path <- system.file("extdata", "plantCARE_output.tab", package = "BioVizSeq") 
#' plantcare_file <- read.table(plantcare_path, header = FALSE, sep = '\t', quote="")
#' plantcare_data <- plantcare_classify(plantcare_file)
#' statistic_data2 <- plantcare_statistic2(plantcare_data)

plantcare_statistic2 <- function(plantcare_data){
  df <- plantcare_data

  df_combined2 <- df %>%
    dplyr::group_by(ID, class) %>% 
    dplyr::summarise(frequency = n(), .groups = 'drop')
  return(df_combined2)
}


#' Visualization of cis-element in plantcare result file
#' 
#' @title plantcare_plot
#' @param plantcare_file The path of plantcare result file (.tab).
#' @param promoter_length The promoter length.
#' @param the_order The path of order file. A List of Gene ID , One ID Per Line.
#' @param shape RoundRect or Rect.
#' @param r The radius of rounded corners.
#' @param legend_size The size of legend.
#' @param element_color The color set of cis-element.
#' @export
#' @author Shiqi Zhao
#' @return p
#' @examples
#' plantcare_path <- system.file("extdata", "plantCARE_output.tab", package = "BioVizSeq") 
#' plantcare_plot(plantcare_path, promoter_length = 2000)

plantcare_plot <- function(plantcare_file, promoter_length = 2000, 
                           the_order = NULL, shape = "Rect", 
                           r = 6, legend_size= 15, element_color = NULL){
  
  plantcare_data <- read.table(plantcare_file, header = F, sep = '\t', quote="")
  
  plantcare_classify_data <- plantcare_classify(plantcare_data)
  
  element_loc <- plantcare_to_loc(plantcare_classify_data)
  
  length_data <- data.frame(ID = unique(element_loc$ID), length=promoter_length)
  
  if (is.null(the_order)) {
    the_order=NULL
  }else{
    the_order <- readLines(the_order)
  }
  
  if (is.null(element_color)) {
    element_color = NULL
  }
  
  motif_plot(element_loc, length_data, the_order = the_order, 
             shape = shape, r=r, legend_size= legend_size, 
             motif_color=element_color) +
    labs(x="", y="")
}

