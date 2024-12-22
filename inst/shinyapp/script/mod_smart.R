mod_smart_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      width = 4,
      h3(strong("The main options:")),
      fileInput(ns("filename"),"Choose Protein File to Upload(.fa/.fasta):", accept = NULL),
      
      actionButton(ns("file_submit"), strong("Submit"), styleclass = "success"),
      br(),
      br(),
      h3(strong("Download options:")),
      downloadButton(ns("downloadlocifile"),label = "Download SMART loci file"),
      br(),
      br(),
      downloadButton(ns("downloadlengthfile"),label = "Download protein length file")
    ),
    
    mainPanel(
      textOutput(ns("status_text")),
      h3("Loci result:"),
      withSpinner(DTOutput(ns("loci_result"))),
      h3("Gene length result:"),
      withSpinner(DTOutput(ns("aa_result")))
    )
  )
}


mod_smart_server <- function(input, output, session){
  ns <- session$ns
  
  motif_loc <- eventReactive(input$file_submit,{
    infile <- input$filename
    if (is.null(infile)){
      return(NULL)
    }else{
      seqs <- seqinr::read.fasta(infile$datapath, seqtype = "AA", as.string = TRUE, whole.header = TRUE)
      time_sub <- length(seqs)*8
      
      withProgress(message = paste('Please be patient. It may take', time_sub, "seconds!"), value = 0, {
        
        for (i in 1:5) {
          incProgress(1/5)
          Sys.sleep(0.1)
        }
        suppressMessages(smart_to_loc(infile$datapath))
      })
      
    }
  })

  
  
  output$loci_result <- renderDT({
    motif_loc()$table_loc
  })
  output$aa_result <- renderDT({
    motif_loc()$gene_length
  })
  
  output$downloadlocifile <- downloadHandler(
    filename = "smart_loci.txt",
    content = function(file) {
      write.table(motif_loc()$table_loc, file, sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
    }
  )
  output$downloadlengthfile <- downloadHandler(
    filename = "protein_length.txt",
    content = function(file) {
      write.table(motif_loc()$gene_length, file, sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
    }
  )
  
}
