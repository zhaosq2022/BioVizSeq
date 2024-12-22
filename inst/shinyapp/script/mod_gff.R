mod_gff_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      width = 4,
      h3(strong("The main options:")),
      fileInput(ns("filename"),"Choose Annotation File to Upload(.gff/.gtf):", accept = NULL),
      
      actionButton(ns("file_submit"), strong("Submit"), styleclass = "success"),
      br(),
      br(),
      h3(strong("Download options:")),
      downloadButton(ns("downloadlocifile"),label = "Download feature loci file"),
      br(),
      br(),
      downloadButton(ns("downloadlengthfile"),label = "Download gene length file")
    ),
    
    mainPanel(
      h3("Loci result:"),
      withSpinner(DTOutput(ns("loci_result"))),
      h3("Gene length result:"),
      withSpinner(DTOutput(ns("aa_result")))
    )
  )
}


mod_gff_server <- function(input, output, session){
  ns <- session$ns
  
  filedata <- eventReactive(input$file_submit,{
    infile <- input$filename
    if (is.null(infile)){
      return(NULL)
    }else{
      read.table(infile$datapath, sep='\t', header = FALSE)
    }
  })
  
  motif_loc <- eventReactive(input$file_submit,{
    df <- filedata()
    gff_to_loc(df)
  })
  
  output$loci_result <- renderDT({
    motif_loc()$table_loc
  })
  output$aa_result <- renderDT({
    motif_loc()$gene_length
  })
  
  output$downloadlocifile <- downloadHandler(
    filename = "feature_loci.txt",
    content = function(file) {
      write.table(motif_loc()$table_loc, file, sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
    }
  )
  output$downloadlengthfile <- downloadHandler(
    filename = "gene_length.txt",
    content = function(file) {
      write.table(motif_loc()$gene_length, file, sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
    }
  )
  
}
