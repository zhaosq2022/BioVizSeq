mod_cdd_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      width = 4,
      h3(strong("The main options:")),
      fileInput(ns("filename"),"Choose CDD File to Upload(.txt):", accept = NULL),
      fileInput(ns("filenamefa"),"Choose Protein File to Upload(.fa/.fasta):", accept = NULL),
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


mod_cdd_server <- function(input, output, session){
  ns <- session$ns
  
  filedata <- eventReactive(input$file_submit,{
    infile <- input$filename
    if (is.null(infile)){
      return(NULL)
    }
    readLines(infile$datapath)
  })
  
  filedatafa <- eventReactive(input$file_submit,{
    infilefa <- input$filenamefa
    if (is.null(infilefa)){
      return(NULL)
    }
    fastaleng(infilefa$datapath)
  })
  
  motif_loc <- eventReactive(input$file_submit,{
    df <- filedata()
    cdd_to_loc(df)
  })
  
  output$loci_result <- renderDT({
    motif_loc()
  })
  output$aa_result <- renderDT({
    filedatafa()
  })
  
  output$downloadlocifile <- downloadHandler(
    filename = "cdd_loci.txt",
    content = function(file) {
      write.table(motif_loc(), file, sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
    }
  )
  output$downloadlengthfile <- downloadHandler(
    filename = "protein_length.txt",
    content = function(file) {
      write.table(filedatafa(), file, sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
    }
  )
  
}
