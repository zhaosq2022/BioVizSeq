mod_info_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      width = 4,
      h3(strong("The main options:")),
      fileInput(ns("filename"),"Choose File to Upload(.gff3/.gtf):", accept = NULL),
      actionButton(ns("file_submit"), strong("Submit"), styleclass = "success"),
      br(),
      br(),
      h3(strong("Download options:")),
      downloadButton(ns("downloadinfofile"),label = "Download file")
    ),
    
    mainPanel(
      h3("Gene information:"),
      withSpinner(DTOutput(ns("info_result"))),
    )
  )
}


mod_info_server <- function(input, output, session){
  ns <- session$ns
  
  filedata <- eventReactive(input$file_submit,{
    infile <- input$filename
    if (is.null(infile)){
      return(NULL)
    }
    read.table(infile$datapath, header = FALSE, sep = '\t')
  })
  
  gene_info <- eventReactive(input$file_submit,{
    df <- filedata()
    gff_statistics(df)
  })
  
  output$info_result <- renderDT({
    gene_info()
  })
  
  output$downloadlocifile <- downloadHandler(
    filename = "gene_info.txt",
    content = function(file) {
      write.table(gene_info(), file, sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
    }
  )
}
