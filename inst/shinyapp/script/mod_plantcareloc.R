mod_plantcareloc_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      width = 4,
      h3(strong("The main options:")),
      fileInput(ns("filename"),"Choose plantcare File to Upload(.tab):", accept = NULL),
      
      actionButton(ns("file_submit"), strong("Submit"), styleclass = "success"),
      br(),
      br(),
      h3(strong("Download options:")),
      downloadButton(ns("downloadlocifile"),label = "Download SMART loci file"),
      
    ),
    
    mainPanel(
      textOutput(ns("status_text")),
      h3("Loci result:"),
      withSpinner(DTOutput(ns("loci_result")))
    )
  )
}


mod_plantcareloc_server <- function(input, output, session){
  ns <- session$ns
  
  filedata <- eventReactive(input$file_submit,{
    infile <- input$filename
    if (is.null(infile)){
      return(NULL)
    }else{
      read.table(infile$datapath, header = FALSE, sep = '\t', quote="")
    }
  })
  
  motif_loc <- eventReactive(input$file_submit,{
    plantcare_data <- plantcare_classify(filedata())
    plantcare_to_loc(plantcare_data)
  })
  
  
  output$loci_result <- renderDT({
    motif_loc()
  })
  
  output$downloadlocifile <- downloadHandler(
    filename = "plantcare_loci.txt",
    content = function(file) {
      write.table(motif_loc(), file, sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
    }
  )

}
