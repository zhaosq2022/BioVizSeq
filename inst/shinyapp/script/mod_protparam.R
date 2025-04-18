mod_protparam_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      width = 4,
      h3(strong("The main options:")),
      fileInput(ns("filenamefa"),"Choose Protein File to Upload(.fa/.fasta):", accept = NULL),
      actionButton(ns("file_submit"), strong("Submit"), styleclass = "success"),
      br(),
      br(),
      h3(strong("Download options:")),
      downloadButton(ns("downloadfile"),label = "Download feature loci file")
      
    ),
    
    mainPanel(
      h3("Protein Paramter Calc result:"),
      withSpinner(DTOutput(ns("calc_result")))
    )
  )
}


mod_protparam_server <- function(input, output, session){
  ns <- session$ns

  protein_calc <- eventReactive(input$file_submit,{
    suppressMessages(ProtParam_calc(input$filenamefa$datapath))
  })
  
  output$calc_result <- renderDT({
    protein_calc()
  })
  
  output$downloadlocifile <- downloadHandler(
    filename = "protparam_calc.txt",
    content = function(file) {
      write.table(protein_calc(), file, sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
    }
  )
  
}
