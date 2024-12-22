mod_smartplot_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      width = 4,
      h3(strong("The main options:")),
      fileInput(ns("filename1"),"Choose Protein File to Upload(.fa/.fasta):", accept = NULL),
      
      fileInput(ns("filename2"),"Choose Order File to Upload(.csv/.txt/.xlsx/.xls):", accept = NULL),
      textInput(ns("motifselect"),"Element select:", value = "NULL"),
      selectInput(ns("shapemotif"), label = "Element shape:",
                  c("RoundRect", "Rect")),
      
      numericInput(ns("roundr"),label = "RoundRect r value",value = 0.3),
      numericInput(ns("legendsize"),label = "Legend size",value = 15),
      
      actionButton(ns("file_submit"), strong("Submit All Data"), styleclass = "success"),
      br(),
      br(),
      h3(strong("Download options:")),
      numericInput(ns("picheigh"),label = "Graph heigh value",value = 9.5),
      numericInput(ns("picwidth"),label = "Graph width value",value = 9.5),
      downloadButton(ns("downloadpic"),label = "Download Picture!"),
      
    ),
    
    mainPanel(
      h3("Domain (SMART):"),
      withSpinner(plotOutput(ns("plot_result"), width='80%', height='800px'))
    )
  )
}


mod_smartplot_server <- function(input, output, session){
  ns <- session$ns
  
  filedata1 <- eventReactive(input$file_submit,{
    infile1 <- input$filename1
    if (is.null(infile1)){
      return(NULL)
    }else{
      infile1$datapath
    }
  })
  
  filedata2 <- eventReactive(input$file_submit,{
    infile2 <- input$filename2
    if (is.null(infile2)){
      return(NULL)
    }else{
      infile2$datapath
    }
  })
  
  show_motifselect <- eventReactive(input$file_submit,{
    if(input$motifselect == "NULL")
      return(NULL)
    else
      return(input$motifselect)
  })
  

  
  element_plot <- eventReactive(input$file_submit,{
    seqs <- seqinr::read.fasta(filedata1(), seqtype = "AA", as.string = TRUE, whole.header = TRUE)
    time_sub <- length(seqs)*8
    
    withProgress(message = paste('Please be patient. It may take', time_sub, "seconds!"), value = 0, {
      
      for (i in 1:5) {
        incProgress(1/5)
        Sys.sleep(0.1)
      }
      suppressMessages(smart_plot(filedata1(), the_order = filedata2(),
                                  domain_select = show_motifselect(), shape = input$shapemotif, 
                                  r = input$roundr, 
                                  legend_size= input$legendsize
      ))
    })
    
  })
  
  output$plot_result <- renderPlot({
    element_plot()
  })
  
  
  output$downloadpic <- downloadHandler(
    filename = function() { 
      paste0("Smart_plot", '.pdf')
    },
    contentType = "image/pdf",
    content = function(file) {
      pdf(file, width = input$picwidth, height = input$picheigh)
      print(element_plot())
      dev.off()
    }
  )
  
}
