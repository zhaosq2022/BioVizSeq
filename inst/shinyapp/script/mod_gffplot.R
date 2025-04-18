mod_gffplot_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      width = 4,
      h3(strong("The main options:")),
      fileInput(ns("filename"),"Choose Annotation File to Upload(.gff/.gtf):", accept = NULL),
      
      fileInput(ns("filename2"),"Choose Order File to Upload(.csv/.txt/.xlsx/.xls):", accept = NULL),
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
      radioButtons(ns("format"), "Figure type:", 
                   choices = c("PNG" = "png", "PDF" = "pdf"), 
                   selected = "png")
    ),
    
    mainPanel(
      h3("Gene structure:"),
      withSpinner(plotOutput(ns("plot_result"), width='80%', height='800px'))
    )
  )
}


mod_gffplot_server <- function(input, output, session){
  ns <- session$ns
  
  filedata <- eventReactive(input$file_submit,{
    infile <- input$filename
    if (is.null(infile)){
      return(NULL)
    }else{
      infile$datapath
    }
  })
  

  element_plot <- eventReactive(input$file_submit,{
    gff_plot(filedata(), the_order = input$filename2$datapath,
              shape = input$shapemotif, r = input$roundr, 
              legend_size= input$legendsize
    )
  })
  
  output$plot_result <- renderPlot({
    element_plot()
  })
  
  
  output$downloadpic <- downloadHandler(
    filename = function() { 
      paste0("Gene_structure_plot", '.pdf')
    },
    content = function(file) {
      ggsave(
        filename = file,
        plot = element_plot(),
        device = input$format,
        path = NULL,
        scale = 1,
        width = input$picwidth,
        height = input$picheigh,
        units = "in",
        dpi = 400,
        limitsize = TRUE
      )
    }
  )
  
}
