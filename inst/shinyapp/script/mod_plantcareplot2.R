mod_plantcareplot2_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      width = 4,
      h3(strong("The main options:")),
      fileInput(ns("filename"),"Choose plantcare File to Upload(.tab):", accept = NULL),
      
      fileInput(ns("filename2"),"Choose Order File to Upload(.csv/.txt/.xlsx/.xls):", accept = NULL),
      selectInput(ns("plottype"), label = "Plot type:",
                  c("Heatmap", "Bar1", "Bar2")),
      
      actionButton(ns("file_submit"), strong("Run"), styleclass = "success"),
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
      h3("Regulation elements (plantcare):"),
      withSpinner(plotOutput(ns("plot_result"), width='80%', height='800px'))
    )
  )
}


mod_plantcareplot2_server <- function(input, output, session){
  ns <- session$ns
  
  filedata <- eventReactive(input$file_submit,{
    infile <- input$filename
    if (is.null(infile)){
      return(NULL)
    }else{
      infile$datapath
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

  element_plot <- eventReactive(input$file_submit,{
    if(input$plottype == "Heatmap"){
      plantcare_plot1(filedata(), the_order = filedata2())
    }else if (input$plottype == "Bar1") {
      plantcare_plot2(filedata(), the_order = filedata2())
    }else{
      plantcare_plot3(filedata(), the_order = filedata2())
    }
  })
  
  output$plot_result <- renderPlot({
    element_plot()
  })
  
  
  output$downloadpic <- downloadHandler(
    filename = function() { 
      paste0("Plantcare_plot2", '.pdf')
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
