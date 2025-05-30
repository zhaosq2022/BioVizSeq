mod_advplot2_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      width = 4,
      h3(strong("The main options:")),
      radioButtons(
        inputId = ns("ordertype"),
        label = "File",
        choices = c("Tree (.newk)" = "tree", "Order (.csv/.txt)" = "orderfile"),
        selected = "tree",
        inline = TRUE 
      ),
      fileInput(ns("filetree"), label = NULL, accept = NULL),
      fileInput(ns("filegroup"),"Choose Group File to Upload(.csv/.txt/.xlsx/.xls):", accept = NULL),
      
      fileInput(ns("fileplantcare"),"Choose Plantcare File to Upload(.tab):", accept = NULL),

      fileInput(ns("filename"),"Choose Renamed File to Upload(.csv/.txt/.xlsx/.xls):", accept = NULL),

      selectInput(ns("plottype"), label = "Picture Type:",
                  c("heatmap + bar1", "heatmap + bar2")),

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
      h3("Advanced Plot:"),
      withSpinner(plotOutput(ns("plot_result"), width='110%', height='600px'))
    )
  )
}


mod_advplot2_server <- function(input, output, session){
  ns <- session$ns

  filedata <- eventReactive(input$file_submit,{
    infile <- input$filename
    if (is.null(infile)){
      return(NULL)
    }else{
      if (grepl(".xls$|.xlsx$", infile$datapath, ignore.case = TRUE)) {
        read.xlsx(infile$datapath,1, header=T)
      } else if (grepl(".csv$", infile$datapath, ignore.case = TRUE)) {
        read.csv(infile$datapath,sep=',', header=T)
      } else if(grepl(".txt$", infile$datapath, ignore.case = TRUE)){
        read.table(infile$datapath,sep = "\t", header = T,)
      }
    }
  })
  
  groupdata <- eventReactive(input$file_submit,{
    infile <- input$filegroup
    if (is.null(infile)){
      return(NULL)
    }else{
      if (grepl(".xls$|.xlsx$", infile$datapath, ignore.case = TRUE)) {
        read.xlsx(infile$datapath,1, header=T)
      } else if (grepl(".csv$", infile$datapath, ignore.case = TRUE)) {
        read.csv(infile$datapath,sep=',', header=T)
      } else if(grepl(".txt$", infile$datapath, ignore.case = TRUE)){
        read.table(infile$datapath,sep = "\t", header = T,)
      }
    }
  })
  
  doplot <- eventReactive(input$file_submit,{
    plot_do <- input$plottype
    if (plot_do == "heatmap + bar1")
      return( FALSE )
    else
      return( TRUE )
  })
  
  element_plot <- eventReactive(input$file_submit,{
    
    if(input$ordertype == "tree"){
      plot_file <- combi_p(tree_path=input$filetree$datapath, 
                           plantcare_path = input$fileplantcare$datapath,
                           renamefile = filedata(), groupfile = groupdata(),
                           promoter_length = 2000
      )
      
      p_tree <- plot_file$p_tree
      p_plantcare1 <- plot_file$p_plantcare1
      p_plantcare2 <- plot_file$p_plantcare2
      p_plantcare3 <- plot_file$p_plantcare3
      if(doplot() == FALSE){
        p_tree + p_plantcare1 + p_plantcare2 + 
          plot_layout(ncol = 3, guides = 'collect', widths = c(1, 3, 1)) + 
          plot_annotation(tag_levels = 'A')
      }else{
        p_tree + p_plantcare1 + p_plantcare3 + 
          plot_layout(ncol = 3, guides = 'collect', widths = c(1, 3, 1)) + 
          plot_annotation(tag_levels = 'A')
      }
    }else{
      plot_file <- combi_p2(order_path=input$filetree$datapath, 
                           plantcare_path = input$fileplantcare$datapath,
                           renamefile = filedata(),
                           promoter_length = 2000
      )

      p_plantcare1 <- plot_file$p_plantcare1
      p_plantcare2 <- plot_file$p_plantcare2
      p_plantcare3 <- plot_file$p_plantcare3
      p_plantcare1 <- p_plantcare1 + 
        theme(axis.text.y = element_text(color = "black", size = 10, angle = 0 ))
      
      if(doplot() == FALSE){
        p_plantcare1 + p_plantcare2 + 
          plot_layout(ncol = 2, guides = 'collect', widths = c(3, 1)) + 
          plot_annotation(tag_levels = 'A')
      }else{
        p_plantcare1 + p_plantcare3 + 
          plot_layout(ncol = 2, guides = 'collect', widths = c(3, 1)) + 
          plot_annotation(tag_levels = 'A')
      }
    }
    
  })
  
  output$plot_result <- renderPlot({
    element_plot()
  })
  
  
  output$downloadpic <- downloadHandler(
    filename = function() { 
      paste0("Advanced_plantcare", '.pdf')
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
