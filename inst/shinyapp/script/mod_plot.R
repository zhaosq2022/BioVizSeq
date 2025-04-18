mod_plot_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      width = 4,
      h3(strong("The main options:")),
      fileInput(ns("filename"),"Choose Loci File to Upload(.csv/.txt/.xlsx/.xls):", accept = NULL, placeholder = "GeneID\tDomainName\tStartPos\tEndPos"),
      radioButtons(ns("uploadSource"), "Choose Length Information to Upload:",
                   choices = c("Length File(.csv/.txt/.xlsx/.xls)" = "local", "Length Number of Sequence" = "lengthnum"),
                   selected = "local"),
      uiOutput(ns("fileInputUI")),
      
      fileInput(ns("filename2"),"Choose Order File to Upload(.csv/.txt/.xlsx/.xls):", accept = NULL),
      textInput(ns("motifselect"),"Element Select:", value = "NULL"),
      selectInput(ns("shapemotif"), label = "Element Shape:",
                  c("RoundRect", "Rect")),
      selectInput(ns("motifid"), label = "Element ID:",
                  c("FALSE", "TRUE")),
      numericInput(ns("roundr"),label = "RoundRect r Value",value = 0.3),
      numericInput(ns("legendsize"),label = "Legend Size",value = 15),
      
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
      h3("Elements plot:"),
      withSpinner(plotOutput(ns("plot_result"), width='80%', height='800px'))
    )
  )
}


mod_plot_server <- function(input, output, session){
  ns <- session$ns
  output$fileInputUI <- renderUI({
    if (input$uploadSource == "lengthnum") {
      fluidRow(
        column(width = 12, numericInput(ns("lengnum"),label = NULL,value = 2000))
      )
    } else {
      fluidRow(
        column(width = 12, fileInput(ns("filename1"), label = NULL, accept = NULL))
      )
    }
  })
  
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
  
  filedata1 <- eventReactive(input$file_submit,{
    if(input$uploadSource == "local" && !is.null(input$filename1$name)){
      infile1 <- input$filename1
      if (grepl(".xls$|.xlsx$", infile1$datapath, ignore.case = TRUE)) {
        read.xlsx(infile1$datapath,1, header=T)
      } else if (grepl(".csv$", infile1$datapath, ignore.case = TRUE)) {
        read.csv(infile1$datapath,sep=',', header=T)
      } else if(grepl(".txt$", infile1$datapath, ignore.case = TRUE)){
        read.table(infile1$datapath,sep = "\t", header = T)
      }
    }else if(input$uploadSource == "lengthnum" && nzchar(input$lengnum)){
      motif_data <- filedata()
      colnames(motif_data)[1] <- "ID"
      data.frame(ID = unique(motif_data$ID), length=input$lengnum)
    }
    
  })
  
  filedata2 <- eventReactive(input$file_submit,{
    infile2 <- input$filename2
    if (is.null(infile2)){
      return(NULL)
    }else{
      readLines(infile2$datapath)
    }
  })
  
  show_motifselect <- eventReactive(input$file_submit,{
    if(input$motifselect == "NULL")
      return(NULL)
    else
      return(input$motifselect)
  })
  
  show_motifid <- eventReactive(input$file_submit,{
    if(input$motifid == "FALSE")
      return(FALSE)
    else
      return(TRUE)
  })
  
  
  element_plot <- eventReactive(input$file_submit,{
    motif_plot(filedata(), filedata1(), the_order = filedata2(),
               motif_select = show_motifselect(), shape = input$shapemotif, 
               show_motif_id = show_motifid(), r = input$roundr, 
               legend_size= input$legendsize
               )
  })

  output$plot_result <- renderPlot({
    element_plot()
  })

  
  output$downloadpic <- downloadHandler(
    filename = function() { 
      paste0("BioVizSeq_plot", '.pdf')
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
