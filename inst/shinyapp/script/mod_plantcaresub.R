mod_plantcaresub_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      width = 4,
      h3(strong("The main options:")),
      fileInput(ns("filename"),"Choose Promoter File to Upload(.fa/.fasta):", accept = NULL),
      textInput(ns("emailname"),"Your email address:", value = ""),
      actionButton(ns("file_submit"), strong("Submit"), styleclass = "success"),

    ),
    
    mainPanel(
      withSpinner(textOutput(ns("sub_info")))
    )
  )
}


mod_plantcaresub_server <- function(input, output, session){
  ns <- session$ns
  
  submit_seq <- eventReactive(input$file_submit,{
    infile <- input$filename
    emailfile <- input$emailname
    if (!is.null(infile) && !is.null(emailfile)){
      
      withProgress(message = "Please be patient. It may take some time.", value = 0, {
        
        for (i in 1:5) {
          incProgress(1/5)
          Sys.sleep(0.1)
        }
        a <- suppressMessages(upload_fa_to_plantcare(infile$datapath, emailfile))
      })
    }
    "Your email will receive the result file."
  })
  
  output$sub_info <- renderText({
    submit_seq()
  })
  
}