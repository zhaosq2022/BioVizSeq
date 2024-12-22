source("global.R")

ui <- fluidPage(
  tags$script(src="www/css/addhash.js"),
  includeCSS("www/css/custom.css"),
  disconnectMessage(
    text = "Your session timed out, reload the application!",
    refresh = "Reload now",
    background = "#f89f43",
    colour = "white",
    overlayColour = "grey",
    overlayOpacity = 0.75,
    top = 250,
    refreshColour = "brown"
  ),
  navbarPage(
    title = "",
    windowTitle = "BioVizSeq",
    theme = shinytheme("flatly"),
    tabPanel("Home", homepage, icon = icon("home")),
    
    navbarMenu(
      title = "One Step plot",  icon = icon("sliders-h"),
      tabPanel("GFF/GTF Plot", mod_gffplot_ui("gffplot"), icon = icon("r-project")),
      tabPanel("MEME Plot", mod_memeplot_ui("memeplot"), icon = icon("r-project")),
      tabPanel("PFAM Plot", mod_pfamplot_ui("pfamplot"), icon = icon("r-project")),
      tabPanel("CDD Plot", mod_cddplot_ui("cddplot"), icon = icon("r-project")),
      tabPanel("SMART Plot", mod_smartplot_ui("smartplot"), icon = icon("r-project")),
      tabPanel("Plantcare Plot", mod_plantcareplot_ui("plantcareplot"), icon = icon("r-project"))
    ),
    
    navbarMenu(
      title = "Pre-processing",  icon = icon("file-lines"),
      tabPanel("GFF/GTF", mod_gff_ui("gff"),icon = icon("r-project")),
      tabPanel("MEME", mod_meme_ui("meme"), icon = icon("r-project")),
      tabPanel("PFAM", mod_pfam_ui("pfam"), icon = icon("r-project"),),
      tabPanel("CDD", mod_cdd_ui("cdd"), icon = icon("r-project"),),
      tabPanel("SMART", mod_smart_ui("smart"), icon = icon("r-project"),),
      tabPanel("Plantcare", icon = icon("r-project"),
               tabsetPanel(
                 tabPanel("1.Upload sequence", mod_plantcaresub_ui("plantcare1")),
                 tabPanel("2.Get loci", mod_plantcareloc_ui("plantcare2")),
                            )
               )
    ),
    
    tabPanel("Basic Plot", mod_plot_ui("plot"), icon = icon("bars-staggered")),
    
    tabPanel("Advance Plot", mod_advplot_ui("advplot"), icon = icon("hands-asl-interpreting")),
    
  )
)

server <- function(input, output, session) {
  callModule(mod_gff_server, "gff")
  callModule(mod_meme_server, "meme")
  callModule(mod_pfam_server, "pfam")
  callModule(mod_cdd_server, "cdd")
  callModule(mod_smart_server, "smart")
  callModule(mod_plantcaresub_server, "plantcare1")
  callModule(mod_plantcareloc_server, "plantcare2")
  callModule(mod_plot_server, "plot")
  callModule(mod_gffplot_server, "gffplot")
  callModule(mod_memeplot_server, "memeplot")
  callModule(mod_pfamplot_server, "pfamplot")
  callModule(mod_cddplot_server, "cddplot")
  callModule(mod_smartplot_server, "smartplot")
  callModule(mod_plantcareplot_server, "plantcareplot")
  callModule(mod_advplot_server, "advplot")
  
  observeEvent(input$disconnect, {
    session$close()
  })
}

shinyApp(ui, server)