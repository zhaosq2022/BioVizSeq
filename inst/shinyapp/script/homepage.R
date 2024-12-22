homepage <- dashboardPage(
  dashboardHeader(disable = T),
  dashboardSidebar(disable = T),
  dashboardBody(
    fluidRow(
      column(
        width = 10,
        offset = 1,
        tags$iframe(src = "html/README.html",
                    width = "100%",
                    height = "780px",
                    style = "border-radius: 10px; border-width: 0px;")
        )
      )
    
    )
  )
