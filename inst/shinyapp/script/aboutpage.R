aboutpage <- dashboardPage(
  dashboardHeader(disable = T),
  dashboardSidebar(disable = T),
  dashboardBody(
    column(
      width = 10,
      offset = 1,
      box(
        title = tagList(icon("r-project", 
                             lib = "font-awesome", class = "white-icon"), 
                        tags$span(
                          style = "color: white;", "| 1. Project Information"
                          )
                        ),
        solidHeader = TRUE, status = "danger", collapsible = TRUE,
        width = 10,
        h3(strong("BioVizSeq: an R package for visualization the element on bio-sequences.")),
        p("GitHub Repository: ",
          a("https://github.com/zhaosq2022/BioVizSeq/",
            href = "https://github.com/zhaosq2022/BioVizSeq/"),
          style = "font-weight:bold"),
        p("Website and API Documents: ",
          a("https://zhaosq2022.github.io/BioVizSeq/",
            href = "https://zhaosq2022.github.io/BioVizSeq/"),
          style = "font-weight:bold"),
        p("Maintain on CRAN: ",
          a("https://cran.r-project.org/package=BioVizSeq/",
            href = "https://cran.r-project.org/package=BioVizSeq/"),
          style = "font-weight:bold"),
      ),
      tags$style(HTML("
      .white-icon {
        color: white !important;
      }
      .box .box-title {
        color: white !important;
      }
    ")),
      box(
        title = tagList(icon("github", 
                             lib = "font-awesome", class = "white-icon"), 
                        tags$span(
                          style = "color: white;", "| 2. Author information"
                          )
                        ),
        solidHeader = TRUE, status = "danger", collapsible = TRUE,
        width = 10,
        p("Author 1: Shiqi Zhao (Zhejiang Ocean University)",
          style = "font-weight:bold"),
        p("GitHub Profile: ",
          a("https://github.com/zhaosq2022/",
            href = "https://github.com/zhaosq2022/"),
          style = "font-weight:bold"),
        p("Email: zhaosq@zjou.edu.cn",
          style = "font-weight:bold"),
      )
    )
  )
)
