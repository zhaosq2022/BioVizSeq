#' BioVizSeq shiny app start function.
#' @title BioVizSeq shiny app start function.
#' 
#' @importFrom shiny runApp
#' @export
#' @author Shiqi Zhao
#' @return Shinyapp: BioVizSeq shiny app.
#' @examples
#' # 1. Library BioVizSeq package
#' library(BioVizSeq)
#'
biovizseq <- function(){
  
  app_path <- system.file("shinyapp", package = "BioVizSeq")
  
  if (is.na(app_path)) {
    stop("Shiny app folder not found in the package directory.")
  }
  
  suppressMessages(shiny::runApp(appDir = app_path))
}