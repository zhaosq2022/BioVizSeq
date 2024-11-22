##' @importFrom utils install.packages
.onAttach <- function(libname, pkgname) {
  
  required_packages <- c("ggplot2", "dplyr", "RColorBrewer", 
                         "httr", "seqinr", "stringr", "tidyr")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is not installed."))
    }
  }

  packageStartupMessage("Package ", pkgname, " loaded successfully!")
}