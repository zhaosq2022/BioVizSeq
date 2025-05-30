% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gff_plot.R
\name{gff_statistics}
\alias{gff_statistics}
\title{gff_statistics}
\usage{
gff_statistics(gff_data, mRNA_ID = NULL)
}
\arguments{
\item{gff_data}{gff file.}

\item{mRNA_ID}{The mRNA you selected. If NULL, it means selecting all mRNAs.}
}
\value{
data.frame
}
\description{
Extract the information of element from gff or gtf file
}
\examples{

gff_path <- system.file("extdata", "idpro.gff3", package = "BioVizSeq")
gff_data <- read.table(gff_path, header = FALSE, sep = '\t')
gff_statistics_data <- gff_statistics(gff_data)

}
\author{
Shiqi Zhao
}
