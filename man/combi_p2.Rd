% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combi_plot2.R
\name{combi_p2}
\alias{combi_p2}
\title{combi_p2}
\usage{
combi_p2(
  order_path,
  gff_path = NULL,
  meme_path = NULL,
  pfam_path = NULL,
  cdd_path = NULL,
  fa_path = NULL,
  smart_path = FALSE,
  plantcare_path = NULL,
  promoter_length = NULL,
  renamefile = NULL,
  shape = "RoundRect",
  r = 0.3,
  legend_size = 6
)
}
\arguments{
\item{order_path}{The path of order file (.txt/.csv).}

\item{gff_path}{The path of .gff/gtf file.}

\item{meme_path}{The path of .meme/mast file.}

\item{pfam_path}{The path of pfam result file (.tsv).}

\item{cdd_path}{The path of cdd result file (.txt).}

\item{fa_path}{The path of protein file (.fa/fasta).}

\item{smart_path}{Do SMART or not. (TRUE or FALSE)}

\item{plantcare_path}{The path of plantcare file (.tab).}

\item{promoter_length}{The length of promoter.}

\item{renamefile}{Rename file. Two cols: new_name and old_name.}

\item{shape}{RoundRect or Rect.}

\item{r}{The radius of rounded corners.}

\item{legend_size}{The size of legend.}
}
\value{
list
}
\description{
Get ggplot2 files to facilitate free combination in patchwork
}
\examples{
order_path <- system.file("extdata", "order.csv", package = "BioVizSeq")
gff_path <- system.file("extdata", "idpro.gff3", package = "BioVizSeq")
plot_file <- combi_p2(order_path, gff_path = gff_path)
}
\author{
Shiqi Zhao
}
