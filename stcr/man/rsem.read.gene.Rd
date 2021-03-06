% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsem.read.gene.R
\name{rsem.read.gene}
\alias{rsem.read.gene}
\title{Reading in RSEM output - Isoform}
\usage{
rsem.read.gene(dir, dirs2, tissue, N)
}
\arguments{
\item{dir}{Global directory. Often points to project start folder. This will be the Parent folder for the folder containing lists of file names for each rsem batch output}

\item{dirs2}{Directory pointing parent folder where RSEM output is}

\item{tissue}{Character value defining name of tissue in file path}

\item{N}{Numerical vector of batch numbers}
}
\description{
This function reads in the RSEM output for each sample and merges the data into one count matrix with rows as gene IDs and columns as sample names
}
\examples{
rsem.read.gene()
}
\keyword{rsem}
