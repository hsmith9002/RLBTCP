% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sep.panel.R
\name{sep.panel}
\alias{sep.panel}
\title{Transcript Categorization}
\usage{
sep.panel(data)
}
\arguments{
\item{data}{The Present absent tabel created by pres.abs.table() function}
}
\description{
This function assigns codes to transcripts based on wether they are in 1 of 4 categories.
If in at least 1 strain in RI (Known and novel datasets
analyzed separatly), then return 1, else return 5.
If in at least 1 strain in IB (Known and novel datasets
analyzed separatly), then return 1, else return 10. This
will cretae 2 columns to be summed. The sum of the columns
has meaning:
2 = present in both data sets
6 = present in IB and NOT in RI
11 = present in RI and NOT in IB
15 = not in either data set
and for the in all subgroup
2 = present in all strains
6 = present in all IB
11 = present in all RI
15 = not no strains
}
\examples{
sep.panel(PA.table)
}
