% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FunDbSearch.R
\name{get_ms2}
\alias{get_ms2}
\title{get_ms2}
\usage{
get_ms2(lab_id, ce = c("10", "20", "40"), polarity = c("positive", "negative"))
}
\arguments{
\item{lab_id}{Dodd lab ID}

\item{ce}{collision enery. "10", "20", "40"}

\item{polarity}{ionization polarity. "positive", "negative"}
}
\description{
search compound ms2 spec via lab_id, ce, and polarity
}
\examples{
\dontrun{
get_ms2(lab_id = 'S00001', ce = '10', polarity = 'positive')
get_ms2(lab_id = 'S00001', ce = '20', polarity = 'positive')
get_ms2(lab_id = 'S00001', ce = '40', polarity = 'positive')
get_ms2(lab_id = 'S00100', ce = '20', polarity = 'negative')
get_ms2(lab_id = 'I00101', ce = '20', polarity = 'negative')
get_ms2(lab_id = 'H00101', ce = '20', polarity = 'negative')
}
}
\author{
Zhiwei Zhou
}
