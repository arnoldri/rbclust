% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{mat2df.rbclust}
\alias{mat2df.rbclust}
\title{Convert an n x m array y in wide format
into a nm x 1 data frame in long format}
\usage{
mat2df.rbclust(
  y,
  ntrials = NULL,
  xr.df = NULL,
  xc.df = NULL,
  responsename = "Y",
  factorname1 = "A",
  factorname2 = "B"
)
}
\arguments{
\item{y}{\code{n} x \code{m} array of binary data}

\item{ntrials}{number of trials - not used}

\item{xr.df}{optional data frame of row covariates
- must have \code{n} rows}

\item{xc.df}{optional data frame of column covariates
- must have \code{m} rows}

\item{responsename}{name for the column containing the \code{y}
values in the output data frame}

\item{factorname1}{name for the column containing the row indices
in the output data frame}

\item{factorname2}{name for the column containing the column indices
in the output data frame}
}
\value{
A long format data frame with \code{n}\code{m} rows
  of the data and covariates
}
\description{
Convert an n x m array y in wide format
into a nm x 1 data frame in long format
}
\details{
Convert an n x m array y in wide format
into a nm x 1 data frame in long format
with additional columns from the n-row data frame of row
covariates \code{xr.df}, and the m-row data frame of column
covariates \code{xc.df}.  The names from \code{xr.df} and
\code{xc.df} are copied into the new data frame.
The column containing the data in \code{y} is named
by the parameter \code{responsename}, and there are
columns indexing the row number, named \code{factorname1}
and the column \code{factorname2}.
}
