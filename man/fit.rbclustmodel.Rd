% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funcs.R
\name{fit.rbclustmodel}
\alias{fit.rbclustmodel}
\title{Maximum likelihood estimation for simple row clustering
with covariates - direct maximisation of the incomplete
data log likelihood}
\usage{
fit.rbclustmodel(
  formula,
  R,
  data,
  xr.df = NULL,
  xc.df = NULL,
  responsename = "Y",
  factorname1 = "Row",
  factorname2 = "Col",
  rowclustername = "RowClust",
  use.gradient = TRUE,
  use.Cpp = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{formula}{Formula object specifying the covariates to use}

\item{R}{Number of clusters to fit}

\item{data}{\code{n} x \code{m} matrix of binary data}

\item{xr.df}{Optional set of row covariates}

\item{xc.df}{Optional set of column covariates}

\item{responsename}{Name for the outcome variable}

\item{factorname1}{Name for the row indices}

\item{factorname2}{Name for the column indices}

\item{rowclustername}{Name for the row clustering factor}

\item{use.gradient}{(logical) Use the gradient of the log likelihood?}

\item{verbose}{(logical) Print out evaluations of the log likelihood
and its derivatives?}
}
\value{
Fitted model object
}
\description{
Maximum likelihood estimation for simple row clustering
with covariates - direct maximisation of the incomplete
data log likelihood
}
