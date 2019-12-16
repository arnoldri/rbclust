#############################################################################
#' logit function
#'
#' @param p A number between 0 and 1
#' @return logit(p) = \code{log(p/(1-p))}
#' @examples
#' logit(0.5)
#' logit(seq(from=0.1,to=0.9,length=8))
#' @export
logit <- function(p) log(0+p/(1-p))

#' inverse logit function
#'
#' @param x A real number
#' @return expit(x) = \code{1/(1+exp(-x))}
#' @examples
#' expit(0.5)
#' logit(seq(from=-10,to=+10,length=9))
#' @export
expit <- function(x) 1/(1+exp(-x))
#############################################################################
# from clustglm
#' Convert an n x m array y in wide format
#' into a nm x 1 data frame in long format
#'
#' @param y \code{n} x \code{m} array of binary data
#' @param ntrials number of trials - not used
#' @param xr.df optional data frame of row covariates
#'   - must have \code{n} rows
#' @param xc.df optional data frame of column covariates
#'   - must have \code{m} rows
#' @param responsename name for the column containing the \code{y}
#'   values in the output data frame
#' @param factorname1 name for the column containing the row indices
#'   in the output data frame
#' @param factorname2 name for the column containing the column indices
#'   in the output data frame
#'
#' @return A long format data frame with \code{n}\code{m} rows
#'   of the data and covariates
#'
#' @details Convert an n x m array y in wide format
#' into a nm x 1 data frame in long format
#' with additional columns from the n-row data frame of row
#' covariates \code{xr.df}, and the m-row data frame of column
#' covariates \code{xc.df}.  The names from \code{xr.df} and
#' \code{xc.df} are copied into the new data frame.
#' The column containing the data in \code{y} is named
#' by the parameter \code{responsename}, and there are
#' columns indexing the row number, named \code{factorname1}
#' and the column \code{factorname2}.
#'
#' @export
mat2df.rbclust <- function(y, ntrials = NULL,
                           xr.df = NULL, xc.df = NULL,
                           responsename = "Y",
                           factorname1 = "A", factorname2 = "B"){
    ## Used to stretch out an I*J data frame or matrix y and up to two
    ## covariate data frames into a data frame with I*J observations.
    ## If binomial data, ntrials is matrix or data frame matching
    ## y with the number of trials in each case.
    ## xr.df is a data frame of row-based covariates (eg. site factors)
    ## xc.df is a data frame of column-based covariates (eg. species traits)
    ## Lengthen the row covariates data frame:

  y <- as.matrix(y)  ## coerce to matrix, in case of dataframe

  if (!(is.data.frame(xr.df))&(!is.null(xr.df)))
      stop("Please supply the row covariates as a data frame")

  if(is.data.frame(xr.df)){
    xr.df2 <- xr.df
    for (i in 1:(ncol(y)-1))
      xr.df2 <- rbind(xr.df2,xr.df)
  }

  if (!(is.data.frame(xc.df))&(!is.null(xc.df)))
      stop("Please supply the column covariates as a data frame")

  ## For xc, construct in wrong order, then rearrange:
  if(is.data.frame(xc.df)){
    xc.df2 <- xc.df
    for (i in 1:(nrow(y)-1))
      xc.df2 <- rbind(xc.df2,xc.df)
    ## Rearrangement preserves factors, numeric
    for (j in 1:ncol(xc.df2))
      xc.df2[,j] <- rep(xc.df[,j],each=nrow(y))
  }
  ## Build new data frame:
  if (is.null(rownames(y))) rownames(y) <- as.character(1:nrow(y))
  if (is.null(colnames(y))) colnames(y) <- as.character(1:ncol(y))
  my.df <- data.frame(Y = c(y),
                      ROW = gl(nrow(y),1,nrow(y)*ncol(y),
                          labels = rownames(y)),
                      COL = gl(ncol(y),nrow(y),nrow(y)*ncol(y),
                          labels = colnames(y)))
  names(my.df)[1] <- responsename
  names(my.df)[2] <- factorname1
  names(my.df)[3] <- factorname2
  if(is.data.frame(xr.df))
    my.df <- cbind(my.df,xr.df2)
  if(is.data.frame(xc.df))
    my.df <- cbind(my.df,xc.df2)
  rownames(my.df) <- as.character(1:nrow(my.df))
  ## If binomial data:
  if (!is.null(ntrials)){
      if (is.data.frame(ntrials)) ntrials <- as.matrix(ntrials)
      if ((nrow(ntrials) != nrow(y))|(ncol(ntrials) != ncol(y)))
          stop("ntrials must be entered as a data frame or matrix matching the dimensions of y")
      my.df$ntrials <- c(ntrials)
      my.df$nfail <- my.df$ntrials - my.df[,1]
  }
  ## Return data frame:
  return(my.df)
}
#############################################################################
