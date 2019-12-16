#############################################################################
#' Simulate a dataset
#'
#' Generate n independent rows of a matrix specified by the
#' object \code{model}
#' @param n number of rows
#' @param model model specification: see details
#' @return a list object with two entries
#'    \code{kvec} a vector of length \code{n} containing the row
#'    group label for each row
#'    \code{ymat} an \code{n}x\code{m} matrix of binary outcomes
#'    generated in the simulation
#' @details The model specification object \code{model} contains
#'    \code{R} the number of cluster components
#'    \code{pi} the mixing probabilities -
#'      a vector of length R, adding to 1, with non-negative entries
#'    \code{theta} the probabilities of success in each row group
#'      - a vector of length R of probabilities (values between 0 and 1)
#'    \code{m} the number of columns
#'
#' @export
simdata <- function(n, model) {
    kvec <- sample(model$R, n, replace=TRUE, prob=model$pi)
    ymat <- t(sapply(model$theta[kvec],
                     function(theta) rbinom(model$m, 1, theta)))
    return(list(kvec=kvec, ymat=ymat))
}
#############################################################################
#' Log likelihood function - simple row clustering
#' without covariates
#'
#' @export
llikef <- function(ymat, model) {
    n <- nrow(ymat)
    m <- model$m
    if(m!=ncol(ymat)) {
        stop("number of columns of ymat does not match model$m")
    }
    llval <- apply(ymat,1,
                   function(yvec) {
                       vmat <- outer(1:model$R, yvec,
                                     function(rvec, yy) {
                                         dbinom(yy, size=1,
                                                prob=model$theta[rvec],
                                                log=TRUE)
                                     })
                       vvec <- apply(vmat,1,sum)
                       return( log(sum(model$pi%*%exp(vvec))))
                   })
    return(sum(llval))
}
#' Gradient of the log likelihood function - simple row clustering
#' without covariates
#'
#' @export
grad.llikef <- function(ymat, model) {
    n <- nrow(ymat)
    m <- model$m
    R <- model$R
    y.i <- apply(ymat,1,sum)
    if(m!=ncol(ymat)) {
        stop("number of columns of ymat does not match model$m")
    }
    # gradient with respect to vr (parvec[1:(R-1)])
    a.ir <- t(apply(ymat,1,
                   function(yvec) {
                       vmat <- outer(1:model$R, yvec,
                                     function(rvec, yy) {
                                         dbinom(yy, size=1,
                                                prob=model$theta[rvec],
                                                log=TRUE)
                                     })
                       vvec <- apply(vmat,1,sum)
                       return(exp(vvec))
                   }))
    f.i <- a.ir%*%model$pi
    #print(c(sum(log(f.i)), llikef(ymat, model)))##!!!==
    grad.pi.r <- as.vector(t(a.ir)%*%(1/f.i))
    grad.llval1 <- (model$pi*(grad.pi.r - sum(model$pi*grad.pi.r)))[-R]
    grad.llval2 <- as.vector(model$pi*( t(a.ir)%*%(y.i/f.i)
                                       -m*model$theta*(t(a.ir)%*%(1/f.i)) ))
    grad.llval <- c(grad.llval1, grad.llval2)
    return(grad.llval)
}

#' Load the parameters of the model with simple row clustering
#' without covariates into a vector
#'
make.parvec <- function(model) {
    # pack the parameters into a vector
    #vr <- log((model$pi[-model$R])/(model$pi[-model$R]+model$pi[2:model$R]))
    vr <- log(model$pi[-model$R]/model$pi[model$R])
    parvec <- c(vr, logit(model$theta))
    return(parvec)
}
#' Interpret the parameter vector for the model with
#' simple row clustering without covariates
#'
make.model <- function(parvec, model) {
    # unpack the parameters into the model object
    #cr <- cumprod(exp(-vr)-1)
    #model$pi <- c(1, cr)/(1+sum(cr))
    expvr <- exp(parvec[1:(model$R-1)])
    piR <- 1/(1+sum(expvr))
    model$pi <- piR*c(expvr,1)
    model$theta <- expit(parvec[model$R-1+1:model$R])
    return(model)
}

#' Maximum likelihood estimation for simple row clustering
#' without covariates - direct maximisation of the incomplete
#' data log likelihood
#'
#' @param ymat Matrix of binary data
#' @param R Number of clusters to fit
#' @param use.gradient (logical) Use the gradient of the log likelihood?
#' @param verbose (logical) Print out evaluations of the log likelihood
#'    and its derivatives?
#'
#' @return Fitted model object
#'
#' @export
fit.model <- function(ymat, R, use.gradient=TRUE, verbose=FALSE) {
    # Direct ML estimation
    if(use.gradient) {
        optim.method <- "L-BFGS-B"
    } else {
        optim.method <- "Nelder-Mead"
    }
    start.time <- Sys.time()
    n <- nrow(ymat)
    m <- ncol(ymat)
    model <- list(R=R, pi=rep(1,R)/R, theta=runif(R), m=m)
    parvec <- make.parvec(model)
    ofunc <- function(parvec, model, ymat, verbose=FALSE) {
        model <- make.model(parvec,model)
        retval <- llikef(ymat,model)
        if(verbose) {
            cat("V*"); print(c(parvec,retval),digits=4)
        }
        return(retval)
    }
    grad.ofunc <- function(parvec, model, ymat, verbose=FALSE) {
        model <- make.model(parvec,model)
        retval <- grad.llikef(ymat,model)
        if(verbose) {
            cat("G*"); print(c(parvec,retval),digits=4)
        }
        return(retval)
    }
    o1 <- optim(parvec, ofunc,
                gr=if(use.gradient) grad.ofunc else NULL,
                method=optim.method,
                control=list(fnscale=-1), hessian=TRUE,
                model=model, ymat=ymat, verbose=verbose)
    model <- make.model(o1$par, model)
    model$llike <- o1$value
    model$npar <- model$np+model$R
    model$aic <- -2*o1$value + 2*model$npar
    model$convergence <- o1$convergence
    model$se.coef <- try(sqrt(-diag(solve(o1$hessian))))
    model$time <- Sys.time()-start.time
    return(model)
}

#############################################################################
llikef.complete <- function(ymat, zmat, model) {
    n <- nrow(ymat)
    m <- model$m
    if(m!=ncol(ymat)) {
        stop("number of columns of ymat does not match model$m")
    }
    llval <- sum(zmat%*%log(model$pi))
    llval <- (llval +
              sum(zmat*
              t(apply(ymat,1,
                   function(yvec) {
                       vmat <- outer(1:model$R, yvec,
                                     function(rvec, yy) {
                                         dbinom(yy, size=1,
                                                prob=model$theta[rvec],
                                                log=TRUE)
                                     })
                       vvec <- apply(vmat,1,sum)
                       return(vvec)
                   }))))
    return(llval)
}

grad.llikef.complete <- function(ymat, zmat, model) {
    n <- nrow(ymat)
    m <- model$m
    if(m!=ncol(ymat)) {
        stop("number of columns of ymat does not match model$m")
    }
    grad.llval <- apply(ymat,1,sum)%*%zmat - m*n*model$pi*model$theta
    return(grad.llval)
}

make.parvec.em <- function(model) {
    # pack the parameters into a vector
    parvec <- logit(model$theta)
    return(parvec)
}
make.model.em <- function(parvec, model) {
    # unpack the parameters into the model object
    model$theta <- expit(parvec)
    return(model)
}

estep <- function(ymat, model) {
    zmat <- t(apply(ymat,1,
                   function(yvec) {
                       vmat <- outer(1:model$R, yvec,
                                     function(rvec, yy) {
                                         dbinom(yy, size=1,
                                                prob=model$theta[rvec],
                                                log=TRUE)
                                     })
                       vvec <- exp(apply(vmat,1,sum))
                       return(vvec)
                   }))
    zmat <- t(model$pi*t(zmat))
    zmat <- t(apply(zmat,1,function(x) x/sum(x)))
    return(zmat)
}

mstep <- function(ymat, model, zmat, use.gradient=TRUE, verbose=FALSE) {
    parvec <- make.parvec.em(model)
    if(use.gradient) {
        optim.method <- "L-BFGS-B"
    } else {
        optim.method <- "Nelder-Mead"
    }
    ofunc <- function(parvec, model, ymat, zmat, verbose=FALSE) {
        model <- make.model.em(parvec,model)
        retval <- llikef.complete(ymat,zmat,model)
        if(verbose) print(c(parvec,retval))
        return(retval)
    }
    grad.ofunc <- function(parvec, model, ymat, zmat, verbose=FALSE) {
        model <- make.model.em(parvec,model)
        retval <- grad.llikef.complete(ymat,zmat,model)
        if(verbose) print(c(parvec,retval))
        return(retval)
    }
    o1 <- optim(parvec, ofunc,
                gr=if(use.gradient) grad.ofunc else NULL,
                method=optim.method,
                control=list(fnscale=-1), hessian=FALSE,
                model=model, ymat=ymat, zmat=zmat, verbose=verbose)
    model <- make.model.em(o1$par, model)
    return(model)
}

#' Maximum likelihood estimation for simple row clustering
#' without covariates - maximisation using the EM algorithm
#'
#' @param ymat Matrix of binary data
#' @param R Number of clusters to fit
#' @param itermax Maximum number of EM iterations
#' @param aeps Absolute tolerance for convergence of parameters
#' @param use.gradient (logical) Use the gradient of the log likelihood?
#' @param verbose (logical) Print out evaluations of the log likelihood
#'    and its derivatives?
#'
#' @return Fitted model object
#'
#' @export
fit.model.em <- function(ymat, R,
                         itermin=10, itermax=100, aeps=1e-6,
                         use.gradient=TRUE, verbose=FALSE) {
    # EM estimation
    start.time <- Sys.time()
    n <- nrow(ymat)
    m <- ncol(ymat)
    model <- list(R=R, pi=rep(1,R)/R, theta=runif(R), m=m)
    iter <- 0
    finished <- FALSE
    converged <- FALSE
    while(!finished) {
       iter <- iter + 1
       old.model <- model
       # E step
       zmat <- estep(ymat, model)
       # M step
       model$pi <- apply(zmat,2,sum)/n
       model <- mstep(ymat, model, zmat, use.gradient=use.gradient)
       if(iter>=itermin) {
          if(iter>=itermax) finished <- TRUE
          if( max(abs(c(model$pi-old.model$pi,
                        model$theta-old.model$theta))) <aeps ) {
              converged <- TRUE
              finished <- TRUE
          }
       }
    }
    model$llike <- llikef(ymat, model)
    model$npar <- 2*R-1
    model$aic <- -2*model$llike + 2*model$npar
    model$time <- Sys.time()-start.time
    model$converged <- converged
    model$iter <- iter
    return(model)
}

#############################################################################
#############################################################################
#' Maximum likelihood estimation for simple row clustering
#' with covariates - using the EM algorithm
#'
#' @param formula Formula object specifying the covariates to use
#' @param R Number of clusters to fit
#' @param data \code{n} x \code{m} matrix of binary data
#' @param xr.df Optional set of row covariates
#' @param xc.df Optional set of column covariates
#' @param responsename Name for the outcome variable
#' @param factorname1 Name for the row indices
#' @param factorname2 Name for the column indices
#' @param rowclustername Name for the row clustering factor
#' @param method "DM" for direct maximisation, "EM" for the EM algorithm
#' @param use.gradient (logical) Use the gradient of the log likelihood?
#' @param use.Cpp use C++ code to evaluate the log likelihood?
#' @param verbose (logical) Print out evaluations of the log likelihood
#'    and its derivatives?
#'
#' @return Fitted model object
#'
#' @export
rbclust <- function(formula, R, data,
                    xr.df=NULL, xc.df=NULL,
                    responsename="Y",
                    factorname1="Row", factorname2="Col",
                    rowclustername="RowClust",
                    method="DM",
                    em.control=list(),
                    use.gradient=TRUE, use.Cpp=TRUE,
                    verbose=FALSE) {
    if(method=="DM") {
        model <- fit.rbclustmodel(formula=formula, R=R, data=data,
                                  xr.df=xr.df, xc.df=xc.df,
                                  responsename=responsename,
                                  factorname1=factorname1,
                                  factorname2=factorname2,
                                  rowclustername=rowclustername,
                                  use.gradient=use.gradient,
                                  use.Cpp=use.Cpp,
                                  verbose=verbose)
    } else if(method=="EM") {
        set.em.control <- list(itermin=10, itermax=100, aeps=1e-6)
        set.em.control[names(em.control)] <- em.control
        model <- fit.rbclustmodel.em(formula=formula, R=R, data=data,
                            xr.df=xr.df, xc.df=xc.df,
                            responsename=responsename,
                            factorname1=factorname1,
                            factorname2=factorname2,
                            rowclustername=rowclustername,
                            itermin=set.em.control$itermin,
                            itermax=set.em.control$itermax,
                            aeps=set.em.control$aeps,
                            use.gradient=use.gradient,
                            use.Cpp=use.Cpp,
                            verbose=verbose)
    } else {
        stop("method not recognised")
    }
    model$use.gradient <- use.gradient
    model$use.Cpp <- use.Cpp
    model$method <- method
    return(model)
}

#' @export
llikef.rbclust <- function(model, n, long.df, xmat.df) {
    # log likelihood
    m <- model$m
    R <- model$R
    np <- model$np
    eta.ijr <- as.vector(xmat.df%*%model$coef) # coef is beta.k
    theta.ijr <- expit(eta.ijr)
    logf.ijr <- long.df[,model$responsename]*eta.ijr - log(1+exp(eta.ijr))
    loga.ir <- apply(array(logf.ijr,dim=c(n,m,R)),c(1,3),sum)
    a.ir <- exp(loga.ir)
    f.i <- as.vector(a.ir%*%model$pi) # pi is pi.r
    llval <- sum(log(f.i))
    return(llval)
}

#' @export
grad.llikef.rbclust <- function(model, n, long.df, xmat.df) {
    # gradient of the log likelihood
    m <- model$m
    R <- model$R
    np <- model$np
    eta.ijr <- as.vector(xmat.df%*%model$coef) # coef is beta.k
    theta.ijr <- expit(eta.ijr)
    logf.ijr <- long.df[,model$responsename]*eta.ijr - log(1+exp(eta.ijr))
    loga.ir <- apply(array(logf.ijr,dim=c(n,m,R)),c(1,3),sum)
    a.ir <- exp(loga.ir)
    f.i <- as.vector(a.ir%*%model$pi) # pi is pi.r
    llval <- sum(log(f.i))

    g.r <- as.vector(t(a.ir)%*%(1/f.i))

    grad.llval.vr.r <- model$pi*(g.r - sum(model$pi*g.r))
    gg.irk <- apply( array((long.df[,model$responsename]-theta.ijr),
                           dim=c(n,m,R,np))
                    *array(xmat.df,dim=c(n,m,R,np)),
                    c(1,3,4), sum)
    b.ik <- apply(
       array(array(rep(model$pi,each=n),dim=c(n,R))*a.ir,dim=c(n,R,np))*gg.irk,
       c(1,3), sum)
    grad.llval.coef.k <- as.vector(t(b.ik)%*%(1/f.i))

    grad.llval <- c(grad.llval.vr.r[-R], grad.llval.coef.k)

    return(grad.llval)
}
make.parvec.rbclust <- function(model) {
    # pack the parameters into a vector
    vr <- log(model$pi[-model$R]/model$pi[model$R])
    parvec <- c(vr, model$coef)
    return(parvec)
}
make.model.rbclust <- function(parvec, model) {
    # unpack the parameters into the model object
    expvr <- exp(parvec[1:(model$R-1)])
    piR <- 1/(1+sum(expvr))
    model$pi <- piR*c(expvr,1)
    model$coef <- parvec[model$R:length(parvec)]
    return(model)
}

#' Long form data of the row clustered model with covariates
#'
#' @export
make.long.rbclustmodel <- function(formula, R, data, xr.df=NULL, xc.df=NULL,
                                   responsename="Y",
                                   factorname1="Row", factorname2="Col",
                                   rowclustername="RowClust") {
    # Make the deisgn matrix for a row clustered model for a binary data matrix
    n <- nrow(data)
    m <- ncol(data)
    ymat.df <- mat2df.rbclust(data, xr.df=xr.df, xc.df=xc.df,
                              responsename=responsename,
                              factorname1=factorname1,
                              factorname2=factorname2)
    long.df <- do.call(rbind,
                       lapply(1:R, function(r) {
                           x <- ymat.df
                           x[,rowclustername] <- r
                           return(x)
                       }))
    long.df[,rowclustername] <- factor(long.df[,rowclustername])
    return(long.df)
}

#' Design matrix of the row clustered model with covariates
#'
#' @export
make.xmat.rbclustmodel <- function(formula, R, data, xr.df=NULL, xc.df=NULL,
                                   responsename="Y",
                                   factorname1="Row", factorname2="Col",
                                   rowclustername="RowClust") {
    # Make the deisgn matrix for a row clustered model for a binary data matrix
    n <- nrow(data)
    m <- ncol(data)
    ymat.df <- mat2df.rbclust(data, xr.df=xr.df, xc.df=xc.df,
                              responsename=responsename,
                              factorname1=factorname1,
                              factorname2=factorname2)
    long.df <- do.call(rbind,
                       lapply(1:R, function(r) {
                          x <- ymat.df
                          x[,rowclustername] <- r
                          return(x)
                       }))
    long.df[,rowclustername] <- factor(long.df[,rowclustername])
    xmat.df <- model.matrix(formula, data=long.df)
    return(xmat.df)
}

#' Prediction for the row clustered model with covariates
#'
#' @export
predictmodel.rbclust <- function(model, nnew, newxr.df=NULL,
                                 responsename="Y",
                                 factorname1="Row", factorname2="Col",
                                 rowclustername="RowClust") {
    m <- model$m
    R <- model$R
    np <- model$np
    newxmat.df <- make.xmat.rbclustmodel(model$formula, model$R,
                                         data=array(1,dim=c(nnew,m)),
                                         xr.df=newxr.df, xc.df=model$xc.df,
                                         responsename=responsename,
                                         factorname1=factorname1,
                                         factorname2=factorname2,
                                         rowclustername=rowclustername)
    eta.ijr <- as.vector(newxmat.df%*%model$coef) # coef is beta.k
    theta.ijr <- expit(eta.ijr)
    dim(theta.ijr) <- c(nnew,m,R)
    return(theta.ijr)
}

#' Maximum likelihood estimation for simple row clustering
#' with covariates - direct maximisation of the incomplete
#' data log likelihood
#'
#' @param formula Formula object specifying the covariates to use
#' @param R Number of clusters to fit
#' @param data \code{n} x \code{m} matrix of binary data
#' @param xr.df Optional set of row covariates
#' @param xc.df Optional set of column covariates
#' @param responsename Name for the outcome variable
#' @param factorname1 Name for the row indices
#' @param factorname2 Name for the column indices
#' @param rowclustername Name for the row clustering factor
#' @param use.gradient (logical) Use the gradient of the log likelihood?
#' @param verbose (logical) Print out evaluations of the log likelihood
#'    and its derivatives?
#'
#' @return Fitted model object
#'
#' @export
fit.rbclustmodel <- function(formula, R, data, xr.df=NULL, xc.df=NULL,
                             responsename="Y",
                             factorname1="Row", factorname2="Col",
                             rowclustername="RowClust",
                             use.gradient=TRUE, use.Cpp=TRUE,
                             verbose=FALSE) {
    # Fit a row clustered model for a binary data matrix
    if(use.gradient) {
        optim.method <- "L-BFGS-B"
    } else {
        optim.method <- "Nelder-Mead"
    }
    start.time <- Sys.time()
    n <- nrow(data)
    m <- ncol(data)
    ymat.df <- mat2df.rbclust(data, xr.df=xr.df, xc.df=xc.df,
                              responsename=responsename,
                              factorname1=factorname1,
                              factorname2=factorname2)
    long.df <- do.call(rbind,
                       lapply(1:R, function(r) {
                          x <- ymat.df
                          x[,rowclustername] <- r
                          return(x)
                       }))
    long.df[,rowclustername] <- factor(long.df[,rowclustername])
    xmat.df <- model.matrix(formula, data=long.df)
    np <- ncol(xmat.df)

    model <- list()
    model$R <- R
    # random start
    model$pi <- runif(R,min=0.2,max=0.8)
    model$pi <- model$pi/sum(model$pi)
    model$coef <- rnorm(np)
    model$m <- m
    model$np <- np
    model$responsename <- responsename

    #llikef.rbclust(model, n, long.df, xmat.df)
    #grad.llikef.rbclust(model, n, long.df, xmat.df)
    #estep.rbclust(model, n, long.df, xmat.df)
    parvec <- make.parvec.rbclust(model)

    if(use.Cpp) {
        ofunc <- function(parvec, model, n, long.df, xmat.df, verbose=FALSE) {
            model <- make.model.rbclust(parvec,model)
            retval <- llikef_rbclust_c(model, n, long.df, xmat.df)
            if(verbose) {
                cat("V*"); print(c(parvec,retval),digits=4)
            }
            return(retval)
        }
        grad.ofunc <- function(parvec, model, n, long.df, xmat.df, verbose=FALSE) {
            model <- make.model.rbclust(parvec,model)
            retval <- grad_llikef_rbclust_c(model, n, long.df, xmat.df)
            if(verbose) {
                cat("G*"); print(c(parvec,retval),digits=4)
            }
            return(retval)
        }
    } else {
        ofunc <- function(parvec, model, n, long.df, xmat.df, verbose=FALSE) {
            model <- make.model.rbclust(parvec,model)
            retval <- llikef.rbclust(model, n, long.df, xmat.df)
            if(verbose) {
                cat("V*"); print(c(parvec,retval),digits=4)
            }
            return(retval)
        }
        grad.ofunc <- function(parvec, model, n, long.df, xmat.df, verbose=FALSE) {
            model <- make.model.rbclust(parvec,model)
            retval <- grad.llikef.rbclust(model, n, long.df, xmat.df)
            if(verbose) {
                cat("G*"); print(c(parvec,retval),digits=4)
            }
            return(retval)
        }
    }
    o1 <- optim(parvec, ofunc,
                gr=if(use.gradient) grad.ofunc else NULL,
                method=optim.method,
                control=list(fnscale=-1), hessian=TRUE,
                model=model, n=n, long.df=long.df, xmat.df=xmat.df,
                verbose=verbose)
    model <- make.model.rbclust(o1$par, model)
    names(model$coef) <- colnames(xmat.df)
    model$llike <- o1$value
    model$npar <- np + R-1
    model$aic <- -2*model$llike + 2*model$npar
    model$convergence <- o1$convergence
    model$se.coef <- try(sqrt(-diag(solve(o1$hessian))))
    model$formula <- formula
    model$xc.df <- xc.df
    model$time <- Sys.time()-start.time
    class(model) <- "rbclust"

    return(model)
}

#############################################################################
#' Complete data log L
#' @export
llikef.complete.rbclust <- function(model, n, long.df,
                                    xmat.df, zmat) {
    m <- model$m
    R <- model$R
    np <- model$np
    eta.ijr <- as.vector(xmat.df%*%model$coef) # coef is beta.k
    logf.ijr <- (long.df[,model$responsename]*eta.ijr
                 - log(1+exp(eta.ijr)))

    llval <- sum(zmat%*%log(model$pi))
    llval <- (llval + sum(zmat*apply(array(logf.ijr,dim=c(n,m,R)),
                                     c(1,3),sum)))
    return(llval)
}

#' Gradient of the complete data log L
#' @export
grad.llikef.complete.rbclust <- function(model, n, long.df,
                                         xmat.df, zmat) {
    m <- model$m
    R <- model$R
    np <- model$np
    eta.ijr <- as.vector(xmat.df%*%model$coef) # coef is beta.k
    theta.ijr <- expit(eta.ijr)
    c.irk <- apply(array(
        (long.df[,model$responsename]-theta.ijr)*xmat.df,
        dim=c(n,m,R,np)),c(1,3,4),sum)
    grad.llval <- apply(array(zmat,dim=c(n,R,np))*c.irk,3,sum)
    return(grad.llval)
}

#' @export
estep.rbclust <- function(model, n, long.df, xmat.df) {
    # posterior probabilities of allocation to each cluster
    m <- model$m
    R <- model$R
    np <- model$np
    eta.ijr <- as.vector(xmat.df%*%model$coef) # coef is beta.k
    theta.ijr <- expit(eta.ijr)
    logf.ijr <- long.df[,model$responsename]*eta.ijr - log(1+exp(eta.ijr))
    loga.ir <- apply(array(logf.ijr,dim=c(n,m,R)),c(1,3),sum)
    a.ir <- exp(loga.ir)
    zmat.ir <- array(rep(model$pi,each=n),dim=c(n,R))*a.ir
    zmat.ir <- t(apply(zmat.ir,1,function(x) x/sum(x)))
    return(zmat.ir)
}

mstep.rbclust <- function(model, n, long.df, xmat.df, zmat,
                          use.gradient=TRUE, use.Cpp=TRUE,
                          verbose=FALSE) {
    parvec <- make.parvec.complete.rbclust(model)
    if(use.gradient) {
        optim.method <- "L-BFGS-B"
    } else {
        optim.method <- "Nelder-Mead"
    }

    if(use.Cpp) {
        ofunc <- function(parvec, model, n,
                          long.df, xmat.df, zmat, verbose=FALSE) {
            model <- make.model.complete.rbclust(parvec,model)
            retval <- llikef_complete_rbclust_c(model, n, long.df,
                                                xmat.df, zmat)
            if(verbose) {
                cat("V*"); print(c(parvec,retval),digits=4)
            }
            return(retval)
        }
        grad.ofunc <- function(parvec, model, n,
                               long.df, xmat.df, zmat, verbose=FALSE) {
            model <- make.model.complete.rbclust(parvec,model)
            retval <- grad_llikef_complete_rbclust_c(model, n,
                                                     long.df,
                                                     xmat.df, zmat)
            if(verbose) {
                cat("G*"); print(c(parvec,retval),digits=4)
            }
            return(retval)
        }
    } else {
        ofunc <- function(parvec, model, n,
                          long.df, xmat.df, zmat, verbose=FALSE) {
            model <- make.model.complete.rbclust(parvec,model)
            retval <- llikef.complete.rbclust(model, n, long.df,
                                              xmat.df, zmat)
            if(verbose) {
                cat("V*"); print(c(parvec,retval),digits=4)
            }
            return(retval)
        }
        grad.ofunc <- function(parvec, model, n,
                               long.df, xmat.df, zmat, verbose=FALSE) {
            model <- make.model.complete.rbclust(parvec,model)
            retval <- grad.llikef.complete.rbclust(model, n, long.df,
                                                   xmat.df, zmat)
            if(verbose) {
                cat("G*"); print(c(parvec,retval),digits=4)
            }
            return(retval)
        }
    }
    o1 <- optim(parvec, ofunc,
                gr=if(use.gradient) grad.ofunc else NULL,
                method=optim.method,
                control=list(fnscale=-1), hessian=TRUE,
                model=model, n=n, long.df=long.df, xmat.df=xmat.df,
                zmat=zmat,
                verbose=verbose)
    model <- make.model.complete.rbclust(o1$par, model)
    names(model$coef) <- colnames(xmat.df)
    return(model)
}

make.parvec.complete.rbclust <- function(model) {
    # pack the parameters into a vector
    parvec <- model$coef
    return(parvec)
}
make.model.complete.rbclust <- function(parvec, model) {
    # unpack the parameters into the model object
    model$coef <- parvec
    return(model)
}

#' Maximum likelihood estimation for simple row clustering
#' with covariates - using the EM algorithm
#'
#' @param formula Formula object specifying the covariates to use
#' @param R Number of clusters to fit
#' @param data \code{n} x \code{m} matrix of binary data
#' @param xr.df Optional set of row covariates
#' @param xc.df Optional set of column covariates
#' @param responsename Name for the outcome variable
#' @param factorname1 Name for the row indices
#' @param factorname2 Name for the column indices
#' @param use.gradient (logical) Use the gradient of the log likelihood?
#' @param verbose (logical) Print out evaluations of the log likelihood
#'    and its derivatives?
#'
#' @return Fitted model object
#'
#' @export
fit.rbclustmodel.em <- function(formula, R, data,
                                xr.df=NULL, xc.df=NULL,
                                responsename="Y",
                                factorname1="Row", factorname2="Col",
                                rowclustername="RowClust",
                                itermin=10, itermax=100, aeps=1e-6,
                                use.gradient=TRUE, use.Cpp=TRUE,
                                verbose=FALSE) {
    # Fit a row clustered model for a binary data matrix
    start.time <- Sys.time()

    ymat.df <- mat2df.rbclust(data, xr.df=xr.df, xc.df=xc.df,
                              responsename=responsename,
                              factorname1=factorname1,
                              factorname2=factorname2)
    long.df <- do.call(rbind,
                       lapply(1:R, function(r) {
                           x <- ymat.df
                           x[,rowclustername] <- r
                           return(x)
                       }))
    long.df[,rowclustername] <- factor(long.df[,rowclustername])
    xmat.df <- model.matrix(formula, data=long.df)
    np <- ncol(xmat.df)

    model <- list()
    model$R <- R
    # random start
    model$pi <- runif(R,min=0.2,max=0.8)
    model$pi <- model$pi/sum(model$pi)
    model$coef <- rnorm(np)
    model$m <- m
    model$np <- np
    model$responsename <- responsename

    iter <- 0
    finished <- FALSE
    converged <- FALSE
    while(!finished) {
        iter <- iter + 1
        old.model <- model
        # E step
        zmat <- estep.rbclust(model, n, long.df, xmat.df)
        # M step
        model$pi <- apply(zmat,2,sum)/n
        model <- mstep.rbclust(model, n, long.df, xmat.df, zmat,
                               use.gradient=use.gradient,
                               use.Cpp=use.Cpp)
        if(iter>=itermin) {
           if(iter>=itermax) finished <- TRUE
           if( max(abs(c(model$pi-old.model$pi,
                         model$coef-old.model$coef))) <aeps ) {
               converged <- TRUE
               finished <- TRUE
           }
        }
        if(verbose) {
            print(c(iter,llikef.rbclust(model, n, long.df, xmat.df)))
        }
    }
    model$llike <- llikef.rbclust(model, n, long.df, xmat.df)
    model$npar <- np+R-1
    model$aic <- -2*model$llike + 2*model$npar
    model$time <- Sys.time()-start.time
    model$converged <- converged
    model$iter <- iter
    class(model) <- "rbclust"

    return(model)
}
