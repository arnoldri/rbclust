#include <Rcpp.h>
using namespace Rcpp;

///////////////////////////////////////////////////////////
//' C++ version of log L of row clustering for binary data
//' @export
// [[Rcpp::export]]
double llikef_rbclust_c(List model,
                        int n,
                        DataFrame long_df,
                        NumericMatrix xmat_df) {
  double llval = 0;
  int m = model["m"];
  int R = model["R"];
  int np = model["np"];
  String responsename = model["responsename"];
  NumericVector y_ijr = as<NumericVector>(long_df[responsename]);
  NumericVector pi_r = as<NumericVector>(model["pi"]);
  NumericVector coef_k = as<NumericVector>(model["coef"]);
  double eta_i_j_r;
  NumericMatrix logf_i_j_r(m,R);
  double loga_i_r;
  double f_i;
  int i,j,r,k,ij,ijr;

  llval = 0;
  for(i=0; i<n; i++) {
    for(j=0; j<m; j++) {
       ij = i + j*n;
       for(r=0; r<R; r++) {
          ijr = ij + r*(m*n);
          eta_i_j_r = 0;
          for(k=0; k<np; k++) {
            eta_i_j_r = eta_i_j_r + xmat_df(ijr,k)*coef_k[k];
          }
          logf_i_j_r(j,r) = ( y_ijr[ijr]*eta_i_j_r
                              - log(1+exp(eta_i_j_r)) );
       }
    }
    f_i = 0;
    for(r=0; r<R; r++) {
       loga_i_r = 0;
       for(j=0; j<m; j++) {
         loga_i_r = loga_i_r + logf_i_j_r(j,r);
       }
       f_i = f_i + pi_r[r]*exp(loga_i_r);
    }
    llval = llval + log(f_i);
  }

  return llval;
}

///////////////////////////////////////////////////////////
//' C++ version of the gradient of log L for
//' row clustering for binary data
//' @export
// [[Rcpp::export]]
NumericVector grad_llikef_rbclust_c(List model,
                                    int n,
                                    DataFrame long_df,
                                    NumericMatrix xmat_df) {
  int m = model["m"];
  int R = model["R"];
  int np = model["np"];
  String responsename = model["responsename"];
  NumericVector y_ijr = as<NumericVector>(long_df[responsename]);
  NumericVector pi_r = as<NumericVector>(model["pi"]);
  NumericVector coef_k = as<NumericVector>(model["coef"]);
  double eta_i_j_r;
  NumericMatrix theta_i_j_r(m,R);
  NumericMatrix logf_i_j_r(m,R);
  NumericMatrix gg_i_r_k(R,np);
  double loga_i_r;
  NumericVector a_i_r(R);
  double f_i;
  double sgr;
  int i,j,r,k,ij,ijr;
  NumericVector g_r(R);
  double b_i_k;
  NumericVector grad_llval_r(R);
  NumericVector grad_llval_k(np);
  NumericVector grad_llval(R-1+np);

  for(r=0; r<R; r++) {
    grad_llval_r[r] = 0;
    g_r[r] = 0;
  }
  for(k=0; k<np; k++) {
    grad_llval_k[k] = 0;
  }
  for(i=0; i<n; i++) {
    for(j=0; j<m; j++) {
      ij = i + j*n;
      for(r=0; r<R; r++) {
        ijr = ij + r*(m*n);
        eta_i_j_r = 0;
        for(k=0; k<np; k++) {
          eta_i_j_r = eta_i_j_r + xmat_df(ijr,k)*coef_k[k];
        }
        theta_i_j_r(j,r) = 1/(1+exp(-eta_i_j_r));
        logf_i_j_r(j,r) = ( y_ijr[ijr]*eta_i_j_r
                              - log(1+exp(eta_i_j_r)) );
      }
    }
    f_i = 0;
    for(r=0; r<R; r++) {
      loga_i_r = 0;
      for(j=0; j<m; j++) {
        loga_i_r = loga_i_r + logf_i_j_r(j,r);
      }
      a_i_r[r] = exp(loga_i_r);
      f_i = f_i + pi_r[r]*a_i_r[r];
    }
    for(r=0; r<R; r++) {
      g_r[r] = g_r[r] + a_i_r[r]/f_i;
      for(k=0; k<np; k++) {
        gg_i_r_k(r,k) = 0;
        for(j=0; j<m; j++) {
           ijr = i + j*n + r*(m*n);
           gg_i_r_k(r,k) = (gg_i_r_k(r,k) +
                      (y_ijr[ijr]-theta_i_j_r(j,r))*xmat_df(ijr,k));
        }
      }
    }
    for(k=0; k<np; k++) {
      b_i_k = 0;
      for(r=0; r<R; r++) {
        b_i_k = b_i_k + pi_r[r]*a_i_r[r]*gg_i_r_k(r,k);
      }
      grad_llval_k[k] = grad_llval_k[k] + b_i_k/f_i;
    }
    //!!==
    // if(i<5) {
    //   Rprintf("i,f_i:%d %.6f\n",i,f_i);
    //   Rprintf("a_i_r:\n");
    //   for(r=0; r<R; r++) {
    //     Rprintf("%f;",a_i_r[r]);
    //   }
    //   Rprintf("\n");
    // }
  }

  //!!==
  // for(r=0; r<R; r++) {
  //    Rprintf("r,pi_r,g_r:%d %.3f %.3f\n",r,pi_r[r],g_r[r]);
  // }

  sgr = 0;
  for(r=0; r<R; r++) {
    sgr = sgr + pi_r[r]*g_r[r];
  }
  for(r=0; r<R; r++) {
    grad_llval_r[r] = pi_r[r]*(g_r[r] - sgr);
  }

  // Pack the two pieces into the output vector
  for(r=0; r<R-1; r++) {
    grad_llval[r] = grad_llval_r[r];
  }
  for(k=0; k<np; k++) {
    grad_llval[R-1+k] = grad_llval_k[k];
  }

  return grad_llval;
}


///////////////////////////////////////////////////////////
//' C++ version of the complete data log L for
//' row clustering for binary data
//' @export
// [[Rcpp::export]]
double llikef_complete_rbclust_c(List model,
                                 int n,
                                 DataFrame long_df,
                                 NumericMatrix xmat_df,
                                 NumericMatrix zmat) {
  int m = model["m"];
  int R = model["R"];
  int np = model["np"];
  String responsename = model["responsename"];
  NumericVector y_ijr = as<NumericVector>(long_df[responsename]);
  NumericVector pi_r = as<NumericVector>(model["pi"]);
  NumericVector coef_k = as<NumericVector>(model["coef"]);
  double eta_i_j_r;
  double logf_i_j_r;
  int i,j,r,k,ijr;
  double llval;

  llval = 0;
  for(i=0; i<n; i++) {
    for(r=0; r<R; r++) {
      for(j=0; j<m; j++) {
        ijr = i + j*n + r*(m*n);
        eta_i_j_r = 0;
        for(k=0; k<np; k++) {
          eta_i_j_r = eta_i_j_r + xmat_df(ijr,k)*coef_k[k];
        }
        logf_i_j_r = ( y_ijr[ijr]*eta_i_j_r
                          - log(1+exp(eta_i_j_r)) );
        llval = llval + zmat(i,r)*logf_i_j_r;
      }
      llval = llval + zmat(i,r)*log(pi_r[r]);
    }
  }
  return llval;
}

//' C++ version of the gradient of the complete data log L
//' @export
// [[Rcpp::export]]
NumericVector grad_llikef_complete_rbclust_c(List model,
                                             int n,
                                             DataFrame long_df,
                                             NumericMatrix xmat_df,
                                             NumericMatrix zmat) {
  // grad.llikef.complete.rbclust <- function(model, n, long.df,
  //                                          xmat.df, zmat) {
  int m = model["m"];
  int R = model["R"];
  int np = model["np"];
  String responsename = model["responsename"];
  NumericVector y_ijr = as<NumericVector>(long_df[responsename]);
  NumericVector pi_r = as<NumericVector>(model["pi"]);
  NumericVector coef_k = as<NumericVector>(model["coef"]);
  double eta_i_j_r;
  double theta_i_j_r;
  NumericVector c_i_r_k(np);
  int i,j,r,k,ijr;
  NumericVector grad_llval(np);

   for(k=0; k<np; k++) {
     grad_llval[k] = 0;
   }
   for(i=0; i<n; i++) {
     for(r=0; r<R; r++) {
       for(k=0; k<np; k++) {
         c_i_r_k[k] = 0;
       }
       for(j=0; j<m; j++) {
         ijr = i + j*n + r*(m*n);
         eta_i_j_r = 0;
         for(k=0; k<np; k++) {
           eta_i_j_r = eta_i_j_r + xmat_df(ijr,k)*coef_k[k];
         }
         theta_i_j_r = 1/(1+exp(-eta_i_j_r));
         for(k=0; k<np; k++) {
           c_i_r_k[k] = (c_i_r_k[k] +
                    (y_ijr[ijr]-theta_i_j_r)*xmat_df(ijr,k));
         }
       }
       for(k=0; k<np; k++) {
         grad_llval[k] = grad_llval[k] + zmat(i,r)*c_i_r_k[k];
       }
     }
  }

  return grad_llval;
}
