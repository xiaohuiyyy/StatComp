#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector MCgnrC(int M, double x0, double sd) {
  NumericVector x(M);
  NumericVector u(M);
  int k = 0;
  double xt = 0;
  double y = 0;
  for(int i = 0; i < M; i++) { 
    u[i] = runif(1,0,1)[0];
  }
  x[0] = x0;
  for(int i = 1; i < M; i++) { 
    xt = x[i-1];
    y = rnorm(1, xt, sd)[0];
    double t1 = exp(-abs(y));
    double t2 = exp(-abs(xt));
    double t = t1/t2;
    if(u[i]<=t)
      x[i]=y;
    else{
      x[i]=x[i-1];
      k++;
    }
  }
  return(x);
}