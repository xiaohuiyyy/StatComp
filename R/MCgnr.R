#' @title a random walk Metropolis sampler 
#' @description generate the standard Laplace distribution
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm runif

#' @param m the length of chain
#' @param x0 the initial value given to the sampler 
#' @param sd the standard deviation of the proposal normal distribution
#' @return a random sample of size \code{m}
#' @examples
#' \dontrun{
#' x <- MCgnr(1e4, 0, 1)
#' par(mfrow=c(1,1));
#' plot(x,type='l')
#' }
#' @export
MCgnr <- function(m, x0, sd){
  f <- function(x){ exp(-abs(x)) }
  x <- numeric(m)
  x[1] <- x0
  k <- 0
  u <- runif(m)
  for (i in 2:m) {
    xt <- x[i-1]
    y <- rnorm(1, mean=xt, sd=sd)
    if (u[i] <= f(y)/f(xt)) x[i] <- y else {
      x[i] <- xt
      k <- k+1 #y is rejected
    }
  }
  x
}