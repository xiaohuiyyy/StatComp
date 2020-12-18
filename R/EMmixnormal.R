#' @title EM algorithm for Gaussian Mixtures with 2 components
#' @description estimate parameters for Gaussian Mixtures

#' @param M the model of initial Gaussian Mixtures
#' @return a list containing at least the following components: 
#' 
#' \code{data} samples to induct EM algorithm on 
#' 
#' \code{lambda} weights for the two normal distributions 
#' 
#' \code{mu} mean of normal distributions 
#' 
#' \code{sigma} standard deviation of normal distributions 
#' @examples
#' \dontrun{
#' data(faithful)
#'mu0 <- list(c(2,60), c(4.5,80))
#'s0 <- list(matrix(c(1,0,0,1),nrow=2), matrix(c(1,0,0,1),nrow=2))
#'M0 <- list(data=faithful, lambda=c(0.5, 0.5), mu=mu0, sigma=s0)
#'M <- myEM(M0) 
#' }
#' @export
myEM <- function(M){
  M0 <- M
  x <- M0$data
  n <- nrow(x)
  lambda0 <- M0$lambda
  mu0 <- M0$mu
  s0 <- M0$sigma
  k <- length(lambda0)
  ## E-step: calculate the weight matrix given priori and samples
  W <- matrix(nrow=n,ncol=k)
  for(i in 1:n){
    for(j in 1:k){
      W[i,j] = mvtnorm::dmvnorm(x[i,], mu0[[j]], s0[[j]])
    }
    W[i,] <- W[i,]/sum(W[i,])
  }
  ## M-step: update the parameters
  x <- as.matrix(x)
  lambda <- apply(W, 2, mean)
  mu <- mu0; s<- s0
  for(j in 1:k){
    mu[[j]] <- as.vector( t(x)%*%W[,j]/sum(W[,j]) )
    xc <- as.matrix(x)-matrix(rep(mu[[j]],n),byrow=TRUE,nrow=n)
    s[[j]] <- t(xc)%*%diag(W[,j])%*%xc/sum(W[,j])
  }
  M <- list(data=M0$data, lambda=lambda, mu=mu, sigma=s)
  dis1 <- dis2 <- numeric(length(M$lambda))
  for (i in 1:length(M$lambda)){
    dis1[i] <- sum( (M$mu[[i]]-M0$mu[[i]])^2 )
    dis2[i] <- sum( (M$sigma[[i]]-M0$sigma[[i]])^2 )
  }
  dis <- sum( (M$lambda-M0$lambda)^2 )+ sum(dis1) + sum(dis2)
  if(dis<1e-8) return(M) # set convergence criterion 1e-8 under mydistance
  else myEM(M)
}