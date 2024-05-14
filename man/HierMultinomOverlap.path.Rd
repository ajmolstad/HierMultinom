\name{HierMultinomOverlap.path}
\alias{HierMultinomOverlap.path}
\title{Compute the solution path when coarse categories overlap}
\description{A function to fit the solution path and assess validation set errors for the multiresolution multinomial logistic regression estimator when groups overlap.}
\usage{
HierMultinomOverlap.path(X, Y, groups, ngamma = 100, delta = 0.01, 
  lambda.vec = 10^seq(-3, 0, length=10), tol = 1e-8, 
  max.iter = 1e4, Xval, Yval, quiet=FALSE)
}



\arguments{
\item{X}{An \eqn{n \times p} design matrix. Do not include a column for an intercept: this is added automatically.}
\item{Y}{An \eqn{n \times K} matrix of response category indicators. Each row should have one a single element equal to one; all others should be zero. See example on the GitHub page. }
\item{groups}{The list of groups, as described in the example on the GitHub page. Note that this function is to be used when groups overlap (though it will still work if they do not). }
\item{Xval}{An \eqn{n_{\rm validation} \times p} matrix of predictors for validation set.}
\item{Yval}{An \eqn{n_{\rm validation} \times K} matrix of responses for validation set. }
\item{ngamma}{Number of candidate \eqn{\gamma} tuning parameter values.}
\item{delta}{The ratio of minimum to maximum candidate \eqn{\gamma} values (i.e., \eqn{\min_k \gamma_k =  \delta \max_k \gamma_k} where \eqn{\gamma_k} is the \eqn{k}th candidate tuning parameters). Check validation errors: if still decreasing at smallest \eqn{\gamma} value, decrease \code{delta} and rerun. }
\item{tol}{Convergence tolerance. The default \eqn{10^{-7}} is stringent. }
\item{lambda.vec}{Candidate \eqn{\lambda} tuning parameter values. Default is  \eqn{(10^{-3}, 10^{-2.67}, 10^{-2.33},\dots, 10^{0})}.} 
\item{max.iter}{Maximum number of iterations for accelerated proximal gradient descent algorithm. }
\item{quiet}{Should progress be printed? }
}

\value{
\item{\code{beta.est}:}{The coefficient matrix (on the standardized predictor scale) for the tuning parameter pair which minimized the validation set deviance.}

\item{\code{beta.array}:}{An array of the estimated coefficient matrices (on the standardized predictor scale) for all tuning parameter pairs.}

\item{\code{val.errs}: }{The validation set deviance for all the tuning parameter pairs. }
\item{\code{X.train.mean}:}{The columnwise mean of the training set predictors: used for other functions.}
\item{\code{X.train.sd}:}{The columnwise standard deviation of the training set predictors: used for other functions.}
}

\example{
# HierMultinom Example
n <- 500
p <- 100
K <- 12
SigmaX <- matrix(0, nrow=p, ncol=p)
for(jj in 1:p){
  for(kk in 1:p){
    SigmaX[jj,kk] <- 0.7^abs(jj-kk)
  }
}
eo <- eigen(SigmaX) 
SigmaXsqrt <- eo$vec%*%diag(eo$val^.5)%*%t(eo$vec)
X <- matrix(rnorm(n*p), nrow=n)%*%SigmaXsqrt
Xval <- matrix(rnorm(n*p), nrow=n)%*%SigmaXsqrt
Xtest <- matrix(rnorm(1e4*p), nrow=1e4)%*%SigmaXsqrt
groups <- list(
  1:3, 4:6, 7:9, 10:12
)

s <- 6
beta <- matrix(0, nrow=p, K)
coarse.inds <- sample(1:p, s)
for(kk in 1:4){
  beta[coarse.inds, groups[[kk]]] <- rnorm(length(coarse.inds), sd = sqrt(5))%*%t(rep(1, length(groups[[kk]])))
}
fine.inds <- sample(c(1:p)[-coarse.inds], 18 - s)
for(kk in 1:4){
  for(ll in 1:length(fine.inds)){
    beta[fine.inds[ll], groups[[kk]]] <- rnorm(length(groups[[kk]]), sd = sqrt(5))
  }
}

t0 <- exp(X%*%beta)
probs <- t0/rowSums(t0)%*%t(rep(1, K))
Y <- matrix(0, nrow = n, ncol = K)
for(jj in 1:n){
  Y[jj,sample(1:K,1, prob = t0[jj,])] <- 1
}

t0 <- exp(Xval%*%beta)
probs <- t0/rowSums(t0)%*%t(rep(1, K))
Yval <- matrix(0, nrow = n, ncol = K)
for(jj in 1:n){
  Yval[jj,sample(1:K,1, prob = t0[jj,])] <- 1
}

\dontrun{

# Running the following will take a few minutes
t0 <- HierMultinom.path(X = X, Y = Y, 
                    groups = groups, 
                    ngamma = 20, delta = 0.005, 
                    tol = 1e-8, lambda.vec = 10^seq(-5,-1, length=10), 
                    max.iter = 1000, 
                    Xval = Xval, Yval = Yval, 
                    quiet = FALSE)

# get coefficients, predict at 2nd gamma, 3rd lambda
fitted.coefs.23 <- HierMultinom.coef(t0, ind1 = 2, ind2 = 3) 
Ypred.23 <- HierMultinom.predict(t0)

# get coefficients, predict  tuning parameters minimizing validation error
best.coefs <- HierMultinom.coef(t0) 
best.Ypred <- HierMultinom.predict(t0)
}
}

