\name{HierMultinom.predict}
\alias{HierMultinom.predict}
\title{Obtain predictions from a fitted model object}
\description{A function to use a fitted multiresolution multinomial logistic regression model for prediction.}
\usage{
HierMultinom.predict(mod.fit, Xtest, ind1 = NULL, ind2 = NULL)
}

\arguments{
\item{mod.fit}{A fitted model from obtained from HierMultinom.path or HierMultinomOverlap.path.}
\item{Xtest}{An \eqn{n_{\rm test} \times p} matrix with new samples' organized by row. }
\item{ind1}{An index for the candidate \eqn{\lambda} at which you want to make a prediction.}
\item{ind2}{An index for the candidate \eqn{\lambda} at which you want to make a prediction.}
}

\value{
\item{\code{probs}:}{The estimated probabilities for each new sample included in \code{Xtest}.  }
\item{\code{preds}: }{The predicted class for each new sample included in \code{Xtest}.}
}

\examples{
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
SigmaXsqrt <- eo$vec\%*\%diag(eo$val^.5)\%*\%t(eo$vec)
X <- matrix(rnorm(n*p), nrow=n)\%*\%SigmaXsqrt
Xval <- matrix(rnorm(n*p), nrow=n)\%*\%SigmaXsqrt
groups <- list(
  1:3, 4:6, 7:9, 10:12
)

s <- 6
beta <- matrix(0, nrow=p, K)
coarse.inds <- sample(1:p, s)
for(kk in 1:4){
  tmp0 <- t(rep(1, length(groups[[kk]])))
  tmp <- rnorm(length(coarse.inds), sd=sqrt(5))\%*\%tmp0
  beta[coarse.inds, groups[[kk]]] <- tmp
}
fine.inds <- sample(c(1:p)[-coarse.inds], 18 - s)
for(kk in 1:4){
  for(ll in 1:length(fine.inds)){
    beta[fine.inds[ll], groups[[kk]]] <- rnorm(length(groups[[kk]]), sd = sqrt(5))
  }
}

t0 <- exp(X\%*\%beta)
probs <- t0/rowSums(t0)\%*\%t(rep(1, K))
Y <- matrix(0, nrow = n, ncol = K)
for(jj in 1:n){
  Y[jj,sample(1:K,1, prob = t0[jj,])] <- 1
}

t0 <- exp(Xval\%*\%beta)
probs <- t0/rowSums(t0)\%*\%t(rep(1, K))
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
