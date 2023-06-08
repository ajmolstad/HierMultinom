\name{HierMultinom.path}
\alias{HierMultinom.path}
\title{HierMultinom.path}
\description{A function to fit the solution path and assess validation set errors for the multiresolution multinomial logistic regression estimator when groups do not overlap.}
\usage{
HierMultinom.path(X, Y, groups, ngamma = 100, delta = 0.01, 
  lambda.vec = 10^seq(-3, 0, length=10), tol = 1e-8, 
  max.iter = 1e4, Xval, Yval, quiet=FALSE)
}



\arguments{
\item{X}{An \eqn{n \times p} design matrix. Do not include a column for an intercept: this is added automatically.}
\item{Y}{An \eqn{n \times K} matrix of response category indicators. Each row should have one a single element equal to one; all others should be zero. See example on the GitHub page. }
\item{groups}{The list of groups, as described in the example on the GitHub page. Note that groups must not overlap to use this function. If any groups overlap, please use HierMultinomOverlap.path. }
\item{Xval}{An \eqn{n_{\rm validation} \times p} matrix of predictors for validation set.}
\item{Yval}{An \eqn{n_{\rm validation} \times K} matrix of responses for validation set. }
For all others, see example on the GitHub page. 
}

\value{
\item{\code{beta.est}:}{The coefficient matrix (on the standardized predictor scale) for the tuning parameter pair which minimized the validation set deviance.}

\item{\code{beta.array}:}{An array of the estimated coefficient matrices (on the standardized predictor scale) for all tuning parameter pairs.}

\item{\code{val.errs}: }{The validation set deviance for all the tuning parameter pairs. }
\item{\code{X.train.mean}:}{The columnwise mean of the training set predictors: used for other functions.}
\item{\code{X.train.sd}:}{The columnwise standard deviation of the training set predictors: used for other functions.}
}

