\name{HierMultinom.coef}
\alias{HierMultinom.coef}
\title{HierMultinom.coef}
\description{Extract coefficient from fitted model obtained from HierMultinom.path or HierMultinomOverlap.path. }\usage{
HierMultinom.coef(mod.fit, ind1 = NULL, ind2 = NULL)
}

\arguments{
\item{mod.fit}{A fitted model from obtained from HierMultinom.path or HierMultinomOverlap.path.}
\item{ind1}{An index for the candidate \eqn{\lambda} at which you want to extract a coefficient.}
\item{ind2}{An index for the candidate \eqn{\gamma} at which you want to extract a coefficient.}}

\value{
\item{beta.hat}{The estimate of the regression coefficient matrix on the original scale of the predictors.}
\item{intercept.hat}{The estimate of the intercept on the original scale of the   predictors.}}
