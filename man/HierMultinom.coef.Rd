\name{HierMultinom.coef}
\alias{HierMultinom.coef}
\title{Extract regression coefficients from fitted model object}
\description{A function used to extract coefficient from fitted model obtained from \code{HierMultinom.path} or \code{HierMultinomOverlap.path{}. }\usage{
HierMultinom.coef(mod.fit, ind1 = NULL, ind2 = NULL)
}

\arguments{
\item{mod.fit}{A fitted model from obtained from HierMultinom.path or HierMultinomOverlap.path.}
\item{ind1}{An index for the candidate \eqn{\lambda} at which you want to extract a coefficient.}
\item{ind2}{An index for the candidate \eqn{\gamma} at which you want to extract a coefficient.}}

\value{
\item{beta.hat}{The estimate of the regression coefficient matrix on the original scale of the predictors.}
\item{intercept.hat}{The estimate of the intercept on the original scale of the   predictors.}}

