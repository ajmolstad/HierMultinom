\name{HierMultinom.predict}
\alias{HierMultinom.predict}
\title{HierMultinom.predict}
\description{A function to use a fitted multiresolution multinomial logistic regression model for prediction. }
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

