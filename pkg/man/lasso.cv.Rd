\name{lasso.cv}
\alias{lasso.cv}
\title{Select Predictors via (10-fold) Cross-Validation of the Lasso}
\description{
  Performs (10-fold) cross-validation of the lasso (via
  \code{\link[glmnet]{cv.glmnet}}) and determines the prediction
  optimal set of parameters.
}
\usage{
lasso.cv(x, y, \dots)
}
\arguments{
  \item{x}{numeric design matrix (without intercept) of dimension \eqn{n
      \times p}{n * p}.}
  \item{y}{response vector of length \eqn{n}.}
  \item{\dots}{further arguments to be passed to
    \code{\link{cv.glmnet}}.}
}
\details{
  The function basically only calls \code{\link{cv.glmnet}}, see source
  code.
}
\value{Vector of selected predictors.}
\author{Lukas Meier}

\seealso{\code{\link{hdi}} which uses \code{lasso.cv()} by default;
  \code{\link{cv.glmnet}}.
  An alternative for \code{hdi()}: \code{\link{lasso.firstq}}.
}
\examples{
x <- matrix(rnorm(100*1000), nrow = 100, ncol = 1000)
y <- x[,1] * 2 + x[,2] * 2.5 + rnorm(100)
sel <- lasso.cv(x, y)
sel
}
\keyword{models}
\keyword{regression}
