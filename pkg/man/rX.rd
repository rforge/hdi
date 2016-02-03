\name{rX}
\title{Generate a design matrix \eqn{X}}
\alias{rX}
\concept{generate datasets}
\description{
  Generate a random design matrix \eqn{X} 
  useful for simulations of (high dimensional) linear models.
  This function is used in rXb to generate reference linear model
  datasets like those used in the hdi paper.
}
\usage{
rX(n, p,
   xtype,
   permuted, do2S = TRUE,
   par = switch(xtype,
                "toeplitz"  = 0.9,
                "equi.corr" = 0.8,
                "exp.decay" = c(0.4, 5)))
}
\arguments{
  \item{n}{integer; the sample size \eqn{n} (paper had always \code{n = 100}).}
  \item{p}{integer; the number of coefficients in the linear
    model. (paper had always \code{p = 500}).}
  \item{xtype}{a \code{\link{character}} string specifying the type of design matrix
    one wants to generate.  Must be one of \code{"toeplitz"},
    \code{"equi.corr"} or \code{"exp.decay"}.}
  \item{permuted}{logical specifying if the columns of the design matrix
    should be permuted.}
  \item{do2S}{logical indicating if in the case of \code{xtype}s
    \code{"toeplitz"} or \code{"equi.corr"}, the \eqn{p \times p}{p * p}
    covariance matrix should be inverted twice.  Must be true, to
    regenerate the \eqn{X} matrices from the hdi paper exactly
    \dQuote{to the last bit}.}
  \item{par}{the parameters to be used for the design matrix.  Must be
    a numeric vector of length one or two.  The default uses the
    parameters also used in the hdi paper.}
}
\details{
  \bold{Generation of the design matrix \eqn{X}:}
  \cr
  For all \code{xtype}'s, the \eqn{X} matrix will be multivariate
  normal, with mean zero and (co)variance matrix \eqn{\Sigma = }\code{C}
  determined from \code{xtype}, \code{x.par} and \eqn{p} as follows:
  \describe{
    \item{\code{xtype = "toeplitz"}:}{\code{C <- par ^ abs(toeplitz(0:(p-1)))}}
    \item{\code{xtype = "equi.corr"}:}{\eqn{\Sigma_{i,j} = \code{par}}
      for \eqn{i \ne j}{i != j}, and \eqn{= 1}
      for \eqn{i = j}, i.e., on the diagonal.}
    \item{\code{xtype = "exp.decay"}:}{\code{C <- solve(par[1] ^
	abs(toeplitz(0:(p-1)) / par[2]))}}
  }
}
\value{
  \item{x}{the generated \eqn{n \times p}{n * p} design matrix \eqn{X}.}
}
\references{
  \dQuote{The} hdi paper:\cr

}
\author{Ruben Dezeure \email{dezeure@stat.math.ethz.ch}}
\examples{
## Generate a design matrix of type "toeplitz"

x <- rX(n = 800, p = 500, xtype = "toeplitz", permuted = FALSE)

##make it reproducible
set.seed(3)
x <- rX(n = 800, p = 500, xtype = "toeplitz", permuted = FALSE)

##permute the columns
set.seed(3)
x <- rX(n = 800, p = 500, xtype = "toeplitz", permuted = TRUE)
}
\keyword{datagen}
\keyword{regression}
