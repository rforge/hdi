\name{generate.reference.dataset}
\title{Generate reference datasets}
\alias{generate.reference.dataset}
\concept{reference datasets}
\description{
  This function generates a design matrix and coefficient vector
  corresponding to the reference linear model datasets of the hdi paper.
}
\usage{
generate.reference.dataset(n, p, s0,
                           xtype = c("toeplitz", "exp.decay", "equi.corr"),
                           btype = "U[-2,2]",
			   permuted = FALSE, iteration = 1, do2S = TRUE)
}
\arguments{
  \item{n}{integer; the sample size \eqn{n} (paper had always \code{n = 100}).}
  \item{p}{integer; the number of coefficients in the linear
    model. (paper had always \code{p = 500}).}
  \item{s0}{integer number of \emph{nonzero} coefficients desired in the
    model; hence at most \code{p}.}
  \item{xtype}{a \code{\link{character}} string specifying the type of design matrix
    one wants to generate.  Must be one of \code{"toeplitz"},
    \code{"equi.corr"} or \code{"exp.decay"}.}
  \item{btype}{a \code{\link{character}} string specifying the type of
    coefficients (\dQuote{beta}) one wants to generate.  In the hdi
    paper, has been one of
    "U[-2,2]", "U[0,2]", "U[0,4]", "bfix1", "bfix2" and "bfix10".  In
    general, any string of the form \code{"U[a,b]"} or \code{"bfix<c>"}
    is allowed, where \code{a}, \code{b}, and \code{<c>} must be numbers
    (with \eqn{a \le b}{a <= b}).}
  \item{permuted}{logical specifying if the columns of the design matrix
    should be permuted.}
  \item{iteration}{integer or \code{NA} specifying which of the 50
    realizations of the design type and coefficients type one wants to
    generate.  If \code{NA}, the current \code{\link{.Random.seed}} is
    taken as usual in \R.}
  \item{do2S}{logical indicating if in the case of \code{xtype}s
    \code{"toeplitz"} or \code{"equi.corr"}, the \eqn{p \times p}{p * p}
    covariance matrix should be inverted twice.  Must be true, to
    regenerate the \eqn{X} matrices from the hdi paper exactly
    \dQuote{to the last bit}.}
}
%\details{}
\value{
  A \code{\link{list}} with components
  \item{x}{the generated \eqn{n \times p}{n * p} design matrix \eqn{X}.}
  \item{beta}{the generated coefficient vector \eqn{\beta} (\sQuote{beta}).}
}
%\references{
%  hdi paper ?
%}
\author{Ruben Dezeure \email{dezeure@stat.math.ethz.ch}}
\examples{
## Generate the first realization of the linear model with design matrix type Toeplitz and
## coefficients type uniform between -2 and 2

dset <- generate.reference.dataset(n = 80, p = 20, s0 = 3,
                                   xtype = "toeplitz", btype = "U[-2,2]")
x <- dset$x
beta <- dset$beta

## generate 100 response vectors of this linear model
y <- as.vector( x \%*\% beta ) + replicate(100, rnorm(nrow(x)))

## Use  'beta_min' fulfilling  beta's  (non standard 'btype'):
str(ds2 <- generate.reference.dataset(n = 50, p = 12, s0 = 3,
                                      xtype = "exp.decay", btype = "U[0.1, 5]"))
}
\keyword{datagen}
\keyword{regression}
