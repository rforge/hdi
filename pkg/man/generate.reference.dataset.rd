\name{generate.reference.dataset}
\alias{generate.reference.dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate reference datasets}
\description{The functions generates a design matrix and coefficient vector corresponding to the reference linear model datasets of the hdi paper.} 
\usage{
generate.reference.dataset(n = 100, p = 500, s0 = 3,
                           xtype, btype,
			   permuted = FALSE,iteration = 1)
}				     
\arguments{
  \item{n}{The sample size}
  \item{p}{The number of coefficients in the linear model}
  \item{s0}{The number of nonzero coefficients desired in the model}
  \item{xtype}{The type of design matrix one wants to generate. This needs to be one of 'toeplitz','equi.corr' and 'exp.decay'}
  \item{btype}{The type of coefficients one wants to generate. This needs to be one of 'U[-2,2]', 'U[0,2]', 'U[0,4]', 'bfix1', 'bfix2' and 'bfix10'.}
  \item{permuted}{If the columns of the design matrix have to be permuted or not. (logical)}
  \item{iteration}{Which of the 50 realizations of the design type and coefficients type one wants to generate.}
}
%\details{}
\value{
\item{x}{The generated design matrix}
\item{beta}{The generated coefficient vector}
}
%\references{
%}
\author{Ruben Dezeure dezeure@stat.math.ethz.ch}
%\note{}

\examples{
## Generate the first realization of the linear model with design matrix type Toeplitz and
## coefficients type uniform between -2 and 2

dataset <- generate.reference.dataset(xtype = "toeplitz", btype = "U[-2,2]", iteration = 1)
x <- dataset$x
beta <- dataset$beta

##generate 100 response vectors of this linear model
y <- x %*% beta + replicate(100,rnorm(nrow(x)))
}
\keyword{reference datasets}
\keyword{linear model}
