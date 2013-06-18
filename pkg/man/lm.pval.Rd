\name{lm.pval}
\alias{lm.pval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to calculate p-values for ordinary multiple linear regression.}
\description{Function to calculate p-values for ordinary multiple linear
  regression.}
\usage{
lm.pval(x, y, exact = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Design matrix (without intercept).}
  \item{y}{Response vector.}
  \item{exact}{Logical. TRUE if p-values based on t-distribution should
    be calculated. FALSE if normal distribution should be used.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Vector of p-values.}
\references{
%% ~put references to the literature/web site here ~
}
\author{Lukas Meier}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{hdi}}}
\examples{
x <- matrix(rnorm(100*5), nrow = 100, ncol = 5)
y <- x[,1] * 2 + x[,2] * 2.5 + rnorm(100)
pval <- lm.pval(x, y)
pval
}

\keyword{models}
\keyword{regression}% __ONLY ONE__ keyword per line
