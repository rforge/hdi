\name{lasso.firstq}
\alias{lasso.firstq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to determine the first q predictors in the lasso path.}
\description{Determines the q predictors that enter the lasso path first.}
\usage{lasso.firstq(x, y, q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Design matrix (without intercept).}
  \item{y}{Response vector.}
  \item{q}{Number of predictors that should be selected.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Vector of selected predictors.}
\references{
%% ~put references to the literature/web site here ~
}
\author{Lukas Meier}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{hdi}}.}
\examples{
x <- matrix(rnorm(100*1000), nrow = 100, ncol = 1000)
y <- x[,1] * 2 + x[,2] * 2.5 + rnorm(100)
sel <- lasso.firstq(x, y, q = 5)
sel
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}% __ONLY ONE__ keyword per line
