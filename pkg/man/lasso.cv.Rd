\name{lasso.cv}
\alias{lasso.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to select predictors based on 10-fold cross-validation
  of the lasso estimator.}
\description{Function to select predictors based on 10-fold cross-validation
  of the lasso estimator.}
\usage{
lasso.cv(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Design matrix (without intercept).}
  \item{y}{Number of predictors that should be selected.}
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
sel <- lasso.cv(x, y)
sel
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}% __ONLY ONE__ keyword per line
