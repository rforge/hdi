\name{multi.split}
\alias{multi.split}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to calculate p-values based on multi-splitting approach}

\description{Function to calculate p-values based on multi-splitting approach}
\usage{
multi.split(x, y, B = 50, fraction = 0.5, model.selector = lasso.cv,
            classical.fit = lm.pval, gamma = seq(0.05, 0.99, by = 0.01),
            args.model.selector = NULL, args.classical.fit = NULL, trace = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Design matrix (without intercept).}
  \item{y}{Response vector.}
  \item{B}{Number of sample-splits.}
  \item{fraction}{Fraction of data used at each of the B iterations.}
  \item{model.selector}{Function to perform model selection. Default is
    \code{\link{lasso.cv}}. Function must have at
    least two arguments: x (the design matrix) and y (the response
    vector). Return value is the index vector of selected columns. See
    \code{\link{lasso.cv}} and \code{\link{lasso.firstq}} for an 
    example. Additional arguments can be passed through
    \code{args.model.selector}.}
  \item{classical.fit}{Function to calculate (classical)
    p-values. Default is \code{\link{lm.pval}}. Function must 
    have at least two arguments: x (the design matrix) and y (the
    response vector). Return value is the vector of p-values. See
    \code{\link{lm.pval}} for an example. Additional arguments can be
    passed through \code{args.classical.fit}.}
  \item{gamma}{Vector of gamma-values.}
  \item{args.model.selector}{Named list of further arguments for
    function \code{model.selector}.}
  \item{args.classical.fit}{Named list of further arguments for function
    \code{classical.fit}.}
  \item{trace}{Should information be printed out while computing (logical).}
}
\details{...}
\value{
  \item{pval}{(only for "multi-split"). Vector of p-values.}
  \item{gamma.min}{(only for "multi-split"). Value of gamma where
    minimal p-values was attained.}
}
\references{
  Meinshausen, N., Meier, L. and B\"uhlmann, P. (2009),
  \emph{P-values for high-dimensional regression}, Journal of the
  American Statistical Association 104, 1671-1681.
}
\author{Lukas Meier}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{...}

\examples{
x <- matrix(rnorm(100*1000), nrow = 100, ncol = 1000)
y <- x[,1] * 2 + x[,2] * 2.5 + rnorm(100)

## Multi-splitting with lasso.firstq as model selector function
fit.multi <- multi.split(x, y, model.selector =lasso.firstq,
                         args.model.selector = list(q = 10))
fit.multi
fit.multi$pval[1:10] ## the first 10 p-values
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression} % __ONLY ONE__ keyword per line
