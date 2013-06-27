\name{stability}
\alias{stability}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to perform stability selection}
\description{Function to perform stability selection}
\usage{
stability(x, y, EV, q, B = 100, fraction = 0.5,
          model.selector = lasso.firstq, args.model.selector = NULL,
          trace = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Design matrix (without intercept).}
  \item{y}{Response vector.}
  \item{EV}{Bound(s) for expected number of false positives. Can
    be a vector.}
  \item{q}{Bound for model size. Must be a scalar.}
  \item{B}{Number of sub-sample iterations.}
  \item{fraction}{Fraction of data used at each of the B sub-samples.}
  \item{model.selector}{Function to perform model selection. Default is
    \code{\link{lasso.firstq}}. User supplied function must have at
    least three arguments: x (the design matrix), y (the response
    vector) and q (the maximal model size). Return value is the index
    vector of selected columns. See \code{\link{lasso.firstq}} for an
    example. Additional arguments can be passed through
    \code{args.model.selector}.}
  \item{args.model.selector}{Named list of further arguments for
    function \code{model.selector}.}
  \item{trace}{Should information be printed out while computing (logical).}
}
\details{...}

\value{
  \item{select}{List with selected predictors
    for the supplied values of EV.}
  \item{EV}{Vector of corresponding values of EV.}
  \item{thresholds}{Thresholds corresponding to different values of EV.}
  \item{freq}{Vector of selection frequencies.}
}

\references{
  Meinshausen, N. and B\"uhlmann, P. (2010), \emph{Stability selection
  (with discussion)}, Journal of the Royal Statistical Society: Series B,
  72, 417-473.}

\author{Lukas Meier}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{...}

\examples{
x <- matrix(rnorm(100*1000), nrow = 100, ncol = 1000)
y <- x[,1] * 2 + x[,2] * 2.5 + rnorm(100)
fit.stab <- stability(x, y, EV = c(1, 2), q = 10)
fit.stab
fit.stab$freq[1:10] ## frequency of the first 10 predictors
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression} % __ONLY ONE__ keyword per line
