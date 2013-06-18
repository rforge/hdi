\name{hdi}
\alias{hdi}
\alias{multi.split}
\alias{stability}
\title{Function to perform inference in high-dimensional (generalized) linear models}
\description{Perform inference in high-dimensional (generalized) linear
  models using various approaches.}
\usage{
hdi(x, y, method = "multi-split", B = NULL, fraction = 0.5,
    model.selector = NULL, EV = 1, q = NULL,
    gamma = seq(0.05, 0.99, by = 0.01),
    classical.fit = NULL,
    args.model.selector = NULL, args.classical.fit = NULL,
    trace = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Design matrix (without intercept).}
  \item{y}{Response vector.}
  \item{method}{Multi-splitting ("multi-split") or Stability-Selection
    ("stability").} 
  \item{B}{Number of sample-splits ("multi-split") or sub-sample
    iterations ("stability"). Default is 50 ("multi-split")or 500
    ("stability").}  
  \item{fraction}{Fraction of data used at each of the B iterations.}
  \item{model.selector}{Function to perform model selection. Default is
    \code{\link{lasso.cv}} ("multi-split") and
    \code{\link{lasso.firstq}} ("stability"). Function must have at
    least two arguments: x (the design matrix) and y (the response
    vector). Return value is the index vector of selected columns. See
    \code{\link{lasso.cv}} and \code{\link{lasso.firstq}} for
    examples. Additional arguments can be passed through
    \code{args.model.selector}.}
  \item{EV}{Bound for expected number of false positives (only for "stability").}
  \item{q}{Bound for model size (only for "stability").}
  \item{gamma}{Sequence of gamma-values (only for "multi-split")}
  \item{classical.fit}{Function to calculate (classical) p-values (only
    for "multi-split"). Default is \code{\link{lm.pval}}. Function must
    have at least two arguments: x (the design matrix) and y (the
    response vector). Return value is the vector of p-values. See
    \code{\link{lm.pval}} for an example. Additional arguments can be
    passed through \code{args.classical.fit}.}
  \item{args.model.selector}{Named list of further arguments for
    function \code{model.selector}.} 
  \item{args.classical.fit}{Named list of further arguments for function
    \code{classical.fit}.}
  \item{trace}{Should information be printed out while computing (logical).}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A list containing the p-values or the selected predictors.}
\references{
  Meinshausen, N., Meier, L. and B\"uhlmann, P. (2009),
  \emph{P-values for high-dimensional regression}, Journal of the
  American Statistical Association 104, 1671-1681.

  Meinshausen, N. and B\"uhlmann, P. (2010), \emph{Stability selection
  (with discussion)}, Journal of the Royal Statistical Society: Series B,
  72, 417-473.}
\author{Lukas Meier}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
x <- matrix(rnorm(100*1000), nrow = 100, ncol = 1000)
y <- x[,1] * 2 + x[,2] * 2.5 + rnorm(100)
fit <- hdi(x, y, method = "stability", EV = 1, q = 10)
fit
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}% __ONLY ONE__ keyword per line
