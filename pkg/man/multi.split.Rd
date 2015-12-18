\name{multi.split}
\alias{multi.split}
\title{Calculate P-values Based on Multi-Splitting Approach}
\description{
  Calculate p-values and confidence intervals  based on the
  multi-splitting approach
}
\usage{
multi.split(x, y, B = 100, fraction = 0.5, ci = TRUE, ci.level = 0.95,
            model.selector = lasso.cv,
            classical.fit = lm.pval, classical.ci = lm.ci,
            parallel = FALSE, ncores = 4,
            gamma = seq(ceiling(0.05*B)/B,1-1/B,by=1/B),
            args.model.selector = NULL, args.classical.fit = NULL,
            args.classical.ci = NULL,
            return.nonaggr = FALSE, return.selmodels = FALSE,
            repeat.max = 20,
            trace = FALSE)
}
\arguments{
  \item{x}{numeric design matrix (without intercept).}
  \item{y}{numeric response vector.}
  \item{B}{the number of sample-splits, a positive integer.}
  \item{fraction}{a number in \eqn{(0,1)}, the fraction of data used at
    each sample split for the model selection process.  The remaining
    data is used for calculating the p-values.}
  \item{ci}{logical indicating if a confidence interval should be
    calculated for each parameter.}
  \item{ci.level}{(if \code{ci} is true:) a number in \eqn{(0,1)},
    typically close to 1, the desired coverage level of the confidence
    intervals.}
  \item{model.selector}{a \code{\link{function}} to perform model
    selection, with default \code{\link{lasso.cv}}.  The function must have at
    least two arguments, \code{x} (the design matrix) and \code{y} (the
    response vector).  Return value is the index vector of selected columns.  See
    \code{\link{lasso.cv}} and \code{\link{lasso.firstq}} for an
    example.  Additional arguments can be passed via \code{args.model.selector}.}
  \item{classical.fit}{a \code{\link{function}} to calculate (classical)
    p-values.  Default is \code{\link{lm.pval}}.  Function must
    have at least two arguments: x (the design matrix) and y (the
    response vector).  Return value is the vector of p-values.  See
    \code{\link{lm.pval}} for an example.  Additional arguments can be
    passed through \code{args.classical.fit}.}
  \item{classical.ci}{a \code{\link{function}} to calculate (classical)
    confidence intervals.  Default is \code{\link{lm.ci}}.  Function must
    have at least 3 arguments: x (the design matrix), y (the
    response vector) and level (the coverage level).  Return value is the
    matrix of confidence intervals.  See \code{\link{lm.ci}} for an
    example.  Additional arguments can be passed through
    \code{args.classical.ci}.}
  \item{parallel}{Should parallelization be used? (logical)}
  \item{ncores}{Number of cores used for parallelization.}
  \item{gamma}{Vector of gamma-values.  In case gamma is a scalar, the
    value Q_j instead of P_j is being calculated (see reference below).}
  \item{args.model.selector}{named \code{\link{list}} of further arguments for
    function \code{model.selector}.}
  \item{args.classical.fit}{named \code{\link{list}} of further
    arguments for function \code{classical.fit}.}
  \item{args.classical.ci}{name \code{\link{list}} of further arguments
    for function \code{classical.ci}.}
  \item{return.nonaggr}{Should the unadjusted p-values be reported? (logical).}
  \item{return.selmodels}{Should the selected models (at each split) be
    reported? (logical).}
  \item{repeat.max}{positive integer indicating the maximal number of
    split trials.   Should not matter in regular cases, but necessary to
    prevent infinite loops in borderline cases.}
  \item{trace}{should information be printed out while computing? (logical).}
}
\details{...}
\value{
  \item{pval.corr}{Vector of multiple testing corrected p-values.}
  \item{gamma.min}{Value of gamma where minimal p-values was attained.}
  \item{clusterGroupTest}{Function to perform groupwise tests based on
    hierarchical clustering.  You can either provide a distance matrix
    and clustering method or the output of hierarchical clustering from
    the function \code{\link{hclust}} as for
    \code{\link{clusterGroupBound}}.  P-values are adjusted for multiple
    testing.}
}
\references{
  Meinshausen, N., Meier, L. and B\"uhlmann, P. (2009).
  P-values for high-dimensional regression.
  \emph{Journal of the American Statistical Association} \bold{104}, 1671--1681.

  %%TODO Jacopo erwaehnen
}
\author{Lukas Meier}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{lasso.cv}}, \code{\link{lasso.firstq}};
  \code{\link{lm.pval}}, \code{\link{lm.ci}}.
}

\examples{
x <- matrix(rnorm(100*1000), nrow = 100, ncol = 1000)
y <- x[,1] * 2 + x[,2] * 2.5 + rnorm(100)

## Multi-splitting with lasso.firstq as model selector function
## 'q' must be specified
fit.multi <- multi.split(x, y, model.selector =lasso.firstq,
                         args.model.selector = list(q = 10))
fit.multi
fit.multi$pval.corr[1:10] ## the first 10 p-values
}
\keyword{models}
\keyword{regression}
