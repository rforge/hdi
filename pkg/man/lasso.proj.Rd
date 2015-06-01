\name{lasso.proj}
\alias{lasso.proj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{P-values based on lasso projection method}
\description{P-values based on lasso projection method}
\usage{
lasso.proj(x, y, family = "gaussian", standardize = TRUE,
           multiplecorr.method = "holm", N = 10000,
           parallel = FALSE, ncores = 4,
           sigma = NULL, Z = NULL, verbose = FALSE, return.Z = FALSE,
           robust = FALSE)
}
\arguments{
  \item{x}{Design matrix (without intercept).}
  \item{y}{Response vector.}
  \item{family}{family}
  \item{standardize}{Should design matrix be standardized to unit column
    standard deviation.}
  \item{multiplecorr.method}{Either "WY" or any of
    \code{\link{p.adjust.methods}}.}
  \item{N}{Number of empirical samples (only used if multiplecorr.method
  == "WY")}
  \item{parallel}{Should parallelization be used? (logical)}
  \item{ncores}{Number of cores used for parallelization.}
  \item{sigma}{Estimate of standard deviation of error term.}
  \item{Z}{user input}
  \item{verbose}{A boolean to enable reporting on the progress of the
    computations. (Only prints out information when Z is not provided by
    the user)}
  \item{return.Z}{An option to return the intermediate result which only
    depends on the design matrix x. This intermediate results
    can be used when calling the function again and the design matrix
    is the same as before.}
  \item{robust}{Uses a robust variance estimation procedure to be able
    to deal with model misspecification.}
}
\value{
  \item{pval}{Individual p-values for each parameter.}
  \item{pval.corr}{Multiple testing corrected p-values for each
    parameter.}
  \item{groupTest}{Function to perform groupwise tests. Groups are
    indicated using an index vector with entries in {1,...,p} or a list thereof.}
  \item{clusterGroupTest}{Function to perform groupwise tests based on
    hierarchical clustering. You can either provide a distance matrix
    and clustering method or the output of hierarchical clustering from
    the function \code{\link{hclust}} as for
    \code{\link{clusterGroupBound}}. P-values are adjusted for multiple testing.}
  %\item{betahat}{initial estimate by the scaled lasso of \eqn{\beta^0}}
  %\item{bhat}{de-sparsified \eqn{\beta^0} estimate used for p-value calculation}
  \item{sigmahat}{\eqn{\widehat{\sigma}} coming from the scaled lasso.}
  \item{Z}{Only different from NULL if the option return.Z is on. This
    is an intermediate result from the computation which only depends on
    the design matrix x. These are the residuals of the nodewise regressions.}
}
\references{van de Geer, S., B\"uhlmann, P., Ritov, Y. and Dezeure,
  R. (2013). On asymptotically optimal confidence regions and tests for
  high-dimensional models. Preprint.

  Zhang, C.,  Zhang, S. (2014), Confidence intervals for low dimensional
  parameters in high dimensional linear models, Journal of the Royal
  Statistical Society: Series B (Statistical Methodology).
  
  B\"uhlmann, P. and van de Geer, S. (2015). High-dimensional
  inference in misspecified linear models. Preprint arXiv:1503.06426.
}
\author{Ruben Dezeure}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
x <- matrix(rnorm(100*20), nrow = 100, ncol = 20)
y <- x[,1] + x[,2] + rnorm(100)
fit.lasso <- lasso.proj(x, y)
which(fit.lasso$pval.corr < 0.05)

fit.lasso.robust <- lasso.proj(x, y, robust = TRUE)
which(fit.lasso.robust$pval.corr < 0.05)

## Group-wise testing of the first two coefficients
fit.lasso$groupTest(1:2)

## Hierarchical testing using distance matrix based on
## correlation matrix
out.clust <- fit.lasso$clusterGroupTest()
plot(out.clust)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
