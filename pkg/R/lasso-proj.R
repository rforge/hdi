lasso.proj <- function(x, y, ci.level = 0.95,
                       multiplecorr.method = "holm",
                       parallel = FALSE, ncores = 4,
                       sigma = NULL, ## sigma estimate provided by the user
                       Z = NULL)##Z or \Thetahat potentially provided by the user
{
  ## Purpose:
  ## An implementation of the LDPE method http://arxiv.org/abs/1110.2563
  ## which is identical to
  ## http://arxiv.org/abs/1303.0518 
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Return values:
  ## individual: p-values for every parameter (individual tests)
  ## corrected:  multiple testing corrected p-values for every parameter
  ## betahat:    initial estimate by the scaled lasso of \beta^0
  ## bhat:       de-sparsified \beta^0 estimate used for p-value calculation
  ## sigmahat:   \sigma estimate coming from the scaled lasso
  ## ci:         confidence intervals calculated for each parameter
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 18 Oct 2013 (initial version),
  ## in part based on an implementation of the ridge projection method
  ## ridge-proj.R by P. Buehlmann + adaptations by L. Meier.

  out <- lasso.proj.Z(x = x,
                      y = y,
                      ci.level = ci.level,
                      multiplecorr.method = multiplecorr.method,
                      parallel = parallel,
                      ncores = ncores,
                      sigma = sigma,
                      Z = Z)
  return(out)
}

lasso.proj.Z <- function(x,
                         y,
                         ci.level,
                         multiplecorr.method,
                         N = 10000,
                         parallel,
                         ncores,
                         sigma,
                         Z)
{
  ## Purpose:
  ## An implementation of the LDPE method http://arxiv.org/abs/1110.2563
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 29 Nov 2012 (initial version),
  ## in part based on an implementation
  ## of the ridge projection method ridge-proj.R
  ## by P. Buehlmann + adaptations by L. Meier.

  n <- nrow(x)
  p <- ncol(x)

  ## Calculate Z
  if(is.null(Z)){
    nodewiselasso.out <- score.nodewiselasso(x = x,
                                             wantTheta = FALSE,
                                             parallel = parallel,
                                             ncores = ncores)
    Z <- nodewiselasso.out$out
  }else{
    ##we assume
    if(all.equal(rep(1,p), colSums(Z*x)/n, tolerance = 10^-8))
      {
        ##all is well
      }else{
        print("Z was not properly normalized with x, we assume diag(Z^T x)/n == 1 in the code")
        print("rescaling Z ourselves")
        Z <- score.rescale(Z=Z,x=x)
      }
  }
  
  
  ## Projection estimator
  bproj <- t(Z) %*% y / n
  
  ## Bias estimate
  y <- as.numeric(y)
  
  scaledlassofit <- scalreg(X = x, y = y, lam0 = "univ")
  betalasso      <- scaledlassofit$coefficients
  
  if(is.null(sigma))
    sigmahat <- scaledlassofit$hsigma
  else
    sigmahat <- sigma
  
  bias <- numeric(p)
  for(j in 1:p){
    bias[j] <- bias[j] + (t(Z[,j]) %*% x[,-j]) %*% betalasso[-j] / n
  }
  
  ## Subtracting bias
  bproj <- bproj - bias

  ## Normalise
  scaleb        <- n / (sigmahat * sqrt(diag(crossprod(Z))))
  bprojrescaled <- bproj*scaleb
  
  ## Calculate pvalue
  pval <- 2 * pnorm(abs(bprojrescaled), lower.tail = FALSE)
  
  ## Multiple testing correction
  if(multiplecorr.method == "WY"){
    ## Westfall-Young like procedure as in ridge projection method,
    ## P.Buhlmann & L.Meier
    ##method related to the Westfall - Young procedure
    cov2 <- crossprod(Z)
    ## constants left out since we'll rescale anyway
    ## otherwise cov2 <- crossprod(Z)/n
    pcorr <- p.adjust.wy(cov = cov2, pval = pval, N = N)
    }else{
      if(multiplecorr.method %in% p.adjust.methods){
        pcorr <- p.adjust(pval,method=multiplecorr.method)
      }else{
        stop("Unknown multiple correction method specified")
      }
    }## end multiple testing correction

  ## Also return the confidence intervals
  se   <- 1 / scaleb
  myci <- calc.ci(bj = bproj,
                  level = ci.level,
                  se = se)

  ##############################################
  ## function to calculate p-value for groups ##
  ##############################################

  group.testing.function <- function(group)
    {
      calculate.pvalue.for.group(brescaled=bprojrescaled,
                                 group=group,
                                 individual=pval,
                                 cov=cov2,
                                 N=N,
                                 correct=TRUE,
                                 alt=TRUE)
    }
  
  ##return list
  list(individual = as.vector(pval),
       corrected  = pcorr,
       bhat       = bproj,
       betahat    = betalasso,
       sigmahat   = sigmahat,
       ci         = cbind(myci$lci,myci$rci),
       group.testing.function=group.testing.function)
}

