##library(MASS)
##library(plus)
source("helpers.R")
source("helpers.nodewise.R")

lasso.proj <- function(x,
                       y,
                       signif.level=0.05,
                       useThetahat=FALSE,
                       multiplecorr.method="holm",
                       parallel=FALSE,
                       ncores=4,
                       sigma=NULL,##sigma estimate potentially provided by the user
                       Z=NULL)##Z Or \Thetahat potentially provided by the user
{
  ## Purpose:
  ## Depending on the argument useThetahat:
  ## An implementation of the LDPE method http://arxiv.org/abs/1110.2563 (useThetahat=FALSE)
  ## or the bias corrected projection estimator
  ## using Thetahat from eq 8 'Asymptotically optimal testing
  ## and confidence regions for high-dimensional models' http://arxiv.org/abs/1303.0518
  ## (useThetahat=TRUE)
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Return values:
  ## individual: the individual testing p-values for each parameter
  ## corrected: the multiple testing corrected p-values for each parameter
  ## betahat: the initial estimate by the scaled lasso of \beta^0
  ## bhat: the de-sparsified \beta^0 estimate used for p-value calculation
  ## sigmahat: the \sigma estimate coming from the scaled lasso
  ## ci: the confidence intervals calculated for each parameter
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 18 Oct 2013 (initial version),
  ## in part based on an implementation of the ridge projection method
  ##ridge-proj.R by P.Buehlmann + adaptations by L.Meier.

  if(useThetahat)
    {
      out <- lasso.proj.thetahat(x=x,
                                 y=y,
                                 signif.level=signif.level,
                                 multiplecorr.method=multiplecorr.method,
                                 parallel=parallel,
                                 ncores=ncores,
                                 sigma=sigma,
                                 Thetahat=Z)
    }else{
      out <- lasso.proj.Z(x=x,
                          y=y,
                          signif.level=signif.level,
                          multiplecorr.method=multiplecorr.method,
                          parallel=parallel,
                          ncores=ncores,
                          sigma=sigma,
                          Z=Z)
    }
  
  return(out)
}

lasso.proj.Z <- function(x,
                         y,
                         signif.level,
                         multiplecorr.method,
                         N=1000,
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
  ## by P.Buehlmann + adaptations by L.Meier.

  if(is.null(Z))
    {
      ##calculate Z
      nodewiselasso.out <- score.nodewiselasso(x=x,
                                               wantTheta=FALSE,
                                               parallel=parallel,
                                               ncores=ncores)
      Z <- nodewiselasso.out$out
    }
  
  n <- nrow(x)
  p <- ncol(x)
  
  ##projection estimator
  bproj <- t(Z) %*% y/n
  
  ##bias estimate
  y <- as.numeric(y)
  pluslasso <- plus(x,y,method="lasso",intercept=FALSE)
  ##take no intercept?
  scaledlassofit <- scale.lasso(pluslasso)
  betalasso <- scaledlassofit$beta
  if(is.null(sigma))
    sigmahat <- scaledlassofit$sigma
  else
    sigmahat <- sigma
  
  bias <- numeric(p)
  for(j in 1:p)
    {
      bias[j] <- bias[j] + (t(Z[,j]) %*% x[,-j]) %*% betalasso[-j]/n
    }
  
  ##substracting bias
  bproj <- bproj - bias

  ##normalise
  scaleb <- n/(sigmahat*sqrt(diag(crossprod(Z))))
  bprojrescaled <- bproj*scaleb
  
  ##calculate pvalue
  pval <- 2*pnorm(abs(bprojrescaled),lower.tail=FALSE)
  
  ##multiple testing correction
  if(multiplecorr.method == "WY")
    {##westfall-young like procedure as in ridge projection method, P.Buhlmann & L.Meier
      ##method related to the Westfall - Young procedure
      cov2 <- crossprod(Z)
      ##constants left out since we'll rescale anyway
      ##otherwise cov2 <- crossprod(Z)/n

      
      pcorr <- p.adjust.wy(cov=cov2,
                           pval=pval,
                           N=N)
    }else{
      if(multiplecorr.method %in% p.adjust.methods)
        {
          pcorr <- p.adjust(pval,method=multiplecorr.method)
        }else{
          stop("Unknown multiple correction method specified")
        }
    }
  ##end multiple testing correction

  ##Also return the confidence intervals
  se <- 1/scaleb
  myci <- calc.ci(bj=bproj,
                  signif.level=signif.level,
                  se=se)
                           
  ##return list
  list(individual = as.vector(pval),
       corrected = pcorr,
       bhat = bproj,
       betahat = betalasso,
       sigmahat=sigmahat,
       ci=cbind(myci$lci,myci$rci))
}

lasso.proj.thetahat <- function(x,
                                y,
                                signif.level,
                                multiplecorr.method,
                                N=1000,
                                parallel,
                                ncores,
                                sigma,
                                Thetahat)
{
  ## Purpose:
  ## An implementation of the bias corrected projection estimator
  ## from 'Asymptotically optimal testing and confidence regions for
  ## high-dimensional models' http://arxiv.org/abs/1303.0518
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 10 Oct 2013,
  ## in part based on an implementation of the ridge projection method
  ## ridge-proj.R by P.Buehlmann + adaptations by L.Meier.

  if(is.null(Thetahat))
     {
       ##calculate Thetahat
       nodewiselasso.out <- score.nodewiselasso(x=x,
                                                wantTheta=TRUE,
                                                parallel=parallel,
                                                ncores=ncores)
       Thetahat <- nodewiselasso.out$out
     }
  n <- nrow(x)
  p <- ncol(x)
  
  ##initial estimate
  y <- as.numeric(y)
  pluslasso <- plus(x,y,method="lasso",intercept=FALSE)
  ##take no intercept?
  scaledlassofit <- scale.lasso(pluslasso)
  betalasso <- scaledlassofit$beta
  if(is.null(sigma))
     sigmahat <- scaledlassofit$sigma
  else
     sigmahat <- sigma
  
  bproj <- betalasso + Thetahat %*% t(x) %*% (y - x %*% betalasso)/n
  
  
  ##normalise
  scaleb <- n/(sigmahat*sqrt(diag(Thetahat %*% crossprod(x) %*% t(Thetahat))))
  bprojrescaled <- bproj*scaleb
  
  ##calculate pvalue
  pval <- 2*pnorm(abs(bprojrescaled),lower.tail=FALSE)
  
  ##multiple testing correction
  if(multiplecorr.method == "WY")
    {##westfall-young like procedure as in ridge projection method, P.Buhlmann & L.Meier
      ##method related to the Westfall - Young procedure
      cov2 <- crossprod(Thetahat %*% crossprod(x) %*% t(Thetahat))
      ##constants left out since we'll rescale anyway
      ##otherwise cov2 <- crossprod(Z)/n

      pcorr <- p.adjust.wy(cov=cov2,
                           pval=pval,
                           N=N)
                           
    }else{
      if(multiplecorr.method %in% p.adjust.methods)
        {
          pcorr <- p.adjust(pval,method=multiplecorr.method)
        }else{
          stop("Unknown multiple correction method specified")
        }
    }
  ##end multiple testing correction
  
  ##Also return the confidence intervals
  se <- 1/scaleb
  myci <- calc.ci(bj=bproj,
                  signif.level=signif.level,
                  se=se)
  
  ##return list
  list(individual = as.vector(pval),
       corrected = pcorr,
       bhat = bproj,
       betahat = betalasso,
       sigmahat=sigmahat,
       ci=cbind(myci$lci,myci$rci))
}
