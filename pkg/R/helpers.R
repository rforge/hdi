lasso.cv <- function(x, y, ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 25 Mar 2013, 17:08
  
  fit.cv <- cv.glmnet(x, y, ...)
  ## Use default value of "lambda1.se" in cv.glmnet optimal lambda sel.
  sel <- predict(fit.cv, type = "nonzero") ## Intercept??? Exceptions???
  sel[[1]] ## ugly...
}

lasso.firstq <- function(x, y, q, ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  2 Apr 2013, 13:42

  ## Use glmnet (dfmax = q+1 because of Intercept)
  fit <- glmnet(x, y, dfmax = q, ...)   ## only need partial path
  m   <- predict(fit, type = "nonzero") ## determine non-zero coefs

  ## determine largest model that is <= q
  delta <- q - unlist(lapply(m, length)) ## deviation from desired model size
  delta[delta < 0] <- Inf ## overshooting not allowed

  take <- which.min(delta) ## takes first occurrence
  m[[take]]
}

lm.pval <- function(x, y, exact = TRUE, ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  2 Apr 2013, 11:34

  fit.lm <- lm(y ~ x, ...) ## Intercept??? Exceptions???
  fit.summary <- summary(fit.lm)

  tstat <- coef(fit.summary)[-1, "t value"] ## Intercept??? Exceptions???

  if(exact){ ## Use appropriate t-dist
    pval.sel <- 2 * pt(abs(tstat), df = fit.lm$df.residual,
                       lower.tail = FALSE)
  }else{ ## p-values based on *normal* distribution
    pval.sel <- 2 * pnorm(abs(tstat), lower.tail = FALSE)
  }
  
  names(pval.sel) <- colnames(x)
  pval.sel
}

lm.ci <- function(x, y, level = 0.95, ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 20 Feb 2014, 13:43

  fit.lm <- lm(y ~ x, ...) ## Intercept??? Exceptions???

  ci.sel <- confint(fit.lm, level = level)[-1,, drop = FALSE]
  ci.sel
}

glm.pval <- function(x, y, family = "binomial", trace = FALSE, ...)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Ruben Dezeure based on lm.pval, Date:  30 Sept 2013, 18:04
    
    fit.glm <- glm(y ~ x, family = family, ...) ## Intercept??? Exceptions???
    fit.summary <- summary(fit.glm)

    if(!fit.glm$converged & trace){ ## should be consistent with lm.pval?
      print(fit.summary)
    }
    
    pval.sel <- coef(fit.summary)[-1,4] ## dangerous with [,4]???
    
##-     if(family %in% c("poisson", "binomial")){
##-       zstat <- fit.summary$coefficients[-1, "z value"]
##-         ## Intercept??? Exceptions???
##-             
##-       ## p-values based on *normal* distribution
##-       pval.sel <- 2 * pnorm(abs(zstat), lower.tail = FALSE)
##-     }else{
##-       tstat <- fit.summary$coefficients[-1, "t value"]
##-       ## Intercept??? Exceptions???
##-       pval.sel <- 2 * pt(abs(tstat), df = fit.lm$df.residual,
##-                          lower.tail = FALSE)
##-     }
    names(pval.sel) <- colnames(x)
    pval.sel
}

fdr.adjust <- function(p)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 17 Jul 2013, 16:42

  p.fin <- p
  use <- (p < 1)
  if(any(use)){
    p.use <- p[use]

    lp <- length(p.use) ## as in p.adjust
    i <- lp:1L 
    o <- order(p.use, decreasing = TRUE)
    ro <- order(o)
    p.fin[use] <- pmin(1, cummin(p.use[o] / i))[ro]
  }
  p.fin
}

calc.ci <- function(bj, se, level = 0.95)
{
  ## Purpose:
  ## calculating confidence intervals with the given coefficients, standard
  ## errors and significance level.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 6 Feb 2014, 14:27
  
  quant <- qnorm(1 - (1 - level) / 2)
  
  lci <- bj - se * quant
  rci <- bj + se * quant
  
  return(list(lci = lci, rci = rci))
}


p.adjust.wy <- function(cov,
                        pval,
                        N=10000)
{
  ## Purpose:
  ## multiple testing correction with a Westfall young-like procedure as
  ## in ridge projection method, http://arxiv.org/abs/1202.1377 P.Buehlmann 
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## cov: covariance matrix of your estimator
  ## pval: the single testing pvalues
  ## N: the number of samples to take for the empirical distribution
  ##    which is used to correct the pvalues
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 6 Feb 2014, 14:27

  ## Simulate distribution  
  zz  <- mvrnorm(N, rep(0, ncol(cov)), cov)
  zz2 <- scale(zz, center = FALSE, scale = sqrt(diag(cov)))
  Gz  <- apply(2 * pnorm(abs(zz2),lower.tail = FALSE), 1, min)
  
  ## Corrected p-values
  pcorr <- ecdf(Gz)(pval) ## is this efficient???
  return(pcorr)
}

calculate.pvalue.for.group <- function(brescaled,
                                       group,
                                       individual,
                                       cov,
                                       N = 10000,
                                       Delta = NULL,
                                       correct = TRUE,
                                       alt = TRUE){
  ## Purpose:
  ## calculation of p-values for groups
  ## using the maximum as test statistic
  ## http://arxiv.org/abs/1202.1377 P.Buehlmann
  ## Author: Ruben Dezeure 2 May 2014

  stopifnot(is.logical(group))
  stopifnot(length(group) == length(brescaled))

  p <- length(brescaled)
  
  if(alt){ ## alt proposed by Nicolai
    pvalue <- min(individual[group])
  }else{
    ## max test statistics according to http://arxiv.org/abs/1202.1377
    ## P.Buehlmann

    ## Simulate distribution
    set.seed(3) ## always give the same result
    zz  <- mvrnorm(N, rep(0, ncol(cov)), cov)
    zz2 <- scale(zz, center = FALSE, scale = sqrt(diag(cov)))
    if(sum(group) > 1){
      group.coefficients <- abs(brescaled[group])
      max.coefficient <- max(group.coefficients)
      group.zz <- abs(zz2[,group])
        
      if(!is.null(Delta)){ ## special treatment of ridge method
        group.Delta <- Delta[group]
        group.zz <- sweep(group.zz, 2, group.Delta, "+")
      }
        
      Gz <- apply(group.zz,1,max)
      pvalue <- 1-ecdf(Gz)(max.coefficient)
    }else{ ## We don't have a group,only a single variable
      pvalue <- individual[group]
    }
  }
  
  if(correct){
    pvalue <- min(1, p * pvalue) ## Bonferroni correction
    ## We cannot use p.adjust since we need to correct for all variables p
    ## and pvalue might be a single value!
  }
  return(pvalue)
}
                                       
