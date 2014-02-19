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

scale.lasso <- function(obj, epsilon = 1e-10, sigma.hat = 1)
{
  ## Purpose:
  ## This is the scaled lasso procedure
  ## Can modify this to accept different values for lambda0
  ## Presumes a = 0 (in Scaled Sparse Linear Regression paper, eq(3))
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## obj: any plus object, which has beta.hat calculated along lambda.path
  ## eps: stopping algorithm when change in sigma.hat < epsilon
  ## sigma.hat: start value for the estimate of sigma
  ## ----------------------------------------------------------------------
  ## Output:
  ## beta: estimated beta.hat from algorithm
  ## sigma: estimated sigma
  ## lambda: the lambda used to estimate beta, = sigma*lambda0
  ## ----------------------------------------------------------------------
  ## Authors: Tingni Sun, Cun-Hui Zhang
  ## http://arxiv.org/abs/1104.4595
  
  lambda0 = sqrt(2/nrow(obj$x)*log(ncol(obj$x)))
  while (TRUE)
    {
      beta.hat <- coef(obj,min(sigma.hat*lambda0,obj$lam.path[1]-0.0001))
      old.sigmahat <- sigma.hat
      sigma.hat <- sqrt(crossprod(obj$y-obj$x%*%beta.hat)/
                        nrow(obj$x))[1,1]
      if (abs(sigma.hat-old.sigmahat)<epsilon) break
    }
  beta.hat2 <- coef(obj,sigma.hat*lambda0/4)
  return(list("beta" = beta.hat, "beta2" = beta.hat2, "sigma"=sigma.hat,
              "lambda"=lambda0*sigma.hat))
}

lm.pval <- function(x, y, exact = TRUE, ci=FALSE,signif.level=0.05, ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  2 Apr 2013, 11:34
  ##    Updated with confidence interval calculation, Ruben Dezeure, Date: 5 Feb 2014

  fit.lm <- lm(y ~ x, ...) ## Intercept??? Exceptions???
  fit.summary <- summary(fit.lm)

  tstat <- coef(fit.summary)[-1, "t value"] ## Intercept??? Exceptions???

  if(exact){ ## Use appropriate t-dist
    pval.sel <- 2 * pt(abs(tstat), df = fit.lm$df.residual,
                       lower.tail = FALSE)
  }else{ ## p-values based on *normal* distribution
    pval.sel <- 2 * pnorm(abs(tstat), lower.tail = FALSE)
  }

  ci.sel <- confint(fit.lm)[-1,,drop=FALSE]

  names(pval.sel) <- colnames(x)
  if(ci)
    {
      ci.sel
    }else{
      pval.sel
    }
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

calc.ci <- function(bj,
                    se,
                    signif.level=0.05)
{
  ## Purpose:
  ## calculating confidence intervals with the given coefficients, standard
  ## errors and significance level.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 6 Feb 2014, 14:27
  numbse <- qnorm(signif.level/2,lower.tail=FALSE)
  
  lci <- bj - se*numbse
  rci <- bj + se*numbse
  return(list(lci=lci,rci=rci))
}


p.adjust.wy <- function(cov,
                        pval,
                        N=10000)
{
  ## Purpose:
  ## multiple testing correction with a westfall young-like procedure as
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
  Gz  <- apply(2 * pnorm(abs(zz2),lower.tail=FALSE), 1, min)
  
  ## Corrected p-values
  pcorr <- ecdf(Gz)(pval) ## is this efficient???
  return(pcorr)
}
