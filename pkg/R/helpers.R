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

lm.pval <- function(x, y, exact = TRUE)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  2 Apr 2013, 11:34

  fit.lm <- lm(y ~ x) ## Intercept??? Exceptions???
  fit.summary <- summary(fit.lm)

  tstat <- fit.summary$coefficients[-1, "t value"] ## Intercept??? Exceptions???

  if(exact){ ## Use appropriate t-dist
    pval.sel <- 2 * pt(abs(tstat), df = fit.lm$df.residual,
                       lower.tail = FALSE)
  }else{ ## p-values based on *normal* distribution
    pval.sel <- 2 * pnorm(abs(tstat), lower.tail = FALSE)
  }
  
  names(pval.sel) <- colnames(x)
  pval.sel
}

glm.pval <- function(x,y,family="binomial",trace=FALSE,...)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Ruben Dezeure based on lm.pval, Date:  30 Sept 2013, 18:04
    
    fit.glm <- glm(y ~ x,family=family,...) ## Intercept??? Exceptions???
    fit.summary <- summary(fit.glm)
    
    if(!fit.glm$converged && trace)
        {
            print(fit.summary)
        }
    
    if(family %in% c("poisson","binomial"))
        {
            zstat <- fit.summary$coefficients[-1, "z value"] ## Intercept??? Exceptions???
            
            ## p-values based on *normal* distribution
            pval.sel <- 2 * pnorm(abs(zstat), lower.tail = FALSE)
        }else{
            tstat <- fit.summary$coefficients[-1, "t value"] ## Intercept??? Exceptions???
            
            pval.sel <- 2 * pt(abs(tstat), df = fit.lm$df.residual,
                               lower.tail = FALSE)
        }
    
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


