ridge.proj <- function(x, y, ci.level = 0.95,
                       multiplecorr.method = "holm",
                       N = 10000, 
                       group.testing = FALSE,
                       groups = NULL,
                       activeset = NULL)
{
  ## Purpose:
  ## calculation of the ridge projection method proposed in
  ## http://arxiv.org/abs/1202.1377 P.Buehlmann
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## x: the design matrix
  ## y: the response vector
  ## ci.level: the significance level to be used for the confidence
  ##           interval calculation
  ## ----------------------------------------------------------------------
  ## Return values:
  ## individual: the individual testing p-values for each parameter
  ## corrected:  the multiple testing corrected p-values for each parameter
  ## betahat:    the initial estimate by the scaled lasso of \beta^0
  ## bhat:       the de-sparsified \beta^0 estimate used for p-value calculation
  ## sigmahat:   the \sigma estimate coming from the scaled lasso
  ## ci:         the confidence intervals calculated for each parameter
  ## ----------------------------------------------------------------------
  ## Author: Peter Buehlmann (initial version),
  ##         adaptations by L. Meier and R. Dezeure

  ## these are some old arguments that are still used in the code below
  ridge.unprojected <- TRUE
  
  x <- scale(x, scale = FALSE) ## *center* the columns
  y <- scale(y, scale = FALSE)

  n  <- nrow(x)
  p  <- ncol(x)
  y  <- as.numeric(y)

  biascorr <- Delta <- numeric(p)
  
  lambda <- 1  ## other choices?

  h1 <- svd(x) ## x must be *centered*, see above
  
  Px               <- tcrossprod(h1$v)
  Px.offdiag       <- Px
  diag(Px.offdiag) <- 0 ## set all diagonal elements to zero

  ## Use svd for getting the inverse for the ridge problem
  hh <- h1$v %*% ((h1$d / (h1$d^2 + lambda)) * t(h1$u))
  
  ## Ruben Note: here the h1$d^2/n 1/n factor has moved out, the value of lambda
  ## used is 1/n! See also comparing to my version
  
  cov2      <- tcrossprod(hh) 
  diag.cov2 <- diag(cov2)

  ## Estimate regression parameters (initial estimator) and noise level sigma
  fit.scaleL <- scalreg(X = x, y = y, lam0 = "univ")
  beta.lasso <- fit.scaleL$coefficients
  hat.sigma2 <- fit.scaleL$hsigma^2 ## added ^2 to fix bug!

  ## bias correction
  biascorr <- crossprod(Px.offdiag, beta.lasso)

  ## ridge estimator
  hat.beta <- hh %*% y

  hat.betacorr <- hat.beta - biascorr

  if(ridge.unprojected){ ## bring it back to the original scale
    hat.betacorr <- hat.betacorr / diag(Px)
  }
  
  ## Ruben Note: a_n = 1 / scale.vec, there is no factor sqrt(n) because this
  ## falls away with the way diag.cov2 is calculated see paper
  scale.vec  <- sqrt(hat.sigma2 * diag.cov2)
  
  ##1, is coupled with 2!! don't shut this off and not the other or vice versa
  if(ridge.unprojected){
    scale.vec <- scale.vec / abs(diag(Px))
  }
  
  hat.betast <- hat.betacorr / scale.vec ## sqrt(hat.sigma2 * diag(cov2))
  Delta      <- apply(abs(Px.offdiag), 2, max) * (log(p) / n)^(0.45) / scale.vec
  ## new
  if(ridge.unprojected ){ ## 2
    Delta <- Delta / abs(diag(Px))
    ## to put it on the same scale when scale.vec is different!
  }
     
  hat.gamma1 <- abs(hat.betast)

  #########################
  ## Individual p-values ##
  #########################
  
  res.pval <- c(pmin(2 * pnorm(hat.gamma1 - Delta, lower.tail = FALSE), 1))

  #########################################
  ## Multiple testing corrected p-values ##
  #########################################
  
  if(multiplecorr.method == "WY"){
    ## Westfall-young like procedure 
    pcorr <- p.adjust.wy(cov = cov2, pval = res.pval, N = N)
  }else{
    if(multiplecorr.method %in% p.adjust.methods){
      pcorr <- p.adjust(res.pval, method = multiplecorr.method)
    }else{
      stop("Unknown multiple correction method specified")
    }
  }

  ##########################
  ## Confidence Intervals ##
  ##########################

  scaleb <- 1 / scale.vec
  se     <- 1 / scaleb
  
  myci <- calc.ci(bj = hat.betacorr, se = se, level = ci.level)
  
  myci$lci <- myci$lci - Delta * se
  myci$rci <- myci$rci + Delta * se

  ## need to multiply Delta with se because it is on the scale of
  ## standard normal dist and we want to bring it to the distribution of
  ## \hat{\beta}_j 

  list(individual    = res.pval,
       corrected     = pcorr,
       ci            = cbind(myci$lci, myci$rci),
       bhatuncorr    = hat.beta,
       biascorr      = biascorr,
       normalisation = 1/scale.vec,
       bhat          = hat.betacorr,
       betahat       = beta.lasso,
       sigmahat      = sqrt(hat.sigma2))
}
