boot.lasso.proj <- function(x, y, family = "gaussian",
                            standardize = TRUE,
                            multiplecorr.method = "WY",
                            parallel = FALSE, ncores = getOption("mc.cores", 2L),
                            betainit = "cv lasso",
                            sigma = NULL, ## sigma estimate provided by the user
                            Z = NULL,     ## Z or Thetahat provided by the user
                            verbose = FALSE,
                            return.Z = FALSE,
                            suppress.grouptesting = FALSE,
                            robust = FALSE,
                            B = 1000,
                            boot.shortcut = FALSE,
                            return.bootdist = FALSE)
{
  ## Purpose:
  ## An implementation of the LDPE method http://arxiv.org/abs/1110.2563
  ## (which is identical to http://arxiv.org/abs/1303.0518)
  ## performing inference using the bootstrap techniques from
  ## http://arxiv.org/abs/1606.03940
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Return values:
  ## pval: p-values for every parameter (individual tests)
  ## pval.corr:  multiple testing corrected p-values for every parameter
  ## betahat:    initial estimate by the scaled lasso of \beta^0
  ## bhat:       de-sparsified \beta^0 estimate used for p-value calculation
  ## sigmahat:   \sigma estimate coming from the lasso
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 15 June 2016 (initial version)
  ## in part based on the implementation of the de-sparsified Lasso
  ## lasso-proj.R

  if(!identical(family,"gaussian"))
    stop("The function boot.lasso.proj is currently not supporting families other than 'gaussian'")
  ## The current code assumes the family is gaussian
  ## residual based resampling would be incorrect for other family
  
  ####################
  ## Get data ready ##
  ####################
  n <- nrow(x)
  p <- ncol(x)

  if(standardize)
    sds <- apply(x, 2, sd)
  else
    sds <- rep(1, p)

  pdata <- prepare.data(x = x,
                        y = y,
                        standardize = standardize,
                        family = family)
  x <- pdata$x
  y <- pdata$y
  
  ## force sigmahat to 1 when doing glm!
  if(family == "binomial")
    sigma <- 1
  
  ######################################
  ## Calculate Z using nodewise lasso ##
  ######################################
  
  Zout <- calculate.Z(x = x,
                      parallel = parallel,
                      ncores = ncores,
                      verbose = verbose,
                      Z = Z)
  Z <- Zout$Z
  scaleZ <- Zout$scaleZ
  
  ###################################
  ## de-sparsified Lasso estimator ##
  ###################################
  
  ## Initial Lasso estimate
  initial.estimate <- initial.estimator(betainit = betainit, sigma = sigma,
                                        x = x, y = y)
  betalasso <- initial.estimate$beta.lasso
  sigmahat <- initial.estimate$sigmahat
  
  ## de-sparsified Lasso
  bproj <- despars.lasso.est(x = x,
                             y = y,
                             Z = Z,
                             betalasso = betalasso)
  se <- est.stderr.despars.lasso(x = x,
                                 y = y,
                                 Z = Z,
                                 betalasso = betalasso,
                                 sigmahat = sigmahat,
                                 robust = robust)
  scaleb <- 1/se
  
  bprojrescaled <- bproj * scaleb
  
  #########################
  ## p-Value calculation ##
  #########################

  ## Compute the correct centered bootstrap distribution

  lambda <- NULL
  if(boot.shortcut)
  {
    stop("TODO implement the bootstrap shortcut technique!")
  }

  ##stop("TODO implement bootstrap")
  cboot.dist <- replicate(B,rnorm(ncol(x),sd=se))## Centered bootstrap distribution

  
  ## Compute p-values based on that distribution
  dist <- bproj-cboot.dist
  counts <- apply(dist >=0,1,sum)
  if(any(counts >= B/2))
    counts[counts >= B/2] <- B-counts[counts >= B/2]
  counts <- 2*counts
  
  pval <- (counts+1)/(B+1)

  ## Compute confidence intervals
  ## How to deal with this? just have in 

  #################################
  ## Multiple testing correction ##
  #################################
  
  cboot.dist.underH0c <- NULL
  pcorr <- if(multiplecorr.method == "WY") {
             ###########################################################################
             ## Bootstrap under H0 complete another B bootstrap samples to perform WY ##
             ###########################################################################
             
             ##stop("TODO implement bootstrap under H0 complete")
             ##Compute the maxT distribution under H0c
             cboot.dist.underH0c <- replicate(B,rnorm(ncol(x),sd=se))
             max.t.dist <- apply(abs(cboot.dist.underH0c),2,max)

             counts.matrix <- sapply(max.t.dist,FUN=">=",abs(bproj))
             counts <- apply(counts.matrix,1,sum)
             
             (counts+1)/(B+1)  
           } else if(multiplecorr.method %in% p.adjust.methods) {
             ####################################################################
             ## Check and warn for accuracy of p-values when using traditional ##
             ## multiple testing correction methods                            ##
             ####################################################################

             ## P-values that are so far in the tail that the number of bootstrap
             ## samples didn't suffice to capture
             inaccurate.pval <- pval <= (4+1)/(B+1)
             
             ## Is the lower bound on p-values too small for multiple testing correction?
             B.toosmall.formc <- (1+4)/(B+1)*ncol(x) > 0.05
             ## we don't know what rejection level the user wants to use, but 0.05 is a good reference?
             
             if(any(inaccurate.pval) & B.toosmall.formc)
             {
               warning(paste("Lacking accuracy for multiple testing:\n",
                             "The provided number of bootstrap samples B was not high enough to ",
                             "accurately compute\nthe smallest multiple testing corrected p-values.\n",
                             "The lowest attainable p-value for the chosen B value is 1/(B+1)=",
                             signif(1/(B+1),3),".\n",
                             "Reasonable accuracy is only attained starting from (4+1)/(B+1)=",
                             signif((4+1)/(B+1),3),".\n",
                             "This issue is easily solved by using 'WY' as multiple testing method.",
                             sep=""))
             }
             
             p.adjust(pval,method = multiplecorr.method)
           } else
             stop("Unknown multiple correction method specified")


  ##############################################
  ## Function to calculate p-value for groups ##
  ##############################################

  ##TODO implement group testing for boot lasso proj
  group.testing.function <- NULL
  cluster.group.testing.function <- NULL

  ############################
  ## Return all information ##
  ############################
  
  out <- list(pval        = as.vector(pval),
              pval.corr   = pcorr,
              groupTest   = group.testing.function,
              clusterGroupTest = cluster.group.testing.function,
              sigmahat    = sigmahat,
              standardize = standardize,
              sds         = sds,
              bhat        = bproj / sds,
              se          = se / sds,
              betahat     = betalasso / sds,
              family      = family,
              method      = "boot.lasso.proj",
              B           = B,
              call        = match.call())
  if(return.Z)
    out <- c(out,
             list(Z = scale(Z,center=FALSE,scale=1/scaleZ)))##unrescale the Zs
  names(out$pval) <- names(out$pval.corr) <- names(out$bhat) <-
    names(out$sds) <- names(out$se) <- names(out$betahat) <-
    colnames(x)
  if(return.bootdist)
  {
    out <- c(out,
             list(cboot.dist = cboot.dist))
    rownames(out$cboot.dist) <- names(out$bhat)
    if(!is.null(cboot.dist.underH0c))
    {
      out <- c(out,
               list(cboot.dist.underH0c = cboot.dist.underH0c))      
      rownames(out$cboot.dist.underH0c) <- names(out$bhat)
    }
  }
  class(out) <- "hdi"
  return(out)
}

