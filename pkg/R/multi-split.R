multi.split <- function(x, y, B, fraction,
                        model.selector,
                        classical.fit,
                        gamma,
                        args.model.selector,
                        args.classical.fit,
                        trace = TRUE)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  2 Apr 2013, 11:52

  n <- nrow(x)
  p <- ncol(x)
  
  ## matrix of bootstrap p-values
  pvals <- matrix(1, nrow = B, ncol = ncol(x))

  n.left <- floor(n * fraction)

  for(b in 1:B){
    if(trace)
      cat("...Split", b, "\n")

    ## Perform sample-splitting; sample *without* replacement
    split <- sample(1:n, n.left, replace = FALSE)
    
    ## x.left is used to perform model selection
    x.left <- x[split,]
    y.left <- y[split]
    
    ## x.right is used to calculate p-values
    x.right <- x[-split,]
    y.right <- y[-split]

    sel.model <- do.call(model.selector,
                         args = c(list(x=x.left, y=y.left),
                           args.model.selector))
    
    if(length(sel.model) > 0){
      sel.pval <- do.call(classical.fit,
                          args = c(list(x = x.right[,sel.model], y = y.right),
                            args.classical.fit))
      pvals[b, sel.model] <- pmin(sel.pval * length(sel.model), 1)
    }
    if(length(sel.model) == 0)
      ## Do nothing in that case. Matrix already filled with 1's
      if(trace)
        cat("......Empty model selected. That's ok...\n")
  }
  ## For loop is not very innovative, but it does it's job...
  pvals.current <- which.gamma <- numeric(p)
  for(j in 1:p){
    quant.gamma <- quantile(pvals[,j], gamma) / gamma
    
    pvals.pre        <- min(quant.gamma) * (1 + log(1/min(gamma)))
    pvals.current[j] <- pmin(pvals.pre, 1)
    
    which.gamma[j] <- which.min(quant.gamma)
  }

  names(pvals.current) <- names(which.gamma) <- colnames(x)
  
  return(list(pval = pvals.current,
              gamma.min = gamma[which.gamma]))
}





