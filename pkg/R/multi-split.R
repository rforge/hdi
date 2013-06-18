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
  n.right <- n - n.left
      
  for(b in 1:B){
    try.again <- TRUE
    while(try.again){
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

      p.sel <- length(sel.model)

      ## Classical situation
      if(p.sel > 0 & p.sel < nrow(x.right) - 1){ ## intercept!
        sel.pval <- do.call(classical.fit,
                            args = c(list(x = x.right[,sel.model],
                                y = y.right), args.classical.fit))
        pvals[b, sel.model] <- pmin(sel.pval * p.sel, 1)
        try.again <- FALSE
      }
      ## Empty model selected
      if(p.sel == 0){
        ## Do nothing in that case. Matrix already filled with 1's
        if(trace)
          cat("......Empty model selected. That's ok...\n")

        try.again <- FALSE
      }
      ## Too large model selected for classical fitter
      if(p.sel >= n.right - 1){ ## p.sel + 1 < n.right for p-val calculation
        try.again <- TRUE ## re-do sample splitting
        warning("Too large model selected in a sample split")
      }
    } ## end while(try.again)
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





