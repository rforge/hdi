stability <- function(x, y, B, fraction, model.selector, EV, q,
                      args.model.selector, trace = TRUE)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 28 Mar 2013, 12:49

  if(length(q) > 1) ## q has to be a scalar (EV can be a vector)
      stop("q must be a scalar")
  
  n <- nrow(x)
  p <- ncol(x)
  
  col.nam <- colnames(x)

  thresholds <- 0.5 * (1 + q^2 / (p * EV)) ## vector of thresholds

  if(any(thresholds > 1))
    warning("Some thresholds larger than 1. Decrease q or increase EV.")

  ## Matrix of selected models:
  ## rows = subsamples
  ## cols = predictors
  sel.mat <- matrix(FALSE, nrow = B, ncol = p) 

  sel.n <- floor(fraction * n)
  
  ## Subsampling
  for(b in 1:B){
    if(trace)
      cat("...Subsample", b, "\n")
    sel <- sample(1:n, sel.n, replace = FALSE)

    ## Current sub-sampled data
    x.sel <- x[sel,]
    y.sel <- y[sel]

    ## Get selected model
    sel.model <- do.call(model.selector, c(list(x = x.sel, y = y.sel, q = q),
                                           args.model.selector))
    sel.mat[b, sel.model] <- TRUE
  }
  
  ## Get selection frequencies
  freq <- colMeans(sel.mat)
  names(freq) <- col.nam
        
  out <- list(); length(out) <- length(EV)
  
  for(i in 1:length(EV)){
    sel.current        <- which(freq >= thresholds[i])
    names(sel.current) <- col.nam[sel.current]
    out[[i]] <- sel.current
  }
  
  return(list(select = out,
              EV = EV,
              thresholds = thresholds,
              freq = freq))
}
