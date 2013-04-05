print.hdi <- function(x, ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  3 Apr 2013, 15:12

  method <- x$method

  if(method == "multi-split"){
    cat("alpha = 0.01:")
    cat(" Selected predictors:", which(x$pval <= 0.01), "\n")
    cat("alpha = 0.05:")
    cat(" Selected predictors:", which(x$pval <= 0.05), "\n")
    cat("------\n")
    cat("Familywise error rate controlled at level alpha.\n")
  }
  else if(method == "stability"){
    for(i in 1:length(x$select)){
      cat("EV = ", x$EV[[i]], ":", sep = "")
      cat(" Selected predictors:", x$select[[i]], "\n")
    }
    cat("------\n")
    cat("Expected number of false positives controlled at level EV.\n")
  }
}

