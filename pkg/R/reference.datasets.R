generate.reference.dataset <- function(n = 100,
                                       p = 500,
                                       s0 = 3,
                                       xtype = c("toeplitz","exp.decay","equi.corr"),
                                       btype = c("U[-2,2]","U[0,2]","U[0,4]","bfix1","bfix2","bfix10"),
                                       permuted = FALSE,
                                       iteration = 1)
{
  ## Purpose: generate the reference datasets used in the paper accompanying 
  ##          this R package
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 6 Oct 2015, 16:27

  ##Checking arguments
  xtype <- match.arg(xtype)
  btype <- match.arg(btype)
  stopifnot(iteration > 0  && iteration <= 50)##only 50 versions of each were created
  stopifnot(s0 >= 0)
  
  seed <- iteration + 2
  set.seed(seed)
  x <- generate.reference.x(n = n,
                            p = p,
                            xtype = xtype,
                            permuted = permuted)
  
  set.seed(seed)
  beta <- generate.reference.beta(p = p,
                                  s0 = s0,
                                  btype = btype)
  out <- list(x=x,beta=beta)
  return(out)
}

generate.reference.x <- function(n,
                                 p,
                                 xtype,
                                 permuted)
{
  ## Purpose: generate the reference design matrix used in the paper accompanying 
  ##          this R package
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 6 Oct 2015, 16:27
  
  Sigma <- switch(tolower(xtype),##in case capitalization was used
                  "toeplitz"={
                    indices <- toeplitz(0:(p-1))
                    cov <- (0.9)^abs(indices)
                    solve(solve(cov))##This looks really stupid but the tiny numerical differences make it identical to the datasets used in the paper
                  },
                  "equi.corr"={
                    cov <- matrix(rep(0.8,p*p),p)
                    diag(cov) <- 1
                    solve(solve(cov))##This looks really stupid but the tiny numerical differences make it identical to the datasets used in the paper
                  },
                  "exp.decay"={
                    indices <- toeplitz(0:(p-1))
                    icov <- 0.4^(abs(indices)/5)
                    solve(icov)
                  },
                  {
                    stop("Invalid xtype: Please provide any of the following xtypes 'toeplitz', 'equi.corr' or 'exp.decay'")
                  })

  x <- mvrnorm(n,rep(0,p),Sigma)

  if(permuted)
    {
      colindexing <- 1:ncol(x)
      x <- x[,sample(colindexing)]
    }

  return(x)
}

generate.reference.beta <- function(p,
                                    s0,
                                    btype)
{
  ## Purpose: generate the reference coefficient vector used in the paper accompanying 
  ##          this R package
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 6 Oct 2015, 16:27

  invalid.btype.txt <- "Invalid btype: Please provide a btype of the form 'bfix*' for a fixed value or 'U[*,*]', where * is a number."
  
  if(grepl("U\\[",btype) &&
     grepl("\\]",btype) &&
     (sapply(regmatches(btype, gregexpr(",", btype)), length)==1))##there should be one and only one comma in btype
    {
      ##generate random coefficients from a uniform distribution
      split.btype <- strsplit(btype,",")[[1]]

      lower <- suppressWarnings(as.numeric(substr(split.btype[1],
                                                  nchar("U[")+1,
                                                  nchar(split.btype[1]))))##extract the lower bound
      upper <- suppressWarnings(as.numeric(substr(split.btype[2],
                                                  0,
                                                  nchar(split.btype[2])-1)))##extract the upper bound
      if(is.na(lower) || is.na(upper))##if there were any errors in the extraction
        stop(invalid.btype.txt)
      
      beta <- c(runif(s0,lower,upper),rep(0,p-s0))      
    }else{
      if(grepl("bfix",btype)){
        
        b <- suppressWarnings(as.numeric(substr(btype,nchar("bfix")+1,nchar(btype))))##extract the bfix number

        if(is.na(b))##if there were any errors in the extraction
          stop(invalid.btype.txt)
        
        beta <- c(rep(b,s0),rep(0,p-s0))
        
      }else{
        stop(invalid.btype.txt)
      }
    }  
  return(beta)
}
