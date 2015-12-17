generate.reference.dataset <-
    function(n, p, s0,
             xtype = c("toeplitz", "exp.decay", "equi.corr"),
             btype = "U[-2,2]", permuted = FALSE,
             iteration = 1, do2S = TRUE)
{
  ## Purpose: generate the reference datasets used in the paper accompanying
  ##          this R package
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 6 Oct 2015; 'do2S' and tweaks: Martin Maechler

  ## Checking arguments
  xtype <- match.arg(xtype)
  ## rather let generate.reference.beta() check

  stopifnot(is.character(btype), length(btype) == 1,
            n == as.integer(n), length(n) == 1, n >= 1,
            p == as.integer(p), length(p) == 1, p >= 1,
            s0== as.integer(s0),length(s0)== 1, 0 <= s0, s0 <= p)
  do.seed <- !is.na(iteration)
  if(do.seed)
      stopifnot(iteration > 0  && iteration <= 50) ## only 50 versions of each were created
  if(do.seed) {
      seed <- iteration + 2
      set.seed(seed)
  }
  x <- generate.reference.x(n = n, p = p, xtype = xtype,
                            permuted = permuted, do2S = do2S)
  if(do.seed) set.seed(seed)
  beta <- generate.reference.beta(p = p, s0 = s0, btype = btype)
  list(x=x, beta=beta)
}

generate.reference.x <- function(n, p, xtype, permuted, do2S = TRUE)
{
  ## Purpose: generate the reference design matrix used in the paper accompanying
  ##          this R package
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 6 Oct 2015, 16:27

  Sigma <- switch(tolower(xtype),##in case capitalization was used
                  "toeplitz" = {
                    indices <- toeplitz(0:(p-1))
                    cov <- (0.9)^abs(indices)
                    if(do2S) ## This looks really stupid but the tiny numerical differences
                             ## make it identical to the datasets used in the paper
                    solve(solve(cov))
                  },
                  "equi.corr" = {
                    cov <- matrix(rep(0.8,p*p),p)
                    diag(cov) <- 1
                    if(do2S) ## see remark above
                      solve(solve(cov))
                  },
                  "exp.decay" = {
                    indices <- toeplitz(0:(p-1))
                    icov <- 0.4^(abs(indices)/5)
                    solve(icov)
                  },
                    stop("Invalid xtype: Please provide any of the following xtypes 'toeplitz', 'equi.corr' or 'exp.decay'"))

  x <- mvrnorm(n, rep(0,p), Sigma)

  if(permuted)
    x[, sample.int(ncol(x))]
  else
    x
}

generate.reference.beta <- function(p, s0, btype)
{
  ## Purpose: generate the reference coefficient vector used in the paper accompanying
  ##          this R package
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 6 Oct 2015;  tweaks by Martin Maechler

  stopifnot(s0 <= p, p >= 0, length(s0) == 1, length(p) == 1,
            is.character(btype), length(btype) == 1)
  invalid.btype.txt <- "Invalid btype: Please provide a btype of the form 'bfix*' for a fixed value or 'U[*,*]', where * are two numbers, the lower and upper bounds."

  if(grepl("^U\\[",btype) &&
     grepl("\\]$", btype) &&
     sapply(regmatches(btype, gregexpr(",", btype)), length) == 1) {
    ## there should be one and only one comma in btype
    ## generate random coefficients from a uniform distribution
    split.btype <- strsplit(sub("^U\\[", '', sub("\\]$", '', btype)),
                            ",")[[1]]

    ## extract the lower and upper bounds:
    lower <- as.numeric(split.btype[1])
    upper <- as.numeric(split.btype[2])
    if(is.na(lower) || is.na(upper) || length(lower) != 1 || length(upper) != 1)
      ## if there were any errors in the extraction
      stop(invalid.btype.txt)

    b <- runif(s0,lower,upper)

  } else if(grepl("^bfix", btype)) {

    b <- as.numeric(sub("^bfix","", btype)) ## extract the bfix number
    if(is.na(b))##if there were any errors in the extraction
      stop(invalid.btype.txt)

    b <- rep(b, s0)
  }
  else
    stop(invalid.btype.txt)

  c(b, rep(0, p-s0))
}
