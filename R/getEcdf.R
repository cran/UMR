## like the 'ecdf' function but works with weights
## the xs must be unique! if they are not unique and you have equal weights
## use the standard 'ecdf' function.

getEcdf <- 
function (x, w=1/length(x)) 
{
  x <- sort.int(x, index.return=TRUE)
  w <- w[x$ix]
  x <- x$x
  n <- length(x)
  if (n < 1) 
    stop("'x' must have 1 or more non-missing values")
  ##vals <- unique(x)
  rval <- approxfun(x, cumsum(w), 
                    method = "constant", yleft = 0, yright = 1, f = 0,
                    ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}
