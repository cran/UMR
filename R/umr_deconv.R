## unmatched mono reg via deconvolution
## (As in Carpentier and Schl\"uter 2016)


#' @rdname umr_deconv
#' @title  Carpentier and  Schluter 2016  deconvolution method  for unmatched
#'     monotone regression
#'
#'
#' @export umr_deconv
#' @export quant_deconv
#' @importFrom stats approxfun
#' @importFrom decon DeconCdf
#'
#' @param xx X (covariate or predictor) observation vector
#' @param yy Y (response) observation vector (numeric)
#' @param sig standard deviation of epsilon (passed to DeconCdf)
#' @param error Must be "normal" or "laplacian" or "snormal"; see help("DeconCdf")
#' @param bw Bandwidth choice or method for kernel estimator; see help("DeconCdf")
#' @param adjust See help("DeconCdf")
#' @param n See help("DeconCdf")
#'
#' @details quant_deconv implements Carpentier and Schluter 2016
#'     deconvolution method for unmatched monotone regression, using deconv
#'     package. Note that because the DeconCdf() function computes the CDF
#'     but there is no direct code for computing the quantile function, we
#'     use approxfun to create the quantile function; this may be slow.
#'     quant_deconv() returns a vector of length length(yy).  Then umr_deconv
#'     is a wrapper for quant_deconv.  NOTE: It returns the output of
#'     approxfun, which is may change over time.  The output value is of type
#'     function.  We linearly interpolate between the points i/n.
#'
#'
#' @examples
#'
#'
#' library(distr)
#' mysig <- 1 ##  std dev
#' errdist <- distr::Norm(0, sd=mysig)
#' mm0 <- function(xx){xx}
#'  nn <- 300
#'  xx <- sort(runif(n=nn, 0, 7))
#'  yy <- mm0(xx) + errdist@r(nn)
#'  ## plot(xx,yy)
#'  modeldistname <- truedistname <- "Gauss" ## used for savefile name
#' myScale <- mysig
#'
#' xx <- sort(runif(n=nn, 0, 7))
#' mmtrue <- mm0(xx)
#' yy <- mmtrue + errdist@r(nn)
#' plot(xx,yy)
#' qq <- quant_deconv(yy, sig=1, error="normal")
#' lines(xx, ## already sorted
#'      qq)


umr_deconv <- function(xx, yy, sig, error="normal",
                       bw="dboot1", adjust=1,
                       n=512){
    qq <- quant_deconv(yy, sig, error, bw=bw, adjust=adjust, n=n)
    approxfun(xx, qq)
}



#' @rdname umr_deconv
#' @name quant_deconv
#'
#' @param monotonize is a function taking a numeric vector argument which
#'     returns an increasing numeric vector of the same length.  This is used
#'     to monotonize the output of the CDF from deconvolution, which is not
#'     guaranteed to be a "bona-fide" CDF in the sense that it may not be
#'     monotone.

quant_deconv <- function(yy, sig, error="normal",
                         bw="dboot1", adjust=1,
                         n=512, monotonize=base::cummax){
    nn <- length(yy)
    ## need specify the accuracy still ...
    dCDF <- DeconCdf(y=yy, sig=sig, error=error,
                     bw=bw, adjust=adjust, n=n)
    yy <- monotonize(dCDF$y)    ## monotonization step
    if (length(yy) != length(dCDF$y))
        stop("The monotonize argument should return a vector of the same length as its input.")
    QF <- approxfun(yy, dCDF$x) ## quantile function
    ## This approxfun() call gives warnings() b/c of nonunique y values
    probs <- (1:nn)/nn
    QF(probs)
}
