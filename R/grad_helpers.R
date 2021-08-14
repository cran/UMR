#' @rdname grad_helpers
#' @name AA
#' @title Helper functions for calculating gradient  of least-squares Shuffled Isotonic Regression criterion, for Laplace or for Gaussian errors
#'
#'
#' 
#' @export AA
#' @export BB
#' @export AAfunc_Laplace_generic
#' @export AAfunc_Gauss_generic
#' @export BBfunc_Laplace_generic
#' @export BBfunc_Gauss_generic
#' @export getAAfunc_est_outer
#' @export getBBfunc_est_outer
#' @export BBfunc_mixGauss_generic
#' @export BBpfunc_mixGauss_generic
#' @export BBpfunc_Gauss_generic
#' @export BBpfunc_Laplace_generic
#' 
#'
#' @importFrom stats pnorm
#' @importFrom distr Norm
#' 
#' @param yy Y (response) observation vector (numeric).  Will apply
#'     as.vector() so it may be a matrix or array with all dimensions trivial
#'     except 1.
#' @param mm Current (unsorted) estimate/iterate at which to compute
#'     gradient. (Length equals length of yy).  Will apply as.vector() so it
#'     may be a matrix or array with all dimensions trivial except 1.
#' @param func This is a function; should be the actual "A" or "B" function
#'     from the paper; AA and BB are just wrappers that call outer() with
#'     func().  func() should accept vector or matrix arguments.
#'
#' @details See helper functions "A" and "B" in paper.
#'
#' @examples
#'
#' 
#' 
#'  ## the "!!" de-quote (see ?partial) so e.g., can save mygradSIR for future runs.
#'
#' ####### gradient settings/setup for Gaussian
#'
#' 
#' ## set.seed(501)
#' library(distr)
#' mysig <- 1 ##  std dev
#' errdist <- Norm(0, sd=mysig)
#' mm0 <- function(xx){xx}
#' nn <- 300
#' xx <- sort(runif(n=nn, 0, 7))
#' yy <- mm0(xx) + errdist@r(nn)
#' ## plot(xx,yy)
#' 
#' myScale <- mysig
#' 
#' AAfunc_Gauss <- purrr::partial(AAfunc_Gauss_generic, sig=!!mysig)
#' AA_Gauss <- purrr::partial(AA, func=!!AAfunc_Gauss)
#' BBfunc_Gauss <- purrr::partial(BBfunc_Gauss_generic, sig=!!mysig)
#' BB_Gauss <- purrr::partial(BB, func=!!BBfunc_Gauss)
#' mygradSIR <- 
#'     grad_SIR_Gauss <- ## just for ease of reference
#'         purrr::partial(grad_SIR_generic,
#'                        rescale=TRUE, ## factor of nn/2
#'                        AAfunc=!!AA_Gauss, BBfunc=!!BB_Gauss)
#'
#' ####### gradient settings/setup for Laplace
#'
#'
#'
#' 
#'set.seed(501)
#'library(distr)
#'myLL <- .7 ## (1/"rate") parameter, aka "mean" parameter (except Laplace mean is 0)
#' errdist <- DExp(1/myLL)
#'
#'
#'nn <- 200
#'mm0 <- function(xx){
#'    (xx<=0)*0 + (0<=xx & xx<=2)*1 +
#'        (2<xx & xx<=3)*3 +
#'        (3<xx)*6
#'}
#'xx <- sort(runif(n=nn, 0, 7))
#'yy <- mm0(xx) + errdist@r(nn)
#'
#'myScale <- myLL;
#'
#'## CS settings
#'#'mysig <- sqrt(2) * myLL;
#'#' 
#' AAfunc_Laplace <- purrr::partial(AAfunc_Laplace_generic, LL=!!myLL)
#' AA_Laplace <- purrr::partial(AA, func=!!AAfunc_Laplace)
#' BBfunc_Laplace <- purrr::partial(BBfunc_Laplace_generic, LL=!!myLL)
#' BB_Laplace <- purrr::partial(BB, func=!!BBfunc_Laplace)
#' mygradSIR <-
#'     grad_SIR_Laplace <- purrr::partial(grad_SIR_generic,
#'                                        rescale=TRUE, ## factor of nn/2
#'                                       AAfunc=!!AA_Laplace, BBfunc=!!BB_Laplace)

## no real point in the A generic functions.  they are just the survival
## function.




## The nomenclature ended up being slightly messy, and inconsistent.
## Probably grad_SIR_generic should take as argument not AAfunc but
## 'AAfunc_outer' or something.  Similarly for BB.  And then the variabels
## passed in could be renamed.

## yy and mm are vectors (usually of same length, nn), say of lengths n1 and
## n2.  returns an n1 times n2 matrix.



AA <- function(yy, mm, func){
    dd  <- outer(as.vector(yy), as.vector(mm), "-")
    func(dd)
}


#' @rdname grad_helpers
#' @name BB
#'
## #' @param mmhat Vector of estimates of mm (increasing, so corresponding to $X_{(1)}, \ldots, X_{(n)}$.)  (Actually: mmhat is not truly needed; one really should use 'mm' twice; but the code has not yet been changed.)

##BB <- function(mmhat, mm, func){
BB <- function(mm, func){
    ## dd <- outer(as.vector(-mmhat), as.vector(mm), "+")
    dd <- outer(as.vector(-mm), as.vector(mm), "+") ## what this shoul dbe ;
    ##mmhat dropped i thnk
    func(dd)
}


#' @rdname grad_helpers
#' @name AAfunc_Laplace_generic
#'
#' @param  LL  Double Exponential "mean" parameter: corresponding density is $exp(-|d|/LL) / (2LL)$.
#' @param dd generic argument to the "A" function; usually of the form m - mmhat, where m is just some value of the regression function
AAfunc_Laplace_generic  <-  function(dd, LL){
    vv <- exp(-abs(dd)/LL)/2
    res <- vv * (dd >= 0) +
        (1-vv) * (dd < 0);
}

## sig is *standard deviation* (not variance)

#' @rdname grad_helpers
#' @name AAfunc_Gauss_generic

AAfunc_Gauss_generic <- function(dd, sig){
    1 - pnorm(dd / sig)
}


#' @rdname grad_helpers
#' @name BBfunc_Laplace_generic
## BBfunc_Laplace_generic_old <- function(dd, LL){
##     T1 <- exp(-abs(dd)/LL)/8
##     T2 <- -dd * exp(dd/LL) / (4*LL) ## used if dd <= 0
##     T3  <- 1/2 - exp(-dd/LL) * (2*LL+dd) / (4*LL)
##     (T1 + T3 + (1/2 - exp(-dd/LL)/8)) * (dd>=0) +
##         (T1 + T2 + exp(dd/LL)*3/8) * (dd<0);
## }
BBfunc_Laplace_generic <- function(dd, LL){
    (1/2 - (dd/ (4*LL))) * exp(dd/LL)  * (dd <0) +
        (1 - ((2*LL + dd)/(4*LL)) * exp(-dd/LL) ) * (dd>=0);
}




#### laplace generic function can trivially be simplified 

#' @rdname grad_helpers
#' @name BBfunc_Gauss_generic
#'
#' @param sig is standard deviation of the normal distribution.
BBfunc_Gauss_generic <- function(dd, sig){
    pnorm(dd/ (sig * sqrt(2)))
}



#### The following two functions return functions that are analogous to say
#### AAfunc_Gauss or BBfunc_Gauss.  Based on using residuals
#### to estimate the AA and BB functions

## not 'generic' because there is no unknown parameters involved.


#### This doesn't work because ecdf doesn't take/return matrix arguments properly
## getAAfunc_est <- function(eps){
##     function(x){1-(ecdf(eps))(x)}
## }


#' @rdname grad_helpers
#' @name getAAfunc_est_outer
#'
#' @param eps is a vector of residuals (or estimated residuals).  In current
#'     coding, it should have been preprocessed to be *unique*.  (If there
#'     are repeats this should be encoded in ww).
#' @param ww is vector of weights of same length as eps, and summing to 1.
#'     Default is a weight of 1/length(eps) for each value of eps; if eps has
#'     been pre-binned then ww is the weights from binning.
#'
#' @details For getAAfunc_est:  returns a function(yy,mm) which is analogous to passing in an estimated 'func' argument to the AA function.  (Reason to not do it that way relates to making sure matrix arguments are handled correctly.)


## for the estimated AA function, I combine the 'AA' and
## 'AAfunc_DISTN_generic' into one, and similarly for BB.
getAAfunc_est_outer <- function(eps, ww=1/length(eps)){
    Phihat <- getEcdf(eps, ww)
    function(yy, mm){
        outer(as.vector(yy), as.vector(mm),
              function(y,m){1-Phihat(y-m)}
              )
    }
}

## getAAfunc_est_outer <- function(eps){
##     Phihat <- ecdf(eps)
##     function(yy, mm){
##         outer(as.vector(yy), as.vector(mm),
##               function(y,m){1-Phihat(y-m)}
##               )
##     }
## }



#' @rdname grad_helpers
#' @name getBBfunc_est
#'
#'
#'
#'
#'
#'
#' @details getBBfunc_est returns a function which as of this coding *MUST*
#'     take only a numeric vector of length 1; longer vectors will not
#'     work. Be careful!  Note that ecdf objects are not intended to be
#'     stored permanently so storing functions returned by
#'     getBBfunc_est_outer or getAAfunc_est_outer may cause issues.
#'

## ## does '@details' work when added on in this way ?


## I don't know the precise details of how "outer" works; using it with the
## output of ecdf() did not work as I expected in my original attempts at
## getBBfunc_est_outer.

## getBBfunc_est_outer <- function(eps){
##     Phihat <- ecdf(eps)
##     OO <- length(eps)
##     function(mm){
##         dd <- outer(as.vector(-mm), as.vector(mm), "+")
##         outprod_args <- outer(eps, dd, "+") ## dims: length(eps) by length(mm) by length(mm)
##         ## #######  ## Seems like combining these two 'outer' calls fails b/c of ecdf().
##         outprod_Phi <- apply(outprod_args, c(1,2,3), Phihat) ## This part could be sped up
##         ## #######  ## would be speedier to sort first and do something a bit more complicated
##         apply(outprod_Phi, c(2,3), function(vec){sum(vec)/OO})
##     }
## }


getBBfunc_est_outer <- function(eps, ww=1/length(eps)){
    Phihat <- getEcdf(eps, ww) ##
    ## OO <- length(eps)
    function(mm){
        dd <- outer(as.vector(-mm), as.vector(mm), "+")
        outprod_args <- outer(eps, dd, "+") ## dims: length(eps) by length(mm) by length(mm)
        ## #######  ## Seems like combining these two 'outer' calls fails b/c of ecdf().
        outprod_Phi <- apply(outprod_args, c(1,2,3), Phihat) ## This part could be sped up
        ## #######  ## would be speedier to sort first and do something a bit more complicated
        apply(outprod_Phi, c(2,3), function(vec){sum(vec * ww)})
    }
}




## the above BB estimated implementation can probably be done more cleverly; this is 'brute force'





#' @rdname grad_helpers
#' @name BBfunc_mixGauss_generic
#'
#' @param locs Vector (length LL) of mixture locations
#' @param wws Vector (length LL, sum to 1) of mixture weights
#' @param sigs Vector (length LL, positive) of component standard deviations
#' ## here dd should be a matrix (usually from a call to outer())
BBfunc_mixGauss_generic <- function(dd,
                                    locs,
                                    wws,
                                    sigs){
    if (!is.matrix(dd)) stop("dd should be a matrix.")
    doOne <- function(dd1){
        if (length(dd1) !=1) stop("length dd != 1")
        arg1  <- outer(as.vector(locs), as.vector(locs), "-") + dd1
        arg2 <- outer(as.vector(sigs), as.vector(sigs),
                      function(sig1, sig2){ sqrt(sig1^2 + sig2^2)})
        scale <- outer(as.vector(wws), as.vector(wws), "*")
        sum(scale * pnorm(arg1/arg2))
    }
    apply(dd, c(1,2), doOne)
}

#' @rdname grad_helpers
#' @name BBpfunc_mixGauss_generic
## vectorize not good enough here (so used apply()).
## here dd should be a matrix (usually from a call to outer())
BBpfunc_mixGauss_generic <- function(dd,
                                     locs,
                                     wws,
                                     sigs){
    if (!is.matrix(dd)) stop("dd should be a matrix.")
    doOne <- function(dd1){
        if (length(dd1) !=1) stop("length dd != 1")
        arg1  <- outer(as.vector(locs), as.vector(locs), "-") + dd1
        arg2 <- outer(as.vector(sigs), as.vector(sigs),
                      function(sig1, sig2){ sqrt(sig1^2 + sig2^2)})
        scale <- outer(as.vector(wws), as.vector(wws), "*")
        sum(scale * stats::dnorm(arg1/arg2) / arg2)
    }
    apply(dd, c(1,2), doOne)
}

#' @rdname grad_helpers
#' @name BBpfunc_Gauss_generic
#'
#' @param xx Point at which to evaluate function
BBpfunc_Gauss_generic <- function(xx, sig){
    distr::Norm(0, sd=sig*sqrt(2))@d(xx) 
}

#' @rdname grad_helpers
#' @name BBpfunc_Laplace_generic
#' @param myLL is Laplace parameter.
## derivation in paper 
BBpfunc_Laplace_generic <- function(xx, myLL){
    exp(-abs(xx)/myLL) * (myLL + abs(xx)) / (4 * myLL^2)
}



## ## AA func in general should be very generic and just have errdist passed in
## AAfunc_mixGauss_generic <- function(dd,
##                                     locs,
##                                     wws,
##                                     sigs){
##     Univar
##     1 - 
## }
