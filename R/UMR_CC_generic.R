## Generic functions related to second derivative computations for unlinked monotone regression ("SIR" comes from "shuffled isotonic regression" although this terminology is now outdated).

## These functions serve as argument to the "CC" or "curv" argument in the
## active set algorithms.  (The "old" active set algorithm has named argument
## CC and the new one has curv, but they perform the identical function.  The
## naming seemed important though.)

#'  @title Second derivative computations of least-squares Unlinked  Isotonic Regression criterion ("SIR"  comes from "shuffled isotonic regression" although this terminology is now outdated).

#'
#' 
#' @export UMR_CC_generic
#' @export UMR_curv_generic
#' @export UMR_curv_generic2
#' 
#' 
#'
#' @param yy Y (response) observation vector (numeric)
#' @param mm Current (unsorted) estimate/iterate at which to compute
#'     gradient. (Length is <= than the number of X observations in the problem).
#' @param ww_y Weights (nonnegative, sum to 1) corresponding to yy.  Same length as yy.
#' Default is just 1/length(yy) for each value.
#' @param ww_m Weights (nonnegative, sum to 1) corresponding to mm.  Same length as mm.
#'
#' @param  densfunc This is the error density, a  function object (Balabdaoui, Doss, Durot (2021+).
#' @param DDfunc This is the function "D" defined in the second derivative 
#'     calculations in the paper (Balabdaoui, Doss, Durot (2021+).
#' @param BBpfunc This is the function B', i.e. derivative of "B" function in the paper. 
#'
#'     @details The  "CC" or  "curv" functions  are used to  be passed  in to
#'     UMRactiveSet_trust() (generally  after 'currying'/substituting  in for
#'     the  parameter  arguments).   UMR_CC_generic  returns  a  1xlength(mm)
#'     matrix giving the  C function defined in  the paper.  UMR_curv_generic
#'     is    returning    also    a   1xlength(mm)    matrix    giving    the
#'     (d^2/dtheta^2)(objective function), where "theta" is as defined in the
#'     paper.   [This  is  mathfrak{C}  in the  paper.]   These  are  similar
#'     quantities, the "curv" quantity is just C rescaled by the weight.  See
#'     calculations  in  paper.   The  more substantive  difference  is  that
#'     UMR_CC_generic requires  a closed  form for  the "D"  function whereas
#'     UMR_curv_generic simply  uses the hessian computation  (i.e., requires
#'     B', the derivative of the "B"  function).  (The closed form of the "D"
#'     function can be found  from the closed form of the  hessian, but it is
#'     not necessary.)
#'
#'      UMR_curv_generic2  is analogous  to UMR_curv_generic  but the  latter
#'      relies on UMR_CC_generic.
#'
#'      UMR_CC_generic1 is  analogous to UMR_CC_generic  (aka CC_SIR_generic)
#'      but the former is calculated  in fashion identical to UMR_curv_generic
#'      (i.e., relying on UMRhess).
#'
#'      DDfunc_Gauss_generic is the "D" function that can be passed in (after
#'      substituting  for  sig) for  DDfunc  in  various other  functions  to
#'      compute the "C" function (e.g., UMR_CC_generic).
#'     
#'      Note: "CC" and "DD", etc., refer to the "C" or "D" functions.  Double
#'      lettering is  a convention  often used  in the code  to refer  to the
#'      single letter.
#' 
#'  

#' @rdname UMR_CC_generic
#' @name UMR_curv_generic
UMR_curv_generic <- function(yy, mm,
                             ww_y = rep(1/length(yy), length(yy)),
                             ww_m = rep(1/length(mm), length(mm)),
                             densfunc,
                             BBpfunc
                             ){
    mm2 = rep(mm, each=2) ## I guess Not good if one of mm is a 'singleton'
    ww_m2 = rep(ww_m, each=2) / 2
    nn2  <- length(mm2)    
    myhess <- UMRhess(mm2, ww_m2, yy, ww_y,
                              dens=densfunc,
                              BBp=BBpfunc)
    vv <- c(1,-1) / sqrt(2)
    curv <- matrix(data=NA,nrow=1, ncol=nn2/2)
    for (ii in 1:(nn2/2)){
        i1 = 2*(ii-1) + 1
        i2 = 2*(ii-1) + 2
        lochess = myhess[i1:i2, i1:i2]
        curv[1,ii] =     t(vv) %*% lochess %*% vv;
    }
    curv  ## corre ct output
    ## curv / ww_m ## testing somthing!
}

#' @rdname UMR_CC_generic
#' @name UMR_curv_generic2
UMR_curv_generic2 <- function(yy,
                              mm,
                              ww_y=rep(1/length(yy), length(yy)),
                              ww_m = rep(1/length(mm), length(mm)),
                              densfunc, ## "phi"
                              DDfunc){
    UMR_CC_generic(yy,
                   mm,
                   ww_y,
                   ww_m,
                   densfunc, ## "phi"
                   DDfunc) * ww_m; ## multiplying by ww_m may slow down.
}


#' @rdname UMR_CC_generic
#' @name UMR_CC_generic
## here we take m_i == m_j 
## hess_SIR_generic
## "C" function defined in paper
UMR_CC_generic <- CC_SIR_generic <- function(yy,
                           mm,
                           ww_y=rep(1/length(yy), length(yy)),
                           ww_m = rep(1/length(mm), length(mm)),
                           densfunc, ## "phi"
                           DDfunc
                           ){
    ## error checks
    if (!is.null(ww_y)){
        if (!isTRUE(all.equal(sum(ww_y), 1))) stop("Need sum(ww_y) == 1")
    } else ww_y <- rep(1/length(yy), length(yy));
    if (!isTRUE(all.equal(sum(ww_m), 1))) stop("Need sum(ww_m) == 1")
    if (!all(ww_y>=0 & ww_y <=1))
        stop("Need ww_y entries in [0,1]")
    if (!all(ww_m>=0 & ww_m <=1))
        stop("Need ww_m entries in [0,1]")
    yymm <- outer(as.vector(yy), as.vector(-mm), "+")
    res1 <- matrix(ww_y, nrow=1) %*% densfunc(yymm)
    mmmm <- outer(as.vector(mm), as.vector(mm), "-")
    res2 <-    matrix(ww_m, nrow=1) %*% DDfunc(mmmm)
    res1-res2
}


## analogous to UMR_curv_generic
UMR_CC_generic1 <- function(yy, mm,
                             ww_y = rep(1/length(yy), length(yy)),
                             ww_m = rep(1/length(mm), length(mm)),
                             densfunc,
                             BBpfunc
                             ){
    mm2 = rep(mm, each=2) ## I guess Not good if one of mm is a 'singleton'
    ww_m2 = rep(ww_m, each=2) / 2
    nn2  <- length(mm2)    
    myhess <- UMRhess(mm2, ww_m2, yy, ww_y,
                              dens=densfunc,
                              BBp=BBpfunc)
    vv <- c(1,-1) / sqrt(2)
    curv <- matrix(data=NA,nrow=1, ncol=nn2/2)
    for (ii in 1:(nn2/2)){
        i1 = 2*(ii-1) + 1
        i2 = 2*(ii-1) + 2
        lochess = myhess[i1:i2, i1:i2]
        curv[1,ii] =     t(vv) %*% lochess %*% vv;
    }
    ##curv  ## "curv" output
    curv / ww_m ##  "CC" output
}

DDfunc_Gauss_generic <- function(dd, sig){
    exp(-dd^2 / (4*sig^2)) / (2 * sig * sqrt(pi))
}

