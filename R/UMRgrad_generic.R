## Generic gradient function for isotonic shuffled regression.


#' @title Gradient of least-squares Shuffled Isotonic Regression criterion
#'
#' @export UMRgrad_generic
#' @export grad_SIR_generic
#'
##  grad_SIR_generic is deprecated
#'
#' @param yy Y (response) observation vector (numeric)
#' @param mm Current (unsorted) estimate/iterate at which to compute
#'     gradient. (Length equals length of yy).
#' 
#' @param ww_y Weight vector for yy.
#' @param ww_m Weight vector for mm.
#' 
#' @param AAfunc This is the function "A" defined in the gradient
#'     calculations in the paper (Balabdaoui, Doss, Durot (2020+).
#' @param BBfunc This is the function "B" defined in the gradient
#'     calculations in the paper (Balabdaoui, Doss, Durot (2020+).
#'
#'  @details  Returns gradient as a  column matrix.  See calculations  in the
#'     paper.
#'     
#'
#'  @examples
#' #### See help for gradDesc_PC, gradDesc, or grad_helpers
#' 

#' @rdname UMRgrad_generic
#' @name UMRgrad_generic 
UMRgrad_generic <- function(yy,
                            ww_y = rep(1/length(yy), length(yy)),
                            mm,
                            ww_m = rep(1/length(mm), length(mm)),
                            AAfunc, BBfunc){
    ## nn <- length(yy);
    ## pp <- length(mm)
    res1 <- matrix(ww_y, nrow=1) %*% AAfunc(yy,mm)
    res2 <- matrix(ww_m, nrow=1) %*% BBfunc(mm)
    res <-  ww_m * 2 * t(res1 - res2) ## column matrix
    return(res)
}


#' @rdname UMRgrad_generic
#' @name grad_SIR_generic
## deprecated
#' @param counts If the function that mm represents is piecewise constant, then mm may be passed in as only the unique entries.  In that case counts contains the number of times each element of mm is repeated.  Thus length(counts)==length(mm).  (Default for counts is thus a vector of all 1's.)
#
#' @param rescale Boolean: if False then the final return value is the
#     gradient; if True the final return value is gradient * length(yy) / 2.
grad_SIR_generic <- function(yy, mm,
                                                counts = rep(1, length(mm)),
                                                AAfunc, BBfunc, rescale=FALSE){
    nn <- length(yy);
    res1 <- apply(AAfunc(yy,mm), 2, sum)
    ##    res2 <- matrix(counts, nrow=1) %*% BBfunc(mm, mm)
    res2 <- matrix(counts, nrow=1) %*% BBfunc(mm)
    res <- t(res1 - res2) ## column matrix
    if (rescale) return(res / nn)
    else return(res * 2 / nn^2)
}
