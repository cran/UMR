#' @title Compute Hessian of  Unlinked Monotone Regression objective function
#'     from Balabdaoui, Doss, and Durot
#'
#' @export UMRhess
#' @export UMRhess_generic
#'
#' @importFrom stats ecdf
#' @importFrom stats integrate
#' 
#' 
#'
#' @param yy Y (response) observation vector (numeric)
#' @param  mm  Current  (unsorted)  estimate/iterate  at  which  to  compute
#'     gradient.  (Length is  <= than  the number  of X  observations in  the
#'     problem).
#' @param  ww_y Weights (nonnegative,  sum to  1) corresponding to  yy.  Same
#'     length as yy.  Default is just 1/length(yy) for each value.
#' @param  ww_m Weights (nonnegative,  sum to  1) corresponding to  mm.  Same
#'     length as mm.
#'
#' @param  dens This is the error density, a  function object (Balabdaoui, Doss, Durot (2020+).  Function accepting vector or matrix arguments.
#' 
#' @param  BBp This  is derivative of  "B" function ("B  prime"), where  B is
#'     defined in the paper.  Function accepting vector or matrix arguments.
#' 
#'
#'
#' @details
#'
#'
#'
#' See paper for derivations. 
#'


#' @rdname UMRhess
#' @name UMRhess
UMRhess_generic <-
    UMRhess <- function(mm, ww_m, yy, ww_y=rep(1/length(yy),length(yy)), dens, BBp){
    ddmm <- outer(as.vector(mm), as.vector(mm), "-")
    BBp_mat <- 2 * BBp(ddmm)
    ## BBp_mat <- diag(ww_m) %*% BBp_mat %*% diag(ww_m) ## slow 
    BBp_mat <- sweep(BBp_mat, 1, ww_m, FUN="*")
    BBp_mat <- sweep(BBp_mat, 2, ww_m, FUN="*")
    diag(BBp_mat) <- 0 ## this 0'd out matrix is used to compute diagonal and
                       ## to compute off-diagonal!
    ## #####
    yymm <- outer(as.vector(yy), as.vector(mm), "-")    
    dens_mat <- 2 * dens(yymm);
    dens_mat <- sweep(dens_mat, 1, ww_y, FUN="*")
    dens_mat <- sweep(dens_mat, 2, ww_m, FUN="*") 
    ## dens_mat <- diag(ww_y) %*% dens_mat %*% diag(ww_m) ## slow
    ## #####
    resdiag <- (t(apply(dens_mat, 2, sum)) - apply(BBp_mat, 1, sum));
    diag(BBp_mat) <- resdiag
    BBp_mat
}


## do hessian calculations by integration; different (slower) formula for checking. 
UMRhess_ij <- hess_SIR_ij <- function(yy, mm, ww_m, ii, jj,
                                      Phi, phi, phipr){
    HH <- stats::ecdf(yy)
    GG <- Vectorize(function(zz){
        sum(ww_m*Phi(zz - mm))
    }) ## vector
    integrand <- function(zz){
        ## 2 *
        ((GG(zz) - HH(zz)) *
         (ww_m[ii] * phipr(zz-mm[ii]) + ww_m[jj]*phipr(zz-mm[jj]))
            + (ww_m[ii] * phi(zz-mm[ii]) - ww_m[jj] * phi(zz-mm[jj]))^2)
    }
    
    stats::integrate(integrand,
                     subdivisions=3000L,
                     -Inf, Inf)$value
}
