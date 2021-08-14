
#' @title Compute Unlinked Monotone Regression objective function numerically
#'
#' @export objective_fn_numint
#' 
#' @param  yy   Y   (response)  observation   vector  (numeric   vector).
#'     Alternatively,  yy  may be  an  ecdf,  i.e.  ecdf(yy)  or  getEcdf(yy,
#'     weights).
#' @param  mm  Current  (unsorted) estimate/iterate  at  which  to  compute
#'     gradient.  (Length  is <=  than the  number of  X observations  in the
#'     problem).
#' @param  ww_y Weights (nonnegative,  sum to  1) corresponding to  yy.  Same
#'     length as yy.  Default is just  1/length(yy) for each value.  If yy is
#'     non-numeric i.e. yy is an ecdf() then ww_y is ignored.
#' @param  ww_m Weights (nonnegative,  sum to  1) corresponding to  mm.  Same
#'     length as mm.
#'
#' @param  Phi This is the error (cumulative) distribution function, a  function object (Balabdaoui, Doss, Durot (2020+).  Function accepting vector or matrix arguments.
#'
#' @param subdivisions Passed argument to integrate().
#'
#' @details
#'
#'
#'
#' See paper for derivations. 
#'





objective_fn_numint <- function(mm, ww_m=NULL, yy, ww_y=NULL,  Phi,
                                subdivisions=1000L){
    integrand <- objective_integrand(mm, ww_m, yy, ww_y, Phi)    
    integrate(integrand,
              subdivisions = subdivisions,
              -Inf, Inf)$value
}


## need verify the code for approxfun to know the order of magnitude speed
## but i think can do n log(n).
objective_integrand <- function(mm, ww_m=NULL, yy, ww_y=NULL, Phi){
    stopifnot(length(ww_m)==length(mm) || is.null(ww_m))
    if (is.numeric(yy)){
        if (is.null(ww_y)) HH <- ecdf(yy)
        else HH <- getEcdf(yy,ww_y)
    }
    else HH <- yy
    if (is.null(ww_m))
        ww_m  <- rep(1 / length(mm), length(mm))
    ##GG <- function(zz){Phi(outer(zz, mm, "-")) / length(mm)}
    ##else
    GG <- function(zz){Phi(outer(zz, mm, "-")) %*%
                           as.matrix(ww_m, ncol=1)}
    ## GG <- Vectorize(function(zz){ sum(ww_m * Phi(zz - mm)) })
    ##    sweep(Phi(outer(zz, mm, "-")), 2, ww_m, "*")
    function(zz){(HH(zz) - GG(zz))^2}
}

#### backup
## ## need verify the code for approxfun to know the order of magnitude speed
## ## but i think can do n log(n).
## objective_integrand <- function(mm, ww_m=NULL, yy, ww_y=NULL, Phi){
##     stopifnot(length(ww_m)==length(mm) || is.null(ww_m))
##     if (is.numeric(yy)){
##         if (is.null(ww_y)) HH <- ecdf(yy)
##         else HH <- getEcdf(yy,ww_y)
##     }
##     else HH <- yy
##     if (is.null(ww_m)) ## ww_m  <- rep(1 / length(mm), length(mm))
##         GG <- function(zz){Phi(outer(zz, mm, "-")) / length(mm)}
##     else GG <- function(zz){Phi(outer(zz, mm, "-")) %*% as.matrix(ww_m, ncol=1)}
##     ## GG <- Vectorize(function(zz){ sum(ww_m * Phi(zz - mm)) })
##     ##    sweep(Phi(outer(zz, mm, "-")), 2, ww_m, "*")
##     function(zz){(HH(zz) - GG(zz))^2}
## }

## unneeded
objective_integrand_old <- function(mm, ww_m=NULL, yy, ww_y=NULL, Phi){
    if (is.null(ww_y)) HH <- ecdf(yy)
    else HH <- getEcdf(yy,ww_y)
    if (is.null(ww_m)) ww_m  <- rep(1 / length(mm), length(mm))
     GG <- Vectorize(function(zz){ sum(ww_m * Phi(zz - mm)) })
    function(zz){(HH(zz) - GG(zz))^2}
}


## ## Phi *does* need to take matrix argument 
## get_objective_integrand <- function(ww_m=NULL, yy, ww_y=NULL, Phi){
    
## }
