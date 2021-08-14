## Gradient descent with a fixed "degrees of freedom", i.e. number of constant pieces.
## Right now, does not save output history for memory conservation.


#' @title Gradient Descent with a fixed number of constant pieces (degrees of
#'     freedom)
#'
#' 
#' @export UMRgradDesc_fixed_df 
#' 
#'
#' @param grad a function(mm) where mm
#'     is the previous iterate value (i.e., the estimate vector).
#' @param  init Initial value  of estimate  ('mm').   The output will be of length length(init).
#' @param stepsize Gradient descent stepsize.  Set carefully!
#' @param MM Number of iterations  in which "support reduction" (combining of
#'     approximately equal  values into a  region of constancy) is  done (see
#'     details and paper).  Depending on tol, may not use all MM iterations.
#' @param tol Tolerance: end algorithm  once sum(abs(mm-mmprev)) < tol or you
#'     hit MM iterations.
#' @param  printevery integer  value (generally  << MM).   Every 'printevery'
#'     iterations, a count will be printed and the output saved.
#' @param filename path1/path2/filename to save output to.
#'
#'
#' @details
#'
#' UMRgradDesc_fixed_df does a gradient descent with a fixed (upper bound) on the number of constant segments of the function. 
#'
#' Output of UMRgradDesc_fixed_df  is unsorted.  Note weights for  'mm' are not
#' passed in; rather they will be contained/used in grad().
#'
#' 



## does NOT sort the output; this way 'counts' still aligns with the output. 
UMRgradDesc_fixed_df <-  function(grad,
                                  init,
                                  stepsize,  MM, tol=1e-7,
                                  printevery=Inf, filename){
    ## myord <- order(init)
    ## mm <- init[myord]
    ## ww_m <- ww_m[myord]
    mm <- init;
    mmprev <- rep(Inf, length(mm))
    ii <- 1
    while (ii<=MM && sum(abs(mm-mmprev))>= tol){
        ii <- ii+1
        mmprev <- mm
        if ((ii %% printevery) == 0) {
            print(paste0("Starting ", ii, "th iteration."));
            save(mm,
                 ii,
                 stepsize, MM,  ## algorithm params
                 file=filename)
        }
        mm <- mm - stepsize * grad(mm=mm)
    }
    mm ## don't sort
}


reconstruct_est <- function(mm, counts){
    mm <- sort(mm) ## mm should be already sorted except for last grad steps
    rep.int(mm, times=counts)
}
