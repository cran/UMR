## Gradient descent with a fixed "degrees of freedom", i.e. number of constant pieces.
## Right now, does not save output history for memory conservation.


#' @title Gradient Descent with a fixed number of constant pieces (degrees of
#'     freedom)
#'
## export reconstruct_est
#' @export gradDesc_fixed_df
#' @importFrom stats median
#'
#'
#' @param yy Y (response) observation vector (numeric)
#' @param grad a function(yy, mm) where  mm may be shorter length than yy and
#'     is the previous iterate value (i.e., the estimate vector).
#' @param  init Initial value  of estimate  ('mm').  I.e., numeric  vector of
#'     length <= length(mm).  The output will be of length length(init).
#' @param counts Vector of length length(init); each entry indicates how many
#'     values of yy the corresponding  value of init (and output) corresponds
#'     to.  Alternatively,  can think of  counts as  a vector of  weights for
#'     each estimator value.
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
#' Prefer using UMRgradDesc_fixed_df now; this function deprecated.
#'
#' xxxx Implements a gradient descent.  See paper for details.  Right now
#'     stepsize is fixed.  Right now: init gets sorted in gradDesc_PC so does
#'     not need to be sorted on input.  Roughly, the difference between this
#'     algorithm and gradDesc() (which is just vanilla gradient descent on
#'     this problem) is that: if mm is the current value of the output
#'     estimate, then gradDesc_PC 'collapses' or combines values of mm that
#'     are (roughly, up to tolerance 'eps') equal.  Because the solution is
#'     generally piecewise constant with a relatively small number of
#'     constant regions this enormously speeds up the later stages of the
#'     algorithm.  Note that once points are combined/collapsed they
#'     contribute identically to the objective function, so they will never
#'     be "uncombined".
#' 


## deprecated
gradDesc_fixed_df <- function(yy,
                              grad,
                              init = stats::median(yy), counts = length(yy),
                              stepsize,  MM, tol=1e-7,
                              printevery=Inf, filename){
    mm <- sort(init)
    mmprev <- rep(Inf, length(mm))
    ii <- 1
    while (ii<=MM && sum(abs(mm-mmprev))>= tol){
        ii <- ii+1
        mmprev <- mm
        if ((ii %% printevery) == 0) {
            print(paste0("Starting ", ii, "th iteration."));
            mmhat_curr <- rep.int(mm, times=counts)
            save(yy,
                 ii,
                 stepsize, MM,  ## algorithm params
                 mmhat_curr,
                 file=filename)
        }
        mm <- mm - stepsize * grad(yy=yy, mm=mm, counts=counts)
    }
    mm ## don't sort
}



