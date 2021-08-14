## Gradient descent for Piecewise Constant function
## Right now, does not save output history for memory conservation.


#' @title Gradient Descent implemented for Piecewise Constant functions
#'
#' @export UMRgradDesc_PC
#'
#'
#' @param yy Y (response) observation vector (numeric)
#' @param grad a function(yy, mm) where mm may be shorter length than yy and
#'     is the previous iterate value (i.e., the estimate vector).
#' @param init Initial value of estimate ('mm').  I.e., numeric vector
#'     usually of same length as yy.
#' @param stepsize Gradient descent stepsize.  Set carefully!
#' @param MM Number of iterations in which "support reduction" (combining of approximately equal values into a region of constancy) is done (see details and paper).
#' @param JJ Total number of gradient steps is MM*JJ.  JJ gradient steps are
#'     taken for each of the MM steps.
#' @param eps Roughly, points that are eps apart are considered to be equal
#'     and are thus collapsed into a single region of piecewise constancy of
#'     the output.  (This is not precisely true because one can have a long
#'     sorted-increasing vector of points that are each eps from their two
#'     neighboring points but such that the first and last points are not eps
#'     apart.  See algorithm description in paper for details. )
#' @param printevery integer value (generally << MM).  Every 'printevery'
#'     iterations, a count will be printed and the output saved.
#' @param filename path1/path2/filename to save output to.
#'
#'
#' @details Implements a gradient descent.  See paper for details.  Right now
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
#' @examples
#'
#' #'
#' #### Set up the gradient function 
#'  mysig <- 1 ##  std dev
#' errdist <- distr::Norm(0, sd=mysig)
#' modeldistname <- truedistname <- "Gauss" ## used for savefile name
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
#'  ## Now run the gradient descent 
#' savefilenameUnique <- paste("graddesc_", modeldistname, "_", truedistname,
#'                             "_n", nn,
#'                             "_", format(Sys.time(), "%Y-%m-%d-%T"), ".rsav", sep="")
#' print(paste("The unique save file name for this run is", savefilenameUnique))
#' stepsize <- nn^(1/2) ## Has to be tuned
#' MM <-  200 ## Total number iterations is MM * JJ 
#' JJ <- 2
#' eps <- (max(yy)-min(yy)) / (1000 * nn^(1/5) * myScale)
#' ## print *and* SAVE every 'printevery' iterations;
#' ## here no save occurs, printevery > MM
#' printevery <- 1000 
#' init <- yy
#' 
#' 
#'  mmhat <- UMRgradDesc_PC(yy=yy, grad=mygradSIR, ## from settings file
#'                       init=init,
#'                       stepsize=stepsize, MM=MM,
#'                       JJ=JJ, eps=eps,
#'                       printevery=printevery,
#'                       filename=paste0("../saves/", savefilenameUnique))
#' ####  some classical/matched [oracle] estimators
#' isoreg_std <- Iso::ufit(y=yy, x=xx, lmode=Inf)
#' mmhat_std = isoreg_std$y ## Isotonic regression
#' linreg_std <- lm(yy~xx)
#' 



UMRgradDesc_PC <-
    gradDesc_PC <- function(yy, grad, init,
                        stepsize, MM,
                        eps,
                        JJ=50, ## integer, number grad steps per
                               ## grouping/combining step ; MM*JJ is total
                               ## num iterations
                        printevery, filename){
    mm <- sort(init)
    nn_i <- length(mm)
    counts <- rep(1, nn_i)
    for (ii in 2:MM){
        if ((ii %% printevery) == 0) {
            print(paste0("Starting ", ii, "th iteration."));
            mmhat_curr <- rep.int(mm, times=counts)
            save(yy,
                 ii,
                 stepsize, MM, JJ, eps,  ## algorithm params
                 ## outhist, ## waste of memory?
                 ## mmhat = mm, ##   replace this 'mm' in save by:  rep.int(mm, times=counts)
                 mmhat_curr,
                 file=filename)
        }

        ## ###### Group (approximately) non-unique entries
        begidx <- 1
        newidx <- 1;
        mm_new <- counts_new <- rep(NA, nn_i)
        for (jj in 2:(nn_i+1)){
            if ((jj==nn_i+1) || ((mm[jj] - mm[begidx]) > eps)){
                mm_new[newidx] <- mean(mm[begidx:(jj-1)])
                counts_new[newidx] <- sum(counts[begidx:(jj-1)])
                begidx <- jj
                newidx <- newidx+1
            }
        }
        nn_i  <- sum(!is.na(mm_new));
        mm <- mm_new[1:nn_i]
        counts <- counts_new[1:nn_i];

        ## ######Do gradient step(s)
        for (jj in 1:JJ)
            mm <- mm - stepsize * grad(yy=yy, mm=mm, counts=counts)
    }
    ## ######## reconstruct estimator
    mm <- sort(mm) ## mm should be already sorted except for last grad steps
    rep.int(mm, times=counts)
}



