## Basic gradient descent implementation.
## Right now, does not save output history for memory conservation.


#' @title Basic gradient descent implementation
#'
#' @export
#'
#' @param yy Y (response) observation vector (numeric)
#' @param grad a function(yy, mm) where mm is same length of yy and is the
#'     previous iterate value (i.e., the estimate vector).
#' @param init Initial value of estimate ('mm').  I.e., numeric vector of same length as yy.
#' @param stepsize Gradient descent stepsize.  Set carefully! (I often use
#'     nn^2 / 2 where nn = length(yy), or nn if gradient is 'rescaled'.)
#' @param MM Number of iterations
#' @param printevery integer value (generally << MM).  Every 'printevery' iterations, a count will be printed and the output saved.
#' @param filename filename (path) to save output to.
#'
#'
#' @details Implements a very basic gradient descent. Right now stepsize is fixed.
#'
#'
#' @examples
#'
#'
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
#' MM <-  100 ## Total number iterations is MM * JJ 
#' JJ <- 2
#' eps <- (max(yy)-min(yy)) / (1000 * nn^(1/5) * myScale)
#' ## print *and* SAVE every 'printevery' iterations.
#' ## here no save occurs, printevery > MM
#' printevery <- 1000 
#' init <- yy
#' 
#' mmhat <- gradDesc(yy=yy, grad=mygradSIR, ## from settings file
#'                      init=init,
#'                      stepsize=stepsize, MM=MM,
#'                      printevery=printevery,
#'                      filename=paste0("../saves/", savefilenameUnique))
#' #### some classical/matched [oracle] estimators
#' isoreg_std <- Iso::ufit(y=yy, x=xx, lmode=Inf)
#' mmhat_std = isoreg_std$y ## Isotonic regression
#' linreg_std <- lm(yy~xx)
#' 

gradDesc <- function(yy, grad, init,
                     stepsize, MM,
                     printevery, filename){
    curr <- init
    for (ii in 2:MM){
        if ((ii %% printevery) == 0) {
            print(paste0("Completed ", ii, "th iteration."));
            save(yy,
                 ii,
                 ## mmhat_f, quantile_f,
                 ## errdist, mm0, ## model params
                 stepsize, MM,  ## algorithm params
                 ## outhist, ## waste of memory?
                 mmhat = curr,
                 file=filename)
        }
        curr <- curr - stepsize * grad(yy=yy, mm=curr)
    }
    sort(curr)
}
