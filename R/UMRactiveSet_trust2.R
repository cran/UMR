## "Active Set" approach to finding (local) minima, in unlinked monotone regression.

## active set trust 2: difference between this and original is that '2' moves in *all* directions at once rather than just in one.  Also it checks for whether or not a given direction of curvature is a 'singleton'.  Not really sure how to thnk about this case but I guess for now may not move in that direction.

## This variant seems much slower (than UMRactiveSet_trust() original).
## Haven't actually studied why.  But the clear/obvious reason seems that
## dimensionality probably grows unchecked.



## ## main to-dos: 1) stepsize, modify the scaling? (2) subdivisions ... not
## ## sure how to deal with integrate() failure.  Maybe modifying the objective function scaling will help.  Would have to modify grad and hess though.  (3) Algorithmic question of
## ## doing one cc-step at a time or multiples.  Need to do comparisons.



#' @title An active set approach to minimizing objective in Unlinked Monotone
#'     Regression
#'
#' @importFrom trust trust
#' @importFrom stats ecdf
#' 
#'
#' @param yy Y (response) observation vector (numeric)
#'
#' @param ww_y  weight vector corresponding to  yy of same length  as yy.  If
#'     NULL then all yy entries get weight 1/length(yy).
#' 
#' @param grad a function(yy, mm) where mm  is the
#'     previous iterate value (i.e., the estimate vector).
#' @param UMR_curv  A curvature function object (giving mathfrak(C) in the paper; and related to "C" in the paper).   See UMR_curv_generic() and examples.  This is generally a "curried" version of UMR_curv_generic with densfunc and BBp passed in.
#' @param init Initial value of estimate ('mm').  Vector, length may be different than length(yy). See 'counts' input.
#' @param counts Together 'init' and 'counts' serve as the initialization; the implied initial vector is rep.int(init, counts).
#' @param stepsize Stepsize for moving out of saddle points.
#' 
#'
#' @param MM  A number  of iterations.   May not  use them  all.  MM  is not
#'     exactly the total  number of iterations used in the  sense that within
#'     each of  MM iterations, we  will possibly run another  algorithm which
#'     may take up to MM iterations (but usually takes many fewer).
#' @param tol_end  Used as tolerance at various points .  Generally algorithm (and
#'     some subalgorithms) end once sum(abs(mm-mmprev))  < tol, or you hit MM
#'     iterations.
#'
#' @param tol_collapse Collapsing roughly equal mm values into each other.
#' 
#' @param  printevery integer  value (generally  << MM).   Every 'printevery'
#'     iterations, a count will be printed and the output saved.
#' @param filename filename (path) to save output to.
#'
#' @param CDF This is the error (cumulative) distribution function, a  function object.  Function accepting vector or matrix arguments.
#' 
#'
#'
#'
#' @param grad Is function(mm, ww_m).  (Will be defined based on yy [and maybe ww_y] before being passed in.)  Returns vector of length(mm).  Gradient of objective function.
#' @param hess Is function(mm, ww_m). (Will be defined based on yy [and maybe ww_y] before being passed in.)  Returns matrix of dimensions length(mm) by length(mm).  Hessian of objective function.
#'
#' @param  ww_y Weights (nonnegative,  sum to  1) corresponding to  yy.  Samelength as yy.  Or NULL in which yy are taken as being evenly weighted.
#'
#'
#' @details  Uses first order  (gradient) for optimization, and  uses certain
#'     second   derivative  computations   to  leave   saddle  points.    See
#'     Balabdaoui, Doss, and Durot (2021).  Note that yy and mm (i.e., number
#'     covariates) may have different length.
#' 
#'
#'
#' 




## Need think on stepsize /nny (curvature 'step')

##unclear to me if "counts" or "weights" is right way to think about the
##estimator

## this is somehow much slower than the other version which does one
## direction at a time.


#'
#' #### dens and bbp are deprecated
#' 
#' param  dens This is the error density, a  function object.  Function accepting vector or matrix arguments.
#' 
#' param  BBp This  is derivative of  "B" function ("B  prime"), where  B is
#'     defined in the paper (Balabdaoui, Doss, Durot (2020+)).  Function accepting vector or matrix arguments.

UMRactiveSet_trust2 <- function(yy,
                               ww_y = NULL,
                               grad,
                               hess,
                               ## CC_SIR,
                               UMR_curv,
                               CDF,
                               ## dens,
                               ## BBp,
                               init,
                               counts = rep(1, length(init)),
                               stepsize, MM, tol_end=1e-4, tol_collapse,
                               printevery, filename){


    ## mmprev and mmcurr_full are only used for stopping conditions
    mmprev <- rep(Inf, length(init))
    mmcurr_full <- rep(0, length(init))


    
    mmcurr <- init
    nnx <- sum(counts)
    nny <- length(yy)
    
    stopifnot(length(ww_y)==length(yy) || is.null(ww_y))
    yyord <- order(yy)
    ## am not sure if ecdf() (taking null ww_y) is more efficient than
    ## getEcdf (taking non-null ww_y; am allowing to pass on the null ww_y
    ## for now.
    if (is.null(ww_y)) {
        yy_ecdf <- stats::ecdf(yy)
        ww_y <- rep(1/length(yy), length(yy)) ## wont affect ecdf
    }        
    else {
        stopifnot(length(yy)==length(ww_y))
        yy_ecdf <- getEcdf(yy,ww_y)
    }
    ww_y <- ww_y[yyord] 
    yy <- yy[yyord]
    prevobjval <- objval <-  Inf

    


    
    ii <- 1
    while (ii<= MM && sum(abs(mmcurr_full-mmprev))>= tol_end){
    ##while (ii==1 || (ii<= MM && (prevobjval-objval)>= tol_end)){
        ## can't do Inf-Inf 
        
        mmprev <- rep.int(mmcurr, times=counts)

        prevobjval <- objval

        
        ## mmprev <- mmcurr_full
        if ((ii %% printevery) == 0) {
            print(paste0("Completed ", ii, "th iteration."));
            save(yy,
                 ii,
                 stepsize, MM,  ## algorithm params
                 mmhat = mmcurr,
                 file=filename)
        }

        ww_m <- counts/sum(counts)

        ## ## setup for trust().
        myobj <- function(mm){
            objective_fn_numint(mm, ww_m=ww_m,
                                ## yy=yy, ww_y=ww_y,
                                yy=yy_ecdf,
                                Phi=CDF,
                                subdivisions= max(300, 3* max(nnx,nny))
                                ##subdivisions=2000
                                )
        }
        mygrad <- function(mm){
            grad(mm=mm, ww_m=ww_m)
        }

        myhess <- function(mm){
            hess(mm, ww_m=ww_m)
        }

        myobjfun <- function(mm){
            list(value=myobj(mm), gradient=mygrad(mm), hessian=myhess(mm))
        }

        ## mmcurr <- gradDesc_fixed_df(yy, grad,
        ##                             init=mmcurr,
        ##                             counts=counts,
        ##                             stepsize=stepsize,
        ##                             MM=ceiling(sqrt(MM)),
        ##                             tol=tol_end,
        ##                             printevery=printevery, filename=filename)


        ##        test the below / cmprae to aboev
        
        out <- trust::trust(objfun=myobjfun,
                            parinit=mmcurr,
                            rinit=5, ## NO IDEA need to plot
                            rmax=Inf,
                            ## parscale=c(1,3,6),
                            ## iterlim=ceiling(sqrt(MM)), ##arbitrary
                            iterlim=200, ##arbitrary
                            minimize=TRUE,
                            ##blather=TRUE
                            blather=FALSE
                            )
        objval <- out$value
        mmcurr <- out$argument

        ## ##### Currently have two sets of code for collapsing non-unique
        ## ##### entries.  Think I only need the latter?
        

        ## ## The "collapse" non-unique entries / "activate constraints" step
        {
            ## sort needed for simplifying vector.  Unclear algorithmically if this
            ## (probabilistically) is the best thing to do or if its better to
            ## just let the length grow over time)
            neword <- order(mmcurr)
            mmcurr <- mmcurr[neword]
            counts <- counts[neword]
            mm_active <- rle(mmcurr)
            mmcurr <- mm_active$values
            metacounts <- mm_active$lengths
            inds <- cumsum(metacounts)
            pp <- length(inds)
            indsstart <- c(0, inds[-pp])+1
            countidcs <- mapply(":", indsstart, inds)
            ## accumulate counts
            counts <- sapply(countidcs, function(xx, bb){sum(bb[xx])}, counts)
        }

        
        ## ###### Group (approximately) non-unique entries
        begidx <- 1
        newidx <- 1;
        nn_i <- length(mmcurr)        
        mm_new <- counts_new <- rep(NA, nn_i)
        for (jj in 2:(nn_i+1)){
            if ((jj==nn_i+1) || ((mmcurr[jj] - mmcurr[begidx]) > tol_collapse)){
                mm_new[newidx] <- mean(mmcurr[begidx:(jj-1)])
                counts_new[newidx] <- sum(counts[begidx:(jj-1)])
                begidx <- jj
                newidx <- newidx+1
            }
        }
        nn_i  <- sum(!is.na(mm_new));
        mmcurr <- mm_new[1:nn_i]
        if (sum(counts_new[!is.na(counts_new)]) != nnx){
            print(counts_new)
            print(counts)
        }
        counts <- counts_new[1:nn_i];

        ## if (sum(counts/nnx) != 1){
        ##     print(counts)
        ## }
        ## if (length(counts) != length(mmcurr)){
        ##     print(counts)
        ##     print(mmcurr)
        ##     stop("length(counts) != length(mmcurr)")
        ## }

        
        curv <- UMR_curv(yy=yy, mm=mmcurr,
                         ww_y=ww_y,
                         ww_m=counts/nnx)
        ## disqualify anyone with count of 1; (could change to 'weights' version)
        curv[counts==1]  <- 0;
        minidx <- which.min(curv)
        if (curv[minidx] >= 0)
            break;
        
        {
            negidcs <- curv<0
            posidcs <- !negidcs
            numnegs <- sum(negidcs)
            evens <- (1:numnegs)*2
            odds <- ((1:numnegs)*2)-1

            newmm <- rep(mmcurr[negidcs], each=2)
            newcounts <- rep(counts[negidcs], each=2)
            newmm[evens] <- newmm[evens] + sqrt(stepsize/nnx)
            newmm[odds] <- newmm[odds] - sqrt(stepsize/nnx)
            mmcurr <- (c(mmcurr[posidcs], newmm))
            ord <- order(mmcurr)
            newcounts[evens] <- floor(newcounts[evens] / 2);
            newcounts[odds]  <- ceiling(newcounts[odds] / 2);
            counts <- c(counts[posidcs], newcounts)
          
            ## mmcurr <- (c(mmcurr[posidcs], newmm))
            mmcurr <- mmcurr[ord]
            counts <- counts[ord]

            ## ## need to do another round of sorting for comparison with
            ## ## previous iterate below ...
            ## neword <- order(mmcurr)
            ## mmcurr <- mmcurr[neword]
            ## counts <- counts[neword]

            if (sum(counts) != nnx){
                print(newcounts)
                print(counts)
                stop("sum(counts) != nnx")
            }
            if (length(counts) != length(mmcurr))        {
                print(counts)
                print(mmcurr)
                stop("length(counts) != length(mmcurr)")
            }

        }

        {
        
        ## curvlen <- length(curv)
        ## ## the following code is inefficient
        ## for (kk in 1:curvlen){
        ##     ## ## take one step in " negatively curved directions" and
        ##     ## ## then iterate
        ##     pp <- length(mmcurr)
        ##     if (curv[kk] < 0){
        ##     }
            
        ##     ## double up the minimum index
        ##     mmcurr <- c(mmcurr[1:minidx], mmcurr[minidx:pp])
        ##     counts <- c(counts[1:minidx], counts[minidx:pp])
        ##     ## take step
        ##     counts[minidx] <- floor(counts[minidx] / 2);
        ##     counts[minidx+1] <- ceiling(counts[minidx+1] / 2);
        ##     ##         stopifnot( counts[minidx])
        ##     ## mmcurr[minidx] <- mmcurr[minidx] - stepsize/nnx;
        ##     ## mmcurr[minidx+1] <- mmcurr[minidx+1] + stepsize/nnx;
        ##     mmcurr[minidx] <- mmcurr[minidx] - sqrt(stepsize/nnx);
        ##     mmcurr[minidx+1] <- mmcurr[minidx+1] + sqrt(stepsize/nnx);
            ## }

        }


        ## {
        ## ## ## take one step in "most negatively curved direction" and
        ## ## ## then iterate
        ## pp <- length(mmcurr)
        ## ## double up the minimum index
        ## mmcurr <- c(mmcurr[1:minidx], mmcurr[minidx:pp])
        ## counts <- c(counts[1:minidx], counts[minidx:pp])
        ## ## take step
        ## counts[minidx] <- floor(counts[minidx] / 2);
        ## counts[minidx+1] <- ceiling(counts[minidx+1] / 2);
        ## if (counts[minidx] <= 0) stop("Have counts <=0")
        ## ##         stopifnot( counts[minidx])
        ## ## mmcurr[minidx] <- mmcurr[minidx] - stepsize/nnx;
        ## ## mmcurr[minidx+1] <- mmcurr[minidx+1] + stepsize/nnx;
        ## mmcurr[minidx] <- mmcurr[minidx] - sqrt(stepsize/nnx);
        ## mmcurr[minidx+1] <- mmcurr[minidx+1] + sqrt(stepsize/nnx);
        ## }

        
        ## for comparison with mmprev.  Not sure if this does or
        ## doesn't slow anything down (if nnx equals nny then
        ## shouldn't be dramatic slowdown).
        mmcurr_full <- rep.int(mmcurr, counts)
        ii <- ii+1                
    }

    

    
    mmord <- order(mmcurr)
    mmcurr <- mmcurr[mmord]
    counts <- counts[mmord]
    
    res <-list(mm=mmcurr,
               counts=counts,
               mm_full =rep.int(mmcurr, times=counts))
    return(res);
}

