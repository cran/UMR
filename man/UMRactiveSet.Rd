% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UMRactiveSet.R
\name{UMRactiveSet}
\alias{UMRactiveSet}
\title{An active set approach to minimizing objective in Unlinked Monotone
    Regression}
\usage{
UMRactiveSet(
  yy,
  grad,
  CC_SIR,
  init,
  counts = rep(1, length(init)),
  stepsize,
  MM,
  tol_end = 1e-04,
  tol_collapse,
  printevery,
  filename
)
}
\arguments{
\item{yy}{Y (response) observation vector (numeric)}

\item{grad}{a function(yy, mm) where mm  is the
previous iterate value (i.e., the estimate vector).}

\item{CC_SIR}{A curvature function object (denoted "C" in the paper).   See CC_SIR_generic() and examples.}

\item{init}{Initial value of estimate ('mm').  Vector, length may be different than length(yy). See 'counts' input.}

\item{counts}{Together 'init' and 'counts' serve as the initialization; the implied initial vector is rep.int(init, counts).}

\item{stepsize}{Gradient descent stepsize.}

\item{MM}{A number  of iterations.   May not  use them  all.  MM  is not
exactly the total  number of iterations used in the  sense that within
each of  MM iterations, we  will possibly run another  algorithm which
may take up to MM iterations (but usually takes many fewer).}

\item{tol_end}{Used as tolerance at various points .  Generally algorithm (and
some subalgorithms) end once sum(abs(mm-mmprev))  < tol, or you hit MM
iterations.}

\item{tol_collapse}{Collapsing roughly equal mm values into each other.}

\item{printevery}{integer  value (generally  << MM).   Every 'printevery'
iterations, a count will be printed and the output saved.}

\item{filename}{filename (path) to save output to.


param
ww_y Weights (nonnegative, sum to 1) corresponding to yy.  Same length as yy.}
}
\description{
An active set approach to minimizing objective in Unlinked Monotone
    Regression
}
\details{
Uses first order  (gradient) for optimization, and  uses certain
    second   derivative  computations   to  leave   saddle  points.    See
    Balabdaoui, Doss, and Durot (20xx).  Note that yy and mm (i.e., number
    covariates) may have different length.
}
