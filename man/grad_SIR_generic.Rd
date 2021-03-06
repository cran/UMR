% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grad_SIR_generic.R
\name{grad_SIR_generic}
\alias{grad_SIR_generic}
\title{@title Gradient of least-squares Shuffled Isotonic Regression criterion}
\usage{
grad_SIR_generic(
  yy,
  mm,
  counts = rep(1, length(mm)),
  AAfunc,
  BBfunc,
  rescale = FALSE
)
}
\arguments{
\item{yy}{Y (response) observation vector (numeric)}

\item{mm}{Current (unsorted) estimate/iterate at which to compute
gradient. (Length equals length of yy).}

\item{counts}{If the function that mm represents is piecewise constant, then mm may be passed in as only the unique entries.  In that case counts contains the number of times each element of mm is repeated.  Thus length(counts)==length(mm).  (Default for counts is thus a vector of all 1's.)}

\item{AAfunc}{This is the function "A" defined in the gradient
calculations in the paper (Balabdaoui, Doss, Durot (2020+).}

\item{BBfunc}{This is the function "B" defined in the gradient
    calculations in the paper (Balabdaoui, Doss, Durot (2020+).

 @details Returns gradient.  See calculations in ShuffReg.pdf.
    

 @examples
#### See help for gradDesc_PC, gradDesc, or grad_helpers}

\item{rescale}{Boolean: if False then the final return value is the
gradient; if True the final return value is gradient * length(yy) / 2.}
}
\description{
@title Gradient of least-squares Shuffled Isotonic Regression criterion
}
