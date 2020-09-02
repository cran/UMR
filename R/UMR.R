#' UMR: For computing an estimator in Unlinked Monotone Regression.
#'
#'
#' A package for computing (via first order, gradient descent type
#' algorithms) an estimator in the problem of univariate Unlinked Monotone
#' Regression. See Balabdaoui, Doss, and Durot (2020+).
#'
#' @section UMR functions:
#'
#' The main functions are gradDesc_PC (for Gradient Descent for Piecwise
#' Constant functions) and gradDesc.  The former is faster and recommended.
#' The latter is the more naive vanilla gradient descent method (can be used
#' for instance to double check results from gradDesc_PC).
#'
#' @docType package
#'
#' @name UMR
NULL
