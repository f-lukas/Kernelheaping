#' Kernel Density Estimation for Heaped Data
#'
#' In self-reported data the user often encounters heaped data, i.e. data which are rounded to a different
#' degree of coarseness. While this is mostly a minor problem in parametric density estimation the bias can be very large
#' for non-parametric methods such as kernel density estimation. This package implements a partly Bayesian algorithm 
#' treating the true, unknown values as additional parameters and estimates the rounding parameters 
#' to give a corrected kernel density estimate. It supports various standard bandwidth selection methods.
#' Additionally varying rounding probabilities (with the true value) and asymmetric rounding is estimable as well.
#'
#'
#' The most important function is \code{\link{dheaping}}. See the help and the attached examples on how to use the package.
#' @docType package
#' @name Kernelheaping
#' @import MASS evmix plyr corpcor
NULL