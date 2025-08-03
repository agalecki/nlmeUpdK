### See: https://ebi-forecast.igb.illinois.edu/redmine/projects/pecan-1-2-5/wiki/Roxygen2
### See: https://github.com/yihui/roxygen2
#' Companion package to nlmeU. 
#' 
#' Companion package to nlmeU. It contains auxiliary functions pertaining to pdKronecker class introduced in Galecki and Burzykowski book (2013). 
#' Package under development.
#' 
#' @aliases nlmeUpdK
#'
#' nlmeUpdK: Companion Package to nlmeU for pdKronecker Class
#'
#' The `nlmeUpdK` package is a companion to the `nlmeU` package, providing auxiliary functions
#' for the `pdKronecker` class introduced in Galecki and Burzykowski (2013). It includes
#' functionality for constructing and manipulating positive-definite Kronecker product matrices
#' and related methods for random effects structures. The package is under development.
#'
#' @aliases nlmeUpdK
#' @title Auxiliary Functions for pdKronecker Class
#' @author Andrzej Galecki \email{agalecki@@umich.edu}, Tomasz Burzykowski \email{tomasz.burzykowski@@uhasselt.be}
#' @keywords package
#' @seealso \code{\link[nlmeU]{nlmeU}}
#' @references
#' Galecki, A., & Burzykowski, T. (2013). *Linear Mixed-Effects Models Using R: A Step-by-Step Approach*.
#' Springer.
#' @importFrom nlme coef<- matrix<- Names<- corMatrix isInitialized logDet Names pdConstruct pdFactor pdMatrix VarCorr
#' @importFrom nlme Dim asOneFormula getCovariateFormula splitFormula getResponseFormula pdMat
#' @importFrom stats coef contrasts formula model.frame model.matrix terms
#' @export coef<- matrix<- Names<- corMatrix isInitialized logDet Names pdConstruct pdFactor pdMatrix VarCorr
#' @export Dim asOneFormula getCovariateFormula splitFormula getResponseFormula pdMat
#' @export pdKronecker model.matrix.reStruct.U unlistFun
#' _PACKAGE

