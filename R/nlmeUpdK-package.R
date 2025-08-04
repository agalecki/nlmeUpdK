#' nlmeUpdK: Companion Package to nlmeU for pdKronecker Class
#'
#' The nlmeUpdK package is a companion to the nlmeU package, providing auxiliary functions
#' for the pdKronecker class introduced in Galecki and Burzykowski (2013). It includes
#' functionality for constructing and manipulating positive-definite Kronecker product matrices
#' and related methods for random effects structures. The package is under development.
#'
#' Upon loading, nlmeUpdK modifies the nlme package namespace by replacing the internal
#' function model.matrix.reStruct with model.matrix.reStruct.U to support enhanced handling
#' of random effects structures, particularly for pdKronecker objects. This modification is
#' session-specific and reversible upon unloading the package.
#'
#' @aliases nlmeUpdK nlmeUpdK-package
#' @title Auxiliary Functions for pdKronecker Class
#' @author Andrzej Galecki \email{agalecki@@umich.edu}, Tomasz Burzykowski \email{tomasz.burzykowski@@uhasselt.be}
#' @keywords package
#' @seealso \code{\link[nlmeU]{nlmeU}}
#' @references
#' Galecki, A., & Burzykowski, T. (2013). *Linear Mixed-Effects Models Using R: A Step-by-Step Approach*.
#' Springer.
#' @importFrom nlme coef<- matrix<- Names<- corMatrix isInitialized logDet Names pdConstruct pdFactor pdMatrix VarCorr Dim asOneFormula getCovariateFormula splitFormula getResponseFormula pdMat
#' @importFrom stats coef contrasts formula model.frame model.matrix terms
#' @export pdKronecker model.matrix.reStruct.U unlistFun
#' @name nlmeUpdK-package
"_PACKAGE"