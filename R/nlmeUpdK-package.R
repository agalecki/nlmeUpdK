### See: https://ebi-forecast.igb.illinois.edu/redmine/projects/pecan-1-2-5/wiki/Roxygen2
### See: https://github.com/yihui/roxygen2
#' Companion package to nlmeU. 
#' 
#' Companion package to nlmeU. It contains auxiliary functions pertaining to pdKronecker class introduced in Galecki and Burzykowski book (2013). 
#' Package under development.
#' 
#' @name nlmeUpdK-package
#' @aliases nlmeUpdK
#' @docType package
#' @title Auxiliary functions pertaining to pdKronecker class for Galecki and Burzykowski book 2013.
#' @author \email{agalecki@@umich.edu}, \email{Tomasz Burzykowski <tomasz.burzykowski@@uhasselt.be>}
#' @keywords package
#' @seealso \code{\link{package:nlmeU}}
NULL


#' @importFrom nlme "coef<-" "matrix<-" "Names<-" corMatrix isInitialized logDet Names pdConstruct
#' @importFrom nlme pdFactor pdMatrix VarCorr

#' @export          "coef<-" "matrix<-" "Names<-" corMatrix isInitialized logDet Names pdConstruct 
#' @export          pdFactor pdMatrix VarCorr

#' @importFrom nlme Dim asOneFormula getCovariateFormula splitFormula getResponseFormula pdMat
#' @export          Dim asOneFormula getCovariateFormula splitFormula getResponseFormula pdMat

########
#' @export pdKronecker
#' @export model.matrix.reStruct.U
#' @export unlistFun

###  #S3method "matrix<-", reStruct)
#' @S3method "["  pdKronecker
#' @S3method coef  pdKronecker
#' @S3method "coef<-"  pdKronecker
#' @S3method corMatrix  pdKronecker
#' @S3method formula  pdKronecker
#' @S3method isInitialized  pdKronecker
#' @S3method logDet  pdKronecker
#' @S3method "matrix<-"  pdKronecker
#' @S3method Names  pdKronecker
#' @S3method "Names<-"  pdKronecker
#' @S3method pdConstruct  pdKronecker
#' @S3method pdFactor  pdKronecker
#' @S3method pdMatrix  pdKronecker
#' @S3method solve  pdKronecker
#' @S3method summary  pdKronecker
#' @S3method VarCorr  pdKronecker
NULL