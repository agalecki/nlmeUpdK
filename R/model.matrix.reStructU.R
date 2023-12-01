## -> model.matrix.reStruct.U function
#' This function is used to replace \code{\link{model.matrix.reStruct}} function in namespace of nlme package 
#'
#' This function is used to replace \code{\link{model.matrix.reStruct}} function in namespace of nlme package 
#'
#' @export
#' @param object an object inheriting from class \code{reStruct}, representing a random effects structure and consisting of a list of pdMat objects.
#' @param data a data frame in which to evaluate the variables defined in \code{formula(object)}.
#' @param contrast an optional named list specifying the contrasts to be used for representing the \code{factor} variables in data.
#'  The components names should match the names of the variables in data for which the contrasts are to be specified. The components of this list
#'  will be used as the contrasts attribute of the corresponding factor. If missing, the default contrast specification is used.
#' 
#' @return a matrix obtained by binding together, column-wise, the model matrices for each element of formula(object)..
#' @author Based on documentation for \code{model.matrix.reStruct} in \code{nlme} package by Jose Pinheiro and Douglas Bates.
##' @seealso \code{\link{model.matrix}}
##' @examples
##'  \dontrun{
##'   Pwr (fm1)
#' }

model.matrix.reStruct.U<-
function (object, data, contrast = NULL, ...) 
{
   .functionLabel <- "model.matrix.reStruct.U"           # Function label (recommended)
   .traceR <- attr(options()$traceR, "fun")
   .traceR <-  if (is.null(.traceR)) function(...){} else .traceR      
   
       
    .traceR(1, lbl = "-> model.matrix.reStruct.U STARTS")
    pdDef <- !(length(object) == 1 && inherits(object[[1]],"pdKronecker"))

    if (is.null(form <- formula(object, asList = TRUE))) {
        stop("Cannot extract model matrix without formula")
    }
    form1 <- asOneFormula(form)
    if (length(form1) > 0) {
        data <- model.frame(form1, data = data)
    }
    else {
        data <- data.frame("(Intercept)" = rep(1, nrow(data)))
    }
###
    .traceR(30)
    any2list <- function(object, data, contrast) {       # <===
        .traceR(940, lbl = "- any2list STARTS", store = FALSE)
        form2list <- function(form, data, contrast) {    # <===
        .traceR(950, lbl =  "- form2list STARTS", store = FALSE)

            if (length(asOneFormula(form)) == 0) {
                tt1 <- list("(Intercept)" = rep(1, dim(data)[1]))
                
                .traceR(951, lbl = " - EXIT951")   
                return(tt1)
            }
            tt2 <- as.data.frame(unclass(model.matrix(form, model.frame(form, 
                data), contrast)))  
                 .traceR(952, lbl = "- EXIT952")  
            tt2    
        }
      
        if (inherits(object, "formula")) {
            tt3 <- form2list(object, data, contrast)
            .traceR(953, lbl = "EXIT953")  
            return(tt3)
        }
        if (is.list(object)) {
            return(unlist(lapply(object, form2list, data = data, 
                contrast = contrast), recursive = FALSE))
        }
        .traceR(940, lbl = "ret940")
        return(NULL)
    }
###
    .traceR(50)
    if (pdDef){     # No pdKronecker
        value <- as.list(lapply(form, any2list, data = data, contrast = contrast))
        ncols <- as.vector(unlist(lapply(value, length)))
        nams <- if (length(value) == 1) {
        names(value[[1]])
    }
    else {
        paste(rep(names(value), ncols), unlist(lapply(value, 
            names)), sep = ".")
    }
    val <- matrix(unlist(value), nrow = nrow(data), dimnames = list(row.names(data), 
        nams))
    attr(val, "nams")  <- as.list(lapply(value, names))

}
    
    if (!pdDef){    # pdKroneck

       form2 <- formula(object, asList = FALSE)
       form2name <- names(form2) 
       form2 <- form2[[1]]
       form  <- form2   ####  <- 
       val <- model.matrix(form,data)
       ncols <- ncol(val)
       nams <- dimnames(val)[[2]]

       namsL <- list(nams)
       names(namsL) <- form2name 
       attr(val,"nams") <- namsL
       .traceR(804, lbl = "!pdDef")
    }    

    contr <- as.list(lapply(as.data.frame(data), function(x) if (inherits(x, 
        "factor") && length(levels(x)) > 1) 
        contrasts(x)
    else NULL))
    contr[names(contrast)] <- contrast

    attr(val, "ncols") <- ncols
    attr(val, "contr") <- contr
    attr(val, "names")  <- NULL
    .traceR(1, lbl = "model.matrix.reStruct.U ENDS <-")
    val
}
