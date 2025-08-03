#  Some of the methods in this file need to be tested
### methods(class=pdKronecker)
#  [1] [.              coef.          coef<-.       
#  [4] corMatrix.*     formula.       isInitialized.
#  [7] logDet.*        matrix<-.Needed      Names.        
# [10] Names<-.*       pdConstruct.   pdFactor.     
# [13] pdMatrix.*      solve.         summary.      
# [16] VarCorr.*   
#' Unlist Names for pdKronecker
#'
#' Creates a character vector of names by combining names from a list using a Kronecker product-like structure.
#'
#' @param namesList A list of character vectors containing names.
#' @return A character vector of combined names.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
unlistFun <- function(namesList) {
    lenx <- sapply(namesList, length)
    namesList <- namesList[lenx>1]
    len1 <- length(namesList)
    namesList <- rev(namesList)
    
    dt <- expand.grid(namesList)
    allNames <- dt[[len1]]
    
    for (i in (len1-1):1) {
        aux <- dt[[i]]       
        allNames <- paste(allNames, dt[[i]], sep=":") 
    }
    
    meltx <- function(data, ...) data.frame(value = data)
    allNames <- meltx(allNames)
    allNames <- allNames[,"value"]
    allNames <- as.character(allNames)
    allNames
}

#' Subset pdKronecker Object
#'
#' Subsets a `pdKronecker` object, returning either a matrix or a new `pdKronecker` object
#' based on the indices provided.
#'
#' @param x A `pdKronecker` object.
#' @param i Row indices.
#' @param j Column indices.
#' @param drop Logical; if `TRUE`, drops dimensions when possible.
#' @param ... Additional arguments (currently ignored).
#' @return A matrix or a `pdKronecker` object.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
"[.pdKronecker" <- function(x, i, j, drop = TRUE, ...) {
    xx <- x
    x <- as.matrix(x)
    mCall <- match.call()
    mCall[[1]] <- get("[")
    mCall[["x"]] <- x
    mCall[["drop"]] <- drop
    if (length(i) == length(j) && mode(i) == mode(j) && all(i == j)) {
        mCall[["drop"]] <- FALSE
        val <- eval(mCall)
        vNames <- colnames(val)
        auxNames <- lapply(Names(xx, TRUE), function(el, vN) {
            aux <- match(vN, el)
            if (any(aux1 <- !is.na(aux))) {
                el[aux[aux1]]
            }
        }, vN = vNames)
        auxWhich <- !unlist(lapply(auxNames, is.null))
        if (sum(auxWhich) == 1) {
            return(pdConstruct(as.list(xx)[auxWhich][[1]], val))
        }
        auxNames <- auxNames[auxWhich]
        auxClass <- unlist(lapply(xx, function(el) class(el)[1]))[auxWhich]
        return(pdConstruct(xx, val, nam = auxNames, form = NULL, pdClass = auxClass))
    } else {
        eval(mCall)
    }
}

#' Extract Coefficients from pdKronecker
#'
#' Extracts coefficients from a `pdKronecker` object, optionally in unconstrained form.
#'
#' @param object A `pdKronecker` object.
#' @param unconstrained Logical; if `TRUE`, returns unconstrained coefficients.
#' @param ... Additional arguments (currently ignored).
#' @return A numeric vector of coefficients.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
coef.pdKronecker <- function(object, unconstrained = TRUE, ...) {
    coefSum <- coef(object[[1]], unconstrained)
    coefL <- lapply(object, coef, unconstrained)
    coefs <- unlist(lapply(coefL, FUN = function(x) x[-1]))
    coefAll <- c(coefSum, coefs)
    nms <- names(coefAll)
    dupi <- duplicated(nms)*(1:length(nms))
    idx <- dupi[dupi>0]
    nmsidx <- nms[idx]
    nms[idx] <- paste(nmsidx, idx, sep="")
    names(coefAll) <- nms
    coefAll
}

#' Assign Coefficients to pdKronecker
#'
#' Assigns new coefficients to a `pdKronecker` object.
#'
#' @param object A `pdKronecker` object.
#' @param value Numeric vector of coefficients.
#' @param ... Additional arguments (currently ignored).
#' @return The modified `pdKronecker` object.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
"coef<-.pdKronecker" <- function(object, ..., value) {
    if (is.null(plen <- attr(object, "plen"))) {
        stop(paste("Cannot change the parameter when", "length of parameters is undefined"))
    }
    plen <- plen -1
    plen[1] <-1
    
    if (length(value) != sum(plen)) {
        stop("Cannot change parameter length of initialized pdMat object")
    }
    ends <- cumsum(plen)
    starts <- 1 + c(0, ends[-length(ends)])
    
    attrx <- attributes(object)
    coef(object[[1]]) <- value[1] 
    for (i in 2:length(object)) {      
        coef(object[[i]]) <- c(0,value[(starts[i]):(ends[i])])
    }
    attributes(object) <- attrx
    object
}

#' Extract Correlation Matrix from pdKronecker
#'
#' Returns the correlation matrix of a `pdKronecker` object.
#'
#' @param object A `pdKronecker` object.
#' @param ... Additional arguments (currently ignored).
#' @return A correlation matrix with standard deviations as an attribute.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
corMatrix.pdKronecker <- function(object, ...) {
    if (!isInitialized(object)) {
        stop("Cannot access the matrix of uninitialized objects")
    }
    if (length(Names(object)) == 0) {
        stop("Cannot access the matrix of object without names")
    }
    namesList <- Names(object, TRUE)
    Ncol <- Dim(object)[2]
    value <- array(0, c(Ncol, Ncol), attr(object, "Dimnames"))
    stdDev <- double(Ncol)
    names(stdDev) <- colnames(value)
    for (i in seq_along(object)) {
        aux <- corMatrix(object[[i]])
        value[namesList[[i]], namesList[[i]]] <- as.vector(aux)
        stdDev[namesList[[i]]] <- attr(aux, "stdDev")
    }
    attr(value, "stdDev") <- stdDev
    value
}

#' Extract Formula from pdKronecker
#'
#' Returns the formula associated with a `pdKronecker` object, optionally as a list.
#'
#' @param x A `pdKronecker` object.
#' @param asList Logical; if `TRUE`, returns a list of formulas.
#' @param ... Additional arguments (currently ignored).
#' @return A formula or list of formulas.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
formula.pdKronecker <- function(x, asList = FALSE, ...) {
    val <- lapply(x, formula)
    isNULL <- unlist(lapply(val, is.null))
    if (all(isNULL)) 
        return(NULL)
    if (any(isNULL)) {
        stop("All elements must have formulas, when any has a formula.")
    }
    if (asList) 
        return(val)
    isTwoSided <- unlist(lapply(val, function(el) {
        inherits(el, "listForm")
    }))
    if (all(isTwoSided)) {
        val <- do.call("c", val)
        class(val) <- "listForm"
        return(val)
    }
    if (any(isTwoSided)) {
        stop(paste("All elements of formula must be list of two-sided formulae", 
            "or two-sided formulae"))
    }
    val <- lapply(rev(val), terms)
    aux <- paste(unlist(lapply(val, function(el) attr(el, "term.labels"))), 
        collapse = ":")
    if (!all(unlist(lapply(val, function(el) attr(el, "intercept"))))) {
        aux <- paste(aux, " - 1")
    }
    eval(parse(text = paste("~", aux)))
}

#' Check Initialization of pdKronecker
#'
#' Checks if a `pdKronecker` object is initialized.
#'
#' @param object A `pdKronecker` object.
#' @return Logical; `TRUE` if all components are initialized.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
isInitialized.pdKronecker <- function(object) {
    InitAll <- all(unlist(lapply(object, isInitialized)))
    InitAll
}

#' Compute Log-Determinant of pdKronecker
#'
#' Computes the log-determinant of a `pdKronecker` object.
#'
#' @param object A `pdKronecker` object.
#' @param ... Additional arguments (currently ignored).
#' @return Numeric; the log-determinant.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
logDet.pdKronecker <- function(object, ...) {
    sum(unlist(lapply(object, logDet)))
}

#' Assign Matrix to pdKronecker
#'
#' Assigns a matrix to a `pdKronecker` object, updating its component matrices.
#'
#' @param object A `pdKronecker` object.
#' @param value A matrix to assign.
#' @return The modified `pdKronecker` object.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
"matrix<-.pdKronecker" <- function(object, value) {
    value <- as.matrix(value)
    namesList <- Names(object, TRUE)
    Ncol <- Dim(object)[2]
    dims <- dim(value)
 
    if (!((dims[1] == dims[2]) && (dims[1] == Ncol))) {
        stop("Cannot change the number of columns on an initialized object")
    }   
    if (is.null(vNames <- rownames(value))) {
        vNames <- unlistFun(namesList)
        dimnames(value) <- list(vNames, vNames)
    } else {
        if (!(all(match(unlistFun(namesList), vNames, nomatch = 0)))) {
            stop("Names of object and value must match.")
        }
        attr(object, "Dimnames") <- list(vNames, vNames)
    }
    
    val1 <- as.matrix(value[1,1])
    nmsi <- namesList[[1]]   
    dimnames(val1) <- list(nmsi,nmsi)
    matrix(object[[1]]) <- val1
    
    len <- length(object)
    val <- value/val1[1,]
    
    for (i in len:2){
        nmsi <- namesList[[i]]
        leni <- length(nmsi)        
        nv <- nrow(val)
        idx <- seq(1,to =leni)
        val1 <- as.matrix(val[idx,idx])
        dimnames(val1) <- list(nmsi,nmsi)
        matrix(object[[i]]) <- val1
        idx <- seq(1, to =nv, by=leni)
        val <- as.matrix(val[idx,idx])/val1[1,1]
    }
     
    object
}

#' Extract Names from pdKronecker
#'
#' Returns the names of a `pdKronecker` object, optionally as a list.
#'
#' @param object A `pdKronecker` object.
#' @param asList Logical; if `TRUE`, returns a list of names.
#' @param ... Additional arguments (currently ignored).
#' @return A character vector or list of names.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
Names.pdKronecker <- function(object, asList = FALSE, ...) {
    if (asList) 
        attr(object, "namesList") 
    else 
        attr(object, "Dimnames")[[2]]
}

#' Assign Names to pdKronecker
#'
#' Assigns new names to a `pdKronecker` object.
#'
#' @param object A `pdKronecker` object.
#' @param value A character vector or list of names.
#' @param ... Additional arguments (currently ignored).
#' @return The modified `pdKronecker` object.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
"Names<-.pdKronecker" <- function(object, ..., value) {
    .functionLabel <- "Names<-.pdKronecker"                 # Function label (recommended)
    .traceR <- attr(options()$traceR, "fun")
    .traceR <-  if (is.null(.traceR)) function(...){} else .traceR      

    .traceR(210, lbl = "-> Names<-.pdKronecker starts")
    
    tmp <- Names(object)
  
    if (!is.null(Names(object))) { 
        clss <- class(object)
        class(object) <- "pdMat"
        Names(object) <- value
        class(object) <- clss
        .traceR(211, lbl = "Names<-.pdKronecker EXIT1")
        object
    } else {
        .traceR(212, lbl = "Names<-.pdKronecker EXIT2")
        object
    }
}

#' Extract Factor of pdKronecker
#'
#' Returns the factor of a `pdKronecker` matrix.
#'
#' @param object A `pdKronecker` object.
#' @return A matrix representing the factor.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
pdFactor.pdKronecker <- function(object) {
    pdMatrix(object, factor = TRUE)
}

#' Extract Matrix from pdKronecker
#'
#' Returns the matrix representation of a `pdKronecker` object, optionally as a factor.
#'
#' @param object A `pdKronecker` object.
#' @param factor Logical; if `TRUE`, returns the factor matrix.
#' @return A matrix.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
pdMatrix.pdKronecker <- function(object, factor = FALSE) {
    if (!isInitialized(object)) {
        stop("Cannot access the matrix of uninitialized objects")
    }
    if (length(Names(object)) == 0) {
        stop("Cannot access the matrix of object without names")
    }
    
    namesList <- Names(object, TRUE)
    Ncol <- Dim(object)[2]
    
    if (factor) {
        lD <- 0
    }
    
    value <- matrix(1)
    len <- length(object)
    for (i in 1:len) {
        aux <- pdMatrix(object[[i]], factor)
        value <- value %x% aux
        if (factor) 
            lD <- lD + attr(aux, "logDet")
    }
    if (factor) 
        attr(value, "logDet") <- lD
    
    dimnames(value) <- attr(object, "Dimnames")
    value
}

#' Unlist Names for pdKronecker
#'
#' Creates a character vector of names by combining names from a list using a Kronecker product-like structure.
#'
#' @param namesList A list of character vectors containing names.
#' @return A character vector of combined names.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
unlistFun <- function(namesList) {
    lenx <- sapply(namesList, length)
    namesList <- namesList[lenx>1]
    len1 <- length(namesList)
    namesList <- rev(namesList)
    
    dt <- expand.grid(namesList)
    allNames <- dt[[len1]]
    
    for (i in (len1-1):1) {
        aux <- dt[[i]]       
        allNames <- paste(allNames, dt[[i]], sep=":") 
    }
    
    meltx <- function(data, ...) data.frame(value = data)
    allNames <- meltx(allNames)
    allNames <- allNames[,"value"]
    allNames <- as.character(allNames)
    allNames
}


#' Subset pdKronecker Object
#'
#' Subsets a `pdKronecker` object, returning either a matrix or a new `pdKronecker` object
#' based on the indices provided.
#'
#' @param x A `pdKronecker` object.
#' @param i Row indices.
#' @param j Column indices.
#' @param drop Logical; if `TRUE`, drops dimensions when possible.
#' @param ... Additional arguments (currently ignored).
#' @return A matrix or a `pdKronecker` object.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
"[.pdKronecker" <- function(x, i, j, drop = TRUE, ...) {
    xx <- x
    x <- as.matrix(x)
    mCall <- match.call()
    mCall[[1]] <- get("[")
    mCall[["x"]] <- x
    mCall[["drop"]] <- drop
    if (length(i) == length(j) && mode(i) == mode(j) && all(i == j)) {
        mCall[["drop"]] <- FALSE
        val <- eval(mCall)
        vNames <- colnames(val)
        auxNames <- lapply(Names(xx, TRUE), function(el, vN) {
            aux <- match(vN, el)
            if (any(aux1 <- !is.na(aux))) {
                el[aux[aux1]]
            }
        }, vN = vNames)
        auxWhich <- !unlist(lapply(auxNames, is.null))
        if (sum(auxWhich) == 1) {
            return(pdConstruct(as.list(xx)[auxWhich][[1]], val))
        }
        auxNames <- auxNames[auxWhich]
        auxClass <- unlist(lapply(xx, function(el) class(el)[1]))[auxWhich]
        return(pdConstruct(xx, val, nam = auxNames, form = NULL, pdClass = auxClass))
    } else {
        eval(mCall)
    }
}

#' Extract Coefficients from pdKronecker
#'
#' Extracts coefficients from a `pdKronecker` object, optionally in unconstrained form.
#'
#' @param object A `pdKronecker` object.
#' @param unconstrained Logical; if `TRUE`, returns unconstrained coefficients.
#' @param ... Additional arguments (currently ignored).
#' @return A numeric vector of coefficients.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
coef.pdKronecker <- function(object, unconstrained = TRUE, ...) {
    coefSum <- coef(object[[1]], unconstrained)
    coefL <- lapply(object, coef, unconstrained)
    coefs <- unlist(lapply(coefL, FUN = function(x) x[-1]))
    coefAll <- c(coefSum, coefs)
    nms <- names(coefAll)
    dupi <- duplicated(nms)*(1:length(nms))
    idx <- dupi[dupi>0]
    nmsidx <- nms[idx]
    nms[idx] <- paste(nmsidx, idx, sep="")
    names(coefAll) <- nms
    coefAll
}

#' Assign Coefficients to pdKronecker
#'
#' Assigns new coefficients to a `pdKronecker` object.
#'
#' @param object A `pdKronecker` object.
#' @param value Numeric vector of coefficients.
#' @param ... Additional arguments (currently ignored).
#' @return The modified `pdKronecker` object.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
"coef<-.pdKronecker" <- function(object, ..., value) {
    if (is.null(plen <- attr(object, "plen"))) {
        stop(paste("Cannot change the parameter when", "length of parameters is undefined"))
    }
    plen <- plen -1
    plen[1] <-1
    
    if (length(value) != sum(plen)) {
        stop("Cannot change parameter length of initialized pdMat object")
    }
    ends <- cumsum(plen)
    starts <- 1 + c(0, ends[-length(ends)])
    
    attrx <- attributes(object)
    coef(object[[1]]) <- value[1] 
    for (i in 2:length(object)) {      
        coef(object[[i]]) <- c(0,value[(starts[i]):(ends[i])])
    }
    attributes(object) <- attrx
    object
}

#' Extract Correlation Matrix from pdKronecker
#'
#' Returns the correlation matrix of a `pdKronecker` object.
#'
#' @param object A `pdKronecker` object.
#' @param ... Additional arguments (currently ignored).
#' @return A correlation matrix with standard deviations as an attribute.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
corMatrix.pdKronecker <- function(object, ...) {
    if (!isInitialized(object)) {
        stop("Cannot access the matrix of uninitialized objects")
    }
    if (length(Names(object)) == 0) {
        stop("Cannot access the matrix of object without names")
    }
    namesList <- Names(object, TRUE)
    Ncol <- Dim(object)[2]
    value <- array(0, c(Ncol, Ncol), attr(object, "Dimnames"))
    stdDev <- double(Ncol)
    names(stdDev) <- colnames(value)
    for (i in seq_along(object)) {
        aux <- corMatrix(object[[i]])
        value[namesList[[i]], namesList[[i]]] <- as.vector(aux)
        stdDev[namesList[[i]]] <- attr(aux, "stdDev")
    }
    attr(value, "stdDev") <- stdDev
    value
}

#' Extract Formula from pdKronecker
#'
#' Returns the formula associated with a `pdKronecker` object, optionally as a list.
#'
#' @param x A `pdKronecker` object.
#' @param asList Logical; if `TRUE`, returns a list of formulas.
#' @param ... Additional arguments (currently ignored).
#' @return A formula or list of formulas.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
formula.pdKronecker <- function(x, asList = FALSE, ...) {
    val <- lapply(x, formula)
    isNULL <- unlist(lapply(val, is.null))
    if (all(isNULL)) 
        return(NULL)
    if (any(isNULL)) {
        stop("All elements must have formulas, when any has a formula.")
    }
    if (asList) 
        return(val)
    isTwoSided <- unlist(lapply(val, function(el) {
        inherits(el, "listForm")
    }))
    if (all(isTwoSided)) {
        val <- do.call("c", val)
        class(val) <- "listForm"
        return(val)
    }
    if (any(isTwoSided)) {
        stop(paste("All elements of formula must be list of two-sided formulae", 
            "or two-sided formulae"))
    }
    val <- lapply(rev(val), terms)
    aux <- paste(unlist(lapply(val, function(el) attr(el, "term.labels"))), 
        collapse = ":")
    if (!all(unlist(lapply(val, function(el) attr(el, "intercept"))))) {
        aux <- paste(aux, " - 1")
    }
    eval(parse(text = paste("~", aux)))
}

#' Check Initialization of pdKronecker
#'
#' Checks if a `pdKronecker` object is initialized.
#'
#' @param object A `pdKronecker` object.
#' @return Logical; `TRUE` if all components are initialized.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
isInitialized.pdKronecker <- function(object) {
    InitAll <- all(unlist(lapply(object, isInitialized)))
    InitAll
}

#' Compute Log-Determinant of pdKronecker
#'
#' Computes the log-determinant of a `pdKronecker` object.
#'
#' @param object A `pdKronecker` object.
#' @param ... Additional arguments (currently ignored).
#' @return Numeric; the log-determinant.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
logDet.pdKronecker <- function(object, ...) {
    sum(unlist(lapply(object, logDet)))
}

#' Assign Matrix to pdKronecker
#'
#' Assigns a matrix to a `pdKronecker` object, updating its component matrices.
#'
#' @param object A `pdKronecker` object.
#' @param value A matrix to assign.
#' @return The modified `pdKronecker` object.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
"matrix<-.pdKronecker" <- function(object, value) {
    value <- as.matrix(value)
    namesList <- Names(object, TRUE)
    Ncol <- Dim(object)[2]
    dims <- dim(value)
 
    if (!((dims[1] == dims[2]) && (dims[1] == Ncol))) {
        stop("Cannot change the number of columns on an initialized object")
    }   
    if (is.null(vNames <- rownames(value))) {
        vNames <- unlistFun(namesList)
        dimnames(value) <- list(vNames, vNames)
    } else {
        if (!(all(match(unlistFun(namesList), vNames, nomatch = 0)))) {
            stop("Names of object and value must match.")
        }
        attr(object, "Dimnames") <- list(vNames, vNames)
    }
    
    val1 <- as.matrix(value[1,1])
    nmsi <- namesList[[1]]   
    dimnames(val1) <- list(nmsi,nmsi)
    matrix(object[[1]]) <- val1
    
    len <- length(object)
    val <- value/val1[1,]
    
    for (i in len:2){
        nmsi <- namesList[[i]]
        leni <- length(nmsi)        
        nv <- nrow(val)
        idx <- seq(1,to =leni)
        val1 <- as.matrix(val[idx,idx])
        dimnames(val1) <- list(nmsi,nmsi)
        matrix(object[[i]]) <- val1
        idx <- seq(1, to =nv, by=leni)
        val <- as.matrix(val[idx,idx])/val1[1,1]
    }
     
    object
}

#' Extract Names from pdKronecker
#'
#' Returns the names of a `pdKronecker` object, optionally as a list.
#'
#' @param object A `pdKronecker` object.
#' @param asList Logical; if `TRUE`, returns a list of names.
#' @param ... Additional arguments (currently ignored).
#' @return A character vector or list of names.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
Names.pdKronecker <- function(object, asList = FALSE, ...) {
    if (asList) 
        attr(object, "namesList") 
    else 
        attr(object, "Dimnames")[[2]]
}

#' Assign Names to pdKronecker
#'
#' Assigns new names to a `pdKronecker` object.
#'
#' @param object A `pdKronecker` object.
#' @param value A character vector or list of names.
#' @param ... Additional arguments (currently ignored).
#' @return The modified `pdKronecker` object.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
"Names<-.pdKronecker" <- function(object, ..., value) {
    .functionLabel <- "Names<-.pdKronecker"                 # Function label (recommended)
    .traceR <- attr(options()$traceR, "fun")
    .traceR <-  if (is.null(.traceR)) function(...){} else .traceR      

    .traceR(210, lbl = "-> Names<-.pdKronecker starts")
    
    tmp <- Names(object)
  
    if (!is.null(Names(object))) { 
        clss <- class(object)
        class(object) <- "pdMat"
        Names(object) <- value
        class(object) <- clss
        .traceR(211, lbl = "Names<-.pdKronecker EXIT1")
        object
    } else {
        .traceR(212, lbl = "Names<-.pdKronecker EXIT2")
        object
    }
}

#' Extract Factor of pdKronecker
#'
#' Returns the factor of a `pdKronecker` matrix.
#'
#' @param object A `pdKronecker` object.
#' @return A matrix representing the factor.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
pdFactor.pdKronecker <- function(object) {
    pdMatrix(object, factor = TRUE)
}

#' Extract Matrix from pdKronecker
#'
#' Returns the matrix representation of a `pdKronecker` object, optionally as a factor.
#'
#' @param object A `pdKronecker` object.
#' @param factor Logical; if `TRUE`, returns the factor matrix.
#' @return A matrix.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
pdMatrix.pdKronecker <- function(object, factor = FALSE) {
    if (!isInitialized(object)) {
        stop("Cannot access the matrix of uninitialized objects")
    }
    if (length(Names(object)) == 0) {
        stop("Cannot access the matrix of object without names")
    }
    
    namesList <- Names(object, TRUE)
    Ncol <- Dim(object)[2]
    
    if (factor) {
        lD <- 0
    }
    
    value <- matrix(1)
    len <- length(object)
    for (i in 1:len) {
        aux <- pdMatrix(object[[i]], factor)
        value <- value %x% aux
        if (factor) 
            lD <- lD + attr(aux, "logDet")
    }
    if (factor) 
        attr(value, "logDet") <- lD
    
    dimnames(value) <- attr(object,"Dimnames")
    value
}

#' Solve pdKronecker Object
#'
#' Computes the inverse of a `pdKronecker` object by solving its component matrices.
#'
#' @param a A `pdKronecker` object.
#' @param b Optional; not used.
#' @param ... Additional arguments (currently ignored).
#' @return A `pdKronecker` object representing the inverse.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
solve.pdKronecker <- function(a, b, ...) {
    if (!isInitialized(a)) {
        stop("Cannot get the inverse of an uninitialized object")
    }
    
    pdx <- lapply(a, solve)
    attributes(pdx) <- attributes(a)
    pdx <- KroneckAux(pdx)
    
    pdx
}

#' Summarize pdKronecker Object
#'
#' @param object A `pdKronecker` object.
#' @param structName Character; name of the structure (default: "pdKronecker").
#' @param ... Additional arguments (currently ignored).
#' @return A summary object of class `summary.pdMatX`.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
summary.pdKronecker <- function(object, structName = "pdKronecker", ...) {
    value <- lapply(object, summary)
    names(value) <- unlist(lapply(object, function(el) paste(Names(el), 
        collapse = ", ")))
    attr(value, "structName") <- structName
    attr(value, "elementName") <- "KBlock"
    class(value) <- "summary.pdMatX"
    value
}

#' Variance-Covariance Matrix for pdKronecker
#'
#' Returns the variance-covariance matrix components of a `pdKronecker` object.
#'
#' @param x A `pdKronecker` object.
#' @param sigma Numeric; standard deviation scaling factor (default: 1).
#' @param ... Additional arguments; `rdig` can be specified as an integer for the number of digits for rounding (default: 3).
#' @return A matrix of variance-covariance components.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
VarCorr.pdKronecker <- function(x, sigma = 1, ...) {
    args <- list(...)
    rdig <- if (!is.null(args$rdig)) args$rdig else 3
    m <- lapply(x, VarCorr, sigma = sigma, rdig = rdig)
    bd <- do.call("rbind", m)
    attr(bd, "formStr") <- paste(sapply(m, attr, which = "formStr"), 
        collapse = ", ")
    bd
}

#' Print pdKronecker Object
#'
#' Prints a `pdKronecker` object, displaying its matrix and component covariance profiles.
#'
#' @param x A `pdKronecker` object.
#' @param opt Integer; if 1, prints component matrices (default: 1).
#' @param ... Additional arguments passed to print.
#' @return Invisibly returns the object.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
print.pdKronecker <- function(x, opt = 1, ...) {
    if (isInitialized(x)) {
        cat("Positive definite matrix structure of class", class(x)[1], 
            "representing\n")
        print(as.matrix(x), ...)
        if (opt==1) {
            cat("Matrix is a Kronecker product of the following covariance profiles: \n")
            lapply(x, FUN=function(el) print(as.matrix(el), ...))
        }
    } else {
        cat("Uninitialized positive definite matrix structure of class ", 
            class(x)[1], ".\n", sep = "")
    }
    invisible(x)
}
