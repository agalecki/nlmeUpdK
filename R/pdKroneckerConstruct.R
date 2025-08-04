#' Construct a Positive-Definite Kronecker Product Object 
#'
#' This function constructs a `pdKronecker` object, representing a positive-definite matrix
#' based on the Kronecker product of underlying positive-definite matrices, each represented
#' by individual `pdMat` objects. The first element is of class `pdIdent` and acts as a scale
#' multiplier.
#'
#' @param value A list with elements that are `pdMat` objects. The first element must be of
#'   class `pdIdent`. Defaults to `numeric(0)`.
#' @param form A list of formulas specifying the structure of the matrices. Defaults to `NULL`.
#' @param nam A list of names for the matrices. Defaults to `NULL`.
#' @param data A data frame used by the source function. Defaults to `sys.frame(sys.parent())`.
#' @return A `pdKronecker` object, also inheriting from class `pdMat`, representing a
#'   positive-definite Kronecker product matrix.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
pdKronecker <- function(value = numeric(0), form = NULL, nam = NULL,
                        data = sys.frame(sys.parent())) {
   .functionLabel <- "pdKronecker"                 # Function label (recommended)
   .traceR <- attr(options()$traceR, "fun")
   .traceR <-  if (is.null(.traceR)) function(...){} else .traceR      

   .traceR(1, lbl = "pdKronecker STARTS")
    
    object <- numeric(0)
    class(object) <- c("pdKronecker", "pdMat")
    .traceR(1, lbl = "pdKronecker ENDS", store = FALSE)
    pdConstruct(object, value, form, nam, data)
}

#' Construct pdKronecker Object
#'
#' Constructs a `pdKronecker` object from provided values, formulas, and names, ensuring
#' compatibility of dimensions and proper initialization of component `pdMat` objects.
#'
#' @param object A `pdKronecker` object to be constructed.
#' @param value A list, matrix, or numeric vector specifying the initial values.
#' @param form A list of formulas for the component matrices.
#' @param nam A list of names for the component matrices.
#' @param data A data frame for evaluating formulas.
#' @param pdClass A character vector specifying the classes of component `pdMat` objects.
#' @param ... Additional arguments passed to methods.
#' @return A constructed `pdKronecker` object.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @export
pdConstruct.pdKronecker <- function(object, value = numeric(0),
                                    form = formula(object, TRUE),
                                    nam = Names(object, TRUE),
                                    data = sys.frame(sys.parent()),
                                    pdClass = sapply(object, function(el) class(el)[1]), ...) {
   .functionLabel <- "PredConstruct.pdKronecker"                 # Function label (recommended)
   .traceR <- attr(options()$traceR, "fun")
   .traceR <-  if (is.null(.traceR)) function(...){} else .traceR      

   .traceR(101)
 
    nmsObject<- attr(object,"Xnames")    # Nov20, 2008

    if (inherits(value, "pdMat")) {
        if (inherits(value, "pdKronecker")) {
            if (length(form) == 0) 
                form <- formula(value, TRUE)
            if (length(nam) == 0) 
                nam <- Names(value, TRUE)
            if (missing(pdClass)) 
                pdClass <- unlist(lapply(value, data.class))
        }
        if (isInitialized(value)) {
            retx <- pdConstruct(object, as.matrix(value), form, 
                nam, data, pdClass)
            .traceR(9001, lbl = "EXIT9001")  
            return(retx)
        } else {
            retx <- pdConstruct(object, form = form, nam = nam, 
                data = data, pdClass = pdClass)
            .traceR(9002, lbl = "EXIT9002")  
            return(retx)
        }
    }
    .traceR(102, store = FALSE)
   
    if (!is.null(form)) {
        if (data.class(form) != "list") {
            stop("\"form\" must be a list")
        }
        nF <- length(form)
    }
    else {
        nF <- 0
    }
    if (!is.null(nam)) {
        if (data.class(nam) != "list") {
            stop("\"nam\" must be a list")
        }
        nN <- length(nam)
        if ((nF > 0) && (nN != nF)) {
            stop("\"form\" and \"nam\" have incompatible lengths")
        }
    }
    else {
        nN <- 0
    }
    if (!missing(pdClass)) {
        if (!is.character(pdClass)) {
            stop("\"pdClass\" must be a character vector")
        }
        nP <- length(pdClass)
        if ((nP > 1)) {
            if ((nF > 0) && (nF != nP)) {
                stop("\"form\" and \"pdClass\" have incompatible lengths")
            }
            if ((nN > 0) && (nN != nP)) {
                stop("\"nam\" and \"pdClass\" have incompatible lengths")
            }
        }
    }
    else {
        nP <- 1
    }

    nB <- max(c(nF, nN, nP))
    
    oVal <- value
    .traceR(103)
    
    if (length(value) == 0 || is.matrix(value) || is.numeric(value)) {
        if (nB == 1) {
            stop("None of the arguments specify more than one block")
        }
        value <- lapply(vector("list", nB), function(el) numeric(0))
    }
    else {
        if (data.class(value) != "list") {
            stop(paste("\"object\" must be a list, when not missing,", 
                "not a matrix, and not numeric"))
        }
        nO <- length(value)
        if ((nB > 1) && (nB != nO)) {
            stop("Arguments imply different number of blocks")
        }
        nB <- nO
    }
 
   if (nP == 1) {
        pdClass <- rep(pdClass, nB)
    }
   .traceR(104)
    
    object <- vector("list", nB)
    namInterc <- rep(FALSE, nB)
    namCoef <- vector("list", nB)
    for (i in 1:nB) {
        if (is.null(nm <- nam[[i]])) {
            if (is.null(frm <- form[[i]])) {
                if (inherits(value[[i]], "formula")) {
                  nm <- Names(getCovariateFormula(value[[i]]))
                  if ((length(nm) == 1) && (nm == "(Intercept)") && 
                    length(value[[i]]) == 3) {
                    nm <- sapply(splitFormula(getResponseFormula(value[[i]])[[2]], 
                      sep = "+"), function(el) deparse(el[[2]]))
                  }
                  if (length(value[[i]]) == 3) {
                    namCoef[[i]] <- sapply(splitFormula(getResponseFormula(value[[i]])[[2]], 
                      sep = "+"), function(el) deparse(el[[2]]))
                  }
                }
            }
            else {
                if (inherits(frm, "formula")) {
                  nm <- Names(getCovariateFormula(frm))
                  if ((length(nm) == 1) && (nm == "(Intercept)") && 
                    length(frm) == 3) {
                    nm <- sapply(splitFormula(getResponseFormula(frm)[[2]], 
                      sep = "+"), function(el) deparse(el[[2]]))
                  }
                  if (length(value[[i]]) == 3) {
                    namCoef[[i]] <- sapply(splitFormula(getResponseFormula(value[[i]])[[2]], 
                      sep = "+"), function(el) deparse(el[[2]]))
                  }
                }
                else {
                  nm <- unique(unlist(lapply(frm, function(el) {
                    Names(getCovariateFormula(el))
                  })))
                  if ((length(nm) == 1) && (nm == "(Intercept)") && 
                    length(frm[[1]]) == 3) {
                    nm <- sapply(frm, function(el) {
                      sapply(splitFormula(getResponseFormula(el)[[2]], 
                        sep = "+"), function(el1) deparse(el1[[2]]))
                    })
                  }
                  namCoef[[i]] <- sapply(frm, function(el) {
                    sapply(splitFormula(getResponseFormula(el)[[2]], 
                      sep = "+"), function(el1) deparse(el1[[2]]))
                  })
                }
            }
        }
        if (!is.null(nm)) {
            namInterc[i] <- (length(nm) == 1) && (nm == "(Intercept)")
        }
        object[[i]] <- pdMat(value[[i]], form[[i]], nam[[i]], 
            data, pdClass[i])
        .traceR(9010, lbl = paste("nB=", i))    
    }
    .traceR(106, lbl = "End of nB loop")
    
    names(object) <- nmsObject   # added Nov.20, 2008
    if (!all(unlist(lapply(object, inherits, "pdMat")))) {
        stop("all elements in the argument must inherit from pdMat objects")
    }
 
    InitX  <- all(unlist(lapply(object,isInitialized)))
 
    if (InitX) object <- KroneckAux(object) 
    namesList <- lapply(object, Names)  

    lNam <- unlist(lapply(namesList, length))
    if (!is.null(namCoef[[1]])) {
 
        namCoef <- unlist(namCoef)
        duplCoef <- unique(namCoef[duplicated(namCoef)])
        if (length(duplCoef) > 0) {
            for (i in 1:nB) {
                wchDupl <- !is.na(match(namesList[[i]], duplCoef))
                if (any(wchDupl)) {
                  namesList[[i]][wchDupl] <- paste(namesList[[i]][wchDupl], 
                    "(Intercept)", sep = ".")
                  Names(object[[i]]) <- namesList[[i]]
                }
            }
        }
    }
    if (sum(namInterc) > 1 && (length(unique(lNam[namInterc])) == 
        1)) {
        stop("Cannot have duplicated column names in a pdMat object")
    }
    if ((sum(namInterc) == length(lNam)) || !any(lNam[!namInterc])) {
        class(object) <- c("pdKronecker", "pdMat")
        if (is.null(formula(object))) {
            stop("Must have formula, when no names are given")
        }
        if (length(oVal) && (is.matrix(oVal) || is.numeric(oVal))) {
            stop("Must give names when initializing from matrix or parameter")
        }
  
       .traceR(900, lbl = "coef(object) EXIT3")  
        return(object)
    } else {
        if (!all(lNam)) {
            stop("All elements must have names, when any has names.")
        }
        attr(object, "namesList") <- namesList
        allNames <- unlistFun(namesList)  # :::
        
        if (any(duplicated(allNames))) {
            stop("Cannot have duplicated column names in a pdMat object")
        }

 
        # plen is Number of coefs per component
        plen <- unlist(lapply(object, function(el) {
            if (isInitialized(el)) {
                length(coef(el, TRUE))
            }
            else {
                matrix(el) <- diag(length(Names(el)))
                length(coef(el, TRUE))
            }
        }))
        if (!all(plen)) {
            stop("All elements must have a non-zero size")
        }
        attr(object, "plen") <- plen
  

        attr(object, "Dimnames") <- list(allNames, allNames)
        class(object) <- c("pdKronecker", "pdMat")
        if (length(oVal) > 0) {
            if (is.matrix(oVal)) {
                matrix(object) <- oVal
            }
            else if (is.numeric(oVal)) {
                coef(object) <- oVal
            }
        }
        
        
        names(object) <- nmsObject  # Added Nov.20, 2008  
        .traceR(9001,lbl = "EXIT9001")
        return(object)
    }
}

#' File: pdKroneckerConstruct.txt
#'
#' Rescale Coefficients for pdKronecker Components
#'
#' Rescales coefficients of component `pdMat` objects in a `pdKronecker` object, ensuring
#' the first component (of class `pdIdent`) holds the sum of coefficients.
#'
#' @param pdx A list or object of class `pdKronecker` with initialized components.
#' @return The rescaled `pdKronecker` object.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @keywords internal
KroneckAux <- function(pdx) {
    attrx <- attributes(pdx)

    coefsList <- lapply(pdx, FUN=function(el) coef(el,unconstrained=TRUE))
    coefSum <- sum(sapply(coefsList, FUN=function(el) el[1]))

    for (i in seq_along(pdx)){
        pdxi <- pdx[[i]]
        clss <- class(pdxi)[1]
        if (! (clss %in% c("pdDiag","pdIdent","pdCompSymm","pdLogChol"))) stop("Class:",
            clss, " not allowed") 
        coefs <- coefsList[[i]]
        coefs1 <- coefs
        if (clss=="pdDiag")     coefs1 <- coefs - coefs[1] 
        if (clss=="pdIdent")    coefs1 <- 0 
        if (clss=="pdCompSymm") coefs1[1] <- 0 
        if (clss=="pdLogChol"){
            Ncol <- round((sqrt(8 * length(coefs) + 1) - 1)/2)
            x1 <- exp(-coefs[1])

            coefs1[1:Ncol] <- coefs[1:Ncol] - coefs[1]

            coefs1[(Ncol+1):length(coefs)] <- x1* coefs[(Ncol+1):length(coefs)] 
        }

        if (is.na(coefs1)[1]) coefs1 <- numeric(0)
        coef(pdxi) <- coefs1
        pdx[[i]] <- pdxi
    }
    coef(pdx[[1]]) <- coefSum
    attributes(pdx) <- attrx
    pdx   # Rescaled coefficients
}
