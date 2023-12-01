###  pdKronecker class
###
### Copyright 2011  Andrzej Galecki <agalecki@umich.edu>,
###                 Tomasz Burzykowski
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#


# pdKronecker
# KroneckAux
# pdConstruct.pdKronecker

pdKronecker <- function (value = numeric(0), 
   form = NULL, nam = NULL, data = sys.frame(sys.parent())) 
{
    funNm <- "pdKronecker"

   .traceRinit <- attr(options()$traceR, "init")
   .traceR <-   if (is.null(.traceRinit))
      function(...){} else .traceRinit(funNm) 

   .traceR(1, "pdKronecker STARTS here", funNm, msg = TRUE)
    
    object <- numeric(0)
    class(object) <- c("pdKronecker", "pdMat")
    .traceR(1, "pdKronecker ENDS here", funNm, msg = TRUE)
    pdConstruct(object, value, form, nam, data)
}

KroneckAux  <- function(pdx){
## pdx is pdk list or object 
## All components need to be initialized

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


pdConstruct.pdKronecker <-
function (object, 
    value = numeric(0), 
    form = formula(object, TRUE),
    nam = Names(object, TRUE), 
    data = sys.frame(sys.parent()), 
    pdClass = sapply(object,FUN=function(el) class(el)[1]) , ...) 
{

    funNm <- "pdConstruct.pdKronecker"

   .traceRinit <- attr(options()$traceR, "init")
   .traceR <-   if (is.null(.traceRinit))
      function(...){} else .traceRinit(funNm) 

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
            .traceR(9001, "pdConstruct.pdKronecker EXIT1", funNm, msg = TRUE)
            return(pdConstruct(object, as.matrix(value), form, 
                nam, data, pdClass))
        }
        else {
            .traceR(9001, "pdConstruct.pdKronecker EXIT2", funNm, msg = TRUE)
            return(pdConstruct(object, form = form, nam = nam, 
                data = data, pdClass = pdClass))
        }
    }
    .traceR(9002, "AFTER if value inherits pdMat", funNm, msg = TRUE)
   
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
    .traceR(9002,"oVal created", funNm, msg=TRUE)
    
    oVal <- value
    .traceR(9002, oVal, funNm)
    
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
   .traceR(9002,"object initiated" , funNm, msg= TRUE)
    
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
    }
    .traceR(9002,"End of nB loop" , funNm, msg = TRUE)
    
    names(object) <- nmsObject   # added Nov.20, 2008
    if (!all(unlist(lapply(object, inherits, "pdMat")))) {
        stop("all elements in the argument must inherit from pdMat objects")
    }
 ### Inserted Starts  
 
    InitX  <- all(unlist(lapply(object,isInitialized)))
 
    if (InitX) {
      object <- KroneckAux(object)  # :::
 
    }
### Inserted Ends
   .traceR(101, "Before namesList" , funNm, msg=TRUE)
    namesList <- lapply(object, Names)
    .traceR(101,"After namesList" , funNm, msg = TRUE)  

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
        .traceR(9001, , funNm, "pdConstruct.pdKronecker EXIT3")
       
        return(object)
    }
    else {
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
        .traceR(9001,"pdConstruct.pdKronecker EXIT4" , funNm, msg=TRUE)
        return(object)
    }
}




