# nlmeUpdK

Installation

```
devtools::install_github("agalecki/nlmeUpdK")
```


## Selected scripts 

Selected scripts from Chapters 17 (Panels 17.14 - 17.23) and 20 (Panels 20.1 - 20.5)

```
library(nlmeUpdK)
path <- system.file(package = "nlmeUpdK")
ch17 <- paste0(path, "/scriptsR4.5.1/Ch17a.R")
source(ch17)
ch20pdK <- paste0(path, "/scriptsR4.5.1/Ch20.2pdK1a.R")
source(ch20pdK)
```

## Illustration

Create pdKronecker object

```
# Define matrices D1 and D2
D1 <- matrix(c(3, 9, 9, 30), nrow = 2)  # Matrix for factor f1
D2 <- matrix(c(2, 4, 4, 10), nrow = 2)  # Matrix for factor f2
D1 %x% D2  # Kronecker product of D1 and D2
  
# Create pdKronecker object
library(nlme)
pdId <- pdIdent(as.matrix(1), form = ~1)  # Identity matrix
pd1 <- pdLogChol(D1, form = ~f1-1)        # pdLogChol for D1
pd2 <- pdLogChol(D2, form = ~f2-1)       # pdLogChol for D2
pdL1 <- list(X = pdId, pD1 = pd1, pD2 = pd2)
  
# Create data frame with factors
f1 <- gl(2, 1, labels = c("A", "B"))
f2 <- gl(2, 1, labels = c("a", "b"))
dt <- data.frame(f1, f2)
  
# Create pdKronecker object
library(nlmeUpdK)
pdK <- pdKronecker(pdL1, data = dt)

```

Inspect/modify pdKronecker object 

```
coef(pdK, unconstrained = FALSE)  # Constrained coefficients
round(coef(pdK), 3)               # Unconstrained coefficients (by default), rounded
formula(pdK)                      # formula
isInitialized(pdK)
logDet(pdK)
Names(pdK)
pdFactor(pdK)
pdMatrix(pdK)
solve(pdK)
summary(pdK)
VarCorr(pdK)
```