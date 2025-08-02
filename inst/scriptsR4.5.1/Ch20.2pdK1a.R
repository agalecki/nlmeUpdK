
date()
sessionInfo()

# R20.1

D1 <- c(3, 9, 
        9, 30)
dim(D1) <- c(2, 2)       # D1 for factor  f1

D2 <- c(2, 4, 
        4, 10)     
dim(D2) <- c(2, 2)       # D2 for factor f2
D1 %x% D2                # D1 @ D2

# R20.2

library(nlme)
(pdId <- pdIdent(as.matrix(1), form = ~1)) # Mandatory
(pd1 <- pdLogChol(D1, form = ~f1-1))       # D1
(pd2 <- pdLogChol(D2, form = ~f2-1))       # D2
pdL1 <-                                    # The main argument
   list(X = pdId,
        pD1 = pd1,
        pD2 = pd2)

# R203a
f1 <- gl(2, 1, labels = c("A","B")) 
f2 <- gl(2, 1, labels = c("a","b"))
(dt <- data.frame(f1, f2))



# R203b
library(nlmeUpdKronecker)
methods(class= pdKronecker)
(pdK <- pdKronecker(pdL1, data = dt))      # D1 @ D2


# R20.4
(c0x <- as.matrix(pdK[[1]]))
(D1x <- as.matrix(pdK[[2]]))               #  Proportional to D1
(D2x <- as.matrix(pdK[[3]]))               #  Proportional to D2
Dx <- c0x %x% D1x %x% D2x                  #  D = D1 @ D2



# R20.5a
formula(pdK, asList = TRUE)               # List of component formulae
formula(pdK)                              # One-sided formula (default)

# R20.5b
(pdKform <- formula(~(f2-1):(f1-1)-1))
pdKterms <- terms(pdKform)                # Terms object  
labels(pdKterms)                          # Expanded formula
attr(pdKterms, "intercept")               # Intercept omitted

library(nlmeForTesting)
(nms <- Names(pdK))
dimnames(Dx) <- list(nms, nms)
Dx
methods(class = pdKronecker)   
packageVersion("nlmeUpdKronecker")
sessionInfo()
detach(package:nlmeUpdKronecker)

