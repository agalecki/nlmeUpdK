#### Replace the function in nlme namespace (i.e. replace model.matrix.reStruct with model.matrix.reStruct.U)

fun <- nlmeUpdK:::model.matrix.reStruct.U
environment(fun) <- .GlobalEnv
environment(fun)

unlockBinding("model.matrix.reStruct", asNamespace("nlme"))
assign("model.matrix.reStruct", fun, envir = asNamespace("nlme"))
lockBinding("model.matrix.reStruct", asNamespace("nlme"))

# Verify the replacement
environment(get("model.matrix.reStruct", envir = asNamespace("nlme")))
# Should output: <environment: R_GlobalEnv>

library(nlme)
## require(R.utils)
R.utils:::reassignInPackage("model.matrix.reStruct", "nlme", fun)
#require(testthat)
#as.character(body(getS3method("model.matrix","reStruct")))[2]
#testhat(

###################################################
### code chunk: Chap17a_init (script from Chapter 17 pertaining to pdKronecker class)
###################################################
options(digits = 5, show.signif.stars = FALSE)

library(nlmeUpdK)  #  including pdKronecker class
library(lattice)
library(reshape)
sessionInfo()

data(prt, package = "nlmeU")
lme.spec.form3 <- spec.fo ~ prt.f + occ.f + sex.f + age.f + bmi + fiber.f +
   prt.f:occ.f + sex.f:age.f + prt.f:fiber.f + occ.f:fiber.f

fm17.3 <- 
   lme(lme.spec.form3,                       # (17.6)
       random = ~occ.f:fiber.f - 1|id,       # D:(17.7)
       data = prt) 

###################################################
### code chunk: R17.14
###################################################
pdId  <- pdIdent(~1)                   # c * I_1  in (17.14)
pd1UN <- pdLogChol(~fiber.f - 1)       # A for FiberType  
pd2UN <- pdLogChol(~occ.f - 1)         # B for PrePos
pdL1  <-                               # List of pdMat objects
   list(X = pdId, 
        FiberType = pd1UN,           
        PrePos = pd2UN)          
(pdKnms <- names(pdL1))                # Names saved for later use
fm17.4 <-                              # M17.4
   lme(lme.spec.form3, 
       random = list(id = pdKronecker(pdL1)),
       data = prt)


###################################################
### code chunk: R17.15
###################################################
DmtxKron1 <- getVarCov(fm17.4)                  # D:(17.14)
rownames(DmtxKron1)                             # Long row/col names
nms <- c("T1:Pre", "T1:Pos", "T2:Pre", "T2:Pos")# Short names...
dimnames(DmtxKron1) <- list(nms, nms)           # ... assigned.
DmtxKron1                                       # D displayed.


sgma    <- fm17.4$sigma                         # sigma...
(sgma2  <- sgma^2)                              # ... and sigma^2
reSt <- fm17.4$modelStruct$reStruct      # Random-effects structure
pdKron1 <- reSt[[1]]
names(pdKron1) <- pdKnms                 # Component names reinstated 
(c1 <- as.numeric(as.matrix(pdKron1$X))) # Mandatory multiplier
c1*sgma2                                       # Compare to d_11
(A1 <- as.matrix(pdKron1$FiberType))           # A for fiber type
(B1 <- as.matrix(pdKron1$PrePos))              # B for occasion
# sgma2 * c1 %x% A1 %x% B1               #  D = sigma2 @ A @ B: (17.14)

###################################################
### code chunk: R17.16
###################################################
pdId  <- pdIdent(~1)                     # c I_1 in ( 17.16)
pd1CS <- pdCompSymm(~fiber.f -1)         # A for  FiberType
pd2CS <- pdCompSymm(~occ.f-1)            # B for  PrePos 
pdL2  <- 
    list(X = pdId,           
         FiberType = pd1CS,              # A 
         PrePos = pd2CS)                 # B 
fm17.5 <-                                # M17.5
   lme(lme.spec.form3,               
       random = list(id = pdKronecker(pdL2)),
       data = prt)


###################################################
### code chunk: R17.17
###################################################
DmtxKron2 <- getVarCov(fm17.5)           # D
dimnames(DmtxKron2) <- list(nms,nms)     # Row/col names shortened
DmtxKron2                                # D printed


###################################################
### code chunk: R17.18
###################################################
anova(fm17.3, fm17.4, fm17.5)        # M17.3 > M17.4  > M17.5
anova(fm17.3, fm17.5)                # M17.3 > M17.5

###################################################
### code chunk: R17.19
###################################################
prt.Dep2 <-                              # Raw data transposed
    melt(prt, measure.var = c("iso.fo", "spec.fo"),
         variable_name = "DV")           # DV indicates dep. variable
lme.DV.form <-                           # Formula with interactions 
    formula(value ~ -1 + DV + DV:(-1 + (prt.f + occ.f + fiber.f)^2 
            + sex.f + age.f + sex.f:age.f + bmi)) 
pdUN <- pdLogChol(~DV-1)                 # E for DV  
pdL <- 
  list(X = pdId,
       DV = pdUN,
       FiberType = pd1UN,
       PrePos = pd2UN)
pdKnms   <- names(pdL)                   # Names saved for later use
lmeC     <- lmeControl(msMaxIter = 100)  # Maximum iterations
fm17.6 <- 
   lme(lme.DV.form,                          # Fixed part  
       random = list(id = pdKronecker(pdL)), # 8 random effects
       data = prt.Dep2,
       weights = varIdent(form = ~1|DV),     # DV-specific variance
       control = lmeC)


###################################################
### code chunk: R17.20
###################################################
DmtxDV2Kron <- getVarCov(fm17.6)
rownames(DmtxDV2Kron)                        # Row/col D names
nms <-                                         # ... shortened 
  c("is:T1:Pre", "is:T1:Pos", "is:T2:Pre", "is:T2:Pos", 
    "sp:T1:Pre", "sp:T1:Pos", "sp:T2:Pre", "sp:T2:Pos")
dimnames(DmtxDV2Kron) <- list(nms, nms)      # ... and assigned.
print(DmtxDV2Kron[1:4, 1:4], digits = 2)     # Block for iso.fo
print(DmtxDV2Kron[5:8, 5:8], digits = 5)     # Block for spec.fo
print(DmtxDV2Kron[1:4,5:8], digits = 3)      # Off-diagonal block
print(cov2cor(DmtxDV2Kron)[1:4, 1:4],        # Corr. block for iso.fo
      digits = 3) 

###################################################
### code chunk: R17.21
###################################################
(sgma   <- fm17.6$sigma)                  # sigma_1 (17.19)
(vStDV2 <- fm17.6$modelStruct$varStruct)  # delta
sgma*coef(vStDV2, unconstrained = FALSE)  # sigma_2 (spec.fo)
reStDV2 <- fm17.6$modelStruct$reStruct    # Random-effects structure
DV2pdxKron  <- reStDV2[[1]]
names(DV2pdxKron) <- pdKnms
(c3 <- as.matrix(DV2pdxKron$X))           # Mandatory multiplier
(E3 <- as.matrix(DV2pdxKron$DV))          # E for dep. variables
(A3 <- as.matrix(DV2pdxKron$FiberType))   # A for fiber type (17.21)
(B3 <- as.matrix(DV2pdxKron$PrePos))      # B for occasion


###################################################
### code chunk: R17.22
###################################################
cKron3 <- sgma^2 %x% c3 %x% E3 %x% A3 %x% B3 # Kronecker product ...
rownames(cKron3) <- nms
print(cKron3, digits = 2)                    # ... printed


###################################################
### code chunk: R17.23
###################################################
fixed.DV2.Kron1 <- summary(fm17.6)$tTable
fxdNms3 <- 
   c("iso:1", "spec:1", "iso:Low", "spec:Low", "iso:Pos",
     "spec:Pos", "iso:T2", "spec:T2", "iso:Male", "spec:Male",
     "iso:Old", "spec:Old", "iso:BMI", "spec:BMI",
     "iso:Low:Pos", "spec:Low:Pos", "iso:Low:T2", "spec:Low:T2",
     "iso:Pos:T2", "spec:Pos:T2", "iso:Male:Old", "spec:Male:Old")
rownames(fixed.DV2.Kron1) <- fxdNms3
printCoefmat(fixed.DV2.Kron1, digits = 3, cs.ind = c(1, 2), 
             dig.tst = 2, zap.ind = 5)
anova(fm17.6, type = "marginal")

#detach("package:R.utils")
#detach("package:R.oo")
#detach("package:R.methodsS3")
unloadNamespace("R.utils")
unloadNamespace("R.oo")
unloadNamespace("R.methodsS3")
sessionInfo()
proc.time()
