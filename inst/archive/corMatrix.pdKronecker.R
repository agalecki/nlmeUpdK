#' Extract Correlation Matrix from pdKronecker (needs work)
#'
#' Returns the correlation matrix of a `pdKronecker` object, as an S3 method for the `corMatrix` generic from the `nlme` package.
#'
#' @param object A `pdKronecker` object.
#' @param ... Additional arguments (currently ignored).
#' @return A correlation matrix with standard deviations as an attribute.
#' @author Andrzej Galecki, Tomasz Burzykowski
#' @method corMatrix pdKronecker
#' @export
#' @examples
#' \dontrun{
#' # Define matrices D1 and D2
#' D1 <- matrix(c(3, 9, 9, 30), nrow = 2)  # Matrix for factor f1
#' D2 <- matrix(c(2, 4, 4, 10), nrow = 2)  # Matrix for factor f2
#' D1 %x% D2  # Kronecker product of D1 and D2
#'
#' # Create pdKronecker object
#' library(nlme)
#' pdId <- pdIdent(as.matrix(1), form = ~1)  # Identity matrix
#' pd1 <- pdLogChol(D1, form = ~f1-1)       # pdLogChol for D1
#' pd2 <- pdLogChol(D2, form = ~f2-1)       # pdLogChol for D2
#' pdL1 <- list(X = pdId, pD1 = pd1, pD2 = pd2)
#'
#' # Create data frame with factors
#' f1 <- gl(2, 1, labels = c("A", "B"))
#' f2 <- gl(2, 1, labels = c("a", "b"))
#' dt <- data.frame(f1, f2)
#'
#' # Create and inspect pdKronecker object
#' library(nlmeUpdK)
#' pdK <- pdKronecker(pdL1, data = dt)
#' corMatrix(pdK)  # Extract correlation matrix
#' }
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
