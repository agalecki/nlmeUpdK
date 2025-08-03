#' @noRd
.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("nlme", quietly = TRUE)) {
    stop("The 'nlme' package is required for 'nlmeUpdK' to function properly.")
  }
  
  # Create a package-specific environment to store the original function
  .nlmeUpdK_env <- new.env(parent = emptyenv())
  assign(".nlmeUpdK_env", .nlmeUpdK_env, envir = asNamespace("nlmeUpdK"))
  
  # Get the replacement function
  fun <- get("model.matrix.reStruct.U", envir = asNamespace("nlmeUpdK"))
  environment(fun) <- .GlobalEnv
  
  # Replace model.matrix.reStruct in nlme namespace
  ns <- asNamespace("nlme")
  if (exists("model.matrix.reStruct", envir = ns, inherits = FALSE)) {
    # Store the original function in nlmeUpdK's environment
    assign("model.matrix.reStruct_original", get("model.matrix.reStruct", envir = ns), envir = .nlmeUpdK_env)
    
    # Unlock and replace the binding
    if (bindingIsLocked("model.matrix.reStruct", ns)) {
      unlockBinding("model.matrix.reStruct", ns)
    }
    assign("model.matrix.reStruct", fun, envir = ns)
    lockBinding("model.matrix.reStruct", ns)
  } else {
    warning("Function 'model.matrix.reStruct' not found in 'nlme' namespace. Replacement skipped.")
  }
}

#' @noRd
.onUnload <- function(libpath) {
  ns <- asNamespace("nlme")
  .nlmeUpdK_env <- get(".nlmeUpdK_env", envir = asNamespace("nlmeUpdK"))
  
  if (exists("model.matrix.reStruct_original", envir = .nlmeUpdK_env, inherits = FALSE)) {
    if (bindingIsLocked("model.matrix.reStruct", ns)) {
      unlockBinding("model.matrix.reStruct", ns)
    }
    assign("model.matrix.reStruct", get("model.matrix.reStruct_original", envir = .nlmeUpdK_env), envir = ns)
    lockBinding("model.matrix.reStruct", ns)
    rm("model.matrix.reStruct_original", envir = .nlmeUpdK_env)
  }
}
