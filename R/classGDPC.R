is.gdpc <- function(object, ...) {
  #This function checks whether an object is of class gdpc
  #First check if object is a list.
  #Second check if the list has the correct entries
  #Third check if the entries have the correct classes
  #Fourth check if the entries have the correct dims
  if (any(!inherits(object, "list"), !inherits(object, "gdpc"))) {
    return(FALSE)
  } else if (any(is.null(object$f), is.null(object$initial_f), is.null(object$beta),
                 is.null(object$alpha), is.null(object$mse), is.null(object$crit),
                 is.null(object$k), is.null(object$expart), is.null(object$call), 
                 is.null(object$conv))) {
    return(FALSE)
  } else if (any(!inherits(object$mse, "numeric"), !inherits(object$crit, "numeric"), !inherits(object$alpha, "numeric"),
                 !inherits(object$beta, "matrix"), !inherits(object$call, "call"), !inherits(object$conv, "logical"), 
                 all(!inherits(object$f,"numeric"), !inherits(object$f, "ts"), !inherits(object$f, "xts")),
                 all(!inherits(object$k, "numeric"), !inherits(object$k, "integer")), !inherits(object$expart, "numeric"),
                 all(!inherits(object$initial_f,"numeric"), !inherits(object$initial_f,"ts"), !inherits(object$initial_f,"xts"))
  )) {
    return(FALSE)
  } else if (any(length(object$alpha) != dim(object$beta)[1], dim(object$beta)[2] != object$k + 1)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

construct.gdpc <- function(out, data) {
  #This function constructs an object of class gdpc.
  #INPUT
  # out: the output of gdpc_priv
  # data: the data matrix passed to gdpc_priv
  #OUTPUT
  # An object of class gdpc, that is, a list with entries:
  # f: coordinates of the Principal Component corresponding to the periods 1,…,T
  # initial_f: Coordinates of the Principal Component corresponding to the periods -k+1,…,0.
  # beta: beta matrix corresponding to f
  # alpha: alpha vector corresponding to f
  # mse: mean (in N and m) squared error of the residuals of the fit
  # k: number of lags used
  # crit: the criterion of the fitted model, according to what was specified in crit
  # expart: proportion of the variance explained
  # call: the matched call
  # conv: logical. Did the iterations converge?
  
  k <- ncol(out$beta) - 2  #number of leads
  out$alpha <- as.numeric(out$beta[, k + 2])
  out$beta <- out$beta[, (k + 1):1]
  if (k != 0) {
    out$initial_f <- out$f[1:k]
  } else {
    out$initial_f <- 0
  }
  out$beta <- as.matrix(out$beta)
  rownames(out$beta) <- colnames(data)
  out$f <- out$f[(k + 1):length(out$f)]
  out$res <- NULL
  if (inherits(data, "xts")) {
    out$f <- reclass(out$f, match.to = data)
  } else if (inherits(data, "ts")) {
    out$f <- ts(out$f, start = start(data), end = end(data), frequency = frequency(data))
  }
  class(out) <- append(class(out), "gdpc")
  return(out)
}

fitted.gdpc <- function(object, ...) {
  # Returns the fitted values of a gdpc object
  if (!is.gdpc(object)){
    stop("object should be of class gdpc")
  }
  fitted <- getFitted(object$f, object$initial_f, object$beta, object$alpha, object$k)
  if (inherits(object$f, "xts")) {
    fitted <- reclass(fitted, match.to = object$f)
  } else if (inherits(object$f, "ts")) {
    fitted <- ts(fitted)
    attr(fitted, "tsp") <- attr(object$f, "tsp")
  }
  return(fitted)
}

plot.gdpc <- function(x, which = "Component", which_load = 0, ...) {
  #Plots a gdpc object
  #INPUT
  # x: An object of class gdpc, the result of gdpc or one of the entries 
  # of the result of auto.gdpc
  # which: String. Indicates what to plot, either 'Component' or 'Loadings'
  # Default is 'Component'
  # which_load: Lag number indicating which loadings should be plotted. 
  # Only used if which = 'Loadings'. Default is 0.
  if (!is.gdpc(x)) {
    stop("x should be of class gdpc")
  }
  if (!which %in% c("Component", "Loadings")) {
    stop("which should be either Component or Loadings ")
  }
  if (!inherits(which_load, "numeric")) {
    stop("which_load should be numeric")
  } else if (any(!which_load == floor(which_load), which_load < 0, which_load > ncol(x$beta) - 1)) {
    stop("which_load should be a non-negative integer, at most equal to the number of lags")
  }
  if (which == "Component"){
    plot(x$f, type = "l", main = "Principal Component", ...) 
  } else if (which == "Loadings"){
    plot(x$beta[, which_load + 1], type = "l", main = c(paste(which_load, "lag loadings")), ...)
  }
}

print.gdpc <- function(x, ...) {
  #Print method for the gdpc class
  if (!is.gdpc(x)) {
    stop("x should be of class gdpc")
  }
  y <- list(x)
  class(y) <- append(class(y), 'gdpcs')
  print(y)
}

construct.gdpcs <- function(out, data, fn_call) {
  #This function constructs an object of class gdpcs that is, a list of length equal to 
  #the number of computed components. The i-th entry of this list is an object of class gdpc.
  #INPUT
  # out: the output of auto.gdpc
  # data: the data matrix passed to auto.gdpc
  # fn_call: the original call to auto.gdpc
  #OUTPUT
  # An object of class gdpcs, that is, a list where each entry is an object of class gdpc.
  out <- lapply(out, function(x, fn_call){ x$call <- fn_call; return(x)}, fn_call)
  out <- lapply(out, construct.gdpc, data)
  class(out) <- append(class(out), "gdpcs")
  return(out)
}



is.gdpcs <- function(object, ...) {
  #This function checks whether an object is of class gdpcs,
  #that is, if each of its entries is a list of class gdpc
  if (any(!inherits(object, "gdpcs"), !inherits(object, "list"))) {
    return(FALSE)
  } else {
    return(all(sapply(object, is.gdpc)))
  }
}


components <- function(object, ...){
  # Generic function for getting components out of an object
  UseMethod("components", object)
}

components.gdpcs <- function(object, which_comp = 1) {
  # This function extracts the desired components from a gdpcs object
  #INPUT
  # object: the output of auto.gdpc
  # which_comp: Integer vector. Indicates which components to get
  #OUTPUT
  # A matrix with the desired components as columns
  if (!is.gdpcs(object)) {
    stop("object should be of class gdpcs")
  }
  if (all(!inherits(which_comp, "numeric"), !inherits(which_comp, "integer"))) {
    stop("which_comp should be numeric")
  } else if (any(!which_comp == floor(which_comp), which_comp <= 0, which_comp > length(object))) {
    stop("The entries of which_comp should be positive integers, at most equal to the number of components")
  }
  object <- object[which_comp]
  comps <- sapply(object, function(object){ object$f })
  colnames(comps) <- paste("Component number", which_comp)
  if (inherits(object[[1]]$f, "xts")) {
    comps <- reclass(comps, match.to = object[[1]]$f)
  } else if (inherits(object[[1]]$f, "ts")) {
    comps <- ts(comps, start = start(object[[1]]$f), end = end(object[[1]]$f), frequency = frequency(object[[1]]$f))
  } 
  return(comps)
}

plot.gdpcs <- function(x, which_comp = 1, ...) {
  #Plots a gdpcs object
  #INPUT
  # x: An object of class gdpcs, the result of auto.gdpc
  # which_comp: Integer vector. Indicates which components to plot
  if (!is.gdpcs(x)) {
    stop("x should be of class gdpcs")
  }
  if (all(!inherits(which_comp, "numeric"), !inherits(which_comp, "integer"))) {
    stop("which_comp should be numeric")
  } else if (any(!which_comp == floor(which_comp), which_comp <= 0, which_comp > length(x))) {
    stop("The entries of which_comp should be positive integers, at most equal to the number of components")
  }
  comps <- components(x, which_comp)
  if (inherits(comps, "xts") & length(which_comp)==1) {
    plot(comps, main = "Principal Components",...)
  } else if (inherits(comps, "ts")) {
    plot(comps, main = "Principal Components", plot.type = 'multiple', ...)
  } else {
    comps <- ts(comps)
    plot(comps, main = "Principal Components", plot.type = 'multiple', ...)
  }
  
}

print.gdpcs <- function(x, ...) {
#   Print method for the gdpcs class
  if (!is.gdpcs(x)) {
    stop("x should be of class gdpcs")
  }
  lags <- sapply(x, function(x){ round(x$k, 3) })
  vars <- sapply(x, function(x){ round(x$expart, 3) })
  mses <- sapply(x, function(x){ round(x$mse, 3) })
  crits <- sapply(x, function(x){ round(x$crit, 3) })
  mat <- cbind(lags, crits, mses, vars)
  nums <- paste(1:length(x))
  nums <- sapply(nums, function(x){ paste("Component", x)})
  fn_call <- x[[1]]$call
  crit_name <- fn_call$crit
  colnames(mat) <- c("Number of lags", crit_name, "MSE", "Explained Variance")
  tab <- data.frame(mat, row.names = nums)
  print(tab)
}