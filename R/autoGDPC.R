auto.gdpc <- function(Z, crit = "AIC", normalize = TRUE, auto_comp = TRUE, expl_var = 0.9, num_comp = 5, tol = 1e-04, 
                      k_max = 10, niter_max = 500, ncores = 1) {
  # Computes Generalized Dynamic Principal Components. The number of components can be supplied by the user 
  # or chosen automatically so that a given proportion of variance is explained. The number of lags is chosen
  # automatically using an AIC or a BIC type criterion.
  # All computations are done internally using leads rather than lags. However, the final result is outputted
  # using lags.
  #INPUT
  # Z: data matrix, series by columns 
  # crit: a string, either 'AIC' or 'BIC', indicating the criterion used to chose the number of lags. Default is 'AIC'
  # normalize: logical, if TRUE the data is standardized to zero mean and unit variance. Default is TRUE
  # auto_comp:  logical, if TRUE compute components until the proportion of explained variance is equal to expl_var, other
  # wise use num_comp components. Default is TRUE
  # expl_var: a number between 0 and 1. Desired proportion of explained variance (only if auto_comp==TRUE).
  # default is 0.9
  # num_comp: integer, number of components to be computed (only if auto_comp==FALSE). Default is 5
  # tol: desired accuracy when computing the components. Default is 1e-4
  # k_max: maximum number of lags. Default is 10
  # niter_max : maximum number of iterations. Default is 500
  # ncores: number of cores to be used for parallel computations. Default is 1
  
  #OUTPUT
  # A list of length equal to the number of computed components. The i-th entry of this list is an object of class
  # gdpc, that is, a list with entries:
  # f: Coordinates of the i-th Principal Component corresponding to the periods 1,…,T
  # initial_f: Coordinates of the i-th Principal Component corresponding to the periods -k+1,…,0.
  # Only for the case k>0
  # beta: beta matrix corresponding to f
  # alpha: alpha matrix corresponding to f
  # mse: mean (in T and m) squared error of the residuals of the fit with the first i components 
  # k_opt: number of lags used chosen using the criterion specified in crit
  # crit: the AIC or BIC of the fitted model, according to what was specified in crit
  # res: matrix of residuals of the fit with the first i components 
  # fitted: matrix of fitted values of the fit with the first i components 
  # expart: proportion of the variance explained by the first i components
  
  if (all(!inherits(Z, "matrix"), !inherits(Z, "mts"), !inherits(Z, "xts"), !inherits(Z, "data.frame"))) {
    stop("Z should belong to one of the following classes: matrix, data.frame, mts, xts")
  }
  if (!crit %in% c("BIC", "AIC")) {
    stop("crit should be AIC or BIC ")
  }
  if (!inherits(auto_comp, "logical")) {
    stop("auto_comp should be logical")
  }
  if (!inherits(normalize, "logical")) {
    stop("normalize should be logical")
  }
  if (!inherits(expl_var, "numeric")) {
    stop("expl_var should be numeric")
  } else if (!all(expl_var < 1, expl_var > 0)) {
    stop("expl_var be between 0 and 1")
  }
  if (!inherits(tol, "numeric")) {
    stop("tol should be numeric")
  } else if (!all(tol < 1, tol > 0)) {
    stop("tol be between 0 and 1")
  }
  if (!inherits(num_comp, "numeric")) {
    stop("num_comp should be numeric")
  } else if (any(!num_comp == floor(num_comp), num_comp <= 0)) {
    stop("num_comp should be a positive integer")
  }
  if (!inherits(k_max, "numeric")) {
    stop("k_max should be numeric")
  } else if (any(!k_max == floor(k_max), k_max <= 0)) {
    stop("k_max should be a positive integer")
  }
  if (!inherits(niter_max, "numeric")) {
    stop("niter_max should be numeric")
  } else if (any(!niter_max == floor(niter_max), niter_max <= 0)) {
    stop("niter_max should be a positive integer")
  }
  if (!inherits(ncores, "numeric")) {
    stop("ncores should be numeric")
  } else if (any(!ncores == floor(ncores), ncores <= 0)) {
    stop("ncores should be a positive integer")
  }
  # Pass to matrix form. Scale and transpose data.
  if (normalize) {
    V <- t(scale(as.matrix(Z)))
    mean_var_V <- 1
  } else {
    V <- t(as.matrix(Z))
    mean_var_V <- mean(apply(V, 1, var))  # Mean variance of the data
  }
  vard <- (1 - expl_var) * mean_var_V
  output <- vector("list")
  
  sel <- 1
  if (crit == "BIC") {
    sel <- 2
  }
  
  ### Set-up cluster for parallel computations
  cores <- min(detectCores(), ncores)
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  ### 
  
  comp_ready <- 1
  print(paste("Computing component number", comp_ready))
  out <- my_autodyc(V, k_max, mean_var_V, tol, niter_max, sel)
  mse <- out$mse  #Mean squared error (in N and m)
  output[[comp_ready]] <- out
  V <- out$res
  
  if (auto_comp) {
    while (mse > vard) {
      comp_ready <- comp_ready + 1
      print(paste("Computing component number", comp_ready))
      out <- my_autodyc(V, k_max, mean_var_V, tol, niter_max, sel)
      mse <- out$mse
      output[[comp_ready]] <- out
      V <- out$res
    }
  } else {
    while (comp_ready < num_comp) {
      comp_ready <- comp_ready + 1
      print(paste("Computing component number", comp_ready))
      out <- my_autodyc(V, k_max, mean_var_V, tol, niter_max, sel)
      output[[comp_ready]] <- out
      V <- out$res
    }
  }
  
  on.exit(stopCluster(cl))
  
  if (normalize) {
    Z <- scale(Z)
  }
  
  output <- lapply(output, construct_gdpc, Z)
  
  return(output)
  
}

my_autodyc <- function(V, k_max, mean_var_V, tol = 1e-04, niter_max = 500, sel = 1) {
  #Auxiliary function to choose the optimal number of leads
  #INPUT
  # V : matrix of original data or residuals where each ROW is a different time series
  # k_max : maximum of numbers of leads to be considered
  # mean_var_V : mean variance of original data
  # tol : relative precision
  # niter_max: maximum number of iterations
  # sel: criterion to be used, AIC = 1, BIC = 2
  #OUTPUT
  # A list with entries:
  # k_opt: optimal number of leads
  # f: dynamic component with k_opt leads
  # beta: matrix beta corresponding to f
  # alfa: matrix beta alfa corresponding to f
  # mse: mean squared error
  # crit: criterion
  # res: matrix of residuals
  # expart: proportion of the variance explained
  
  
  
  m <- nrow(V)
  N <- ncol(V)
  
  exports <- c("gdpc.priv")
  k_lag <- NULL
  crits <- rep(0, k_max + 1)
  fits <- vector("list", k_max + 1)
  fits <- foreach(k_lag = 1:(k_max + 1), .export = exports, .packages = "gdpc") %dopar% {
    gdpc.priv(V, k = k_lag - 1, tol = tol, niter_max = niter_max, sel = sel)
  }
  
  crits <- sapply(fits, function(x) {
    x$crit
  })  #Get criterion corresponding to each lead
  k_opt <- which.min(crits) - 1
  out <- fits[[k_opt + 1]]
  expart <- 1 - out$mse/mean_var_V
  out$k_opt <- k_opt
  out$expart <- expart
  return(out)
}


gdpc.priv <- function(V, k, f_ini = NULL, tol = 1e-04, niter_max = 500, sel = 1) {
  # This function computes a single GDPC with a given number of leads.
  #INPUT
  # V: data matrix each ROW is a different time series
  # k: number of leads used
  # tol: relative precision, stopping criterion
  # niter_max: maximum number of iterations
  # f_ini: (optional) initial principal component. If no argument is passed, the standard
  # first principal component with k leads is used
  # sel: AIC (1) or BIC (2)
  #OUTPUT
  # f: final principal component
  # beta: matrix beta corresponding to f_fin
  # alpha: matrix alpha corresponding to f_fin
  # mse:  mean squared error (in N and m)
  # crit: criterion used to evaluate the fit, that is, sel
  # res: matrix of residuals
  
  m <- nrow(V)  #Number of time series
  N <- ncol(V)  #Length of the time series
  niter <- 0  #Number of iterations
  
  if (is.null(f_ini)) {
    f_ini <- t(V) %*% svd(scale(t(V), scale = FALSE), nu = 0, nv = 1)$v[, 1]
    f_ini <- c(f_ini, rep(f_ini[N], k))
  }
  
  if (length(f_ini) != N + k) {
    warning("Length of initial factor does not equal N+k")
  }
  f_ini <- scale(f_ini)
  out <- betaf(V, f_ini, k, sel)
  mse_ini <- out$mse/m
  criter <- 10  #Stopping criterion for the iterations
  while (criter > tol & niter < niter_max) {
    niter <- niter + 1
    f_fin <- scale(matrix_ff(V, out$beta, out$alpha, k))
    out <- betaf(V, f_fin, k, sel)
    mse_fin <- out$mse/m
    criter <- 1 - mse_fin/mse_ini
    mse_ini <- mse_fin
  }
  if (niter >= niter_max) {
    warning("Iterations did not converge. Consider increasing niter_max.")
  }
  out$mse <- out$mse/m
  out$f <- f_fin
  return(out)
  
}


gdpc <- function(Z, k, f_ini = NULL, tol = 1e-04, niter_max = 500, crit = "AIC") {
  # A wrapper function for gdpc.priv.
  #INPUT
  # Z: data matrix each COLUMN is a different time series
  # k: number of lags used to predict
  # tol: relative precision, stopping criterion
  # niter_max: maximum number of iterations
  # f_ini: (optional) initial principal component. If no argument is passed, the standard
  # first principal component with k lags is used
  # crit: AIC or BIC
  #OUTPUT
  # An object of class gdpc, that is, a list with entries:
  # f: the dynamic component 
  # beta: beta matrix corresponding to f
  # alpha: alpha matrix corresponding to f
  # mse: mean (in N and m) squared error of the residuals of the fit with the first i components 
  # k_opt: number of lags, chosen using the criterion specified in crit, used to predict with the i-th component
  # crit: the AIC or BIC of the fitted model, according to what was specified in crit
  # res: matrix of residuals of the fit
  # fitted: matrix of fitted values of the fit
  # expart: proportion of the variance explained
  
  
  if (all(!inherits(Z, "matrix"), !inherits(Z, "mts"), !inherits(Z, "xts"), !inherits(Z, "data.frame"))) {
    stop("Z should belong to one of the following classes: matrix, data.frame, mts, xts")
  }
  if (all(!is.null(f_ini), !inherits(f_ini, "numeric"), !inherits(f_ini, "matrix"), !inherits(f_ini, "ts"), !inherits(f_ini, "xts"))) {
    stop("f_ini should belong to one of the following classes: numeric, matrix, ts, xts")
  }
  if (!crit %in% c("BIC", "AIC")) {
    stop("crit should be AIC or BIC ")
  }
  if (!inherits(tol, "numeric")) {
    stop("tol should be numeric")
  } else if (!all(tol < 1, tol > 0)) {
    stop("tol be between 0 and 1")
  }
  if (!inherits(k, "numeric")) {
    stop("k_max should be numeric")
  } else if (any(!k == floor(k), k <= 0)) {
    stop("k_max should be a positive integer")
  }
  if (!inherits(niter_max, "numeric")) {
    stop("niter_max should be numeric")
  } else if (any(!niter_max == floor(niter_max), niter_max <= 0)) {
    stop("niter_max should be a positive integer")
  }
  
  
  sel <- 1
  if (crit == "BIC") {
    sel <- 2
  }
  out <- gdpc.priv(t(Z), k, f_ini, tol, niter_max, sel)
  out$k_opt <- k
  out$expart <- 1 - out$mse/mean(apply(Z, 2, var))
  out <- construct_gdpc(out, Z)
  return(out)
}

construct_gdpc <- function(out, data) {
  #This function constructs an object of class gdpc.
  #INPUT
  # out: the output of gdpc.priv
  # data: the data matrix passed to gdpc.priv
  #OUTPUT
  # An object of class gdpc, that is, a list with entries:
  # f: the dynamic component 
  # beta: beta matrix corresponding to f
  # alpha: alpha matrix corresponding to f
  # mse: mean (in N and m) squared error of the residuals of the fit with the first i components 
  # k_opt: number of lags used
  # crit: the AIC or BIC of the fitted model, according to what was specified in crit
  # res: matrix of residuals of the fit
  # fitted: matrix of fitted values of the fit
  # expart: proportion of the variance explained
  
  k <- ncol(out$beta) - 1  #number of leads
  out$beta <- out$beta[, (k + 1):1]
  if (k != 0) {
    out$initial_f <- out$f[1:k]
  } else {
    out$initial_f <- 0
  }
  out$f <- out$f[(k + 1):length(out$f)]
  out$res <- t(out$res)
  out$fitted <- data - out$res
  colnames(out$res) <- colnames(data)
  colnames(out$fitted) <- colnames(data)
  if (inherits(data, "xts")) {
    out$f <- xts(out$f, order.by = index(data), frequency = frequency(data))
    out$res <- xts(out$res, order.by = index(data), frequency = frequency(data))
    out$fitted <- xts(out$fitted, order.by = index(data), frequency = frequency(data))
  } else if (inherits(data, "ts")) {
    out$f <- ts(out$f, start = start(data), end = end(data), frequency = frequency(data))
    out$res <- ts(out$res, start = start(data), end = end(data), frequency = frequency(data))
    out$fitted <- ts(out$fitted, start = start(data), end = end(data), frequency = frequency(data))
  }
  class(out) <- append(class(out), "gdpc")
  return(out)
}

fitted.gdpc <- function(object, ...) {
  # Returns the fitted values of a gdpc object
  return(object$fitted)
}

residuals.gdpc <- function(object, ...) {
  # Returns the residuals of a gdpc object
  return(object$res)
}

plot.gdpc <- function(x, which = "Loadings", which_load = 0, ...) {
  #Plots a gdpc object
  #INPUT
  # x: An object of class gdpc, the result of gdpc or one of the entries 
  # of the result of auto.gdpc
  # which: String. Indicates what to plot, either 'Component' or 'Loadings'
  # Default is 'Loadings'
  # which_load: Lag number indicating which loadings should be plotted. 
  # Only used if which = 'Loadings'. Default is 0.
  switch(which, Component = plot(x$f, type = "l", main = "Principal Component", ...), Loadings = plot(x$beta[, 
                                                                                                             which_load + 1], type = "l", main = c(paste(which_load, "lag loadings")), ...))
}