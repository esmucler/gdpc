auto.gdpc <- function(Z, crit = 'AIC', normalize = TRUE, auto_comp = TRUE, expl_var = 0.9, num_comp = 5, tol = 1e-4, k_max = 10, niter_max = 500, ncores = 1) {
  # This is the main function to compute the Generalized Dynamic Principal Components, Pena and Yohai (2016) JASA
  # All computations are done internally using leads rather than lags. However, the final result is outputted
  # using lags using the auxiliary leads2lags function.
  #Input
  #Z: data matrix, series by columns 
  #crit: a string, either 'AIC' or 'BIC', indicating the criterion used to chose the number of lags. Default is 'AIC'
  #normalize: logical, if TRUE the data is standardized to zero mean and unit variance
  #auto_comp:  logical, if TRUE compute components until the proportion of explained variance is equal to expl_var, other
  #wise use num_comp components
  #expl_var: a number between 0 and 1. Desired proportion of explained variance (only if auto_comp==TRUE)
  #num_comp: integer, number of components to be computed (only if auto_comp==FALSE)
  #tol: desired accuracy when computing the components
  #k_max: maximum number of lags
  #niter_max : maximum number of iterations
  #ncores: number of cores to be used for parallel computations
    
  #Output
  #A list of length equal to the number of computed components. The i-th entry of this list is a list with entries
  # $f: the i-th dynamic component 
  # $beta: beta matrix corresponding to f_fin
  # $alpha: alpha matrix corresponding to f_fin
  # $mse: mean (in N and m) squared error of the residuals of the fit with the first i components 
  # $k_opt: number of lags, chosen using the criterion specified in crit, used to predict with the i-th component
  # $crit: the AIC or BIC of the fitted model, according to what was specified in crit
  # $res: matrix of residuals of the fit with the first i components 
  # $expart: proportion of the variance explained by the first i components
  
  
  if (normalize) {
    V <- as.matrix(t(scale(Z)))
    mean_var_V <- 1
  } else {
    V <- as.matrix(t(Z))
    mean_var_V <- mean(apply(V, 1, var)) # Mean variance of the data
  }
  vard <- (1 - expl_var) * mean_var_V
  output <- vector("list")
  
  sel <- 1
  if (crit == 'BIC'){
    sel <- 2
  }
  
  ###Set-up cluster for parallel computations
  cores<-min(detectCores(),ncores)
  cl<-makeCluster(cores)
  registerDoParallel(cl)
  ###
  
  comp_ready <- 1
  print(paste('Computing component number',comp_ready))
  out <- my_autodyc(V, k_max, mean_var_V, tol, niter_max, sel)
  mse <- out$mse #Mean squared error (in N and m)
  output[[comp_ready]] <- out
  V <- out$res
  
  if(auto_comp){
     while (mse > vard){
      comp_ready <- comp_ready + 1
      print(paste('Computing component number',comp_ready))
      out <- my_autodyc(V, k_max, mean_var_V, tol, niter_max, sel)
      mse <- out$mse
      output[[comp_ready]] <- out
      V <- out$res
      }
  } else{
    while (comp_ready < num_comp){
      comp_ready <- comp_ready + 1
      print(paste('Computing component number',comp_ready))
      out <- my_autodyc(V, k_max, mean_var_V, tol, niter_max, sel)
      output[[comp_ready]] <- out
      V <- out$res
    }
  }
  
  on.exit(stopCluster(cl))
  
  output <- lapply(output, leads2lags)
  
  return(output)
  
}

my_autodyc <- function(V, k_max, mean_var_V, tol = 1e-04, niter_max = 500, sel = 1) {
  #Auxiliary function to choose the optimal number of leads
  
  
  #INPUT
  #V : matrix of original data or residuals where each row is a different time series
  #k_max : maximum of numbers of leads to be considered
  #mean_var_V : mean variance of original data
  #tol : relative accuracy to stop the search of the DPC ( 0<ep<1)
  #niter_max: maximum number of iterations
  #sel: criterion to be used, AIC = 1, BIC = 2
  
  #OUTPUT
  #A list with entries
  #k_opt: optimal number of leads
  #f: dynamic component with k_opt leads
  #beta: matrix beta corresponding to f
  #alfa: matrix beta alfa corresponding to f
  #mse: mean squared error
  #crit: criterion
  #res: matrix of residuals
  
  
  
  m <- nrow(V)
  N <- ncol(V)
  
  exports <- c("gdpc.priv")
  k_lag <- NULL
  crits <- rep(0, k_max + 1)
  fits <- vector('list', k_max + 1)
  fits <- foreach(k_lag = 1:(k_max + 1), .export = exports, .packages='gdpc') %dopar% {
    gdpc.priv(V, k = k_lag - 1, tol = tol, niter_max = niter_max, sel = sel)
  }
  
  crits <- sapply(fits, function(x){ x$crit })
  k_opt <- which.min(crits) - 1
  out <- fits[[k_opt + 1]]
  expart <- 1 - out$mse / mean_var_V
  out$k_opt <- k_opt
  out$expart <- expart
  return(out)
} 


gdpc.priv <- function(V, k, f_ini = NULL, tol = 1e-4, niter_max = 500, sel = 1) {
  #INPUT
  #V: data matrix each row is a different time series
  #k: number of leads used to predict
  #tol: relative precision, stopping criterion
  #niter_max: maximum number of iterations
  #f_fin: (optional) initial principal component. If no argument is passed, the standard
  #first principal component with k leads is used
  #sel: AIC (1) or BIC (2)
  #OUTPUT
  #f: final principal component
  #beta: matrix beta corresponding to f_fin
  #alpha: matrix alpha corresponding to f_fin
  #mse:  mean squared error (in N and m)
  #crit: criterion
  #res: matrix of residuals
  
  m <- nrow(V)  #Number of time series
  N <- ncol(V)  #Length of the time series
  niter <- 0  #Number of iterations
  
  if (is.null(f_ini)) {
    f_ini <- t(V) %*% svd(scale(t(V),scale=FALSE), nu = 0, nv = 1)$v[,1]
    f_ini <- c(f_ini, rep(f_ini[N], k))
  }
  
  if(length(f_ini)!=N+k){
    warning('Length of initial factor does not equal N+k')
  }
  f_ini <- scale(f_ini)
  out <- betaf(V, f_ini, k, sel)
  mse_ini <- out$mse / m
  criter <- 10  #Stopping criterion for the iterations
  while (criter > tol & niter < niter_max) {
    niter <- niter + 1
    f_fin <- scale(matrix_ff(V, out$beta, out$alpha, k))
    out <- betaf(V, f_fin, k, sel)
    mse_fin <- out$mse / m
    criter <- 1 - mse_fin/mse_ini
    mse_ini <- mse_fin
  }
  if (niter >= niter_max){
    warning('Iterations did not converge. Consider increasing niter_max.')
  }
  out$mse <- out$mse / m
  out$f <- f_fin
  return(out)
  
}


gdpc <- function(Z, k, f_ini = NULL, tol = 1e-4, niter_max = 500, crit = 'AIC') {
  # All computations are done using leads rather than lags. However, the final result is outputted
  # using lags using the auxiliary leads2lags function.
  #INPUT
  #Z: data matrix each column is a different time series
  #k: number of lags used to predict
  #tol: relative precision, stopping criterion
  #niter_max: maximum number of iterations
  #f_ini: (optional) initial principal component. If no argument is passed, the standard
  #first principal component with k lags is used
  #crit: AIC or BIC
  #OUTPUT
  #f: final principal component
  #beta: matrix beta corresponding to f
  #alpha: matrix alpha corresponding to f
  #mse:  mean squared error (in N and m)
  #crit: criterion
  #res: matrix of residuals
  
  sel = 1
  if (crit=='BIC'){
    sel = 2
  }
  out <- gdpc.priv(t(Z), k, f_ini, tol, niter_max, sel)
  out <- leads2lags(out)
  return(out)
}

leads2lags <- function(out){
  # Auxiliary function used to pass from leads to lags form.
  k <- ncol(out$beta) - 1 #number of leads
  out$beta <- out$beta[,(k+1):1]
  if (k !=0 ){
    out$initial_f <- out$f[1:k]
  }
  out$f <- out$f[(k+1):length(out$f)]
  out$res <- t(out$res)
  return(out)
}