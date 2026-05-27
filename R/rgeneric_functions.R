
#' Rgeneric MIDAS model with Beta polynomial constraint where gamma_1 is fixed to 1
#'
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using a Beta polynomial structure.
#'

#' @param cmd Character string indicating the INLA command.
#'   One of \code{"graph"}, \code{"Q"}, \code{"mu"}, \code{"initial"},
#'   \code{"log.norm.const"}, \code{"log.prior"}, or \code{"quit"}.
#' @param theta Numeric vector of hyperparameters controlling the MIDAS weights.
#'
#' @return Depends on \code{cmd}:
#' \itemize{
#'   \item \code{graph}: Sparse precision structure
#'   \item \code{Q}: Precision matrix
#'   \item \code{mu}: Mean vector
#'   \item \code{initial}: Initial values for \code{theta}
#'   \item \code{log.norm.const}: Normalizing constant (numeric(0))
#'   \item \code{log.prior}: Log prior density
#' }
#'
#' @details
#' The MIDAS lag weights are constructed using a normalized Beta polynomial
#' transformation of lag indices. The parameter \code{theta[1]} controls the
#' shape of the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
#' @export
rgeneric.Beta.midas = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                       "log.prior", "quit"),
                               theta = NULL){

  envir = parent.env(environment())
  x <- envir$x

  ## artificial high precision to be added to the mean-model
  prec.high = exp(15)

  interpret.theta = function() {

    lag_k <- ncol(x)-2
    gamma2 <- exp(theta[1L]) + 1
    beta1 <- theta[2L]

    psi <- numeric(lag_k + 1)
    for(lag in 0:lag_k){
      x_temp <- 0.0001 + (1-0.0001)*((lag-1)/(lag_k-1))
      psi[lag + 1] <- 1 + (1-x_temp)^(gamma2-1)
    }

    w <- psi / sum(psi)

    out_list <- list()
    for(lag in 0:lag_k){
      out_list[[paste0("w", lag)]] <- w[lag + 1]
    }

    out_list[["beta1"]] <- beta1

    return(out_list)
  }
  graph = function() {
    G = Matrix::Diagonal(n = length(x$lag0), x=1)
    return(G)
  }
  Q = function() {
    Q = prec.high * graph()
    return(Q)
  }
  mu = function() {
    par = interpret.theta()

    lag_k = ncol(x)-2
    compile_lag_label <- c()
    for(lag in 0:lag_k){
      compile_lag_label <- c(compile_lag_label, paste0("lag",lag))
    }
    compile_w_label <- c()
    for(lag in 0:lag_k){
      compile_w_label <- c(compile_w_label, paste0("w",lag))
    }

    agg <- 0
    for(lag in 0:lag_k){
      agg <- agg + par[[compile_w_label[[lag+1]]]] * x[,which(names(x) == compile_lag_label[[lag+1]])]
    }

    return(par$beta1 * agg)
  }
  log.norm.const = function() {
    return(numeric(0))
  }
  log.prior = function() {
    val = (stats::dnorm(theta[1L], mean=0, sd=1, log=TRUE) +
             stats::dnorm(theta[2L], mean=0, sd=1, log=TRUE))
    return(val)
  }
  initial = function() {
    return(rep(1, 2))
  }
  quit = function() {
    return(invisible())
  }

  val = do.call(match.arg(cmd), args = list())
  return(val)
}


#' Rgeneric MIDAS model with Beta polynomial constraint
#'
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using a Beta polynomial structure.
#'

#' @param cmd Character string indicating the INLA command.
#'   One of \code{"graph"}, \code{"Q"}, \code{"mu"}, \code{"initial"},
#'   \code{"log.norm.const"}, \code{"log.prior"}, or \code{"quit"}.
#' @param theta Numeric vector of hyperparameters controlling the MIDAS weights.
#'
#' @return Depends on \code{cmd}:
#' \itemize{
#'   \item \code{graph}: Sparse precision structure
#'   \item \code{Q}: Precision matrix
#'   \item \code{mu}: Mean vector
#'   \item \code{initial}: Initial values for \code{theta}
#'   \item \code{log.norm.const}: Normalizing constant (numeric(0))
#'   \item \code{log.prior}: Log prior density
#' }
#'
#' @details
#' The MIDAS lag weights are constructed using a normalized Beta polynomial
#' transformation of lag indices. The parameters \code{theta[1]} and \code{theta[2]}
#' control the shape of the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
#' @export
rgeneric.Beta2.midas = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                        "log.prior", "quit"),
                                theta = NULL){

  envir = parent.env(environment())
  x <- envir$x

  ## artificial high precision to be added to the mean-model
  prec.high = exp(15)

  interpret.theta = function() {

    lag_k <- ncol(x)-2
    gamma1 <- exp(theta[1L]) + 1
    gamma2 <- exp(theta[2L]) + 1
    beta1 <- theta[3L]
    for(lag in 0:lag_k){
      x_temp <- 0.0001 + (1-0.0001)*((lag-1)/(lag_k-1))
      psi[lag + 1] <- (x_temp^(gamma1-1))*((1-x_temp)^(gamma2-1))
    }

    w = psi / sum(psi)

    out_list <- list()

    for(lag in 0:lag_k){
      out_list[[paste0("w", lag)]] <- w[lag + 1]
    }

    out_list[["beta1"]] <- beta1

    return(out_list)
  }
  graph = function() {
    G = Matrix::Diagonal(n = length(x$lag0), x=1)
    return(G)
  }
  Q = function() {
    Q = prec.high * graph()
    return(Q)
  }
  mu = function() {
    par = interpret.theta()

    lag_k = ncol(x)-2
    compile_lag_label <- c()
    for(lag in 0:lag_k){
      compile_lag_label <- c(compile_lag_label, paste0("lag",lag))
    }
    compile_w_label <- c()
    for(lag in 0:lag_k){
      compile_w_label <- c(compile_w_label, paste0("w",lag))
    }

    agg <- 0
    for(lag in 0:lag_k){
      agg <- agg + par[[compile_w_label[[lag+1]]]] * x[,which(names(x) == compile_lag_label[[lag+1]])]
    }

    return(par$beta1 * agg)
  }
  log.norm.const = function() {
    return(numeric(0))
  }
  log.prior = function() {
    val = (stats::dnorm(theta[1L], mean=0, sd=1, log=TRUE) +
             stats::dnorm(theta[2L], mean=0, sd=1, log=TRUE) +
             stats::dnorm(theta[3L], mean=0, sd=1, log=TRUE))
    return(val)
  }
  initial = function() {
    return(rep(1, 3))
  }
  quit = function() {
    return(invisible())
  }

  val = do.call(match.arg(cmd), args = list())
  return(val)
}



#' Rgeneric MIDAS model with Almon polynomial constraint, d = 2
#'
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using an Almon polynomial structure.
#'

#' @param cmd Character string indicating the INLA command.
#'   One of \code{"graph"}, \code{"Q"}, \code{"mu"}, \code{"initial"},
#'   \code{"log.norm.const"}, \code{"log.prior"}, or \code{"quit"}.
#' @param theta Numeric vector of hyperparameters controlling the MIDAS weights.
#'
#' @return Depends on \code{cmd}:
#' \itemize{
#'   \item \code{graph}: Sparse precision structure
#'   \item \code{Q}: Precision matrix
#'   \item \code{mu}: Mean vector
#'   \item \code{initial}: Initial values for \code{theta}
#'   \item \code{log.norm.const}: Normalizing constant (numeric(0))
#'   \item \code{log.prior}: Log prior density
#' }
#'
#' @details
#' The MIDAS lag weights are constructed using a normalized Almon polynomial
#' transformation of lag indices. The parameters \code{theta[2]} and
#' \code{theta[3]} control the shape of the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
#' @export
rgeneric.Almon2.midas = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                         "log.prior", "quit"),
                                 theta = NULL){

  envir = parent.env(environment())
  x <- envir$x

  ## artificial high precision to be added to the mean-model
  prec.high = exp(15)

  interpret.theta = function() {

    lag_k <- ncol(x)-2
    gamma1 <- 0.01*sin(theta[2L])
    gamma2 <- 0.01*sin(theta[3L])
    for(lag in 0:lag_k){
      temp <- exp(gamma1*(lag^1) + gamma2*(lag^2))
      assign(paste0("psi",lag), temp)
    }

    compile_sum <- 0
    for(lag in 0:lag_k){
      compile_sum <- compile_sum + get(paste0("psi", lag))
    }

    for(lag in 0:lag_k){
      temp <- get(paste0("psi",lag))
      assign(paste0("w",lag),temp/compile_sum)
    }

    out_list <- vector(mode = "list", length = 1 + lag_k + 1)
    out_list[[1]] <- theta[1L]
    for(lag in 0:lag_k){
      out_list[[lag+2]] <- get(paste0("w",lag))
    }

    compile_names <- c()
    for(lag in 0:lag_k){
      compile_names <- c(compile_names, paste0("w",lag))
    }
    names_vec <- c("beta1", compile_names)
    names(out_list) <- names_vec
    return(out_list)

  }
  graph = function() {
    G = Matrix::Diagonal(n = length(x$lag0), x=1)
    return(G)
  }
  Q = function() {
    Q = prec.high * graph()
    return(Q)
  }
  mu = function() {
    par = interpret.theta()

    lag_k = ncol(x)-2
    compile_lag_label <- c()
    for(lag in 0:lag_k){
      compile_lag_label <- c(compile_lag_label, paste0("lag",lag))
    }
    compile_w_label <- c()
    for(lag in 0:lag_k){
      compile_w_label <- c(compile_w_label, paste0("w",lag))
    }

    agg <- 0
    for(lag in 0:lag_k){
      agg <- agg + par[[compile_w_label[[lag+1]]]] * x[,which(names(x) == compile_lag_label[[lag+1]])]
    }

    return(par$beta1 * agg)
  }
  log.norm.const = function() {
    return(numeric(0))
  }
  # log.prior = function() {
  #   par = interpret.theta()
  #   n.params <- length(par)
  #   val <- 0
  #   for(i in 1:n.params){
  #     temp <- dnorm(par[[i]], mean=0, sd=1, log=TRUE)
  #     val <- val + temp
  #   }
  # }
  log.prior = function() {
    val = (stats::dnorm(theta[1L], mean=0, sd=1, log=TRUE) +
             stats::dnorm(theta[2L], mean=0, sd=1, log=TRUE) +
             stats::dnorm(theta[3L], mean=0, sd=1, log=TRUE))
    return(val)
  }
  initial = function() {
    return(rep(1, 3))
  }
  quit = function() {
    return(invisible())
  }

  val = do.call(match.arg(cmd), args = list())
  return(val)
}



#' Rgeneric MIDAS model with Almon polynomial constraint, d = 3
#'
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using an Almon polynomial structure.
#'

#' @param cmd Character string indicating the INLA command.
#'   One of \code{"graph"}, \code{"Q"}, \code{"mu"}, \code{"initial"},
#'   \code{"log.norm.const"}, \code{"log.prior"}, or \code{"quit"}.
#' @param theta Numeric vector of hyperparameters controlling the MIDAS weights.
#'
#' @return Depends on \code{cmd}:
#' \itemize{
#'   \item \code{graph}: Sparse precision structure
#'   \item \code{Q}: Precision matrix
#'   \item \code{mu}: Mean vector
#'   \item \code{initial}: Initial values for \code{theta}
#'   \item \code{log.norm.const}: Normalizing constant (numeric(0))
#'   \item \code{log.prior}: Log prior density
#' }
#'
#' @details
#' The MIDAS lag weights are constructed using a normalized Almon polynomial
#' transformation of lag indices. The parameters \code{theta[2]}, \code{theta[3]},
#' and \code{theta[4]} control the shape of the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
#' @export
rgeneric.Almon3.midas = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                         "log.prior", "quit"),
                                 theta = NULL){

  envir = parent.env(environment())
  x <- envir$x

  ## artificial high precision to be added to the mean-model
  prec.high = exp(15)

  interpret.theta = function() {

    lag_k <- ncol(x)-2
    gamma1 <- 0.01*sin(theta[2L])
    gamma2 <- 0.01*sin(theta[3L])
    gamma3 <- 0.01*sin(theta[4L])
    for(lag in 0:lag_k){
      temp <- exp(gamma1*(lag^1) + gamma2*(lag^2) + gamma3*(lag^3))
      assign(paste0("psi",lag), temp)
    }

    compile_sum <- 0
    for(lag in 0:lag_k){
      compile_sum <- compile_sum + get(paste0("psi", lag))
    }

    for(lag in 0:lag_k){
      temp <- get(paste0("psi",lag))
      assign(paste0("w",lag),temp/compile_sum)
    }

    out_list <- vector(mode = "list", length = 1 + lag_k + 1)
    out_list[[1]] <- theta[1L]
    for(lag in 0:lag_k){
      out_list[[lag+2]] <- get(paste0("w",lag))
    }

    compile_names <- c()
    for(lag in 0:lag_k){
      compile_names <- c(compile_names, paste0("w",lag))
    }
    names_vec <- c("beta1", compile_names)
    names(out_list) <- names_vec
    return(out_list)

  }
  graph = function() {
    G = Matrix::Diagonal(n = length(x$lag0), x=1)
    return(G)
  }
  Q = function() {
    Q = prec.high * graph()
    return(Q)
  }
  mu = function() {
    par = interpret.theta()

    lag_k = ncol(x)-2
    compile_lag_label <- c()
    for(lag in 0:lag_k){
      compile_lag_label <- c(compile_lag_label, paste0("lag",lag))
    }
    compile_w_label <- c()
    for(lag in 0:lag_k){
      compile_w_label <- c(compile_w_label, paste0("w",lag))
    }

    agg <- 0
    for(lag in 0:lag_k){
      agg <- agg + par[[compile_w_label[[lag+1]]]] * x[,which(names(x) == compile_lag_label[[lag+1]])]
    }

    return(par$beta1 * agg)
  }
  log.norm.const = function() {
    return(numeric(0))
  }
  # log.prior = function() {
  #   par = interpret.theta()
  #   n.params <- length(par)
  #   val <- 0
  #   for(i in 1:n.params){
  #     temp <- dnorm(par[[i]], mean=0, sd=1, log=TRUE)
  #     val <- val + temp
  #   }
  # }
  log.prior = function() {
    val = (stats::dnorm(theta[1L], mean=0, sd=1, log=TRUE) +
             stats::dnorm(theta[2L], mean=0, sd=1, log=TRUE) +
             stats::dnorm(theta[3L], mean=0, sd=1, log=TRUE) +
             stats::dnorm(theta[4L], mean=0, sd=1, log=TRUE))
    return(val)
  }
  initial = function() {
    return(rep(1, 4))
  }
  quit = function() {
    return(invisible())
  }

  val = do.call(match.arg(cmd), args = list())
  return(val)
}


#' Rgeneric MIDAS model with Hyperbolic scheme constraint
#'
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using a Hyperbolic scheme structure.
#'

#' @param cmd Character string indicating the INLA command.
#'   One of \code{"graph"}, \code{"Q"}, \code{"mu"}, \code{"initial"},
#'   \code{"log.norm.const"}, \code{"log.prior"}, or \code{"quit"}.
#' @param theta Numeric vector of hyperparameters controlling the MIDAS weights.
#'
#' @return Depends on \code{cmd}:
#' \itemize{
#'   \item \code{graph}: Sparse precision structure
#'   \item \code{Q}: Precision matrix
#'   \item \code{mu}: Mean vector
#'   \item \code{initial}: Initial values for \code{theta}
#'   \item \code{log.norm.const}: Normalizing constant (numeric(0))
#'   \item \code{log.prior}: Log prior density
#' }
#'
#' @details
#' The MIDAS lag weights are constructed using a normalized Hyperbolic scheme
#' polynomial transformation of lag indices. The parameter \code{theta[2]}
#' controls the shape of the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
#' @export
rgeneric.Hyperbolic.midas = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                             "log.prior", "quit"),
                                     theta = NULL){

  envir = parent.env(environment())
  x <- envir$x

  ## artificial high precision to be added to the mean-model
  prec.high = exp(15)

  interpret.theta = function() {

    lag_k = ncol(x)-2
    gamma_val = exp(theta[2L])/(1+exp(theta[2L]))
    for(lag in 0:lag_k){
      temp <- gamma(lag+gamma_val) / (gamma(lag+1)*gamma(gamma_val))
      assign(paste0("psi",lag), temp)
    }

    compile_sum <- 0
    for(lag in 0:lag_k){
      compile_sum <- compile_sum + get(paste0("psi", lag))
    }

    for(lag in 0:lag_k){
      temp <- get(paste0("psi",lag))
      assign(paste0("w",lag),temp/compile_sum)
    }

    out_list <- vector(mode = "list", length = 1 + lag_k + 1)
    out_list[[1]] <- theta[1L]
    for(lag in 0:lag_k){
      out_list[[lag+2]] <- get(paste0("w",lag))
    }

    compile_names <- c()
    for(lag in 0:lag_k){
      compile_names <- c(compile_names, paste0("w",lag))
    }
    names_vec <- c("beta1", compile_names)
    names(out_list) <- names_vec
    return(out_list)

  }
  graph = function() {
    G = Matrix::Diagonal(n = length(x$lag0), x=1)
    return(G)
  }
  Q = function() {
    Q = prec.high * graph()
    return(Q)
  }
  mu = function() {
    par = interpret.theta()

    lag_k = ncol(x)-2
    compile_lag_label <- c()
    for(lag in 0:lag_k){
      compile_lag_label <- c(compile_lag_label, paste0("lag",lag))
    }
    compile_w_label <- c()
    for(lag in 0:lag_k){
      compile_w_label <- c(compile_w_label, paste0("w",lag))
    }

    agg <- 0
    for(lag in 0:lag_k){
      agg <- agg + par[[compile_w_label[[lag+1]]]] * x[,which(names(x) == compile_lag_label[[lag+1]])]
    }

    return(par$beta1 * agg)
  }
  log.norm.const = function() {
    return(numeric(0))
  }
  log.prior = function() {
    val = (stats::dnorm(theta[1L], mean=0, sd=1, log=TRUE) +
             stats::dnorm(theta[2L], mean=0, sd=1, log=TRUE))
    return(val)
  }
  initial = function() {
    return(rep(1, 2))
  }
  quit = function() {
    return(invisible())
  }
  10
  val = do.call(match.arg(cmd), args = list())
  return(val)
}



#' Rgeneric MIDAS model with Gaussian constraint
#'
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using a Gaussian structure.
#'

#' @param cmd Character string indicating the INLA command.
#'   One of \code{"graph"}, \code{"Q"}, \code{"mu"}, \code{"initial"},
#'   \code{"log.norm.const"}, \code{"log.prior"}, or \code{"quit"}.
#' @param theta Numeric vector of hyperparameters controlling the MIDAS weights.
#'
#' @return Depends on \code{cmd}:
#' \itemize{
#'   \item \code{graph}: Sparse precision structure
#'   \item \code{Q}: Precision matrix
#'   \item \code{mu}: Mean vector
#'   \item \code{initial}: Initial values for \code{theta}
#'   \item \code{log.norm.const}: Normalizing constant (numeric(0))
#'   \item \code{log.prior}: Log prior density
#' }
#'
#' @details
#' The MIDAS lag weights are constructed using a normalized Hyperbolic scheme
#' polynomial transformation of lag indices. The parametesr \code{theta[2]} and
#' \code{theta[3]} control the shape of the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
#' @export
rgeneric.Gaussian.midas = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                           "log.prior", "quit"),
                                   theta = NULL){

  envir = parent.env(environment())
  x <- envir$x

  ## artificial high precision to be added to the mean-model
  prec.high = exp(15)

  interpret.theta = function() {

    lag_k = ncol(x)-2
    mu_val <- lag_k * (1 / (1 + exp(-theta[2L])))
    sigma_val = exp(theta[3L])
    for(lag in 0:lag_k){
      temp <- exp( -(lag-mu_val)^2 / (2 * (sigma_val ^ 2)) )
      assign(paste0("psi",lag), temp)
    }

    compile_sum <- 0
    for(lag in 0:lag_k){
      compile_sum <- compile_sum + get(paste0("psi", lag))
    }

    for(lag in 0:lag_k){
      temp <- get(paste0("psi",lag))
      assign(paste0("w",lag),temp/compile_sum)
    }

    out_list <- vector(mode = "list", length = 1 + lag_k + 1)
    out_list[[1]] <- theta[1L]
    for(lag in 0:lag_k){
      out_list[[lag+2]] <- get(paste0("w",lag))
    }

    compile_names <- c()
    for(lag in 0:lag_k){
      compile_names <- c(compile_names, paste0("w",lag))
    }
    names_vec <- c("beta1", compile_names)
    names(out_list) <- names_vec
    return(out_list)

  }
  graph = function() {
    G = Matrix::Diagonal(n = length(x$lag0), x=1)
    return(G)
  }
  Q = function() {
    Q = prec.high * graph()
    return(Q)
  }
  mu = function() {
    par = interpret.theta()

    lag_k = ncol(x)-2
    compile_lag_label <- c()
    for(lag in 0:lag_k){
      compile_lag_label <- c(compile_lag_label, paste0("lag",lag))
    }
    compile_w_label <- c()
    for(lag in 0:lag_k){
      compile_w_label <- c(compile_w_label, paste0("w",lag))
    }

    agg <- 0
    for(lag in 0:lag_k){
      agg <- agg + par[[compile_w_label[[lag+1]]]] * x[,which(names(x) == compile_lag_label[[lag+1]])]
    }

    return(par$beta1 * agg)
  }
  log.norm.const = function() {
    return(numeric(0))
  }
  log.prior = function() {
    val = (stats::dnorm(theta[1L], mean=0, sd=1, log=TRUE) +
             stats::dnorm(theta[2L], mean=0, sd=1, log=TRUE) +
             stats::dnorm(theta[3L], mean=0, sd=1, log=TRUE))
    return(val)
  }
  initial = function() {
    return(rep(1, 3))
  }
  quit = function() {
    return(invisible())
  }

  val = do.call(match.arg(cmd), args = list())
  return(val)
}



