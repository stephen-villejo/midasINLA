



#' Rgeneric MIDAS spatial model with Gaussian constraint
#'
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using a Gaussian structure,
#' and where \code{beta[i]} varies for each location
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
#' polynomial transformation of lag indices. The parameters \code{theta[1]} and
#' \code{theta[2]} control the shape of the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
#' @export
rgeneric.Gaussian.varybeta.midas = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                                    "log.prior", "quit"),
                                            theta = NULL){

  envir = parent.env(environment())
  x <- envir$x

  regions = sort(unique(x$region))
  n_regions = length(regions)

  ## artificial high precision to be added to the mean-model
  prec.high = exp(15)

  interpret.theta = function() {

    lag_k = ncol(x)-3

    beta = theta[3:(n_regions + 2)]

    mu_val <- lag_k * (1 / (1 + exp(-theta[1L])))
    sigma_val = exp(theta[2L])

    psi = numeric(lag_k + 1)
    for(lag in 0:lag_k){
      psi[lag + 1] = exp( -(lag-mu_val)^2 / (2 * (sigma_val ^ 2)) )
    }
    w = psi / sum(psi)

    out_list = list()
    for(lag in 0:lag_k){
      out_list[[paste0("w", lag)]] = w[lag + 1]
    }

    for(b in 1:n_regions){
      out_list[[paste0("beta", b)]] = beta[b]
    }

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

    lag_k = ncol(x)-3
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

    beta_vec <- numeric(length(x$region))

    for(i in seq_along(x$region)){
      beta_vec[i] <- par[[paste0("beta",x$region[i])]]
    }

    return(beta_vec * agg)
  }
  log.norm.const = function() {
    return(numeric(0))
  }
  log.prior = function() {
    lp_gauss= stats::dnorm(theta[1L], mean = 0, sd = 1, log = TRUE) +
      stats::dnorm(theta[2L], mean = 0, sd = 1, log = TRUE)
    lp_beta = sum(
      stats::dnorm(theta[3:(n_regions + 2)], mean = 0, sd = 1, log = TRUE)
    )
    return(lp_gauss + lp_beta)
  }
  initial = function() {
    return(rep(0, n_regions+2))
  }
  quit = function() {
    return(invisible())
  }

  val = do.call(match.arg(cmd), args = list())
  return(val)
}




#' Rgeneric MIDAS spatial model with Hyperbolic scheme constraint
#'
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using a Hyperbolic scheme structure,
#' and where \code{beta[i]} varies for each location
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
#' polynomial transformation of lag indices. The parameter \code{theta[1]}
#' controls the shape of the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
#' @export
rgeneric.Hyperbolic.varybeta.midas = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                                      "log.prior", "quit"),
                                              theta = NULL){

  envir = parent.env(environment())
  x <- envir$x

  regions = sort(unique(x$region))
  n_regions = length(regions)

  ## artificial high precision to be added to the mean-model
  prec.high = exp(15)

  interpret.theta = function() {

    lag_k = ncol(x)-3

    beta = theta[2:(n_regions + 1)]

    gamma_val = exp(theta[1L])/(1+exp(theta[1L]))
    psi = numeric(lag_k + 1)
    for(lag in 0:lag_k){
      psi[lag + 1] = gamma(lag+gamma_val) / (gamma(lag+1)*gamma(gamma_val))
    }
    w = psi / sum(psi)

    out_list = list()
    for(lag in 0:lag_k){
      out_list[[paste0("w", lag)]] = w[lag + 1]
    }

    for(b in 1:n_regions){
      out_list[[paste0("beta", b)]] = beta[b]
    }

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

    lag_k = ncol(x)-3
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

    beta_vec <- numeric(length(x$region))

    for(i in seq_along(x$region)){
      beta_vec[i] <- par[[paste0("beta",x$region[i])]]
    }

    return(beta_vec * agg)
  }
  log.norm.const = function() {
    return(numeric(0))
  }
  log.prior = function() {
    lp_gamma = stats::dnorm(theta[1L], mean = 0, sd = 1, log = TRUE)
    lp_beta = sum(
      stats::dnorm(theta[2:(n_regions + 1)], mean = 0, sd = 1, log = TRUE)
    )
    return(lp_gamma + lp_beta)
  }
  initial = function() {
    return(rep(0, n_regions+1))
  }
  quit = function() {
    return(invisible())
  }

  val = do.call(match.arg(cmd), args = list())
  return(val)
}






#' Rgeneric MIDAS spatial model with Hyperbolic scheme constraint
#'
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using a Hyperbolic scheme structure,
#' and where \code{beta[i]} and the lag structure varies for each location
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
#' polynomial transformation of lag indices. The parameters \code{theta[i]}
#' control the shape of the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
#' @export
rgeneric.Hyperbolic.varylagstr.midas = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                                        "log.prior", "quit"),
                                                theta = NULL){

  envir = parent.env(environment())
  x <- envir$x
  lag_k_region = envir$args$lag_k_region

  regions = sort(unique(x$region))
  n_regions = length(regions)

  ## artificial high precision to be added to the mean-model
  prec.high = exp(15)

  interpret.theta = function() {

    gamma_raw <- theta[1:n_regions]
    beta <- theta[(n_regions + 1):(2*n_regions)]
    gamma_val <- exp(gamma_raw)/(1+exp(gamma_raw))

    w_region <- list()

    for(r in 1:n_regions){

      K_r <- lag_k_region[r]
      psi <- numeric(K_r + 1)

      for(lag in 0:K_r){
        psi[lag + 1] <- gamma(lag+gamma_val[r]) / (gamma(lag+1)*gamma(gamma_val[r]))
      }

      w_region[[r]] <- psi / sum(psi)

    }
    return(list(beta = beta,
                gamma = gamma_val,
                w = w_region))
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

    par <- interpret.theta()

    eta <- numeric(nrow(x))

    for(i in 1:nrow(x)){

      r <- as.integer(x$region[i])

      K_r <- lag_k_region[r]

      agg <- 0

      for(lag in 0:K_r){
        agg <- agg +
          par$w[[r]][lag + 1] *
          x[i, paste0("lag", lag)]
      }

      eta[i] <- par$beta[r] * agg
    }

    return(eta)
  }
  log.norm.const = function() {
    return(numeric(0))
  }
  log.prior = function() {
    lp_gamma <- sum(
      stats::dnorm(theta[1:n_regions],
                   mean = 0,
                   sd = 1,
                   log = TRUE))
    lp_beta <- sum(stats::dnorm(theta[(n_regions + 1):(2*n_regions)],
                                mean = 0,
                                sd = 1,
                                log = TRUE))
    return(lp_gamma + lp_beta)
  }
  initial = function() {
    return(rep(0, 2*n_regions))
  }
  quit = function() {
    return(invisible())
  }

  val = do.call(match.arg(cmd), args = list())
  return(val)
}



#' Rgeneric MIDAS spatial model with Gaussian constraint
#'
#' Defines a custom \code{rgeneric} model for use with the \code{INLA} framework,
#' implementing MIDAS-type lag weights using a Gaussian structure,
#' and where \code{beta[i]} and the lag structure varies for each location
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
#' polynomial transformation of lag indices. The parameters \code{theta[i]}
#' control the shape of the lag weighting function.
#'
#' @importFrom Matrix Diagonal
#'
#' @export
rgeneric.Gaussian.varylagstr.midas = function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                                      "log.prior", "quit"),
                                              theta = NULL){

  envir = parent.env(environment())
  x <- envir$x
  lag_k_region <- envir$args$lag_k_region

  regions = sort(unique(x$region))
  n_regions = length(regions)

  ## artificial high precision to be added to the mean-model
  prec.high = exp(15)

  interpret.theta = function() {

    mu_raw <- theta[1:n_regions]
    sigma_raw <- theta[(n_regions + 1):(2*n_regions)]
    beta <- theta[(2*n_regions + 1):(3*n_regions)]

    sigma_val = exp(sigma_raw)

    w_region <- list()

    for(r in 1:n_regions){

      K_r <- lag_k_region[r]

      mu_val <- K_r * (1 / (1 + exp(-mu_raw[r])))

      psi <- numeric(K_r + 1)
      for(lag in 0:K_r){
        psi[lag + 1] <- exp( -(lag-mu_val)^2 / (2 * (sigma_val[r] ^ 2)) )
      }

      w_region[[r]] <- psi / sum(psi)

    }

    return(list(beta = beta,
                w = w_region))

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

    eta <- numeric(nrow(x))

    for(i in 1:nrow(x)){

      r <- as.integer(x$region[i])
      K_r <- lag_k_region[r]
      agg <- 0

      for(lag in 0:K_r){
        agg <- agg + par$w[[r]][lag+1] * x[i, paste0("lag",lag)]
      }

      eta[i] <- par$beta[r] * agg

    }
    return(eta)
  }
  log.norm.const = function() {
    return(numeric(0))
  }
  log.prior = function() {

    lp_mu <- sum(
      stats::dnorm(theta[1:n_regions], 0, 1, log = TRUE)
    )
    lp_sigma <- sum(
      stats::dnorm(theta[(n_regions + 1):(2*n_regions)], 0, 1, log = TRUE)
    )
    lp_beta <- sum(
      stats::dnorm(theta[(2*n_regions + 1):(3*n_regions)], 0, 1, log = TRUE)
    )
    lp_mu + lp_sigma + lp_beta
  }
  initial = function() {
    return(rep(0, n_regions*3))
  }
  quit = function() {
    return(invisible())
  }

  val = do.call(match.arg(cmd), args = list())
  return(val)
}
