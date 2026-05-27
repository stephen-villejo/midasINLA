






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


compute_weights <- function(model,
                            constraint,
                            lag_k){

  n.samples = 200

  if(constraint == "hyperbolic"){

    theta2_samples <- INLA::inla.rmarginal(n.samples, model$marginals.hyperpar$`Theta2 for idx`)
    gamma_samples <- exp(theta2_samples)/(1+exp(theta2_samples))
    compile_sum <- vector(length = n.samples)
    for(lag in 0:lag_k){
      temp <- gamma(lag+gamma_samples)/(gamma(lag+1)*gamma(gamma_samples))
      assign(paste0("psi_sample_",lag),temp)
      compile_sum <- compile_sum + temp
    }
    for(lag in 0:lag_k){
      temp <- get(paste0("psi_sample_", lag))
      assign(paste0("w",lag),temp/compile_sum)
    }
    compile_w_list <- vector("list", length = lag_k + 1)
    for(i in 0:lag_k){
      temp <- get(paste0("w",i))
      compile_w_list[[i+1]] <- temp
    }
    mean <- lapply(1:(lag_k+1), function(x) mean(compile_w_list[[x]]))
    q2.5 <- lapply(1:(lag_k+1), function(x) stats::quantile(compile_w_list[[x]], probs = 0.025))
    q97.5 <- lapply(1:(lag_k+1), function(x) stats::quantile(compile_w_list[[x]], probs = 0.975))

    weights_ests_df <- data.frame(mean = unlist(mean),
                                  q2.5 = unlist(q2.5),
                                  q97.5 = unlist(q97.5),
                                  lag = 0:lag_k)
  }else if(constraint == "gaussian"){
    theta1_samples <- INLA::inla.rmarginal(200, model$marginals.hyperpar$`Theta1 for idx`)
    theta2_samples <- INLA::inla.rmarginal(200, model$marginals.hyperpar$`Theta2 for idx`)

    mu_val_samples <- lag_k * (1 / (1 + exp(-theta1_samples)))
    sigma_val_samples <- exp(theta2_samples)
    compile_sum <- vector(length = 200)
    for(lag in 0:lag_k){
      temp <- exp( -(lag-mu_val_samples)^2 / (2 * (sigma_val_samples^2)) )
      assign(paste0("psi_sample_",lag),temp)
      compile_sum <- compile_sum + temp
    }
    for(lag in 0:lag_k){
      temp <- get(paste0("psi_sample_", lag))
      assign(paste0("w",lag),temp/compile_sum)
    }
    compile_w_list <- vector("list", length = lag_k + 1)
    for(i in 0:lag_k){
      temp <- get(paste0("w",i))
      compile_w_list[[i+1]] <- temp
    }
    mean <- lapply(1:(lag_k+1), function(x) mean(compile_w_list[[x]]))
    q2.5 <- lapply(1:(lag_k+1), function(x) stats::quantile(compile_w_list[[x]],probs=0.025))
    q97.5 <- lapply(1:(lag_k+1), function(x) stats::quantile(compile_w_list[[x]],probs=0.975))

    weights_ests_df <- data.frame(mean = unlist(mean),
                                  q2.5 = unlist(q2.5),
                                  q97.5 = unlist(q97.5),
                                  lag = 0:lag_k)

  }


  return(out = weights_ests_df)

}
