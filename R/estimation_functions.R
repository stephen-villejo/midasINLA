


#' Create lag matrix
#' @param tsdata High-frequency covariate data
#' @param lags Lags to be considered
#' @param frequency Number of covariates values associated with each response
#' @return Matrix of lagged values
#' @export
create_lag_Xmatrix <- function(tsdata,
                               lags,
                               frequency){
  X_matrix <- midasr::mls(tsdata, lags, frequency)
  compile_names <- c()
  for(i in 0:(ncol(X_matrix)-1)){
    compile_names <- c(compile_names,paste0("lag",i))
  }
  colnames(X_matrix) <- compile_names

  return(X_matrix)
}


#' Prepare MIDAS objects for INLA estimation
#' @param xdata High-frequency covariate data
#' @param ydata Low-frequency response data
#' @param constraint constraint function for lag-response association
#' @param K Lags to be considered
#' @param m Number of covariates values associated with each response
#' @param lagY Number of lagged response values to be considered
#' @param family Likelihood family for response data
#' @param Ntrials Number of trials for a binomial family response
#' @return A list of objects that will be used to run inla
#' @export
fit_Minla <- function(xdata,
                      ydata,
                      constraint,
                      K,
                      m,
                      lagY = 0,
                      family = "gaussian",
                      Ntrials){

  X_matrix <- create_lag_Xmatrix(tsdata = xdata,
                                 lags = K,
                                 frequency = m)

  temp_data <- as.data.frame(X_matrix)

  rm.row <- which(stats::complete.cases(temp_data) == FALSE)

  if(family == "gaussian"){
    if(length(rm.row > 0) > 0){
      temp_data$y <- as.vector(ydata)
      temp_data <- temp_data[-rm.row,]
    }else{
      temp_data$y <- as.vector(ydata)
    }
  }else if(family == "binomial"){
    if(length(rm.row > 0) > 0){
      temp_data$y <- as.vector(ydata)
      #temp_data$Ntrials <- as.vector(Ntrials)
      temp_data <- temp_data[-rm.row,]
    }else{
      temp_data$y <- as.vector(ydata)
    }
  }else if(family == "poisson"){
    if(length(rm.row > 0) > 0){
      temp_data$y <- as.vector(ydata)
      temp_data <- temp_data[-rm.row,]
    }else{
      temp_data$y <- as.vector(ydata)
    }
  }


  if(constraint == "beta"){
    rgen = INLA::inla.rgeneric.define(model = rgeneric.Beta.midas,
                                      x = temp_data)
  }else if(constraint == "beta2"){
    rgen = INLA::inla.rgeneric.define(model = rgeneric.Beta2.midas,
                                      x = temp_data)
  }else if(constraint == "almon2"){
    rgen = INLA::inla.rgeneric.define(model = rgeneric.Almon2.midas,
                                      x = temp_data)
  }else if(constraint == "almon3"){
    rgen = INLA::inla.rgeneric.define(model = rgeneric.Almon3.midas,
                                      x = temp_data)
  }else if(constraint == "hyperbolic"){
    rgen = INLA::inla.rgeneric.define(model = rgeneric.Hyperbolic.midas,
                                      x = temp_data)
  }else if(constraint == "gaussian"){
    rgen = INLA::inla.rgeneric.define(model = rgeneric.Gaussian.midas,
                                      x = temp_data)
  }

  if(lagY == 0){
    data = temp_data
  }else{
    # lagYdata <- matrix(NA, nrow = nrow(temp_data), ncol = lagY)
    # for(i in 1:lagY){
    #   lagYdata[,i] <- as.vector(mls(temp_data$y, i, 1))
    # }
    # lagYdata <- as.data.frame(lagYdata)
    # compile_names <- c()
    # for(i in 1:lagY){
    #   compile_names <- c(compile_names,paste0("lagy",i))
    # }
    # names(lagYdata) <- compile_names
    # data <- cbind(temp_data, lagYdata)
  }

  if(family == "binomial"){
    if(length(rm.row > 0) > 0){
      return(out = list(data = data,
                        X_matrix = X_matrix,
                        rgen = rgen,
                        rm.row = rm.row,
                        Ntrials = Ntrials[-rm.row]))
    }else{
      return(out = list(data = data,
                        X_matrix = X_matrix,
                        rgen = rgen,
                        Ntrials = Ntrials))
    }

  }else{
    return(out = list(data = data,
                      X_matrix = X_matrix,
                      rgen = rgen,
                      rm.row = rm.row))
  }


}


#' Prepare spatial MIDAS objects for INLA estimation
#' @param xdata High-frequency covariate data
#' @param ydata Low-frequency response data
#' @param loc_x Spatial index for the covariate data
#' @param loc_y Spatial index for the response data
#' @param constraint constraint function for lag-response association
#' @param K Lags to be considered
#' @param m Number of covariates values associated with each response
#' @param lagY Number of lagged response values to be considered
#' @param family Likelihood family for response data
#' @param Ntrials Number of trials for a binomial family response
#' @return A list of objects that will be used to run inla
#' @export
fit_Minla_spatial <- function(xdata,
                              ydata,
                              loc_x,
                              loc_y,
                              constraint,
                              K,
                              m,
                              lagY = 0,
                              family = "gaussian",
                              Ntrials){

  compile_X_matrix <- NULL
  counter_time <- c()
  counter_loc <- c()
  unique_loc_x <- unique(loc_x)
  for(i in unique_loc_x){
    temp_xdata <- xdata[which(loc_x == i)]
    X_matrix <- create_lag_Xmatrix(tsdata = temp_xdata,
                                   lags = K,
                                   frequency = m)

    compile_X_matrix <- rbind(compile_X_matrix,
                              X_matrix)
    counter_time <- c(counter_time, 1:length(which(stats::complete.cases(X_matrix) == TRUE)))
    counter_loc <- c(counter_loc, rep(i, times = length(which(stats::complete.cases(X_matrix) == TRUE))))
  }

  temp_data <- as.data.frame(compile_X_matrix)

  rm.row <- which(stats::complete.cases(temp_data) == FALSE)

  if(family == "gaussian"){
    if(length(rm.row > 0) > 0){
      temp_data$y <- as.vector(ydata)
      temp_data <- temp_data[-rm.row,]
    }else{
      temp_data$y <- as.vector(ydata)
    }
  }else if(family == "binomial"){
    if(length(rm.row > 0) > 0){
      temp_data$y <- as.vector(ydata)
      temp_data <- temp_data[-rm.row,]
    }else{
      temp_data$y <- as.vector(ydata)
    }
  }else if(family == "poisson"){
    if(length(rm.row > 0) > 0){
      temp_data$y <- as.vector(ydata)
      temp_data <- temp_data[-rm.row,]
    }else{
      temp_data$y <- as.vector(ydata)
    }
  }


  if(constraint == "beta"){
    rgen = INLA::inla.rgeneric.define(model = rgeneric.Beta.midas,
                                      x = temp_data)
  }else if(constraint == "beta2"){
    rgen = INLA::inla.rgeneric.define(model = rgeneric.Beta2.midas,
                                      x = temp_data)
  }else if(constraint == "almon2"){
    rgen = INLA::inla.rgeneric.define(model = rgeneric.Almon2.midas,
                                      x = temp_data)
  }else if(constraint == "almon3"){
    rgen = INLA::inla.rgeneric.define(model = rgeneric.Almon3.midas,
                                      x = temp_data)
  }else if(constraint == "hyperbolic"){
    rgen = INLA::inla.rgeneric.define(model = rgeneric.Hyperbolic.midas,
                                      x = temp_data)
  }

  if(lagY == 0){
    data = temp_data
  }else{

  }

  if(family == "binomial"){
    if(length(rm.row > 0) > 0){
      return(out = list(data = data,
                        X_matrix = compile_X_matrix,
                        rgen = rgen,
                        rm.row = rm.row,
                        Ntrials = Ntrials[-rm.row],
                        counter_time = counter_time,
                        counter_loc = counter_loc))
    }else{
      return(out = list(data = data,
                        X_matrix = compile_X_matrix,
                        rgen = rgen,
                        Ntrials = Ntrials,
                        counter_time = counter_time,
                        counter_loc = counter_loc))
    }

  }else{
    return(out = list(data = data,
                      X_matrix = X_matrix,
                      rgen = rgen,
                      rm.row = rm.row))
  }


}


#' Make predictions form a MIDAS model output
#' @param model MIDAS model output
#' @param data Data object to fit MIDAS
#' @param family  Likelihood family for response data
#' @param Ntrials Number of trials for a binomial family response
#' @param nsamples Number of posterior samples
#' @return A list containing the predictions and the samples
#' @export
predict_midas <- function(model,
                          data,
                          family = "gaussian",
                          Ntrials,
                          nsamples){

  if(family == "gaussian"){

    # generate posterior samples
    temp <- INLA::inla.posterior.sample(n = nsamples, model)

    # samples of linear predictor and gaussian variance
    latent_predictor <- sapply(1:nsamples, function(i) temp[[i]]$latent[1:nrow(data)])
    sigma2 <- sapply(1:nsamples, function(i) 1/unname(temp[[i]]$hyperpar[which(names(temp[[i]]$hyperpar) == "Precision for the Gaussian observations")]))
    samples <- vector(mode = "list", length = 2)
    samples[[1]] <- latent_predictor
    samples[[2]] <- sigma2
    names(samples) <- c("latent_predictor","sigma2")

    # samples of y
    sample_y <- sapply(1:nsamples,
                       function(i) temp[[i]]$latent[1:nrow(data)] +
                         stats::rnorm(nrow(data),
                                      mean = 0,
                                      sd = sqrt(1/unname(temp[[i]]$hyperpar[which(names(temp[[i]]$hyperpar) == "Precision for the Gaussian observations")]))))
    # summarise samples of y
    computed_y <- list(mean = rowMeans(sample_y),
                       sd = apply(sample_y, 1, stats::sd),
                       q2.5 = matrixStats::rowQuantiles(sample_y, probs = 0.025),
                       q97.5 = matrixStats::rowQuantiles(sample_y, probs = 0.975))

  }else if(family == "poisson"){

    # generate posterior samples
    temp <- INLA::inla.posterior.sample(n = nsamples, model)

    # samples of linear predictor
    latent_predictor <- sapply(1:nsamples, function(i) temp[[i]]$latent[1:nrow(data)])
    samples <- vector(mode = "list", length = 1)
    samples[[1]] <- latent_predictor
    names(samples) <- c("latent_predictor")

    # samples of y
    sample_y <- sapply(1:nsamples,
                       function(i) stats::rpois(n = nrow(data),
                                                lambda = exp(temp[[i]]$latent[1:nrow(data)])))

    # summarise samples of y
    computed_y <- list(mean = rowMeans(sample_y),
                       sd = apply(sample_y, 1, stats::sd),
                       q2.5 = matrixStats::rowQuantiles(sample_y, probs = 0.025),
                       q97.5 = matrixStats::rowQuantiles(sample_y, probs = 0.975))

  }else if(family == "binomial"){

    # generate posterior samples
    temp <- INLA::inla.posterior.sample(n = nsamples, model)

    # samples of linear predictor
    latent_predictor <- sapply(1:nsamples, function(i) temp[[i]]$latent[1:nrow(data)])
    samples <- vector(mode = "list", length = 1)
    samples[[1]] <- latent_predictor
    names(samples) <- c("latent_predictor")

    # samples of y
    sample_y <- sapply(1:nsamples,
                       function(i) stats::rbinom(n = nrow(data),
                                                 size = Ntrials,
                                                 prob = exp(temp[[i]]$latent[1:nrow(data)])/(1+exp(temp[[i]]$latent[1:nrow(data)]))))
    # summarise samples of y
    computed_y <- list(mean = rowMeans(sample_y),
                       sd = apply(sample_y, 1, stats::sd),
                       q2.5 = matrixStats::rowQuantiles(sample_y, probs = 0.025),
                       q97.5 = matrixStats::rowQuantiles(sample_y, probs = 0.975))

  }

  return(out = list(computed_y = computed_y,
                    samples = samples))

}

#' Compute scores for model assessment
#' @param y Response vector
#' @param Midas_object Midas object
#' @param pred_res Predictions from the MIDAS model
#' @param family  Likelihood family for response data
#' @return Computed scores
#' @export
compute_forecast_scores <- function(y,
                                    Midas_object,
                                    pred_res,
                                    family = "binomial"){


  if(family == "gaussian"){

    y_forecast <- y[-c(1:(length(y)-length(Midas_object$data$y)))][-c(1:max(which(!is.na(Midas_object$data$y))))]

    SE = (y_forecast - pred_res$computed_y$mean[-c(1:max(which(!is.na(Midas_object$data$y))))])^2

    compile_logscore <- matrix(NA,nrow=length(y_forecast),ncol=ncol(pred_res$samples$latent_predictor))
    for(i in 1:ncol(pred_res$samples$latent_predictor)){
      score <- stats::dnorm(y_forecast,
                            mean = pred_res$samples$latent_predictor[-c(1:max(which(!is.na(Midas_object$data$y)))),i],
                            sd = sqrt(pred_res$samples$sigma2[i]))
      compile_logscore[,i] <- log(score)
    }

    post_E <- rowMeans(pred_res$samples$latent_predictor[-c(1:max(which(!is.na(Midas_object$data$y)))),])
    post_Var <- mean(pred_res$samples$sigma2) + apply(pred_res$samples$latent_predictor, 1, stats::var)[-c(1:max(which(!is.na(Midas_object$data$y))))]

    scores <- data.frame(
      SE = SE,
      DS = (y_forecast - post_E)^2 / post_Var + log(post_Var),
      NLS = -rowMeans(compile_logscore)
    )

  }else if(family == "poisson"){

    y_forecast <- y[-c(1:(length(y)-length(Midas_object$data$y)))][-c(1:max(which(!is.na(Midas_object$data$y))))]

    SE = (y_forecast - pred_res$computed_y$mean[-c(1:max(which(!is.na(Midas_object$data$y))))])^2

    compile_logscore <- matrix(NA,nrow=length(y_forecast),ncol=ncol(pred_res$samples$latent_predictor))
    for(i in 1:ncol(pred_res$samples$latent_predictor)){
      score <- stats::dpois(y_forecast,
                            lambda = exp(pred_res$samples$latent_predictor[-c(1:max(which(!is.na(Midas_object$data$y)))),i]))
      compile_logscore[,i] <- log(score)
    }

    post_E <- rowMeans(exp(pred_res$samples$latent_predictor[-c(1:max(which(!is.na(Midas_object$data$y)))),]))
    post_Var <- rowMeans(exp(pred_res$samples$latent_predictor[-c(1:max(which(!is.na(Midas_object$data$y)))),])) +
      apply(exp(pred_res$samples$latent_predictor), 1, stats::var)[-c(1:max(which(!is.na(Midas_object$data$y))))]

    scores <- data.frame(
      SE = SE,
      DS = (y_forecast - post_E)^2 / post_Var + log(post_Var),
      NLS = -rowMeans(compile_logscore)
    )

  }else if(family == "binomial"){

    y_forecast <- y[-c(1:(length(y)-length(Midas_object$data$y)))][-c(1:max(which(!is.na(Midas_object$data$y))))]

    SE = (y_forecast - pred_res$computed_y$mean[-c(1:max(which(!is.na(Midas_object$data$y))))])^2

    compile_logscore <- matrix(NA,nrow=length(y_forecast),ncol=ncol(pred_res$samples$latent_predictor))
    for(i in 1:ncol(pred_res$samples$latent_predictor)){
      score <- stats::dbinom(y_forecast,
                             size = Midas_object$Ntrials[-c(1:max(which(!is.na(Midas_object$data$y))))],
                             prob = exp(pred_res$samples$latent_predictor[-c(1:max(which(!is.na(Midas_object$data$y)))),i])/(1+exp(pred_res$samples$latent_predictor[-c(1:max(which(!is.na(Midas_object$data$y)))),i])))
      compile_logscore[,i] <- log(score)
    }

    p_temp <- exp(pred_res$samples$latent_predictor[-c(1:max(which(!is.na(Midas_object$data$y)))),])/(1+exp(pred_res$samples$latent_predictor[-c(1:max(which(!is.na(Midas_object$data$y)))),]))
    post_E <- rowMeans(p_temp) * Midas_object$Ntrials[-c(1:max(which(!is.na(Midas_object$data$y))))]
    post_Var <- (Midas_object$Ntrials[-c(1:max(which(!is.na(Midas_object$data$y))))] * rowMeans(p_temp * (1 - p_temp))) +
      ((Midas_object$Ntrials[-c(1:max(which(!is.na(Midas_object$data$y))))]^2) * apply(p_temp, 1, stats::var))

    scores <- data.frame(
      SE = SE,
      DS = (y_forecast - post_E)^2 / post_Var + log(post_Var),
      NLS = -rowMeans(compile_logscore)
    )

  }

  return(list(scores = scores))

}



#' Compute weights estimates
#' @param model midas model result
#' @param constraint constraint function
#' @param lag_k maximum lag considered
#' @return estimates of the weights
#' @export
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
  }


  return(out = weights_ests_df)

}
