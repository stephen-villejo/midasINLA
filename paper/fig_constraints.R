


##### Almon #####

vals <- vector(mode = "list", length = 2)
vals[[1]]$gamma1 = 0.002 # rapidly declining weights
vals[[1]]$gamma2 = -.005
vals[[2]]$gamma1 = 0.001 # slowly declining weights
vals[[2]]$gamma2 = -.0005

compile_weights <- vector(mode = "list", length = 2)
for(j in 1:2){
  for(i in 0:100){
    compile_weights[[j]] <- c(compile_weights[[j]],
                              exp((vals[[j]]$gamma1*(i^1)) +
                                    vals[[j]]$gamma2*(i^2)))
  }
  compile_weights[[j]] <- compile_weights[[j]]/sum(compile_weights[[j]])
}

png("paper/figures/Almon_weights_illustration.png", width=25, height=10, units = 'cm', res = 300)
par(mfrow=c(1,2))
plot(compile_weights[[1]], type="l", col="blue", lwd = 2,
     ylab="weights", xlab="lag",
     main=expression(paste(gamma[1]==0.002,",", " ",gamma[2]==-0.005)))
plot(compile_weights[[2]], type="l", col="blue", lwd = 2,
     ylab="weights", xlab="lag",
     main=expression(paste(gamma[1]==0.001,",", " ",gamma[2]==-0.0005)))
dev.off()



##### Beta #####

vals <- vector(mode = "list", length = 2)
vals[[1]]$gamma1 = 1
vals[[1]]$gamma2 = 4
vals[[2]]$gamma1 = 1
vals[[2]]$gamma2 = 8

compile_weights <- vector(mode = "list", length = 2)
for(j in 1:2){
  for(i in 0:100){
    x <- 0.0001 + (1-0.0001)*((i-1)/(100-1))
    compile_weights[[j]] <- c(compile_weights[[j]],
                              (x^(vals[[j]]$gamma1-1))*((1-x))^(vals[[j]]$gamma2-1))
  }
  compile_weights[[j]] <- compile_weights[[j]]/sum(compile_weights[[j]])

}

png("paper/figures/Beta_weights_illustration.png", width=25, height=10, units = 'cm', res = 300)
par(mfrow=c(1,2))
plot(compile_weights[[1]], type="l", col="blue", lwd = 2,
     ylab="weights", xlab="lag",
     main=expression(paste(gamma[1]==1,",", " ",gamma[2]==4)))
plot(compile_weights[[2]], type="l", col="blue", lwd = 2,
     ylab="weights", xlab="lag",
     main=expression(paste(gamma[1]==1,",", " ",gamma[2]==8)))
dev.off()



##### Hyperbolic scheme #####

vals <- vector(mode = "list", length = 2)
vals[[1]] = 0.5
vals[[2]] = 0.9

compile_weights <- vector(mode = "list", length = 2)
for(j in 1:2){
  for(lag in 0:100){
    compile_weights[[j]] <- c(compile_weights[[j]],
                              gamma(lag+vals[[j]]) / (gamma(lag+1)*gamma(vals[[j]])))
  }
  compile_weights[[j]] <- compile_weights[[j]]/sum(compile_weights[[j]])
}

png("paper/figures/Hyperbolic_scheme_weights_illustration.png", width=25, height=10, units = 'cm', res = 300)
par(mfrow=c(1,2))
plot(compile_weights[[1]], type="l", col="blue", lwd = 2,
     ylab="weights", xlab="lag",
     main=expression(paste(gamma==0.5)))
plot(compile_weights[[2]], type="l", col="blue", lwd = 2,
     ylab="weights", xlab="lag",
     main=expression(paste(gamma==0.9)))
dev.off()





##### Gaussian constraint #####

vals <- vector(mode = "list", length = 2)
vals[[1]]$gamma1 = 20
vals[[1]]$gamma2 = 15
vals[[2]]$gamma1 = 20
vals[[2]]$gamma2 = 25


compile_weights <- vector(mode = "list", length = 2)
for(j in 1:2){
  for(lag in 0:100){
    compile_weights[[j]] <- c(compile_weights[[j]],
                              exp( -(lag-vals[[j]]$gamma1)^2 / (2 * (vals[[j]]$gamma2^2)) ))
  }
  compile_weights[[j]] <- compile_weights[[j]]/sum(compile_weights[[j]])
}

png("paper/figures/Gaussian_weights_illustration.png", width=25, height=10, units = 'cm', res = 300)
par(mfrow=c(1,2))
plot(compile_weights[[1]], type="l", col="blue", lwd = 2,
     ylab="weights", xlab="lag",
     main=expression(paste(gamma[1]==20,",", " ",gamma[2]==15)))
plot(compile_weights[[2]], type="l", col="blue", lwd = 2,
     ylab="weights", xlab="lag",
     main=expression(paste(gamma[1]==20,",", " ",gamma[2]==25)))
dev.off()
