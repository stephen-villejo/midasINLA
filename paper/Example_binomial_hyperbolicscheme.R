library(midasINLA)
library(ggplot2)
library(midasr)

data("data_binom_hyperbolic")
str(data_binom_hyperbolic)

y_ts <- ts(data_binom_hyperbolic$y, frequency = 1)
x_ts <- ts(data_binom_hyperbolic$x, frequency = 1)

png("paper/figures/Example_binomial_hyperbolic_covariatedata.png", width=17, height=10, units = 'cm', res = 300)
autoplot(x_ts, color = "brown") + theme_bw() + geom_point(color = "red", size = .01) +
  theme(axis.text.x=element_text(size=16),
        axis.title.x=element_text(size=16,face="bold"),
        axis.text.y=element_text(size=16),
        axis.title.y=element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=20, face = "plain"),
        legend.title=element_blank())
dev.off()
png("paper/figures/Example_binomial_hyperbolic_response.png", width=17, height=10, units = 'cm', res = 300)
autoplot(y_ts, color = "black") + theme_bw() + geom_point(color = "gray", size = .01) +
  theme(axis.text.x=element_text(size=16),
        axis.title.x=element_text(size=16,face="bold"),
        axis.text.y=element_text(size=16),
        axis.title.y=element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=20, face = "plain"),
        legend.title=element_blank())
dev.off()


m <- 7
lag_k <- 13
x <- data_binom_hyperbolic$x
y <- data_binom_hyperbolic$y
Ntrials <- data_binom_hyperbolic$Ntrials
y[(length(y)-15):length(y)] <- NA # remove last 16 time points as forecast test set

#### Model fitting ####

Midas_objects <- fit_Minla(xdata = x,
                           ydata = y,
                           constraint = "hyperbolic",
                           K = 0:lag_k,
                           m = m,
                           lagY = 0,
                           Ntrials = Ntrials,
                           family = "binomial")

res = inla(y ~ 1 + f(idx, model = Midas_objects$rgen, n = nrow(Midas_objects$data)) +
             f(trend, model = "linear"),
           data = data.frame(y = Midas_objects$data$y,
                             trend = 1:nrow(Midas_objects$data),
                             idx = 1:nrow(Midas_objects$data),
                             Ntrials = Midas_objects$Ntrials),
           verbose = TRUE,
           family = "binomial",
           control.family = list(link = "logit"),
           Ntrials = Midas_objects$Ntrials,
           control.compute=list(config = TRUE))
