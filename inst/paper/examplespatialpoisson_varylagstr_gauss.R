

rm(list = ls())

library(midasINLA)
library(ggplot2)
library(midasr)
library(INLA)
library(ggplot2)
library(tidyr)

data("data_spatialpoisson_varylagstr_gauss", package = "midasINLA")
str(data_spatialpoisson_varylagstr_gauss)

for_plot <- data_spatialpoisson_varylagstr_gauss$data_y
for_plot <- for_plot[which(for_plot$loc %in% c(1:6)),]
for_plot$loc <- as.factor(for_plot$loc)

ggplot(for_plot, aes(y = y, x= Time, group = loc, color = loc)) +
  geom_line() +
  geom_point() +
  theme_minimal()


data_y <- data_spatialpoisson_varylagstr_gauss[["data_y"]]
data_y[which(data_y[["Time"]] %in% 55:60),"y"] <- NA

Midas_objects <- fit_Minla_spatial_varylagstr(xdata = data_spatialpoisson_varylagstr_gauss[["data_x"]][["x"]],
                                   ydata = data_y[["y"]],
                                   loc_x = data_spatialpoisson_varylagstr_gauss[["data_x"]][["loc"]],
                                   loc_y = data_y[["loc"]],
                                   constraint = "hyperbolic",
                                   K = rep(29,5),
                                   m = 30,
                                   family = "poisson")


res = inla(y ~ 1 +
             f(idx, model = Midas_objects[["rgen"]],
               n = nrow(Midas_objects[["data"]])) +
             f(iid_loc, model = "iid"),
           data = data.frame(y = Midas_objects[["data"]][["y"]],
                             idx = 1:nrow(Midas_objects[["data"]]),
                             iid_loc = Midas_objects[["idx_loc"]]),
           verbose = TRUE,
           family = "poisson",
           #control.inla = list(strategy = "eb")
           control.compute=list(config = TRUE))


summary(res)
