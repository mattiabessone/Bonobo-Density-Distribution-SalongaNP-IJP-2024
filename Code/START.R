#always make sure you are in the R_library in "C"
library(rstan)
library(rethinking)
rstan_options(auto_write = TRUE)
rstan_options(javascript = FALSE)
options(mc.cores = parallel::detectCores())
library(maptools)
library(sf)
library(INLA)
library(Matrix)
library(spdep)
library(rgdal)

