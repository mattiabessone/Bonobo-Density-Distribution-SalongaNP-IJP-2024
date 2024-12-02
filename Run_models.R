# 1) SET WORKING DIRECTORY AND LOAD LIBRARIES #####
setwd("~/GitHub/Bonobo-Density-Distribution-SalongaNP-IJP-2024")
source("Code/START.R")

# 2) RUN occupancy trend analysis in SNP, entire Park #####
source("Code/list_data.R") # List data
detach(package:rethinking,unload=TRUE)
fit_O<-stan("Models/bonobo_occupancy.stan", data = data, chains = 1, cores = 1, warmup = 1000, iter = 2000, save_warmup = FALSE,control=list(adapt_delta=0.8)) # Run model

# CHECK results
library(rethinking)
precis(fit_O,pars=c("alpha","eta","alpha1_1","alpha1_2","beta1_1","beta1_2","beta2_1","beta2_2","beta3_1","beta3_2",
                    "beta4_1","beta4_2","beta5_1","beta5_2"),depth=3,prob=0.95)

dens(extract(fit_O,pars="beta5_1")$beta5_1[,1]-extract(fit_O,pars="beta5_1")$beta5_1[,2])
dens(extract(fit_O,pars="beta5_2")$beta5_2[,1]-extract(fit_O,pars="beta5_2")$beta5_2[,2])

precis(fit_O,pars=c("occupied","trend_O"),depth=3,prob=0.95)

post_O<-extract(fit_O)
library(bayestestR)
p_direction(post_O$trend_O - 1)
p_direction(fit_O,parameters=c("beta1_1","beta1_2","beta2_1","beta2_2","beta3_1","beta3_2",
                    "beta4_1","beta4_2"))

# 2) RUN assessment of factors effecting on bonobo density #####
source("Code/START.R")
source("Code/list_data_density.R") # List data
detach(package:rethinking,unload=TRUE)
fit_D<-stan("Models/divers_bonobo_density.stan", data = data_D, chains = 1, cores = 1, warmup = 1000, iter = 2000,control=list(max_treedepth=10),init="0") # Run model

# CHECK results
library(rethinking)
library(bayestestR)
precis(fit_D,pars=c("delta1","delta2","delta3","delta4","delta5","delta6","delta7","delta8","delta9"),depth=3,prob=0.95)
p_direction(fit_D,parameters = c("delta1","delta2","delta3","delta4","delta5","delta6","delta7","delta8"))
dens(extract(fit_D,pars="delta5")$delta5[,2])

