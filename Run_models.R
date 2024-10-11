# 1) SET WORKING DIRECTORY AND LOAD LIBRARIES #####
setwd("~/R/bonobo_occupancy")
source("Code/START.R")

# 2) RUN occupancy trend analysis in SNP, entire Park #####
source("Code/list_data.R") # List data
detach(package:rethinking,unload=TRUE)
fit_O<-stan("Models/bonobo_occupancy.stan", data = data, chains = 1, cores = 1, warmup = 1000, iter = 2000, save_warmup = FALSE,control=list(adapt_delta=0.8)) # Run model
saveRDS(fit_O,"Results/results_O.rds")

# 3) CHECK results
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

# Make predictions
predicted_psi<-matrix(nrow=1069,ncol=T+1)
predicted_SD_psi<-matrix(nrow=1069,ncol=T+1)

for (r in 1:Ro){
  for (t in 1:T){
    predicted_psi[r,t] <- mean(post_O$O_out[,r,t])
    predicted_SD_psi[r,t] <- sd(post_O$O_out[,r,t])
  }
} 

for (r in 1:1069){                                               
  for (t in 1:T){
    predicted_psi[r,t] = mean(post_O$pred_psi[,(42*r-41):(42*r),t])
    predicted_SD_psi[r,t] = sd(post_O$pred_psi[,(42*r-41):(42*r),t])
  }
}
    
pred_data <- read.csv("Data/prediction_M2.csv")
predicted_psi[,3]<-pred_data$cell_ID
predicted_SD_psi[,3]<-pred_data$cell_ID

write.csv(predicted_psi,"Results/predicted_psi.csv")
write.csv(predicted_SD_psi,"Results/predicted_SD_psi.csv")

# 4) assess variables effect on bonobo nest density #####
source("Code/START.R")
source("Code/list_data_density.R") # List data
detach(package:rethinking,unload=TRUE)
fit_D<-stan("Models/divers_bonobo_density.stan", data = data_D, chains = 1, cores = 1, warmup = 1000, iter = 2000,control=list(max_treedepth=10),init="0") # Run model
#saveRDS(fit_D,"Results/results_drivers_D.rds")
#fit_D_SNPs<-stan("Models/divers_bonobo_density.stan", data = data_D_SNPs, chains = 1, cores = 1, warmup = 1000, iter = 2000,control=list(max_treedepth=10),init="0") # Run model
#Check results
library(rethinking)
library(bayestestR)
precis(fit_D,pars=c("delta1","delta2","delta3","delta4","delta5","delta6","delta7","delta8","delta9"),depth=3,prob=0.95)
p_direction(fit_D,parameters = c("delta1","delta2","delta3","delta4","delta5","delta6","delta7","delta8"))
dens(extract(fit_D,pars="delta5")$delta5[,2])

