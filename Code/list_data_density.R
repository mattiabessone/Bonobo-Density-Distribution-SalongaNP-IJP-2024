#BRING DATA IN#######
d_dat<-read.csv("Data/drivers.csv")

#Define data dimensions
Rd<-length(d_dat$Transect)#number of transects
Jd<-2L#number of methods used = 2
n_sect<-9L#number of sectors

#Datasets containing presence/absence and density data 
#occupancy covariates
#density

d_dat$for_2<-as.vector(standardize(d_dat$for_2))#forest
d_dat$cit_2<-as.vector(standardize(d_dat$cit_2)) #distance to cities
d_dat$vill_2<-as.vector(standardize(d_dat$vill_2))#distance to villages
d_dat$riv_2<-as.vector(standardize(d_dat$riv_2))#distance to rivers
d_dat$hum_100m<-as.vector(standardize(d_dat$hum_100m))#human encounter rate per 100m
d_dat$bft<-as.vector(standardize(d_dat$bft))#proportion bonobo feeding trees
d_dat$mar<-as.vector(standardize(d_dat$mar))#proportion marantacea
d_dat$bm.D<-as.vector(standardize(d_dat$bm.D))#back mangabey group size

miss_H<-which(is.na(d_dat$hum_100m))
miss_T<-which(is.na(d_dat$bft))
miss_M<-which(is.na(d_dat$mar))
miss_B<-which(is.na(d_dat$bm.D))

d_dat[is.na(d_dat)] <- -10

d<-as.matrix(cbind(d_dat$D_LT,d_dat$D_CT)) #density matrix for estimating number of bonobos
z<-as.matrix(cbind(d_dat$obs_LT,d_dat$obs_CT))
#Discrete covariates
sect_d<-d_dat$sect_2#sectors for density array
pp_d<-as.integer(d_dat$PP_15_2)#patrol post within 15km (yes/no)  for density array

#########LIST DATA######

data_D<-list(Rd=Rd,Jd=Jd,n_sect=n_sect,#data dimnsions
               d=d,z=z,#main datasets
               F2=d_dat$for_2,C2=d_dat$cit_2,V2=d_dat$vill_2,R2=d_dat$riv_2,H2=d_dat$hum_100m,M=d_dat$mar,T=d_dat$bft,B=d_dat$bm.D,#variables for density model
               N_miss=length(miss_B),
               miss_H=miss_H,miss_T=miss_T,miss_M=miss_M,miss_B=miss_B,
               sect_d=sect_d,K_d=pp_d+1)#discrete variables
               
str(data_D)

####### SNPs only ##############
#BRING DATA IN#######
library(rethinking)
d_dat<-read.csv("Data/drivers.csv")

#Define data dimensions
d_dat$for_2<-as.vector(standardize(d_dat$for_2))#forest
d_dat$cit_2<-as.vector(standardize(d_dat$cit_2)) #distance to cities
d_dat$vill_2<-as.vector(standardize(d_dat$vill_2))#distance to villages
d_dat$riv_2<-as.vector(standardize(d_dat$riv_2))#distance to rivers
d_dat$hum_100m<-as.vector(standardize(d_dat$hum_100m))#human encounter rate per 100m
d_dat$bft<-as.vector(standardize(d_dat$bft))#proportion bonobo feeding trees
d_dat$mar<-as.vector(standardize(d_dat$mar))#proportion marantacea
d_dat$bm.D<-as.vector(standardize(d_dat$bm.D))#back mangabey group size

d_dat$mar[is.na(d_dat$mar)] <- 100

d_dat<-subset(d_dat,d_dat$mar<100)

miss_H<-which(is.na(d_dat$hum_100m))
miss_T<-which(is.na(d_dat$bft))
miss_M<-which(is.na(d_dat$mar))
miss_B<-which(is.na(d_dat$bm.D))
d<-as.matrix(cbind(d_dat$D_LT,d_dat$D_CT)) #density matrix for estimating number of bonobos

#Discrete covariates
sect_d<-d_dat$sect_2#sectors for density array
sect_d[sect_d==2]<-1
sect_d[sect_d==3]<-2
sect_d[sect_d==7]<-3
sect_d[sect_d==8]<-4
sect_d[sect_d==5]<-1

pp_d<-as.integer(d_dat$PP_15_2)#patrol post within 15km (yes/no)  for density array
d_dat[is.na(d_dat)] <- -10
Rd<-length(d_dat$Transect)#number of transects
Jd<-2L#number of methods used = 2
n_sect<-4L#number of sectors

#########LIST DATA######

data_D_SNPs<-list(Rd=Rd,Jd=Jd,n_sect=n_sect,#data dimnsions
             d=d,#main datasets
             F2=d_dat$for_2,C2=d_dat$cit_2,V2=d_dat$vill_2,R2=d_dat$riv_2,H2=d_dat$hum_100m,M=d_dat$mar,T=d_dat$bft,B=d_dat$bm.D,#variables for density model
             N_miss=length(miss_B),
             miss_H=miss_H,miss_T=miss_T,miss_M=miss_M,miss_B=miss_B,
             sect_d=sect_d,K_d=pp_d+1)#discrete variables
