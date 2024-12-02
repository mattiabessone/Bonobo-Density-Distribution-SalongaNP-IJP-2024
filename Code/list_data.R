#BRING DATA IN#######
o_dat<-read.csv("Data/occurrence.csv")

# Data dimensions
Ro<-length(o_dat[,1]) # number of cells
T<-2L                 # number of periods
J<-2L                 # number of methods used (transect and recces)
n_sect1<-5L           # number of subsectors surveyed in 2002-2008
n_sect<-9L            # number fo subsectors surveyed in 2012-2018

# Data
o<-array(data=c(as.matrix(cbind(o_dat$o_lt1,o_dat$o_re1)),
                as.matrix(cbind(o_dat$o_lt2,o_dat$o_re2))),
         dim=c(Ro,J,T))
# Continuous covariates
F<-as.vector(standardize(o_dat$for.))#forest
C<-as.matrix(cbind(standardize(o_dat$cit_1),standardize(o_dat$cit_2))) #distance to cities
V<-as.matrix(cbind(standardize(o_dat$vill_1),standardize(o_dat$vill_2)))#distance to villages
R<-as.matrix(cbind(standardize(o_dat$riv_1),standardize(o_dat$riv_2)))#distance to rivers

# Discrete cavariates
sect_o<-o_dat$sect #sectors for occupancy array
pp_o<-as.matrix(cbind(as.integer(o_dat$pp_1),as.integer(o_dat$pp_2))) #patrol post within 15km (yes/no) for occupancy array

# Effort 
Lo<-array(data=c(as.matrix(cbind(standardize(o_dat$e_lt1),standardize(o_dat$e_re1))),
                 as.matrix(cbind(standardize(o_dat$e_lt2),standardize(o_dat$e_re2)))),
          dim=c(Ro,J,T))
####Observed bonobo occurrence - 42km2 grid cells####
o_obs<-read.csv("Data/observed_O.csv")
O_obs<-array(data=c(as.matrix(cbind(o_obs$observed_lt1,o_obs$observed_re1)),
                as.matrix(cbind(o_obs$observed_lt2,o_obs$observed_re2))),
         dim=c(1069L,J,T))

#LIST AND FIT###########
data<-list(Ro=Ro,J=J,T=T,n_sect=n_sect,n_sect1=n_sect1,R_pred=1069L,
              o=o,L=Lo,O_observed=O_obs,
              F=F,C=C,V=V,R=R,sect_o=sect_o,K=pp_o+1)
str(data)
