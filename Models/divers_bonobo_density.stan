// MODEL EVALUATING EFFECT OF COVARIATES ON NEST COUNT AND CAMERA-TRAP DENSITY ESTIMATES

functions{
  
  vector merge_missing(int[] miss_indexes, vector x_obs, vector x_miss){
    int N = dims(x_obs)[1];
    int N_miss = dims(x_miss)[1];
    vector[N] merged;
    merged = x_obs;
    for (i in 1:N_miss)
      merged[miss_indexes[i]] = x_miss[i];
    return merged;
  }
}
// DEFINE DATA AND COVARIATES
data {
 //Data dimensions
  int<lower=0> Rd;                         // Number of transects for density estiamtion
  int<lower=0> Jd;                         // Number of methods used on transects = 2
  int<lower=0> n_sect;                     // Number of sectors = 4

  //Datasets containing presence/absence and density data 
  real<lower=-10> d [Rd,Jd];                // density matrix for estimating number of bonobos
  int<lower=-1,upper=1> z [Rd,Jd];         // 0-1 matrix for estimating process producing zeroes on transects
   
////Continuous covariates
  //Density covariates
  vector [Rd] F2;                          // forest cover
  vector [Rd] C2;                          // distance to cities
  vector [Rd] V2;                          // distance to villages
  vector [Rd] R2;                          // distance to rivers
  vector [Rd] H2;                          // human encounter rate per 100m
  vector [Rd] T;                           // proportion bonobo feeding trees
  vector [Rd] M;                           // proportion Marantaceae
  vector [Rd] B;                           // black mangabey density
  
////Discrete covariates
  int sect_d[Rd];                          // sectors for density array
  int K_d[Rd];                             // patrol post within 15km (yes/no) for density array
////Imputation
  int<lower=0> N_miss;
  int<lower=0> miss_H [N_miss];
  int<lower=0> miss_T [N_miss];
  int<lower=0> miss_M [N_miss];
  int<lower=0> miss_B [N_miss];

}

// DECLARE PARAMETERS TO BE ESTIMATED  
parameters {
  real alpha;
  matrix [n_sect,Jd]  alpha2;                  // sector specific intercept for mean density mu
 
//Slopes density model
  vector [Jd] delta1;                          // forest cover by method
  vector [Jd] delta2;                          // distance to cities by method
  vector [Jd] delta3;                          // distance to vilalges by method
  vector [Jd] delta4;                          // distance to rivers by method
  vector [Jd] delta5;                          // human encounter rate per 100m, by method
  vector [Jd] delta6;                          // proportion of bonobo feeding trees by method
  vector [Jd] delta7;                          // proportion of Marantaceae by method
  vector [Jd] delta8;                          // black mangabey density
  vector [2]  delta9 [Jd];                     // proximity to a PP, by method
  vector [Jd] delta10;
  vector<lower=0,upper=1>[Jd] phi;         // probability (by method) to find a bonobo nest / image on transect
  vector<lower=0>  [Jd] theta;             // overdispersion parameter (by method) of mean density
  
//Imputation
  vector [N_miss] imputed_H;
  vector [N_miss] imputed_T;
  vector [N_miss] imputed_M;
  vector [N_miss] imputed_B;
  
 }

transformed parameters{
  vector [Rd] merged_H;
  vector [Rd] merged_T;
  vector [Rd] merged_M;
  vector [Rd] merged_B;
  
  merged_H = merge_missing(miss_H,to_vector(H2),imputed_H);
  merged_T = merge_missing(miss_T,to_vector(T),imputed_T);
  merged_M = merge_missing(miss_M,to_vector(M),imputed_M);
  merged_B = merge_missing(miss_B,to_vector(B),imputed_B);
}
// DECLARE THE MODEL
model { 
// define perameters
  real mu;                                 // mean density
  
// define priors
   alpha~ normal(0,5);

   for (i in 1:n_sect){
     
     alpha2[i] ~ normal(0,5);

   }
   
   delta1 ~ normal(0,0.5);
   delta2 ~ normal(0,0.5);
   delta3 ~ normal(0,0.5);
   delta4 ~ normal(0,0.5);
   delta5 ~ normal(0,0.5);
   delta6 ~ normal(0,0.5);
   delta7 ~ normal(0,0.5);
   delta8 ~ normal(0,0.5);
 
   
   for (j in 1:Jd){
     
     delta9[j] ~ normal(0,5);

   }
  
  delta10 ~ normal(0,5);
     
   phi ~ beta(2,2);
   theta ~ gamma(0.3,0.3);
   
   imputed_H ~ std_normal();
   imputed_T ~ std_normal();
   imputed_M ~ std_normal();
   imputed_B ~ std_normal();
   
//////////COUNT MODEL

///// Model density using line transects (nest) and camera-traps (bonobos)
    for (r in 1:Rd){
       for (j in 1:Jd){
        if(z[r,j] >-1){                            // if we surveyed the transect with a particular method

          z[r,j] ~ bernoulli(phi[j]);
            
          }
        }
     }

    for (r in 1:Rd){
      for (j in 1:Jd){
        if(d[r,j]>0){
          
          mu = exp(alpha + alpha2[sect_d[r],j] + delta1[j] * F2[r] + delta2[j] * C2[r] + delta3[j] * V2[r] + delta4[j] * R2[r] + delta5[j] * merged_H[r] + delta6[j] * merged_T[r] + delta7[j] * merged_M[r] + delta8[j] * merged_B[r] + delta9[K_d[r],j] +delta10[j]);
          d[r,j] ~ gamma(mu * theta[j], theta[j]);
          
        }
      }
    }
}

