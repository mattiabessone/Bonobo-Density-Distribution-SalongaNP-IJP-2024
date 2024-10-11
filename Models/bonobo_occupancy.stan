/////////// MODEL INTEGRATING RECCES AND LINE TRANSECTS - POPUALTION TREND IN SNP
////////// DEFINE DATA AND COVARIATES
data {
 //Data dimensions
  int<lower=0> Ro;                               // number of cells for occupancy estimation
  int<lower=0> J;                                // number of methods used = 2
  int<lower=0> T;                                // number of periods investigated = 2
  int<lower=0> n_sect;                           // number of sectors = 9
  int<lower=0> n_sect1;                          // number of sectors surveyed in P1, occupancy = 5
  int<lower=0> R_pred;                           // number of cells for prediction
  //Datasets containing presence/absence and density data 
  int<lower=-5,upper=1> o [Ro,J,T];              // presence absence matrix (with NAs) for final estiamtes of occurrence probability
  int<lower=-5,upper=1> O_observed [R_pred,J,T];              // presence absence matrix (with NAs) for final estiamtes of occurrence probability

  ///// Continuous covariates
  // Occupancy
  vector [Ro] F;                                 // forest cover
  matrix [Ro,T] C;                               // distance to cities
  matrix [Ro,T] V;                               // distance to villages
  matrix [Ro,T] R;                               // distance to rivers
  
/////Discrete covariates
  int<lower=1> sect_o[Ro];                       // sectors for occupancy array
  int<lower=1> K[Ro,T];                          // patrol post within 15km (yes/no) for occupancy array
  
/////Effort
  vector [T]  L [Ro,J];                         // path length, i.e. survey effort
  
  }
  
// DECLARE PARAMETERS TO BE ESTIMATED  
parameters {
  vector [J] alpha [T];                          // method-specific intercept for detection probability p
  vector<lower=0> [J] eta [T];                   // method-specific control parameter for detection probability p
  vector[n_sect1] alpha1_1;                      // sector specific intercept for occurrence probability psi in P1
  real<lower=0,upper=1> alpha1_1_unsurveyed;     // intercept for unsurveyd sectors in P1 
  vector[n_sect] alpha1_2;                       // sector specific intercept for occurrence probability psi in P2
  
  /////// Slopes occupancy model
  // P1
  real beta1_1;                                  // forest cover
  vector [2] beta2_1 [n_sect1];                  // distance to cities by sector and proximity to PP (yes or no)
  vector [2] beta3_1 [n_sect1];                  // distance to villages by sector and proximity to PP (yes or no)
  vector [2] beta4_1 [n_sect1];                  // distance to rivers by sector and proximity to PP (yes or no)
  vector [2] beta5_1;                            // proximity to PP
  // P2
  real beta1_2;                                  // forest cover
  vector [2] beta2_2 [n_sect];                   // distance to cities by sector and proximity to PP (yes or no)
  vector [2] beta3_2 [n_sect];                   // distance to villages by sector and proximity to PP (yes or no)
  vector [2] beta4_2 [n_sect];                   // distance to rivers by sector and proximity to PP (yes or no)
  vector [2] beta5_2;                            // proximity to PP
 
 }
 

model {
///// define perameters   
  real p1;                                        // detection probability in P1
  real p2;                                        // detection probability in P2
  real psi1;
  real psi2;

/////priors
  
  for(j in 1:J){
    alpha[j] ~ normal(0,1.4);
    eta[j] ~ normal(0,0.5);
  }

   /////intercepts
  alpha1_1 ~ normal(0,1.4);                      
  alpha1_1_unsurveyed ~ normal(0,1.4);
  alpha1_2 ~ normal(0,1.4);
  
  /////slopes
  //occupancy P1
  beta1_1 ~ normal(0,0.5);
  
  for (i in 1:n_sect1){
   beta2_1[i] ~ normal(0,0.5);
   beta3_1[i] ~ normal(0,0.5);
   beta4_1[i] ~ normal(0,0.5);
  }
  //occupancy P2
  beta1_2~ normal(0,0.5);
  
  for (i in 1:n_sect){
   beta2_2[i] ~ normal(0,0.5);
   beta3_2[i] ~ normal(0,0.5);
   beta4_2[i] ~ normal(0,0.5);
  }
  
  beta5_1 ~ normal(0,1.4);
  beta5_2 ~ normal(0,1.4);
  
for (r in 1:Ro){
      
      if(sect_o[r] < 6){                           // if sector < 6 (i.e. surveyd in both periods), then we use sector specific intercepts
      
        psi1 = inv_logit(alpha1_1[sect_o[r]] + beta1_1 * F[r] + beta2_1[sect_o[r],K[r,1]] * C[r,1] + beta3_1[sect_o[r],K[r,1]] * V[r,1] + beta4_1[sect_o[r],K[r,1]] * R[r,1] + beta5_1[K[r,1]]);
        psi2 = inv_logit(alpha1_2[sect_o[r]] + beta1_2 * F[r] + beta2_2[sect_o[r],K[r,2]] * C[r,2] + beta3_2[sect_o[r],K[r,2]] * V[r,2] + beta4_2[sect_o[r],K[r,2]] * R[r,2] + beta5_2[K[r,2]]);  
      
      }else{                                       // if sector > 6 (i.e. surveyd in in P2 only), then we use sector specific intercept in P2, but a single intercept in P1
      
        psi1 = inv_logit(alpha1_1_unsurveyed);
        psi2 = inv_logit(alpha1_2[sect_o[r]] + beta1_2 * F[r] + beta2_2[sect_o[r],K[r,2]] * C[r,2] + beta3_2[sect_o[r],K[r,2]] * V[r,2] + beta4_2[sect_o[r],K[r,2]] * R[r,2] + beta5_2[K[r,2]]);
      
      }
      for (j in 1:J){
        // declare linear model on  p
        p1 = inv_logit(alpha[j,1] + eta[j,1] * L[r,j,1]); 
        p2 = inv_logit(alpha[j,2] + eta[j,2] * L[r,j,2]);
          
          for (t in 1:T){
            if (t == 1){
                   if(o[r,j,t] == 1){ 
          
                     target += log(psi1) + log(p1);                       // then 1 is observed if cell is occupied (psi) and detected (p1)
          
                     }else{
                       if (o[r,j,t] == 0){                                // site is not apparently occupied, two paths:
          
                       target += log_sum_exp(log(psi1) + log(1-p1),       // either bonobos are there and not detected
                                 log(1-psi1));                            // no bonobos
                       }
                     }
              }else{
                if(o[r,j,t] == 1){
          
                  target += log(psi2) + log(p2);                          // then 1 is observed if cell is occupied (psi) and detected (p1)
          
                }else{
                  if (o[r,j,t] == 0){                                     // site is not apparently occupied, two paths:
          
                     target += log_sum_exp(log(psi2) + log(1-p2),         // either bonobos are there and not detected
                               log(1-psi2))           ;                   // no bonobos
          
                  }
                }
              }
          }
      }
    }
}

///////// MAKE PREDICTIONS
generated quantities{
///// Define generated quantities
  real pred_p [Ro,J,T];
  real pred_psi [Ro,T];                                                      // predicted occupancy by method and period  
  int O_pred [R_pred,J,T];
  int O_out [R_pred,T];                                              // predicted occupancy
  real occupied[T];                                                       // proportion of occupied cells for trend analysis, by period
  real occupied_P2;
  real trend_O;                                                           // trend occupancy in sectors surveyed twice
  
///// Generate true occupancy estimate conditional on method and period
  for ( r in 1:Ro){                                                  // predict detection probability to all cells (1km)
    for (j in 1:J){
    for (t in 1:T){
        
        pred_p[r,j,t] = inv_logit(alpha[j,t] + eta[j,t] * L[r,j,t]); 
        
    }
    }
  }
  
  for (r in 1:Ro){
    if(sect_o[r] < 6){                           // if sector < 6 (i.e. surveyd in both periods), then we use sector specific intercepts
    
    pred_psi[r,1] = inv_logit(alpha1_1[sect_o[r]] + beta1_1 * F[r] + beta2_1[sect_o[r],K[r,1]] * C[r,1] + beta3_1[sect_o[r],K[r,1]] * V[r,1] + beta4_1[sect_o[r],K[r,1]] * R[r,1] + beta5_1[K[r,1]]);
    pred_psi[r,2] = inv_logit(alpha1_2[sect_o[r]] + beta1_2 * F[r] + beta2_2[sect_o[r],K[r,2]] * C[r,2] + beta3_2[sect_o[r],K[r,2]] * V[r,2] + beta4_2[sect_o[r],K[r,2]] * R[r,2] + beta5_2[K[r,2]]);
    
    }else{                                       // if sector > 6 (i.e. surveyd in in P2 only), then we use sector specific intercept in P2, but a single intercept in P1
    
    pred_psi[r,1] = inv_logit(alpha1_1_unsurveyed);
    pred_psi[r,2] = inv_logit(alpha1_2[sect_o[r]] + beta1_2 * F[r] + beta2_2[sect_o[r],K[r,2]] * C[r,2] + beta3_2[sect_o[r],K[r,2]] * V[r,2] + beta4_2[sect_o[r],K[r,2]] * R[r,2] + beta5_2[K[r,2]]);
    
    }
  }
      
for (r in 1:R_pred){                                               //predict occupancy to prediction grid (42km2) by averaging occurence and detection probability estimated above (1 km2)
      for (j in 1:J){
        for (t in 1:T){
          {
          real psi_pred;
          real p_pred;
          psi_pred = mean(pred_psi[(42*r-41):(42*r),t]);
          p_pred = mean(pred_p[(42*r-41):(42*r),j,t]);
          
        // If we detected an animal, then we set the prediction to 1. i.e., we assume no false positives 
        if(O_observed[r,j,t]==1){
          O_pred[r,j,t] = 1;
        }
        
        // If we did not detect an animal, then we set the prediction to the prob of observeing a zero, even if occupancy is truly 1.
        // i.e., we assume that false negatives do occur
        if(O_observed[r,j,t]==0){  
          O_pred[r,j,t] = bernoulli_rng((psi_pred * (1-p_pred)) / (psi_pred * (1-p_pred) + (1-psi_pred)));
        } 

       // If no methods were used to improve predictive accuracy, then occupancy estimate alone is used        
       if(O_observed[r,j,t]<0){ 
          O_pred[r,j,t] = bernoulli_rng(psi_pred);
          }
        }
      }
    }
  }
  
  // Apply threshold assertion that a 1 with any method implies occupation 
   for (r in 1:R_pred){
     for (t in 1:T){
       if(sum(O_pred[r,,t])>0){
         
         O_out[r,t] = 1;
         
         } else{
           
           O_out[r,t] = 0;
         }
     }
   }

///// Final sample estimates

  for (t in 1:T){
   
    occupied[t] = sum(O_out[1:371,t]);
   
  }
 
  occupied_P2 = sum(O_out[1:R_pred,2]);                            //number of occupied cells in P2
  trend_O = occupied[2]/occupied[1];                               //bonobo occupancy trend in sub-sectors surveyed twice for scenario i
    
    
}


