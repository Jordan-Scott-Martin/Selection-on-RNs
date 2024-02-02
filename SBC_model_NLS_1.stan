
data {
  int<lower=1> I; //total individuals
  int<lower=1> npar; //number of RN parameters 
  int<lower=1> N_z; //total number of observations (Z)
  int<lower=1> N_w; //total number of observations (W)
  array[N_z] int<lower=1> ind_z; //index of individual observations (Z)
  array[N_w] int<lower=1> ind_w; //index of individual observations (W)
  vector[N_z] x; //environmental covariate
  vector[N_z] z; //phenotype
  vector[N_w] W; //fitness
}

transformed data {
int nc = (npar * (npar-1))/2;
}

parameters {
  //fixed population effects
  real mu_0z; //z population intercept
  real beta_z; //z population slope
  real theta_0z; //z population dispersion
  
  real mu_0w; //w population intercept
  
  vector[npar] b; //direct selection 
  vector[npar] q; //stabilizing/disruptive selection 
  vector[nc] qc; //correlational selection
  
  //random effects
  vector<lower=0>[npar] sd_I_z; //RN parameter sds
  cholesky_factor_corr[npar] LP; //RN parameter correlations
  matrix[I,npar] std_dev_z; //individual-level RN deviations
  
  vector[I] std_dev_w; //individual-level selection deviations
  real<lower=0> sd_I_w; //unexplained selection sd
  
  real<lower=0> theta_0w; //w dispersion (sigma for Gaussian)
}

transformed parameters {
matrix[I,npar] zp =  std_dev_z * diag_pre_multiply(sd_I_z, LP)'; //individual phenotypic RN parameter values
vector[I] W_0 = std_dev_w * sd_I_w; //individual unexplained selection effects
}

model{
  
  //separate RN parameters
  vector[I] zp_mu = col(zp,1); //intercepts
  vector[I] zp_beta = col(zp,2); //slopes
  vector[I] zp_theta = col(zp,3); //residuals
  
  //initialize vectors for response models
  vector[N_z] z_eta; //linear predictor (trait expectation)
  vector[N_z] z_theta; //linear predictor (trait dispersion)
  vector[N_w] w_eta; //linear predictor (fitness expectation)

  //behavioral RN response model
  z_eta = mu_0z + zp_mu[ind_z] + (beta_z + zp_beta[ind_z]) .* x;
  z_theta = sqrt(exp(theta_0z + zp_theta[ind_z]));
  z ~ normal(z_eta, z_theta);
    
  //fitness response model
  w_eta = mu_0w + W_0[ind_w] +
          
          b[1] * zp_mu[ind_w] + 
          b[2] * zp_beta[ind_w] + 
          b[3] * zp_theta[ind_w] +
          
          q[1] * (zp_mu[ind_w] .* zp_mu[ind_w]) + 
          q[2] * (zp_beta[ind_w] .* zp_beta[ind_w]) +
          q[3] * (zp_theta[ind_w] .* zp_theta[ind_w]) +
                    
          qc[1] * (zp_mu[ind_w] .* zp_beta[ind_w]) + 
          qc[2] * (zp_mu[ind_w] .* zp_theta[ind_w]) + 
          qc[3] * (zp_beta[ind_w] .* zp_theta[ind_w]) ;
  
  W ~ normal(w_eta, theta_0w);

  //model priors
  
  //fixed effects
  mu_0z ~ normal(0,1);
  beta_z ~ normal(0,1);
  theta_0z ~ normal(0,1);
  mu_0w ~ normal(0,1);
  b ~ normal(0,1);
  q ~ normal(0,1);
  qc ~ normal(0,1);
  
  //random effects
  sd_I_z ~ exponential(2);
  LP ~ lkj_corr_cholesky(2);
  to_vector(std_dev_z) ~ std_normal();
  
  sd_I_w ~ exponential(2);
  to_vector(std_dev_w) ~ std_normal();
  theta_0w ~ exponential(2);
  }
