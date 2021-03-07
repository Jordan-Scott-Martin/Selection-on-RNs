data {
  int<lower=1> I; //total individuals
  int<lower=1> N; //total number of observations
  int<lower=1> ind[N]; //index of individual observations
  vector[N] x; //environmental covariate
  vector[N] z; //behavioral measurements
  vector[I] w; //fitness measurements
}
parameters {
  //fixed population effects
  real mu_0z; //z population intercept
  real beta_1z; //z population slope
  real theta_0z; //z population dispersion
  real mu_0; //w population intercept
  vector[3] betas; //fitness regression coefficients
  
  //random effects (standard deviations)
  real<lower=0> sd_P; //RN parameter sds
  vector[I] std_dev; //individual-level RN deviations
  real<lower=0> theta_0; //w dispersion (sigma for Gaussian)
}
transformed parameters {
  vector[I] zp_mu; //individual phenotypic RN parameter values
  zp_mu =  std_dev * sd_P;
}
model{
  //initialize vectors for response models
  vector[N] z_eta; //linear predictor of behavior expectation
  vector[I] w_eta; //linear predictor of fitness expectation

  //behavioral RN response model
  z_eta = mu_0z + zp_mu[ind];
  z ~ normal(z_eta, theta_0z);
    
  //fitness response model
  w_eta = mu_0 + betas[1]*zp_mu + betas[2]*(zp_mu .*zp_mu);
  w ~ normal(w_eta, theta_0);

  //model priors
  
  //fixed effects
  mu_0z ~ normal(0,1);
  beta_1z ~ normal(0,1);
  mu_0 ~ normal(0,1);
  to_vector(betas) ~ normal(0,1); 
  
  //random effects
  sd_P ~ exponential(1);
  std_dev ~ std_normal();
  theta_0 ~ exponential(1);
  theta_0z ~ exponential(1);
}



