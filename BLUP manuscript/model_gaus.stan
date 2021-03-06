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
  vector[9] betas; //fitness regression coefficients
  
  //random effects
  real<lower=0> sigma_0; //w dispersion (sigma for Gaussian)
  vector<lower=0>[3] sd_zp; //RN parameter sds
  matrix[I,3] std_dev; //individual-level RN deviations
  cholesky_factor_corr[3] R_chol; //RN parameter correlations
}
transformed parameters {
  matrix[I,3] zp; //individual phenotypic RN parameter values
  zp =  std_dev * diag_pre_multiply(sd_zp, R_chol)' ;
}
model{
  //separate RN parameters
  vector[I] zp_mu = col(zp,1); //personality
  vector[I] zp_beta = col(zp,2); //plasticity
  vector[I] zp_theta = col(zp,3); //predictability
  
  //initialize vectors for response models
  vector[N] z_mu; //linear predictor of behavior expectation
  vector[N] z_sigma; //linear predictor of behavior dispersion
  vector[I] w_eta; //linear predictor of fitness expectation

  //behavioral RN response model
  z_mu = mu_0z + zp_mu[ind] + (beta_1z + zp_beta[ind]).*x ;
  z_sigma = exp(theta_0z + zp_theta[ind]) ;
  z ~ normal(z_mu, z_sigma);
    
  //fitness response model
  w_eta = mu_0 + betas[1]*zp_mu + betas[2]*zp_beta + betas[3]*zp_theta +
                 betas[4]*(zp_mu .*zp_mu) + betas[5]*(zp_beta .*zp_beta) + betas[6]*(zp_theta .*zp_theta) +   
                 betas[7]*(zp_mu .*zp_beta) + betas[8]*(zp_mu .*zp_theta) + betas[9]*(zp_beta .*zp_theta) ;
  w ~ normal(w_eta, sigma_0);

  //model priors
  
  //fixed effects
  mu_0z ~ normal(0,1);
  beta_1z ~ normal(0,1);
  theta_0z ~ normal(0,1);
  mu_0 ~ normal(1,1);
  betas ~ normal(0,1); 
  
  //random effects
  sd_zp ~ exponential(1);
  R_chol ~ lkj_corr_cholesky(2);
  to_vector(std_dev) ~ std_normal();
  sigma_0 ~ exponential(1);
}

generated quantities{
matrix[3,3] R = R_chol * R_chol'; //RN correlation matrix
matrix[3,3] S = diag_matrix(sd_zp); //RN correlation matrix
matrix[3,3] P = S*R*S; //RN covariance
vector<lower=0>[3] V_P = sd_zp .* sd_zp; //RN variances
}
