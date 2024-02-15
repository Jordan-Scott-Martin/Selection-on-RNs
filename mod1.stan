
data {
  int<lower=1> J; //total individuals
  int<lower=1> N_z; //total number of pheontype meausres (Z)
  int<lower=1> N_W; //total number of fitness measures (W)
  array[N_z] int<lower=1> ind_z; //index of individual measurements (z)
  array[N_W] int<lower=1> ind_W; //index of individual measurements (W)
  vector[N_z] x; //environmental covariate
  vector[N_z] z; //phenotype
  vector[N_W] W; //fitness
}

parameters {
  //fixed population effects
  real mu_0; //z population intercept
  real beta_x; //z population slope
  real sigma_0; //z population dispersion
  
  real W_0; //W population intercept
  vector[3] b; //direct selection 
  vector[3] q; //stabilizing/disruptive selection 
  vector[3] qc; //correlational selection 

  //random effects for z
  vector<lower=0>[3] sd_RN; //RN parameter sds
  matrix[J,3] std_dev_RN; //individual-level RN deviations
  cholesky_factor_corr[3] R_chol; //RN parameter correlations

  //random effects for W
  real<lower=0> sd_W0; //unexplained selection sd
  vector[J] std_dev_W; //individual-level selection deviations
  real<lower=0> delta; //W dispersion (SD of residuals)
}

transformed parameters {
  vector[J] W_0j = std_dev_W * sd_W0; //scaled random intercepts for fitness
  matrix[J,3] RNj = std_dev_RN * diag_pre_multiply(sd_RN, R_chol)' ;
}

model{
  //separate RN parameters
  vector[J] mu_0j = col(RNj,1); //intercepts
  vector[J] beta_xj = col(RNj,2); //slopes
  vector[J] sigma_0j = col(RNj,3); //residuals
  
  //initialize vectors for response models
  vector[N_z] mu; //linear predictor of phenotype expectation
  vector[N_z] sigma; //linear predictor of phenotype dispersion
  vector[N_W] theta; //linear predictor of fitness expectation
  
  //RN model
  mu = mu_0 + mu_0j[ind_z] + (beta_x + beta_xj[ind_z]) .* x;
  sigma = sqrt(exp(sigma_0 + sigma_0j[ind_z]));
  z ~ normal(mu, sigma);
  
  //fitness model
  theta = W_0 + W_0j[ind_W] +
          
          b[1] * mu_0j[ind_W] + 
          b[2] * beta_xj[ind_W] + 
          b[3] * sigma_0j[ind_W] +
          
          q[1] * (mu_0j[ind_W] .* mu_0j[ind_W]) + 
          q[2] * (beta_xj[ind_W] .* beta_xj[ind_W]) +
          q[3] * (sigma_0j[ind_W] .* sigma_0j[ind_W]) +
                    
          qc[1] * (mu_0j[ind_W] .* beta_xj[ind_W]) + 
          qc[2] * (mu_0j[ind_W] .* sigma_0j[ind_W]) + 
          qc[3] * (beta_xj[ind_W] .* sigma_0j[ind_W]) ;
  W ~ normal(theta, delta);
  
  //model priors
  
  //fixed effects
  mu_0 ~ normal(0,1);
  beta_x ~ normal(0,1);
  sigma_0 ~ normal(0,1);
  W_0 ~ normal(0,1);
  
  b ~ normal(0,1);
  q ~ normal(0,1);
  qc ~ normal(0,1);
  
  //random effects
  sd_RN ~ exponential(2);
  R_chol ~ lkj_corr_cholesky(2);
  to_vector(std_dev_RN) ~ std_normal();
  
  sd_W0 ~ exponential(2);
  std_dev_W ~ std_normal();
  delta ~ exponential(2);
}

generated quantities{
  matrix[3,3] R = R_chol * R_chol'; //RN correlation matrix
  matrix[3,3] S = diag_matrix(sd_RN); //RN SD matrix
  matrix[3,3] P = S*R*S; //RN covariance matrix
  vector<lower=0>[3] V_P = sd_RN .* sd_RN; //RN variances
}


