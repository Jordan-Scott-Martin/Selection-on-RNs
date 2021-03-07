functions {
  vector rep_each(vector x, int K) {
    int N = rows(x);
    vector[N * K] y;
    int pos = 1;
    for (n in 1:N) {
      for (k in 1:K) {
        y[pos] = x[n];
        pos += 1;
      }
    }
    return y;
  }
}

data {
  int<lower=1> I; //total individuals
  int<lower=1> N; //total number of observations
  int<lower=1> ind[N]; //index of individual observations
}

transformed data{
  //fixed population effects
  real mu_0z_sim = normal_rng(0,1); //z population intercept
  real beta_1z_sim = normal_rng(0,1); //z population slope
  real theta_0z_sim = normal_rng(0,1); //z population dispersion
  real mu_0_sim = normal_rng(0,1); //w population intercept
  real betas_sim[9] = normal_rng([0,0,0,0,0,0,0,0,0],[1,1,1,1,1,1,1,1,1]); //fitness regression coefficients
  
  //random effects (standard deviations)
  //real<lower=0>sd_P_sim[3] = cauchy_rng([0,0,0],[1,1,1]); //RN parameter sds
  real sd_P_sim[3] = fabs(cauchy_rng([0,0,0],[1,1,1])); //RN parameter sds
  cholesky_factor_corr[3] LP_sim = lkj_corr_cholesky_rng(3,2); //RN parameter correlations
  matrix[I,3] std_dev_sim = append_col(
                              append_col(
                              to_vector(normal_rng(rep_each([0]',I),rep_each([1]',I) ) ),
                              to_vector(normal_rng(rep_each([0]',I),rep_each([1]',I) ) )),
                              to_vector(normal_rng(rep_each([0]',I),rep_each([1]',I) ) )); //individual-level RN deviations
  //real<lower=0>theta_0_sim = normal_rng(0,1); //w dispersion (sigma for Gaussian)
  real theta_0_sim = fabs(cauchy_rng(0,1)); //w dispersion (sigma for Gaussian)
  
  //random effect values
  matrix[I,3] zp_sim =  std_dev_sim * diag_pre_multiply(to_vector(sd_P_sim), LP_sim);
  vector[I] zp_mu_sim = col(zp_sim,1); //personality
  vector[I] zp_beta_sim = col(zp_sim,2); //plasticity
  vector[I] zp_theta_sim = col(zp_sim,3); //predictability
  
  //RN model
  vector[N] x = to_vector( normal_rng(rep_each([0]', N), rep_each([1]', N) ) );
  vector[N] z_eta_sim = mu_0z_sim + zp_mu_sim[ind] + (beta_1z_sim + zp_beta_sim[ind]).*x ;
  vector[N] z_theta_sim = exp(theta_0z_sim + zp_theta_sim[ind]) ;
  vector[N] z = to_vector(normal_rng(z_eta_sim, z_theta_sim));
  
  //fitness model
  vector[I] w_eta_sim = mu_0_sim + betas_sim[1]*zp_mu_sim + betas_sim[2]*zp_beta_sim + betas_sim[3]*zp_theta_sim +
                 betas_sim[4]*(zp_mu_sim .*zp_mu_sim) + betas_sim[5]*(zp_beta_sim .*zp_beta_sim) +
                 betas_sim[6]*(zp_theta_sim .*zp_theta_sim) +   
                 betas_sim[7]*(zp_mu_sim .*zp_beta_sim) + betas_sim[8]*(zp_mu_sim .*zp_theta_sim) +
                 betas_sim[9]*(zp_beta_sim .*zp_theta_sim) ;
  vector[I] w = to_vector(normal_rng(w_eta_sim, theta_0_sim));
  
}

parameters {
  //fixed population effects
  real mu_0z; //z population intercept
  real beta_1z; //z population slope
  real theta_0z; //z population dispersion
  real mu_0; //w population intercept
  vector[9] betas; //fitness regression coefficients
  
  //random effects (standard deviations)
  vector<lower=0>[3] sd_P; //RN parameter sds
  cholesky_factor_corr[3] LP; //RN parameter correlations
  matrix[I,3] std_dev; //individual-level RN deviations
  real<lower=0> theta_0; //w dispersion (sigma for Gaussian)
}
transformed parameters {
  matrix[I,3] zp; //individual phenotypic RN parameter values
  zp =  std_dev * diag_pre_multiply(sd_P, LP);
}
model{
  //separate RN parameters
  vector[I] zp_mu = col(zp,1); //personality
  vector[I] zp_beta = col(zp,2); //plasticity
  vector[I] zp_theta = col(zp,3); //predictability
  
  //initialize vectors for response models
  vector[N] z_eta; //linear predictor of behavior expectation
  vector[N] z_theta; //linear predictor of behavior dispersion
  vector[I] w_eta; //linear predictor of fitness expectation

  //behavioral RN response model
  z_eta = mu_0z + zp_mu[ind] + (beta_1z + zp_beta[ind]).*x ;
  z_theta = exp(theta_0z + zp_theta[ind]) ;
  z ~ normal(z_eta, z_theta);
    
  //fitness response model
  w_eta = mu_0 + betas[1]*zp_mu + betas[2]*zp_beta + betas[3]*zp_theta +
                 betas[4]*(zp_mu .*zp_mu) + betas[5]*(zp_beta .*zp_beta) + betas[6]*(zp_theta .*zp_theta) +   
                 betas[7]*(zp_mu .*zp_beta) + betas[8]*(zp_mu .*zp_theta) + betas[9]*(zp_beta .*zp_theta) ;
  w ~ normal(w_eta, theta_0);

  //model priors
  
  //fixed effects
  mu_0z ~ normal(0,1);
  beta_1z ~ normal(0,1);
  theta_0z ~ normal(0,1);
  mu_0 ~ normal(0,1);
  theta_0 ~ normal(0,1);
  to_vector(betas) ~ normal(0,1); 
  
  //random effects
  to_vector(sd_P) ~ exponential(1);
  LP ~ lkj_corr_cholesky(2);
  to_vector(std_dev) ~ normal(0,1);
  theta_0 ~ exponential(1);
}

generated quantities{
  int<lower = 0, upper = 1> lt_sim[9]
      = { betas[1] < betas_sim[1], betas[2] < betas_sim[2], betas[3] < betas_sim[3],
          betas[4] < betas_sim[4], betas[5] < betas_sim[5], betas[6] < betas_sim[6],
          betas[7] < betas_sim[7], betas[8] < betas_sim[8], betas[9] < betas_sim[9]};
  
}

