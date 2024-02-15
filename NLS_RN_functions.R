list.of.packages = c("reshape2", "ggplot2","tidybayes","plyr","mvtnorm","rethinking", "devtools")
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
if("rethinking" %in% new.packages) devtools::install_github("rmcelreath/rethinking")
library(reshape2)
suppressMessages(library(ggplot2))
library(tidybayes)
library(plyr)
library(mvtnorm)
suppressMessages(library(rethinking))

#this function generates a dataset of Gaussian phenotype and fitness
#measures for a nonlinear selection analysis on reaction norms
#The function assumes balanced sampling and so requires 
#input values for repeated measures
#J = total number of individuals (even number)
#rep_z = repeated phenotype measures per subject (>1)
#rep_W = repeated fitness measures per subject (>0)
#l_es = lower selection effect size (absolute value)
#u_es = upper selection effect size (absolute value)
#stan_list = should the data be returned as a list (T) or data frame (F)
sim_RN_Gaus = function(J, rep_z, rep_W, l_es, u_es,
                       mu_0 = 0, beta_x = 0, sigma_0 = log(0.6), W_0 = 1, delta = 1,
                       V_P = 0.3, V_W0 = 0.3, lkj_cor = 10){
  try(if(rep_z < 2) stop("Estimating individual reaction norms requires repeated measures."))
  try(if(l_es < 0) stop("l_es is an absolute value. Please input a positive number."))
  try(if(u_es < 0) stop("u_es is an absolute value. Please input a positive number."))
  
  #index for individuals
  ind_z = rep(seq(1:J), each = rep_z) 
  ind_W = rep(seq(1:J), each = rep_W)
  
  #total observations
  N_z = J * rep_z
  N_W = J * rep_W
  
  #fixed population effects
  mu_0z = mu_0 #z intercept
  beta_z = beta_x #z slope
  sigma_0z = sigma_0 #z dispersion
  theta_0w = W_0 #w intercept
  delta_0w = delta
  
  min = l_es 
  max = u_es
  b = runif(3,min=min,max=max) * sample(c(-1,1), 3, prob = c(0.5,0.5), replace = T)
  q = runif(3,min=min,max=max) * sample(c(-1,1), 3, prob = c(0.5,0.5), replace = T)
  qc = runif(3,min=min,max=max) * sample(c(-1,1), 3, prob = c(0.5,0.5), replace = T)
  
  #random effects (z)
  sd_P_z = rep(sqrt(V_P), 3)
  LP = t(chol(rethinking::rlkjcorr(1,3,lkj_cor))) #RN parameter correlations
  std_dev_z = cbind(rnorm(J,0,1),rnorm(J,0,1),rnorm(J,0,1)) #individual-level RN deviations
  std_dev_z = matrix(scale(std_dev_z), ncol = 3)
  zp = std_dev_z %*% t( diag(sd_P_z) %*% LP )
  
  #RN parameters
  zp_mu = zp[,1] #intercepts
  zp_beta = zp[,2] #slopes
  zp_sigma = zp[,3] #residuals
  
  #sample environmental states
  x = rnorm(N_z, 0, 1)
  
  #phenotype expectations and dispersions
  z_eta = mu_0z + zp_mu[ind_z] + (beta_z + zp_beta[ind_z])*x
  z_sigma = sqrt(exp(sigma_0z + zp_sigma[ind_z]))
  z = rnorm(N_z, z_eta, z_sigma)
  
  #fitness expectation
  
  if(rep_W > 1){
  #unexplained selection
  sd_I_w = sqrt(V_W0) 
  W_0 = rnorm(J, 0, sd_I_w)
  
  #fitness model
  W_eta = theta_0w + W_0[ind_W]+
    b[1]*zp_mu[ind_W] + 
    b[2]*zp_beta[ind_W] + 
    b[3]*zp_sigma[ind_W] +
    q[1]*(zp_mu[ind_W] * zp_mu[ind_W]) + 
    q[2]*(zp_beta[ind_W] * zp_beta[ind_W]) + 
    q[3]*(zp_sigma[ind_W] * zp_sigma[ind_W]) +   
    qc[1]*(zp_mu[ind_W] * zp_beta[ind_W]) + 
    qc[2]*(zp_mu[ind_W] * zp_sigma[ind_W]) + 
    qc[3]*(zp_beta[ind_W] * zp_sigma[ind_W])
  W = rnorm(N_W, W_eta, delta_0w)
  }
  else{
  #fitness model
  W_eta = theta_0w +
    b[1]*zp_mu[ind_W] + 
    b[2]*zp_beta[ind_W] + 
    b[3]*zp_sigma[ind_W] +
    q[1]*(zp_mu[ind_W] * zp_mu[ind_W]) + 
    q[2]*(zp_beta[ind_W] * zp_beta[ind_W]) + 
    q[3]*(zp_sigma[ind_W] * zp_sigma[ind_W]) +   
    qc[1]*(zp_mu[ind_W] * zp_beta[ind_W]) + 
    qc[2]*(zp_mu[ind_W] * zp_sigma[ind_W]) + 
    qc[3]*(zp_beta[ind_W] * zp_sigma[ind_W])
  W = rnorm(N_W, W_eta, delta_0w) 
  }
  
  #list of variables and generated modeled data
  return(list(J = J, N_z = N_z, N_W = N_W, ind_z = ind_z, ind_W = ind_W, 
              x = x, z = z, W = W, 
              true_b = b, true_q = q, true_qc = qc))
}


#more functions for simulating data and plotting results are in development
#and should be available for general use within a month or two