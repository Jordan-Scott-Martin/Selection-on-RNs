#############################################################
#prep workspace

#set directory
setwd("...")

#load packages
library(SBC)
library(rstan)
library(rethinking)
library(cmdstanr)

#cache for SBC
cache_dir = "./_basic_usage_cmdstan_SBC_cache"
if(!dir.exists(cache_dir)) {dir.create(cache_dir) }

#directory for cmdstan installation
set_cmdstan_path("...")

#############################################################
#load Gaussian NLS selection model
SBC_mod1 = cmdstan_model(stan_file = "SBC_model_NLS_1.stan", 
                         stanc_options = list("O1"))

res = SBC_mod1$sample(data = df1$generated[[1]], iter_sampling = 500, iter_warmup = 500,
                      parallel_chains = 4, init = 0, adapt_delta = 0.9)
shinystan::launch_shinystan(res)

#function for generating datasets and parameter values
generator_function1 = function(I, rep_z, rep_w){
  
  #generate values for the data and model parameters defined in the model
  #see Stan file for details
  
  #index for individuals
  ind_z = rep(seq(1:I), each = rep_z) 
  ind_w = rep(seq(1:I), each = rep_w)
  
  #total observations
  N_z = I * rep_z
  N_w = I * rep_w

  #fixed population effects
  mu_0z = rnorm(1,0,1) #z population intercept
  beta_z = rnorm(1,0,1) #z population slope
  theta_0z = rnorm(1,0,1) #z population dispersion
  mu_0w = 1 #w population intercept
  b = rnorm(3,0,1) #fitness regression coefficients
  q = rnorm(3,0,1) #fitness regression coefficients
  qc = rnorm(3,0,1) #fitness regression coefficients
  
  #random effects
  sd_I_z = rexp(3,2) #RN parameter sds
  LP = t(chol(rethinking::rlkjcorr(1,3,2))) #RN parameter correlations
  std_dev_z = cbind(rnorm(I,0,1),rnorm(I,0,1),rnorm(I,0,1)) #individual-level RN deviations
  std_dev_z = matrix(scale(std_dev_z), ncol = 3)
  
  sd_I_w = rexp(1,2) #unexplained selection sd
  std_dev_w = rnorm(I, 0, 1) #individual unexplained selection deviations
  std_dev_w = matrix(scale(std_dev_w), ncol = 1)
  theta_0w = rexp(1,2) #w dispersion (sigma for Gaussian)
  
  #generate data
  zp = std_dev_z %*% t( diag(sd_I_z) %*% LP )
  
  #RN parameters
  zp_mu = zp[,1] #intercepts
  zp_beta = zp[,2] #slopes
  zp_theta = zp[,3] #residuals
  
  #sample environmental states
  x = rnorm(N_z, 0, 1)

  #phenotype expectations and dispersions
  z_eta = mu_0z + zp_mu[ind_z] + (beta_z + zp_beta[ind_z])*x
  z_theta = sqrt(exp(theta_0z + zp_theta[ind_z]))
  
  #phenotype function
  z = rnorm(N_z, z_eta, z_theta)
    
  #fitness expectation
  
  #unexplained selection
  W_0 = rnorm(I, 0, sd_I_w)
  
  w_eta = mu_0w +  W_0[ind_w] +
          b[1]*zp_mu[ind_w] + 
          b[2]*zp_beta[ind_w] + 
          b[3]*zp_theta[ind_w] +
          q[1]*(zp_mu[ind_w] * zp_mu[ind_w]) + 
          q[2]*(zp_beta[ind_w] * zp_beta[ind_w]) + 
          q[3]*(zp_theta[ind_w] * zp_theta[ind_w]) +   
          qc[1]*(zp_mu[ind_w] * zp_beta[ind_w]) + 
          qc[2]*(zp_mu[ind_w] * zp_theta[ind_w]) + 
          qc[3]*(zp_beta[ind_w] * zp_theta[ind_w])
  W = rnorm(N_w, w_eta, theta_0w)
  
  #list of variables and generated modeled data
  list(variables = 
        list(mu_0z = mu_0z, beta_z = beta_z, theta_0z = theta_0z, 
             mu_0w = mu_0w, theta_0w = theta_0w, 
             b = b, q = q, qc = qc,
             sd_I_z = sd_I_z, sd_I_w = sd_I_w, LP = LP, 
             std_dev_z = std_dev_z, std_dev_w = as.vector(std_dev_w)),
        generated = list(I = I, npar = 3, N_z = N_z, N_w = N_w, 
                         ind_z = ind_z, ind_w = ind_w, 
                          x = x, z = z, W = W))
}

#generate objects for SBC
set.seed(9)
n_sims = 300  # Number of sim per run (x 3 total)
I = 100 #number of subjects
rep_z = 3 #number of repeated phenotype measures
rep_w = 2 #number of repeated fitness measures

SBCf1 = SBC_generator_function(generator_function1, I = I, rep_z = rep_z, rep_w = rep_w)
df1 = generate_datasets(SBCf1, n_sims)
backend1 = SBC_backend_cmdstan_sample(SBC_mod1, iter_sampling = 500, iter_warmup = 1000, 
                                      init = 0, chains = 6,  parallel_chains = 6, 
                                      adapt_delta = 0.95, max_treedepth = 10) 

#conduct SBC 
results1 =  compute_SBC(datasets = df1, backend = backend1,
                        cache_mode = "results", 
                        cache_location = file.path(cache_dir, "results"),
                        keep_fits = FALSE)
saveRDS(results1, "results_SBC_NLS_1.RDS")
results1 = readRDS("results_SBC_NLS.RDS")
#diagnostics
hist(results1$stats$rhat)
hist(results1$stats$ess_tail)

#subset to remove problematic simulations
keep = results1$backend_diagnostics$sim_id[
  results1$backend_diagnostics$n_divergent == 0]
results1nd = results1[keep]

hist(unlist(results1_1ndt$stats[results1_1ndt$stats$variable=="betas[9]","simulated_value"]))

#plot results
bp = plot_ecdf_diff(results1, variables = paste0("b","[",1:3,"]"), prob = 0.95)
qp = plot_ecdf_diff(results1, variables = paste0("q","[",1:3,"]"), prob = 0.95)
qcp = plot_ecdf_diff(results1, variables = paste0("qc","[",1:3,"]"), prob = 0.95)

library(ggplot2)
bdiff = 
  bp +    
    ylab(" \n")+
    xlab("")+
    ggtitle("Directional selection\n ")+
  coord_cartesian(ylim=c(-0.17,0.17))+
  scale_x_continuous(labels=c("0","","0.5","","1"))+
    theme(plot.title =element_text(size=12, face="bold",hjust=0.5),
          axis.ticks.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_text(size=12,face="bold"),
          axis.title.y=element_text(size=12,face="bold", vjust = 0.5),
          axis.text.x=element_text(size=9),
          axis.text.y=element_text(size=9),
          axis.line = element_line(linewidth = 1),
          
          panel.border=element_rect(fill=NA,color="black", linewidth=1, 
                                    linetype="solid"),
          strip.text = element_blank(),
          strip.background = element_blank(),
          panel.background= element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.05,0.5,0.05,0.5), "cm"))+
    guides(fill="none", color ="none")
  
  qdiff = 
  qp +    
    ylab("ECDF difference\n")+
    xlab("")+
    ggtitle("Stabilizing / disruptive selection\n ")+
    coord_cartesian(ylim=c(-0.17,0.17))+
    scale_x_continuous(labels=c("0","","0.5","","1"))+
    theme(plot.title =element_text(size=12, face="bold",hjust=0.5),
          axis.ticks.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_text(size=12,face="bold"),
          axis.title.y=element_text(size=12,face="bold", vjust = 0.5),
          axis.text.x=element_text(size=9),
          axis.text.y=element_text(size=9),
          axis.line = element_line(linewidth = 1),
          
          panel.border=element_rect(fill=NA,color="black", linewidth=1, 
                                    linetype="solid"),
          strip.text = element_blank(),
          strip.background = element_blank(),
          panel.background= element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.05,0.5,0.05,0.5), "cm"))+
    guides(fill="none", color ="none")
  
  qcdiff = 
  qcp +    
    ylab(" \n")+
    xlab("\n Fractional rank")+
    ggtitle("Correlational selection\n ")+
    coord_cartesian(ylim=c(-0.17,0.17))+
    scale_x_continuous(labels=c("0","","0.5","","1"))+
    theme(plot.title =element_text(size=12, face="bold",hjust=0.5),
          axis.ticks.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_text(size=12,face="bold"),
          axis.title.y=element_text(size=12,face="bold", vjust = 0.5),
          axis.text.x=element_text(size=9),
          axis.text.y=element_text(size=9),
          axis.line = element_line(linewidth = 1),
          
          panel.border=element_rect(fill=NA,color="black", linewidth=1, 
                                    linetype="solid"),
          strip.text = element_blank(),
          strip.background = element_blank(),
          panel.background= element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.05,0.5,0.05,0.5), "cm"))+
    guides(fill="none", color ="none")

library(cowplot)
sbcp = plot_grid(bdiff, qdiff, qcdiff, ncol = 1, align = "hv")
save_plot("SBC NLS.png", sbcp, base_height = 6.5, base_width = 5.5, units = "in")

bdiff$facet$params$nrow=4
bdiff$facet$params$ncol=3
bdiff = bdiff + theme(plot.title =element_text(size=12, face="bold",hjust=0.5),
                      axis.ticks.y=element_blank(),
                      axis.ticks.x=element_blank(),
                      axis.title.x=element_text(size=12,face="bold"),
                      axis.title.y=element_text(size=12,face="bold", angle = 0, vjust = 0.5),
                      axis.text.x=element_blank(),
                      axis.text.y=element_blank(),
                      axis.line = element_line(linewidth = 1),
                      
                      panel.border=element_rect(fill=NA,color="black", linewidth=1, 
                                                linetype="solid"),
                      strip.text = element_blank(),
                      strip.background = element_blank(),
                      panel.background= element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"))+
  guides(fill="none", color ="none")

cdiff$facet$params$nrow=4
cdiff$facet$params$ncol=3
cdiff = cdiff + theme(plot.title =element_text(size=12, face="bold",hjust=0.5),
                      axis.ticks.y=element_blank(),
                      axis.ticks.x=element_blank(),
                      axis.title.x=element_text(size=12,face="bold"),
                      axis.title.y=element_text(size=12,face="bold", angle = 0, vjust = 0.5),
                      axis.text.x=element_blank(),
                      axis.text.y=element_blank(),
                      axis.line = element_line(linewidth = 1),
                      
                      panel.border=element_rect(fill=NA,color="black", linewidth=1, 
                                                linetype="solid"),
                      strip.text = element_blank(),
                      strip.background = element_blank(),
                      panel.background= element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"))+
  guides(fill="none", color ="none")

library(cowplot)
cbp = plot_grid(bdiff, cdiff, ncol = 2, rel_widths = c(1,1))
save_plot("fig 2 crn.png", cbp, base_height = 3.5, base_width = 6)


