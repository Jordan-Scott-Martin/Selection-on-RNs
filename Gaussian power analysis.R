#############################################################
#prep workspace

#set directory
setwd("...")

#load packages
library(SBC)
library(rstan)
library(rethinking)
library(cmdstanr)
library(ggplot2)
library(reshape2)

#directory for cmdstan installation
set_cmdstan_path("...")

#############################################################
#simulate data

#sim function
simf = function(){
  
  I = as.integer(runif(1, 200, 2001))
  rep_z = as.integer(runif(1,3,8))
  rep_w = as.integer(runif(1,1,6))
  
  #index for individuals
  ind_z = rep(seq(1:I), each = rep_z) 
  ind_w = rep(seq(1:I), each = rep_w)
  
  #total observations
  N_z = I * rep_z
  N_w = I * rep_w
  
  #fixed population effects
  mu_0z = 0 #z intercept
  beta_z = 0 #z slope
  theta_0z = log(0.6) #z dispersion
  mu_0w = 1 #w intercept
  theta_0w = 1
  
  min = 0.1 
  max = 0.5
  b = runif(3,min=min,max=max)
  q = runif(3,min=min,max=max)
  qc = runif(3,min=min,max=max)
  
  #random effects (z)
  sd_I_z = rep(sqrt(0.3), 3)
  LP = t(chol(rethinking::rlkjcorr(1,3,5))) #RN parameter correlations
  std_dev_z = cbind(rnorm(I,0,1),rnorm(I,0,1),rnorm(I,0,1)) #individual-level RN deviations
  std_dev_z = matrix(scale(std_dev_z), ncol = 3)
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
  z = rnorm(N_z, z_eta, z_theta)
  
  #fitness expectation
  
  #unexplained selection
  sd_I_w = runif(1, min=0.3, max=0.7) 
  W_0 = rnorm(I, 0, sd_I_w)
  
  #fitness model
  w_eta = mu_0w + W_0[ind_w]+
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
  return(list(params = 
                list(b = b, q = q, qc = qc, I = I, rep_z = rep_z, rep_w = rep_w,
                     sds = sd_I_z, cor = cor(zp)[lower.tri(cor(zp))],
                     mean_cor = mean(abs(cor(zp)[lower.tri(cor(zp))]))),
              datal = 
                list(I = I, npar = 3, N_z = N_z, N_w = N_w, 
                     ind_z = ind_z, ind_w = ind_w, rep_z = rep_z, rep_w = rep_w,
                     x = x, z = z, W = W)))
}

#generate datasets
dl_sim2 = list()
for(i in 1:1000){
  dl_sim2[[i]] = simf()
}
saveRDS(dl_sim2, "dl_sim2.RDS")

#############################################################
#estimate models

#compile
NLS_mod_fnr = cmdstan_model(stan_file = "NLS power mod_gg_norepw.stan")
NLS_mod_f = cmdstan_model(stan_file = "NLS power mod_gg.stan")

#matrix to hold values
p_mat = matrix(NA, nrow = length(dl_sim2), ncol = 3*3 + length(unlist(dl_sim2[[1]]$params)), 
               dimnames = list(1:length(dl_sim2), 
                               c(c("b1_p","b2_p","b3_p","q1_p","q2_p","q3_p","qc1_p","qc2_p","qc3_p"),
                               names(unlist(dl_sim2[[1]]$params))) ))

#run simulation
length(dl_sim2)
for(i in 1:length(dl_sim2)){
  
  #if 1 fitness measure
  if(dl_sim2[[i]]$params$rep_w==1){
    
    fit = NLS_mod_fnr$sample(
      data = dl_sim2[[i]]$datal,
      iter_sampling = 500,
      iter_warmup = 1000,
      init = 1e-4,
      chains = 4,  
      parallel_chains = 4,
      adapt_delta = 0.80,
      refresh = 10)
    
    b_s = dl_sim2[[i]]$params$b
    q_s = dl_sim2[[i]]$params$q
    qc_s = dl_sim2[[i]]$params$qc
    
    post = fit$draws(format = "data.frame")
    b_p = apply(cbind(post$`b[1]`, post$`b[2]`, post$`b[3]`), 2, function(x) sum(x>0)/length(x))
    q_p = apply(cbind(post$`q[1]`, post$`q[2]`, post$`q[3]`), 2, function(x) sum(x>0)/length(x))
    qc_p = apply(cbind(post$`qc[1]`, post$`qc[2]`, post$`qc[3]`), 2, function(x) sum(x>0)/length(x))
    
    p_mat[i,1:3] = b_p
    p_mat[i,4:6] = q_p
    p_mat[i,7:9] = qc_p
    p_mat[i,10:ncol(p_mat)] = unlist(dl_sim2[[i]]$params)
  }
  
  #if repeated fitness measures
  else{
    fit = NLS_mod_f$sample(
      data = dl_sim2[[i]]$datal,
      iter_sampling = 500,
      iter_warmup = 1000,
      init = 1e-4,
      chains = 4,  
      parallel_chains = 4,
      adapt_delta = 0.80,
      refresh = 10)
    
    b_s = dl_sim2[[i]]$params$b
    q_s = dl_sim2[[i]]$params$q
    qc_s = dl_sim2[[i]]$params$qc
    
    post = fit$draws(format = "data.frame")
    b_p = apply(cbind(post$`b[1]`, post$`b[2]`, post$`b[3]`), 2, function(x) sum(x>0)/length(x))
    q_p = apply(cbind(post$`q[1]`, post$`q[2]`, post$`q[3]`), 2, function(x) sum(x>0)/length(x))
    qc_p = apply(cbind(post$`qc[1]`, post$`qc[2]`, post$`qc[3]`), 2, function(x) sum(x>0)/length(x))
    
    p_mat[i,1:3] = b_p
    p_mat[i,4:6] = q_p
    p_mat[i,7:9] = qc_p
    p_mat[i,10:ncol(p_mat)] = unlist(dl_sim2[[i]]$params)
  }
  
  saveRDS(p_mat, "p_mat2.RDS")
}



#############################################################
#plot results

#load results
p_mat = readRDS("p_mat2.RDS")

#wide to long
p_mat = data.frame(p_mat)
pml = melt(p_mat, id.vars = colnames(p_mat)[10:ncol(p_mat)])
pml$truev = melt(p_mat[,10:18])[,2]
pml = melt(pml[10:ncol(pml)], id.var = c("variable", "value"), value.name = "parv")
colnames(pml) = c("RNpar", "power", "par", "parv")
pml = pml[! pml$par %in% c("sds1","sds2","sds3", "cor1", "cor2", "cor3"),]
pml = droplevels(pml)
pml$sg = ifelse(grepl("b",pml$RNpar), "beta", ifelse(grepl("qc",pml$RNpar), "corr", "quad"))
pml$col = ifelse(grepl("1",pml$RNpar), "red", ifelse(grepl("2",pml$RNpar), "blue", "yellow"))
pml$sg = factor(pml$sg, levels = c("beta","quad","corr"))
pml$col = factor(pml$col, levels = c("red","blue","yellow"))


scatter.smooth(pml$I[pml$variable=="b1_p"], pml$value[pml$variable=="b1_p"],
               lpars = list(col = 'blue', lwd = 3))
scatter.smooth(pml$I[pml$variable=="b2_p"], pml$value[pml$variable=="b2_p"],
               lpars = list(col = 'blue', lwd = 3))
scatter.smooth(pml$I[pml$variable=="b3_p"], pml$value[pml$variable=="b3_p"],
               lpars = list(col = 'blue', lwd = 3))

scatter.smooth(pml$I[pml$variable=="q1_p"], pml$value[pml$variable=="q1_p"],
               lpars = list(col = 'blue', lwd = 3))
scatter.smooth(pml$I[pml$variable=="q2_p"], pml$value[pml$variable=="q2_p"],
               lpars = list(col = 'blue', lwd = 3))
scatter.smooth(pml$I[pml$variable=="q3_p"], pml$value[pml$variable=="q3_p"],
               lpars = list(col = 'blue', lwd = 3))

scatter.smooth(pml$I[pml$variable=="qc1_p"], pml$value[pml$variable=="qc1_p"],
               lpars = list(col = 'blue', lwd = 3))
scatter.smooth(pml$I[pml$variable=="qc2_p"], pml$value[pml$variable=="qc2_p"],
               lpars = list(col = 'blue', lwd = 3))
scatter.smooth(pml$I[pml$variable=="qc3_p"], pml$value[pml$variable=="qc3_p"],
               lpars = list(col = 'blue', lwd = 3))

scatter.smooth(pml$rep_z[pml$variable=="qc1_p"], pml$value[pml$variable=="qc1_p"],
               lpars = list(col = 'blue', lwd = 3))
scatter.smooth(pml$rep_z[pml$variable=="qc2_p"], pml$value[pml$variable=="qc2_p"],
               lpars = list(col = 'blue', lwd = 3))
scatter.smooth(pml$rep_w[pml$variable=="qc3_p"], pml$value[pml$variable=="qc3_p"],
               lpars = list(col = 'blue', lwd = 3))

pow.f = 
  ggplot(pml, aes(x = parv, y = power, group = RNpar, color = col, fill = col)) +
  facet_wrap(.~ sg*par, scales = "free_x", nrow = 3) + 
  coord_cartesian(ylim = c(0.5, 1))+
  scale_y_continuous(expand=c(0,0))+
  scale_color_manual(values = c("red","blue","#fad416"))+
  scale_fill_manual(values = colorspace::lighten(c("red","blue","#fad416"), amount = 0.2))+
  labs(x = "", y = "Posterior probability (power)\n")+
  geom_line(se = FALSE, stat = "smooth", method = "lm", formula = y ~ poly(x,2), 
           linewidth = 1.5, alpha = 0.5)+
  theme(axis.title.y = element_text(face = "bold", size = 12),
        legend.title = element_blank(),
        legend.position = "top",
        strip.background = element_blank(),
        panel.border=element_rect(fill=NA,color="black", 
                                  linewidth=1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(1, "lines")) +  
  guides(colour = "none")

ggsave("power_plot.tiff", pow.f, dpi = 600, units = "in", width = 8, height = 6)

