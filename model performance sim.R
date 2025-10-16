#############################################################
#prep workspace

#set directory
setwd("...")

#load packages
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
  
  I = as.integer(runif(1, 100, 1000))
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
  theta_0z = log(2) #z dispersion
  mu_0w = 1 #w intercept
  theta_0w = 1
  
  min = 0.1 
  max = 0.5
  b = runif(3,min=min,max=max)
  q = runif(3,min=min,max=max)
  qc = runif(3,min=min,max=max)
  
  #random effects (z)
  sd_I_z = rep(1, 3)
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
  sd_I_w = 1
  W_0 = rnorm(I, 0, sd_I_w)
  
  #fitness model
  w_eta = mu_0w + W_0[ind_w]+
    b[1]*zp_mu[ind_w] + 
    b[2]*zp_beta[ind_w] + 
    b[3]*zp_theta[ind_w] +
    0.5*q[1]*(zp_mu[ind_w] * zp_mu[ind_w]) + 
    0.5*q[2]*(zp_beta[ind_w] * zp_beta[ind_w]) + 
    0.5*q[3]*(zp_theta[ind_w] * zp_theta[ind_w]) +   
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
dl_sim = list()
for(i in 1:500){
  dl_sim[[i]] = simf()
}
saveRDS(dl_sim, "dl_sim.RDS")
dl_sim = readRDS("dl_sim.RDS")

#############################################################
#estimate models

#compile
NLS_mod_fnr = cmdstan_model(stan_file = "nls power mod_gg_norepw.stan")
NLS_mod_f = cmdstan_model(stan_file = "nls power mod_gg.stan")

#matrix to hold values
p_mat = matrix(NA, nrow = length(dl_sim), ncol = 3*3*3 + length(unlist(dl_sim[[1]]$params)), 
               dimnames = list(1:length(dl_sim), 
                               c(c("b1_p","b2_p","b3_p","q1_p","q2_p","q3_p","qc1_p","qc2_p","qc3_p"),
                                 c("b1_b","b2_b","b3_b","q1_b","q2_b","q3_b","qc1_b","qc2_b","qc3_b"),
                                 c("b1_a","b2_a","b3_a","q1_a","q2_a","q3_a","qc1_a","qc2_a","qc3_a"),
                               names(unlist(dl_sim[[1]]$params))) ))

#run simulation
for(i in 1:length(dl_sim)){
  
  #if 1 fitness measure
  if(dl_sim[[i]]$params$rep_w==1){
    
    fit = NLS_mod_fnr$sample(
      data = dl_sim[[i]]$datal,
      iter_sampling = 500,
      iter_warmup = 1000,
      init = 1e-4,
      chains = 4,  
      parallel_chains = 4,
      adapt_delta = 0.80,
      refresh = 10)
    
    b_s = dl_sim[[i]]$params$b
    q_s = dl_sim[[i]]$params$q
    qc_s = dl_sim[[i]]$params$qc
    
    post = fit$draws(format = "data.frame")
    b_p = apply(cbind(post$`b[1]`, post$`b[2]`, post$`b[3]`), 2, function(x) sum(x>0)/length(x))
    q_p = apply(cbind(post$`q[1]`, post$`q[2]`, post$`q[3]`), 2, function(x) sum(x>0)/length(x))
    qc_p = apply(cbind(post$`qc[1]`, post$`qc[2]`, post$`qc[3]`), 2, function(x) sum(x>0)/length(x))
    
    b_b = apply(cbind(post$`b[1]`, post$`b[2]`, post$`b[3]`), 2, median) - b_s
    q_b = apply(cbind(post$`q[1]`, post$`q[2]`, post$`q[3]`)*2 , 2, median) - q_s
    qc_b = apply(cbind(post$`qc[1]`, post$`qc[2]`, post$`qc[3]`), 2, median) - qc_s
    
    b_a = apply(cbind(post$`b[1]`, post$`b[2]`, post$`b[3]`) - b_s, 2, function(x) mean(x*x) )
    q_a = apply(cbind(post$`q[1]`, post$`q[2]`, post$`q[3]`)*2 - q_s, 2, function(x) mean(x*x) )
    qc_a = apply(cbind(post$`qc[1]`, post$`qc[2]`, post$`qc[3]`) - qc_s, 2, function(x) mean(x*x) )
    
    p_mat[i,1:3] = b_p
    p_mat[i,4:6] = q_p
    p_mat[i,7:9] = qc_p
    p_mat[i,10:12] = b_b
    p_mat[i,13:15] = q_b
    p_mat[i,16:18] = qc_b
    p_mat[i,19:21] = b_a
    p_mat[i,22:24] = q_a
    p_mat[i,25:27] = qc_a
    p_mat[i,28:ncol(p_mat)] = unlist(dl_sim[[i]]$params)
  }
  
  #if repeated fitness measures
  else{
    fit = NLS_mod_f$sample(
      data = dl_sim[[i]]$datal,
      iter_sampling = 500,
      iter_warmup = 1000,
      init = 1e-4,
      chains = 4,  
      parallel_chains = 4,
      adapt_delta = 0.80,
      refresh = 10)
    
    b_s = dl_sim[[i]]$params$b
    q_s = dl_sim[[i]]$params$q
    qc_s = dl_sim[[i]]$params$qc
    
    post = fit$draws(format = "data.frame")
    b_p = apply(cbind(post$`b[1]`, post$`b[2]`, post$`b[3]`), 2, function(x) sum(x>0)/length(x))
    q_p = apply(cbind(post$`q[1]`, post$`q[2]`, post$`q[3]`), 2, function(x) sum(x>0)/length(x))
    qc_p = apply(cbind(post$`qc[1]`, post$`qc[2]`, post$`qc[3]`), 2, function(x) sum(x>0)/length(x))
    
    b_b = apply(cbind(post$`b[1]`, post$`b[2]`, post$`b[3]`), 2, median) - b_s
    q_b = apply(cbind(post$`q[1]`, post$`q[2]`, post$`q[3]`)*2 , 2, median) - q_s
    qc_b = apply(cbind(post$`qc[1]`, post$`qc[2]`, post$`qc[3]`), 2, median) - qc_s
    
    b_a = apply(cbind(post$`b[1]`, post$`b[2]`, post$`b[3]`) - b_s, 2, function(x) mean(x*x) )
    q_a = apply(cbind(post$`q[1]`, post$`q[2]`, post$`q[3]`)*2 - q_s, 2, function(x) mean(x*x) )
    qc_a = apply(cbind(post$`qc[1]`, post$`qc[2]`, post$`qc[3]`) - qc_s, 2, function(x) mean(x*x) )
    
    p_mat[i,1:3] = b_p
    p_mat[i,4:6] = q_p
    p_mat[i,7:9] = qc_p
    p_mat[i,10:12] = b_b
    p_mat[i,13:15] = q_b
    p_mat[i,16:18] = qc_b
    p_mat[i,19:21] = b_a
    p_mat[i,22:24] = q_a
    p_mat[i,25:27] = qc_a
    p_mat[i,28:ncol(p_mat)] = unlist(dl_sim[[i]]$params)
  }
  
  saveRDS(p_mat, "p_mat.RDS")
}


#############################################################
#plot results

#load results
p_mat = readRDS("p_mat.RDS")

#wide to long
p_mat$wsd = sqrt(rowSums(p_mat[,c(28:36)]^2)+1)

prob = reshape2::melt(p_mat[,c(1:9)], value.name = "prob")
bias = reshape2::melt(p_mat[,c(10:18)], value.name = "bias")
es = reshape2::melt(p_mat[,c(28:36)], value.name = "es")
pml =  reshape2::melt(p_mat[,c(19:27, 37:39, 46, 47)], value.name = "rmsd", 
                      id.vars = colnames(p_mat[,c(37:39, 46, 47)]))
pml$bias = bias$bias
pml$power = prob$prob
pml$es = es$es
pml$variable = es$variable
pml$bias = pml$bias
pml$rmsd = sqrt(pml$rmsd) 
pml = subset(pml, select = -wsd)
pml = droplevels(pml)
pml$sg = ifelse(grepl("b",pml$variable), "beta", ifelse(grepl("qc",pml$variable), "corr", "quad"))
pml$col = ifelse(grepl("1",pml$variable), "red", ifelse(grepl("2",pml$variable), "blue", "yellow"))
pml = reshape2::melt(pml, id.vars = colnames(pml[,c(5:8,10:11)]), variable.name = "cond", value.name = "condv")
pml = reshape2::melt(pml, id.vars = colnames(pml[,c("variable","sg","col","cond","condv")]), variable.name = "par", value.name = "val")
pml$sg = factor(pml$sg, levels = c("beta","quad","corr"))
pml$par = factor(pml$par, levels = c("bias","rmsd","power"))
pml$variable = factor(pml$variable, levels = c("b1","b2","b3", "q1", "q2", "q3", "qc1", "qc2", "qc3"))

#manual limits and breaks
Ib = c(100, 400, 700, 1000)
cb = c(0.1, 0.3, 0.5)
eb = c(0.1, 0.3, 0.5)

r1y = c(-1, 0, 1)
r2y = c(0.1, 0.5, 0.9)
r3y = c(0.5, 0.75, 1)

r1l = c(-1, 1)
r2l = c(0.1, 1)
r3l = c(0.5, 1)

perform = 
  ggplot(pml, aes(x = condv, y = val, group = variable, color = variable, fill = variable)) +
  facet_wrap(.~ par*cond, scales = "free", ncol = 5) +
  scale_color_manual(breaks = levels(pml$variable),
                     values = colorspace::qualitative_hcl(9, palette = "Dark 3"),
                     labels = c(
                       expression(beta[mu]),
                       expression(beta[beta]),
                       expression(beta[sigma]),
                       expression(gamma[mu]),
                       expression(gamma[beta]),
                       expression(gamma[sigma]),
                       expression(gamma[mu*","*beta]),
                       expression(gamma[mu*","*sigma]),
                       expression(gamma[beta*","*sigma])
                     ))+
  scale_fill_manual(breaks = levels(pml$variable),
                       values = colorspace::qualitative_hcl(9, palette = "Dark 3") )+
  labs(x = "", y = "")+
  geom_point(alpha = 0, pch = 1, size = 0.5)+
  geom_line(se = FALSE, stat = "smooth", method = "lm", formula = y ~ poly(x,1), 
            linewidth = 0.5, alpha = 0.7)+
    
  ggh4x::facetted_pos_scales(
    x = list(
  scale_x_continuous(breaks = Ib), NULL, NULL, scale_x_continuous(breaks = cb), scale_x_continuous(breaks = eb),
  scale_x_continuous(breaks = Ib), NULL, NULL, scale_x_continuous(breaks = cb), scale_x_continuous(breaks = eb),
  scale_x_continuous(breaks = Ib), NULL, NULL, scale_x_continuous(breaks = cb), scale_x_continuous(breaks = eb)),
  y = list(
    scale_y_continuous(limits = r1l, breaks = r1y), scale_y_continuous(limits = r1l, breaks = r1y),
    scale_y_continuous(limits = r1l, breaks = r1y), scale_y_continuous(limits = r1l, breaks = r1y), 
    scale_y_continuous(limits = r1l, breaks = r1y),
    
    scale_y_continuous(limits = r2l, breaks = r2y), scale_y_continuous(limits = r2l, breaks = r2y),
    scale_y_continuous(limits = r2l, breaks = r2y), scale_y_continuous(limits = r2l, breaks = r2y), 
    scale_y_continuous(limits = r2l, breaks = r2y),
    
    scale_y_continuous(limits = r3l, breaks = r3y), scale_y_continuous(limits = r3l, breaks = r3y),
    scale_y_continuous(limits = r3l, breaks = r3y), scale_y_continuous(limits = r3l, breaks = r3y), 
    scale_y_continuous(limits = r3l, breaks = r3y))
  )+
    
  theme(axis.title.y = element_text(face = "bold", size = 12),
        legend.title = element_blank(),
        legend.position = "top",
        #legend.key.width = unit(1.5, "cm"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_text(size = 6),
        panel.border=element_rect(fill=NA,color="black", 
                                  linewidth=1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.25, "lines")) +
  guides(fill = "none", color = guide_legend(nrow = 1))

ggsave("performance_plot.png", perform, dpi = 600, units = "in", width = 6, height = 3.5)
