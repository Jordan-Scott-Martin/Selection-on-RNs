#load packages
library(mvtnorm)

## increase memory to avoid crashing with loop
memory.limit(size=100000)

#set directory
setwd("C:/Users/jormar/Dropbox/Selection-on-RNs/BLUP manuscript")

#Stan settings
options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)

#simulation parameters
cor = 0.3 #RN correlations
sd = sqrt(0.3) #RN standard deviation
beta = 0.3
res = sqrt(0.5) #observation-level residual variance
popint = 0 #population RN intercept
popslope = 0 #population RN slope
popdisp = sqrt(0.5) #population RN dispersion
I = 100 #number of individuals
repm = 3 #repeated behavioral measures

#simulation algorithms
sim_p = function(cor,sd,beta,res,popint,popslope,popdisp,N,repm){

  #RN parameters
  z_p = rnorm(I, mean = 0, sd = sd) 
  personality = z_p
  
  #environmental covariate (z-score)
  x = rnorm(I*repm, 0, 1)
  
  #behavioral response model
  ind = rep(1:I, each = repm) #index of repeated individual measures
  z_mu = popint + personality[ind]
  z = rnorm(I*repm, mean = z_mu, sd = popdisp )
  
  #beta coefficients
  betas = sample(beta, size=2, replace=TRUE)
  
  #fitness response model
  w_mu = 1 + personality * betas[1] +  betas[2]*(personality^2)
  w = rnorm(I, mean = w_mu, sd = res)
  
  data = list(x = x, z = z, w = w, ind = ind, I = I, N = I*repm, betas = betas)
  return(data)
}
sim_pp = function(cor,sd,beta,res,popint,popslope,popdisp,N,repm){

  #RN parameters
  cors = sample(c(cor, -1*cor), size=2*2, replace=TRUE) #random sign
  R = matrix(cors, nrow=2, ncol=2 )
  R[lower.tri(R)] = t(R)[lower.tri(R)] #force symmetric
  diag(R) = 1
  S = matrix( c(sd,0,0,sd), nrow=2, ncol=2)
  P = S %*% R %*% S
  z_p = rmvnorm(I, mean = rep(0,2), sigma = P) 
  personality = z_p[,1]
  plasticity = z_p[,2]
  
  #environmental covariate (z-score)
  x = rnorm(I*repm, 0, 1)
  
  #behavioral response model
  ind = rep(1:I, each = repm) #index of repeated individual measures
  
  z_mu = popint + personality[ind] + (popslope + plasticity[ind])*x
  z = rnorm(I*repm, mean = z_mu, sd = popdisp )
  
  #beta coefficients
  betas = sample(beta, size=5, replace=TRUE)
  
  #fitness response model
  w_mu = 1 + cbind(personality,plasticity) %*% betas[1:2] +
             betas[3]*(personality^2) + betas[4]*(plasticity^2)+
             betas[5]*(personality*plasticity)
  w = rnorm(I, mean = w_mu, sd = res)
  
  data = list(x = x, z = z, w = w, ind = ind, I = I, N = I*repm, R = R, betas = betas)
  
  return(data)
}
sim_ppp = function(cor,sd,beta,res,popint,popslope,popdisp,N,repm){

  #RN parameters
  cors = sample(c(cor, -1*cor), size=3*3, replace=TRUE) #random sign
  R = matrix(cors, nrow=3, ncol=3 )
  R[lower.tri(R)] = t(R)[lower.tri(R)] #force symmetric
  diag(R) = 1
  S = matrix( c(sd,0,0,0,sd,0,0,0,sd), nrow=3, ncol=3 )
  P = S %*% R %*% S
  z_p = rmvnorm(I, mean = rep(0,3), sigma = P) 
  personality = z_p[,1]
  plasticity = z_p[,2]
  predictability = z_p[,3]
  
  #environmental covariate (z-score)
  x = rnorm(I*repm, 0, 1)
  
  #behavioral response model
  ind = rep(1:I, each = repm) #index of repeated individual measures
  
  z_mu = popint + personality[ind] + (popslope + plasticity[ind])*x
  z_sigma = log(popdisp) + predictability[ind]
  z = rnorm(I*repm, mean = z_mu, sd = exp(z_sigma) )
  
  #beta coefficients
  betas = sample(beta, size=9, replace=TRUE)
  
  #fitness response model
  w_mu = 1 + cbind(personality,plasticity,predictability) %*% betas[1:3] +
             betas[4]*(personality^2)+betas[5]*(plasticity^2)+betas[6]*(predictability^2)+
             betas[7]*(personality*plasticity)+betas[8]*(personality*predictability)+
             betas[9]*(plasticity*predictability)
  w = rnorm(I, mean = w_mu, sd = res)
  
  data = list(x = x, z = z, w = w, ind = ind, I = I, N = I*repm, R = R, betas = betas)
  
  return(data)
}

#simulate dataset
sim_p(cor,sd,beta,res,popint,popslope,popdisp,N,repm)
sim_pp(cor,sd,beta,res,popint,popslope,popdisp,N,repm)
sim_ppp(cor,sd,beta,res,popint,popslope,popdisp,N,repm)

##########################################
#Personality only, N 100, 3
##########################################
I = 100; repm = 3; sim = 200

df.p100.3 = list()
  for(i in 1:sim){
    df.p100.3[[i]] = sim_p(cor,sd,beta,res,popint,popslope,popdisp,N,repm) }
  saveRDS(df.p100.3, "df_p100_3.RDS")

#compile Stan model
mod = rstan::stan_model(file="gaus_p.stan")

#MCMC settings
n_iter <- 2000
n_warm <- 1000
n_chains <- 4

#dataframe for results
p100.3_pw = data.frame(iter = seq(1:sim), pb1 = NA, pb2 = NA) 

#estimate models
for(i in 1:sim){
  df = df.p100.3[[i]]
  est_mod <- rstan::sampling(mod, data=df, init="0", iter=n_iter, warmup=n_warm, seed = i, 
                    chains=n_chains, cores=n_chains, control=list(adapt_delta=0.99, max_treedepth=10))
  #posterior samples
  post = rstan::extract(est_mod)

  #selection gradients
  betas = data.frame(post[grepl("betas", names(post)) ] )
  p100.3_pw[i,"pb1"] = sum(sign(df$betas[1])==sign(betas[,1]))/length(betas[,1])
  p100.3_pw[i,"pb2"] = sum(sign(df$betas[2])==sign(betas[,2]))/length(betas[,2])
  saveRDS(p100.3_pw,"p100_3pw.RDS")
  }

##########################################
#Personality only, N 200, 3
##########################################
I = 200; repm = 3; sim = 200

df.p200.3 = list()
for(i in 1:sim){
  df.p200.3[[i]] = sim_p(cor,sd,beta,res,popint,popslope,popdisp,N,repm) }
saveRDS(df.p200.3, "df_p200_3.RDS")

#compile Stan model
mod = rstan::stan_model(file="gaus_p.stan")

#MCMC settings
n_iter <- 2000
n_warm <- 1000
n_chains <- 4

#dataframe for results
p200.3_pw = data.frame(iter = seq(1:sim), pb1 = NA, pb2 = NA) 

#estimate models
for(i in 1:sim){
  df = df.p200.3[[i]]
  est_mod <- rstan::sampling(mod, data=df, init="0", iter=n_iter, warmup=n_warm, seed = i, 
                             chains=n_chains, cores=n_chains, control=list(adapt_delta=0.99, max_treedepth=10))
  #posterior samples
  post = rstan::extract(est_mod)
  
  #selection gradients
  betas = data.frame(post[grepl("betas", names(post)) ] )
  p200.3_pw[i,"pb1"] = sum(sign(df$betas[1])==sign(betas[,1]))/length(betas[,1])
  p200.3_pw[i,"pb2"] = sum(sign(df$betas[2])==sign(betas[,2]))/length(betas[,2])
  saveRDS(p200.3_pw,"p200_3pw.RDS")
}

apply(p200.3_pw, 2, function(x) sum(x>0.95)/length(x)) 

##########################################
#Personality only, N 300, 3
##########################################
I = 300; repm = 3; sim = 200

df.p300.3 = list()
for(i in 1:sim){
  df.p300.3[[i]] = sim_p(cor,sd,beta,res,popint,popslope,popdisp,N,repm) }
saveRDS(df.p300.3, "df_p300_3.RDS")

#compile Stan model
mod = rstan::stan_model(file="gaus_p.stan")

#MCMC settings
n_iter <- 2000
n_warm <- 1000
n_chains <- 4

#dataframe for results
p300.3_pw = data.frame(iter = seq(1:sim), pb1 = NA, pb2 = NA) 

#estimate models
for(i in 1:sim){
  df = df.p300.3[[i]]
  est_mod <- rstan::sampling(mod, data=df, init="0", iter=n_iter, warmup=n_warm, seed = i, 
                             chains=n_chains, cores=n_chains, control=list(adapt_delta=0.99, max_treedepth=10))
  #posterior samples
  post = rstan::extract(est_mod)
  
  #selection gradients
  betas = data.frame(post[grepl("betas", names(post)) ] )
  p300.3_pw[i,"pb1"] = sum(sign(df$betas[1])==sign(betas[,1]))/length(betas[,1])
  p300.3_pw[i,"pb2"] = sum(sign(df$betas[2])==sign(betas[,2]))/length(betas[,2])
  saveRDS(p300.3_pw,"p300_3pw.RDS")
}

apply(p300.3_pw, 2, function(x) sum(x>0.95)/length(x)) 

##########################################
#Personality only, N 400, 3
##########################################
I = 400; repm = 3; sim = 200

df.p400.3 = list()
for(i in 1:sim){
  df.p400.3[[i]] = sim_p(cor,sd,beta,res,popint,popslope,popdisp,N,repm) }
saveRDS(df.p400.3, "df_p400_3.RDS")

#compile Stan model
mod = rstan::stan_model(file="gaus_p.stan")

#MCMC settings
n_iter <- 2000
n_warm <- 1000
n_chains <- 4

#dataframe for results
p400.3_pw = data.frame(iter = seq(1:sim), pb1 = NA, pb2 = NA) 

#estimate models
for(i in 1:sim){
  df = df.p400.3[[i]]
  est_mod <- rstan::sampling(mod, data=df, init="0", iter=n_iter, warmup=n_warm, seed = i, 
                             chains=n_chains, cores=n_chains, control=list(adapt_delta=0.99, max_treedepth=10))
  #posterior samples
  post = rstan::extract(est_mod)
  
  #selection gradients
  betas = data.frame(post[grepl("betas", names(post)) ] )
  p400.3_pw[i,"pb1"] = sum(sign(df$betas[1])==sign(betas[,1]))/length(betas[,1])
  p400.3_pw[i,"pb2"] = sum(sign(df$betas[2])==sign(betas[,2]))/length(betas[,2])
  saveRDS(p400.3_pw,"p400_3pw.RDS")
}

apply(p400.3_pw, 2, function(x) sum(x>0.95)/length(x)) 

##########################################
#Personality only, N 500, 3
##########################################
I = 500; repm = 3; sim = 200

df.p500.3 = list()
for(i in 1:sim){
  df.p500.3[[i]] = sim_p(cor,sd,beta,res,popint,popslope,popdisp,N,repm) }
saveRDS(df.p500.3, "df_p500_3.RDS")

#compile Stan model
mod = rstan::stan_model(file="gaus_p.stan")

#MCMC settings
n_iter <- 2000
n_warm <- 1000
n_chains <- 4

#dataframe for results
p500.3_pw = data.frame(iter = seq(1:sim), pb1 = NA, pb2 = NA) 

#estimate models
for(i in 1:sim){
  df = df.p500.3[[i]]
  est_mod <- rstan::sampling(mod, data=df, init="0", iter=n_iter, warmup=n_warm, seed = i, 
                             chains=n_chains, cores=n_chains, control=list(adapt_delta=0.99, max_treedepth=10))
  #posterior samples
  post = rstan::extract(est_mod)
  
  #selection gradients
  betas = data.frame(post[grepl("betas", names(post)) ] )
  p500.3_pw[i,"pb1"] = sum(sign(df$betas[1])==sign(betas[,1]))/length(betas[,1])
  p500.3_pw[i,"pb2"] = sum(sign(df$betas[2])==sign(betas[,2]))/length(betas[,2])
  saveRDS(p500.3_pw,"p500_3pw.RDS")
}

apply(p500.3_pw, 2, function(x) sum(x>0.95)/length(x)) 


#################################################################
#################################################################
#################################################################
#################################################################


#...

#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#Personality + plasticity, N 100, 3
#################################################################
I = 100; repm = 3; sim = 200

df.pp100.3 = list()
for(i in 1:sim){
  df.pp100.3[[i]] = sim_pp(cor,sd,beta,res,popint,popslope,popdisp,N,repm) }
saveRDS(df.pp100.3, "df_pp100_3.RDS")

#compile Stan model
mod = rstan::stan_model(file="gaus_pp.stan")

#MCMC settings
n_iter <- 2000
n_warm <- 1000
n_chains <- 4

#dataframe for results
pp100.3_pw = data.frame(iter = seq(1:sim), pb1 = NA, pb2 = NA, pb3 = NA, pb4 = NA, pb5 = NA) 

#estimate models
for(i in 1:sim){
  df = df.pp100.3[[i]]
  est_mod <- rstan::sampling(mod, data=df, init="0", iter=n_iter, warmup=n_warm, seed = i, 
                             chains=n_chains, cores=n_chains, control=list(adapt_delta=0.99, max_treedepth=10))
  #posterior samples
  post = rstan::extract(est_mod)
  
  #selection gradients
  betas = data.frame(post[grepl("betas", names(post)) ] )
  pp100.3_pw[i,"pb1"] = sum(sign(df$betas[1])==sign(betas[,1]))/length(betas[,1])
  pp100.3_pw[i,"pb2"] = sum(sign(df$betas[2])==sign(betas[,2]))/length(betas[,2])
  pp100.3_pw[i,"pb3"] = sum(sign(df$betas[3])==sign(betas[,3]))/length(betas[,3])
  pp100.3_pw[i,"pb4"] = sum(sign(df$betas[4])==sign(betas[,4]))/length(betas[,4])
  pp100.3_pw[i,"pb5"] = sum(sign(df$betas[5])==sign(betas[,5]))/length(betas[,5])
  saveRDS(pp100.3_pw,"pp100_3pw.RDS")
}

apply(pp100.3_pw, 2, function(x) sum(x>0.95)/length(x)) 

#################################################################
#Personality + plasticity, N 500, 3
#################################################################
I = 500; repm = 3; sim = 200

df.pp500.3 = list()
for(i in 1:sim){
  df.pp500.3[[i]] = sim_pp(cor,sd,beta,res,popint,popslope,popdisp,N,repm) }
saveRDS(df.pp500.3, "df_pp500_3.RDS")

#compile Stan model
mod = rstan::stan_model(file="gaus_pp.stan")

#MCMC settings
n_iter <- 2000
n_warm <- 1000
n_chains <- 4

#dataframe for results
pp500.3_pw = data.frame(iter = seq(1:sim), pb1 = NA, pb2 = NA, pb3 = NA, pb4 = NA, pb5 = NA) 

#estimate models
for(i in 1:sim){
  df = df.pp500.3[[i]]
  est_mod <- rstan::sampling(mod, data=df, init="0", iter=n_iter, warmup=n_warm, seed = i, 
                             chains=n_chains, cores=n_chains, control=list(adapt_delta=0.99, max_treedepth=10))
  #posterior samples
  post = rstan::extract(est_mod)
  
  #selection gradients
  betas = data.frame(post[grepl("betas", names(post)) ] )
  pp500.3_pw[i,"pb1"] = sum(sign(df$betas[1])==sign(betas[,1]))/length(betas[,1])
  pp500.3_pw[i,"pb2"] = sum(sign(df$betas[2])==sign(betas[,2]))/length(betas[,2])
  pp500.3_pw[i,"pb3"] = sum(sign(df$betas[3])==sign(betas[,3]))/length(betas[,3])
  pp500.3_pw[i,"pb4"] = sum(sign(df$betas[4])==sign(betas[,4]))/length(betas[,4])
  pp500.3_pw[i,"pb5"] = sum(sign(df$betas[5])==sign(betas[,5]))/length(betas[,5])
  saveRDS(pp500.3_pw,"pp500_3pw.RDS")
}

apply(pp500.3_pw, 2, function(x) sum(x>0.95)/length(x)) 

#################################################################
#Personality + plasticity, N 700, 3
#################################################################
I = 700; repm = 3; sim = 200

df.pp700.3 = list()
for(i in 1:sim){
  df.pp700.3[[i]] = sim_pp(cor,sd,beta,res,popint,popslope,popdisp,N,repm) }
saveRDS(df.pp700.3, "df_pp700.3.RDS")

#compile Stan model
mod = rstan::stan_model(file="gaus_pp.stan")

#MCMC settings
n_iter <- 2000
n_warm <- 1000
n_chains <- 4

#dataframe for results
pp700.3_pw = data.frame(iter = seq(1:sim), pb1 = NA, pb2 = NA, pb3 = NA, pb4 = NA, pb5 = NA) 

#estimate models
for(i in 1:sim){
  df = df.pp700.3[[i]]
  est_mod <- rstan::sampling(mod, data=df, init="0", iter=n_iter, warmup=n_warm, seed = i, 
                             chains=n_chains, cores=n_chains, control=list(adapt_delta=0.99, max_treedepth=10))
  #posterior samples
  post = rstan::extract(est_mod)
  
  #selection gradients
  betas = data.frame(post[grepl("betas", names(post)) ] )
  pp700.3_pw[i,"pb1"] = sum(sign(df$betas[1])==sign(betas[,1]))/length(betas[,1])
  pp700.3_pw[i,"pb2"] = sum(sign(df$betas[2])==sign(betas[,2]))/length(betas[,2])
  pp700.3_pw[i,"pb3"] = sum(sign(df$betas[3])==sign(betas[,3]))/length(betas[,3])
  pp700.3_pw[i,"pb4"] = sum(sign(df$betas[4])==sign(betas[,4]))/length(betas[,4])
  pp700.3_pw[i,"pb5"] = sum(sign(df$betas[5])==sign(betas[,5]))/length(betas[,5])
  saveRDS(pp700.3_pw,"pp700_3pw.RDS")
}

apply(pp700.3_pw, 2, function(x) sum(x>0.95)/length(x)) 


#################################################################
#################################################################
#################################################################
#################################################################

#################################################################
#Personality + plasticity, N 500, 6
#################################################################
I = 500; repm = 6; sim = 200

df.pp500.6 = list()
for(i in 1:sim){
  df.pp500.6[[i]] = sim_pp(cor,sd,beta,res,popint,popslope,popdisp,N,repm) }
saveRDS(df.pp500.6, "df_pp500_6.RDS")

#compile Stan model
mod = rstan::stan_model(file="gaus_pp.stan")

#MCMC settings
n_iter <- 2000
n_warm <- 1000
n_chains <- 4

#dataframe for results
pp500.6_pw = data.frame(iter = seq(1:sim), pb1 = NA, pb2 = NA, pb3 = NA, pb4 = NA, pb5 = NA) 

#estimate models
for(i in 1:sim){
  df = df.pp500.6[[i]]
  est_mod <- rstan::sampling(mod, data=df, init="0", iter=n_iter, warmup=n_warm, seed = i, 
                             chains=n_chains, cores=n_chains, control=list(adapt_delta=0.99, max_treedepth=10))
  #posterior samples
  post = rstan::extract(est_mod)
  
  #selection gradients
  betas = data.frame(post[grepl("betas", names(post)) ] )
  pp500.6_pw[i,"pb1"] = sum(sign(df$betas[1])==sign(betas[,1]))/length(betas[,1])
  pp500.6_pw[i,"pb2"] = sum(sign(df$betas[2])==sign(betas[,2]))/length(betas[,2])
  pp500.6_pw[i,"pb3"] = sum(sign(df$betas[3])==sign(betas[,3]))/length(betas[,3])
  pp500.6_pw[i,"pb4"] = sum(sign(df$betas[4])==sign(betas[,4]))/length(betas[,4])
  pp500.6_pw[i,"pb5"] = sum(sign(df$betas[5])==sign(betas[,5]))/length(betas[,5])
  saveRDS(pp500.6_pw,"pp500_6pw.RDS")
}

apply(pp500.6_pw, 2, function(x) sum(x>0.95)/length(x)) 

#################################################################
#Personality + plasticity, N 600, 6
#################################################################
I = 600; repm = 6; sim = 200

df.pp600.6 = list()
for(i in 1:sim){
  df.pp600.6[[i]] = sim_pp(cor,sd,beta,res,popint,popslope,popdisp,N,repm) }
saveRDS(df.pp600.6, "df_pp600_6.RDS")

#compile Stan model
mod = rstan::stan_model(file="gaus_pp.stan")

#MCMC settings
n_iter <- 2000
n_warm <- 1000
n_chains <- 4

#dataframe for results
pp600.6_pw = data.frame(iter = seq(1:sim), pb1 = NA, pb2 = NA, pb3 = NA, pb4 = NA, pb5 = NA) 

#estimate models
for(i in 1:sim){
  df = df.pp600.6[[i]]
  est_mod <- rstan::sampling(mod, data=df, init="0", iter=n_iter, warmup=n_warm, seed = i, 
                             chains=n_chains, cores=n_chains, control=list(adapt_delta=0.99, max_treedepth=10))
  #posterior samples
  post = rstan::extract(est_mod)
  
  #selection gradients
  betas = data.frame(post[grepl("betas", names(post)) ] )
  pp600.6_pw[i,"pb1"] = sum(sign(df$betas[1])==sign(betas[,1]))/length(betas[,1])
  pp600.6_pw[i,"pb2"] = sum(sign(df$betas[2])==sign(betas[,2]))/length(betas[,2])
  pp600.6_pw[i,"pb3"] = sum(sign(df$betas[3])==sign(betas[,3]))/length(betas[,3])
  pp600.6_pw[i,"pb4"] = sum(sign(df$betas[4])==sign(betas[,4]))/length(betas[,4])
  pp600.6_pw[i,"pb5"] = sum(sign(df$betas[5])==sign(betas[,5]))/length(betas[,5])
  saveRDS(pp600.6_pw,"pp600_6pw.RDS")
}

apply(pp600.6_pw, 2, function(x) sum(x>0.95)/length(x)) 

#################################################################
#Personality + plasticity, N 700, 6
#################################################################
I = 700; repm = 6; sim = 200

df.pp700.6 = list()
for(i in 1:sim){
  df.pp700.6[[i]] = sim_pp(cor,sd,beta,res,popint,popslope,popdisp,N,repm) }
saveRDS(df.pp700.6, "df_pp700.6.RDS")

#compile Stan model
mod = rstan::stan_model(file="gaus_pp.stan")

#MCMC settings
n_iter <- 2000
n_warm <- 1000
n_chains <- 4

#dataframe for results
pp700.6_pw = data.frame(iter = seq(1:sim), pb1 = NA, pb2 = NA, pb3 = NA, pb4 = NA, pb5 = NA) 

#estimate models
for(i in 1:sim){
  df = df.pp700.6[[i]]
  est_mod <- rstan::sampling(mod, data=df, init="0", iter=n_iter, warmup=n_warm, seed = i, 
                             chains=n_chains, cores=n_chains, control=list(adapt_delta=0.99, max_treedepth=10))
  #posterior samples
  post = rstan::extract(est_mod)
  
  #selection gradients
  betas = data.frame(post[grepl("betas", names(post)) ] )
  pp700.6_pw[i,"pb1"] = sum(sign(df$betas[1])==sign(betas[,1]))/length(betas[,1])
  pp700.6_pw[i,"pb2"] = sum(sign(df$betas[2])==sign(betas[,2]))/length(betas[,2])
  pp700.6_pw[i,"pb3"] = sum(sign(df$betas[3])==sign(betas[,3]))/length(betas[,3])
  pp700.6_pw[i,"pb4"] = sum(sign(df$betas[4])==sign(betas[,4]))/length(betas[,4])
  pp700.6_pw[i,"pb5"] = sum(sign(df$betas[5])==sign(betas[,5]))/length(betas[,5])
  saveRDS(pp700.6_pw,"pp700_6pw.RDS")
}

apply(pp700.6_pw, 2, function(x) sum(x>0.95)/length(x)) 

#################################################################
#################################################################
#################################################################
#################################################################


#...


#################################################################
#PPP, N 500, 6
#################################################################
I = 500; repm = 6; sim = 200

df.ppp500.6 = list()
for(i in 1:sim){
  df.ppp500.6[[i]] = sim_ppp(cor,sd,beta,res,popint,popslope,popdisp,N,repm) }
saveRDS(df.ppp500.6, "df_ppp500.6.RDS")

#compile Stan model
mod = rstan::stan_model(file="gaus_ppp.stan")

#MCMC settings
n_iter <- 2000
n_warm <- 1000
n_chains <- 4

#dataframe for results
ppp500.6_pw = data.frame(iter = seq(1:sim), pb1 = NA, pb2 = NA, pb3 = NA, pb4 = NA, pb5 = NA) 

#estimate models
for(i in 1:sim){
  df = df.ppp500.6[[i]]
  est_mod <- rstan::sampling(mod, data=df, init="0", iter=n_iter, warmup=n_warm, seed = i, 
                             chains=n_chains, cores=n_chains, control=list(adapt_delta=0.99, max_treedepth=10))
  #posterior samples
  post = rstan::extract(est_mod)
  
  #selection gradients
  betas = data.frame(post[grepl("betas", names(post)) ] )
  ppp500.6_pw[i,"pb1"] = sum(sign(df$betas[1])==sign(betas[,1]))/length(betas[,1])
  ppp500.6_pw[i,"pb2"] = sum(sign(df$betas[2])==sign(betas[,2]))/length(betas[,2])
  ppp500.6_pw[i,"pb3"] = sum(sign(df$betas[3])==sign(betas[,3]))/length(betas[,3])
  ppp500.6_pw[i,"pb4"] = sum(sign(df$betas[4])==sign(betas[,4]))/length(betas[,4])
  ppp500.6_pw[i,"pb5"] = sum(sign(df$betas[5])==sign(betas[,5]))/length(betas[,5])
  saveRDS(ppp500.6_pw,"ppp500.6pw.RDS")
}

apply(ppp500.6_pw, 2, function(x) sum(x>0.95)/length(x)) 

#################################################################
#PPP, N 600, 6
#################################################################
I = 600; repm = 6; sim = 200

df.ppp600.6 = list()
for(i in 1:sim){
  df.ppp600.6[[i]] = sim_ppp(cor,sd,beta,res,popint,popslope,popdisp,N,repm) }
saveRDS(df.ppp600.6, "df_ppp600.6.RDS")

#compile Stan model
mod = rstan::stan_model(file="gaus_ppp.stan")

#MCMC settings
n_iter <- 2000
n_warm <- 1000
n_chains <- 4

#dataframe for results
ppp600.6_pw = data.frame(iter = seq(1:sim), pb1 = NA, pb2 = NA, pb3 = NA, pb4 = NA, pb5 = NA) 

#estimate models
for(i in 1:sim){
  df = df.ppp600.6[[i]]
  est_mod <- rstan::sampling(mod, data=df, init="0", iter=n_iter, warmup=n_warm, seed = i, 
                             chains=n_chains, cores=n_chains, control=list(adapt_delta=0.99, max_treedepth=10))
  #posterior samples
  post = rstan::extract(est_mod)
  
  #selection gradients
  betas = data.frame(post[grepl("betas", names(post)) ] )
  ppp600.6_pw[i,"pb1"] = sum(sign(df$betas[1])==sign(betas[,1]))/length(betas[,1])
  ppp600.6_pw[i,"pb2"] = sum(sign(df$betas[2])==sign(betas[,2]))/length(betas[,2])
  ppp600.6_pw[i,"pb3"] = sum(sign(df$betas[3])==sign(betas[,3]))/length(betas[,3])
  ppp600.6_pw[i,"pb4"] = sum(sign(df$betas[4])==sign(betas[,4]))/length(betas[,4])
  ppp600.6_pw[i,"pb5"] = sum(sign(df$betas[5])==sign(betas[,5]))/length(betas[,5])
  saveRDS(ppp600.6_pw,"ppp600.6pw.RDS")
}

apply(ppp600.6_pw, 2, function(x) sum(x>0.95)/length(x)) 

#################################################################
#################################################################
#################################################################
#################################################################

#plot


p100.3_pw = readRDS("p100_3pw.RDS")
p100.6_pw = readRDS("p100_6pw.RDS")
#p100.9_pw = readRDS("p100_9pw.RDS")

#summary 
  data.frame(I = 100, rep = c(3,6), type = "p",
           pb1 = c(apply(p100.3_pw, 2, function(x) sum(x>0.95)/length(x) )[2],
                   apply(p100.6_pw, 2, function(x) sum(x>0.95)/length(x) )[2]),
                   #apply(p100.9_pw, 2, function(x) sum(x>0.95)/length(x) )[2]),
           pb2 = c(apply(p100.3_pw, 2, function(x) sum(x>0.95)/length(x) )[3],
                   apply(p100.6_pw, 2, function(x) sum(x>0.95)/length(x) )[3]))
                   #apply(p100.9_pw, 2, function(x) sum(x>0.95)/length(x) )[3]))
      
p100_pw =rbind(p100.3_pw,p100.6_pw) #continuous
p100_pw$rep = c(rep(3,nrow(p100.3_pw)), rep(6,nrow(p100.6_pw))) #, rep(9,nrow(p100.9_pw)))


library(tidyr)
lp100 = gather(p100_pw, param, prob, pb1:pb2, factor_key=FALSE)
lp100d = gather(p100_d, param, prob, pb1:pb2, factor_key=FALSE)

library(ggplot2); library(tidybayes)
#pp.plot = 
  ggplot(lp100, aes(x = prob, y = param, group = rep, color = as.factor(rep), fill = as.factor(rep)))+
  stat_dotsinterval(.width=c(0.5,0.9), point_interval = median_qi, show_slab=FALSE,
                    position = position_dodge(-0.1)) +
  xlab("\n Posterior Probability")+
  ylab("\n\n\nParameter")+
  geom_vline(xintercept=c(0.95,1), linetype="dashed")+
  coord_cartesian(xlim=c(0.6,1))+
  scale_x_continuous(breaks = c(seq(0.6,1,by=0.1)))+
  theme(axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_blank(),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=12),
        axis.line = element_line(size = 1),
        panel.border=element_rect(fill=NA,color="black", size=1, 
                                  linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))#+
  #guides(fill=FALSE, color=FALSE)

