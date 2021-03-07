#load packages
library(mvtnorm)
library(rstan)
library(shinystan)

## increase memory to avoid crashing with loop
memory.limit(size=100000)

#set directory
setwd("C:/Temp/BLUPs/BLUP manuscript")

#Stan settings
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

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
  #cors = sample(c(cor, -1*cor), size=3*3, replace=TRUE) #random sign
  #R = matrix(cors, nrow=3, ncol=3 )
  #R[lower.tri(R)] = t(R)[lower.tri(R)] #force symmetric
  #diag(R) = 1
  #S = matrix( c(sd,0,0,0,sd,0,0,0,sd), nrow=3, ncol=3 )
  #P = S %*% R %*% S
  #z_p = rmvnorm(I, mean = rep(0,3), sigma = P) 
  z_p = rnorm(I, mean = 0, sd = sd) 
  personality = z_p
  #plasticity = z_p[,2]
  #predictability = z_p[,3]
  
  #environmental covariate (z-score)
  x = rnorm(I*repm, 0, 1)
  
  #behavioral response model
  ind = rep(1:I, each = repm) #index of repeated individual measures
  popslope = sample(c(popslope, -1*popslope),size=1) #random sign
  
  z_mu = popint + personality[ind] #+ (popslope + plasticity[ind])*x
  #z_sigma = log(popdisp) + predictability[ind]
  #z = rnorm(I*repm, mean = z_mu, sd = exp(z_sigma) )
  z = rnorm(I*repm, mean = z_mu, sd = popdisp )
  
  #beta coefficients
  betas = sample(c(beta, -1*beta), size=3, replace=TRUE)
  
  #fitness response model
  w_mu = 1 + cbind(personality) * betas[1] +
             betas[2]*(personality^2) #betas[5]*(plasticity^2)+betas[6]*(predictability^2)+
             #betas[7]*(personality*plasticity)+betas[8]*(personality*predictability)+
             #betas[9]*(plasticity*predictability)
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
  #predictability = z_p[,3]
  
  #environmental covariate (z-score)
  x = rnorm(I*repm, 0, 1)
  
  #behavioral response model
  ind = rep(1:I, each = repm) #index of repeated individual measures
  popslope = sample(c(popslope, -1*popslope),size=1) #random sign
  
  z_mu = popint + personality[ind] + (popslope + plasticity[ind])*x
  #z_sigma = popdisp + predictability[ind]
  z = rnorm(I*repm, mean = z_mu, sd = popdisp )
  
  #beta coefficients
  betas = sample(c(beta, -1*beta), size=5, replace=TRUE)
  
  #fitness response model
  w_mu = 1 + cbind(personality,plasticity) %*% betas[1:2] +
             betas[3]*(personality^2) + betas[4]*(plasticity^2)#+betas[6]*(predictability^2)+
             betas[5]*(personality*plasticity)#+betas[8]*(personality*predictability)+
             #betas[9]*(plasticity*predictability)
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
  popslope = sample(c(popslope, -1*popslope),size=1) #random sign
  
  z_mu = popint + personality[ind] + (popslope + plasticity[ind])*x
  z_sigma = log(popdisp) + predictability[ind]
  z = rnorm(I*repm, mean = z_mu, sd = exp(z_sigma) )
  
  #beta coefficients
  betas = sample(c(beta, -1*beta), size=9, replace=TRUE)
  
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
  for(i in 1:sim){df.p100.3[[i]] = 
    sim_p(cor,sd,beta,res,popint,popslope,popdisp,N,repm) }
  saveRDS(df.p100.3, "df_p100_3.RDS")

#compile Stan model
mod = stan_model(file="gaus_p.stan")

#MCMC settings
n_iter <- 2000
n_warm <- 1000
n_chains <- 4

#estimate model
est_mod <- sampling(mod, data=stan_data, init="0", iter=n_iter, warmup=n_warm, seed = 9, 
                    chains=n_chains, cores=n_chains, control=list(adapt_delta=0.99, max_treedepth=10))
#posterior samples
post = extract(est_mod)

#selection gradients
betas = data.frame(post[grepl("betas", names(post)) ] )
apply(betas,2,FUN = median)
apply(betas,2,FUN = mad)
(apply(betas,2,FUN = median) - stan_data$betas ) / stan_data$betas 
abs(apply(betas,2,FUN = function(x) mad(x)/ median(x) ))
apply(betas,2,FUN = function(x) sum(sign(median(x))==sign(x))/length(x) )

#GUI overview of model
launch_shinystan(est_mod)






