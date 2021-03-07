#load packages
library(rstan)

## increase memory to avoid crashing with loop
memory.limit(size=100000)

#Stan settings
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#set directory
setwd("C:/Temp/BLUPs/BLUP manuscript")

#compile SBC model
mod = stan_model(file="SBC stan.stan")

#MCMC settings
n_iter <- 1255
n_warm <- 1000
n_thin <- 5
n_chains <- 4

#simulation settings
I = 100
rep = 5
N = I*rep
ind = rep(1:I, each = rep) 

sim = 3 #number of simulations
params = 9 #number of parameters

#initialize matrix for holding selection parameter ranks
ranks = matrix(NA, nrow = sim, ncol = params)

#progress bar
pb = txtProgressBar(min = 0, max = sim, initial = 0) 

#run simulation
for(i in 1:sim){
  #estimate model
  est <- sampling(mod, data=list(I = I, N = N, ind = ind), init="0", iter=n_iter, warmup=n_warm, 
                    chains=n_chains, cores=n_chains, thin = n_thin, control=list(adapt_delta=0.99, max_treedepth=15))
  #posterior samples
  post = extract(est)
  ranks[i,] = colSums(post$lt_sim)
  
  #set progress
  setTxtProgressBar(pb,i)
}

#binning process
M = ((n_iter - n_warm)*n_chains)/n_thin #possible ranks
J = sim/20
bcount = (M + 1)/J #expected count per bin



