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

#sim algorithm
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

