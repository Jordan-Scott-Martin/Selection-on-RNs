library(brms)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#set directory
setwd("C:/Temp/BLUPs/BLUP manuscript")

## increase memory to avoid crashing with loop
memory.limit(size=100000)

n.ind = 50 #sample size
rep = 2 #repeated measurements
ind = rep(1:n.ind, each = rep) #ind index
V.ind = 0.3  #blup variance
V.res = 0.7 #residual (R = 0.3)
beta = 0 #BLUP effect

#structure of sim
blup = rnorm(n.ind,0, sqrt(V.ind))
y = 0 + blup[ind] + rnorm(length(ind),0, sqrt(V.res))
w = 0 + beta*blup + rnorm(n.ind,0,sqrt(V.res))

#fit initial models
ymod = brm(y ~ (1|ind), data = data.frame(y, ind))
blup_pe = ranef(ymod)$ind[,1,1]
wmod = brm(w ~ blup_pe, data = data.frame(w, blup_pe))

#test false positive rate / 100
n.sim = 100

#matrices for results
ci90 = matrix(NA, nrow = n.sim, ncol = 2)
ci95 = matrix(NA, nrow = n.sim, ncol = 2)

#run simulation
for(i in 1:n.sim){

#random data
blup = rnorm(n.ind,0, sqrt(V.ind))
y = 0 + blup[ind] + rnorm(length(ind),0, sqrt(V.res))
w = 0 + beta*blup + rnorm(n.ind,0,sqrt(V.res))

#run model
ymod1 = update(ymod, newdata = data.frame(y, ind), seed = i )
blup_pe = ranef(ymod1)$ind[,1,1]
wmod1 = update(wmod, newdata = data.frame(w, blup_pe), seed = i )
  
#results
ci90[i,1:2] = fixef(wmod1, prob = c(0.05,0.95))["blup_pe",3:4]
ci95[i,1:2] = fixef(wmod1, prob = c(0.025,0.975))["blup_pe",3:4]
print(i)
gc()
}

saveRDS(ci90,"ci90_50n.RDS")
saveRDS(ci95,"ci95_50n.RDS")

#fpr_90_50n
sum(ci90[,1] > 0 & ci90[,2] > 0) +
sum(ci90[,1] < 0 & ci90[,2] < 0) 

#fpr_95_50n
sum(ci95[,1] > 0 & ci95[,2] > 0) +
sum(ci95[,1] < 0 & ci95[,2] < 0) 

#######################################

n.ind = 200 #sample size

#matrices for results
ci90 = matrix(NA, nrow = n.sim, ncol = 2)
ci95 = matrix(NA, nrow = n.sim, ncol = 2)

#run simulation
for(i in 1:n.sim){

#random data
blup = rnorm(n.ind,0, sqrt(V.ind))
y = 0 + blup[ind] + rnorm(length(ind),0, sqrt(V.res))
w = 0 + beta*blup + rnorm(n.ind,0,sqrt(V.res))

#run model
ymod1 = update(ymod, newdata = data.frame(y, ind), seed = i )
blup_pe = ranef(ymod1)$ind[,1,1]
wmod1 = update(wmod, newdata = data.frame(w, blup_pe), seed = i )
  
#results
ci90[i,1:2] = fixef(wmod1, prob = c(0.05,0.95))["blup_pe",3:4]
ci95[i,1:2] = fixef(wmod1, prob = c(0.025,0.975))["blup_pe",3:4]
print(i)
gc()
}

saveRDS(ci90,"ci90_100n.RDS")
saveRDS(ci95,"ci95_100n.RDS")

#fpr_90_50n
sum(ci90[,1] > 0 & ci90[,2] > 0) +
sum(ci90[,1] < 0 & ci90[,2] < 0) 

#fpr_95_50n
sum(ci95[,1] > 0 & ci95[,2] > 0) +
sum(ci95[,1] < 0 & ci95[,2] < 0) 

#######################################

n.ind = 50 #sample size
V.ind = 1  #blup variance
V.res = 1 #residual (R = 0.3)

#matrices for results
ci90 = matrix(NA, nrow = n.sim, ncol = 2)
ci95 = matrix(NA, nrow = n.sim, ncol = 2)

#run simulation
for(i in 1:n.sim){

#random data
blup = rnorm(n.ind,0, sqrt(V.ind))
y = 0 + blup[ind] + rnorm(length(ind),0, sqrt(V.res))
w = 0 + beta*blup + rnorm(n.ind,0,sqrt(V.res))

#run model
ymod1 = update(ymod, newdata = data.frame(y, ind), seed = i )
blup_pe = ranef(ymod1)$ind[,1,1]
wmod1 = update(wmod, newdata = data.frame(w, blup_pe), seed = i )
  
#results
ci90[i,1:2] = fixef(wmod1, prob = c(0.05,0.95))["blup_pe",3:4]
ci95[i,1:2] = fixef(wmod1, prob = c(0.025,0.975))["blup_pe",3:4]
print(i)
gc()
}

saveRDS(ci90,"ci90_100n.RDS")
saveRDS(ci95,"ci95_100n.RDS")

#fpr_90_50n
sum(ci90[,1] > 0 & ci90[,2] > 0) +
sum(ci90[,1] < 0 & ci90[,2] < 0) 

#fpr_95_50n
sum(ci95[,1] > 0 & ci95[,2] > 0) +
sum(ci95[,1] < 0 & ci95[,2] < 0) 

