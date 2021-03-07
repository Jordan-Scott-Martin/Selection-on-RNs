require(arm)

n.inds<-100 ##Number of individuals
n.sim=100    ##Number of simulations
V.res.v<-c(0.01, seq(0.1,1, 0.1)) ##Within individual variances
V.res2<-1 ## residual variance in the blup-trait relationtship
V.ind<-0.3 ## among individual variance
B1<-0.5 ##BLUP-Trait relationtship

est.B1<-matrix(NA, n.sim, 11)
est.B2<-matrix(NA, n.sim, 11)


for(i in 1:11){
  V.res<-0.6
  
  for(j in 1:n.sim){
    I <- rnorm(n.inds, 0, sqrt(V.ind))
    y <- I*B1 + rnorm(n.inds,0, 0.5)
    sim.B1<-coef(lm(y~I))[2]
    
    d<-data.frame(ID = rep(1:n.inds, each=2), I = rep(I, each=2) )
    d$e<-rnorm(nrow(d), 0, sqrt(V.res) )
    d$x <- d$I + d$e
    
    mod1<-lmer(x ~ 1 + (1|ID), data = d)
    summary(mod1)
    
    smod<-sim(mod1, n.sim)
    blups.dist<-smod@ranef$ID
    blupm<-apply(blups.dist, 2, mean)
    blupsd<-apply(blups.dist, 2, sd)
    
    est.B1.tmp<-rep(NA, n.sim)
    
    for(k in 1:n.sim){
      m<-lm(y~blups.dist[k,,])
      est.B1.tmp[k]<-coef(m)[2] - sim.B1 }
    
    est.B1[j,i]<-mean(est.B1.tmp)
    est.B2[j,i]<-coef(lm(y~unlist(coefficients(mod1))))[2]-sim.B1
    
    blupm = unlist(coefficients(mod1))
    blupsd = sqrt( VarCorr(mod1)$ID[1] )
    
    d2 = data.frame(y, blupm, blupsd)
    mod2 = brm(y ~ me(blupm, sdx = blupsd), data = d2,
               prior = c(prior(normal(0,1), class = "Intercept"),
                         prior(normal(0,1), class = "b"),
                         prior(cauchy(0,1), class = "sigma")) )
    #est.B3[j,i]<-
      
    
  }
}

means<-apply(est.B1,2,mean)
upper<-apply(est.B1,2,quantile, 0.025)
lower<-apply(est.B1,2,quantile, 0.975)
Table_S4A <-cbind(V.res.v, means, lower, upper)
Table_S4A


means2<-apply(est.B2,2,mean)
upper2<-apply(est.B2,2,quantile, 0.025)
lower2<-apply(est.B2,2,quantile, 0.975)
Table_S4B<-cbind(V.res.v, means2, lower2, upper2)
Table_S4B


means3<-apply(est.B3,2,mean)
upper3<-apply(est.B3,2,quantile, 0.025)
lower3<-apply(est.B3,2,quantile, 0.975)
Table_S4C<-cbind(V.res.v, means3, lower3, upper3)
Table_S4C
     