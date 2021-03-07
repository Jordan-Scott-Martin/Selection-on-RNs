library(rstan)

setwd("C:/Temp/BLUPs/BLUP manuscript")
ranks = readRDS("ranks3_blup.RDS")
ranks = na.omit(ranks)

hist(ranks)
ranks = ranks[1:250,]
cutranks = as.character(cut(ranks, seq(-1,410, by=51)))
table(cutranks)
chisq.test(table(cutranks))

#binning process
M = 400 #possible ranks
bins = 20
J = 260/bins
bcount = (M + 1)/J #expected count per bin

hist(ranks, breaks = 10)

#plot
library(reshape2); library(ggplot2)

colnames(ranks) = paste0("beta",seq(1:9))
lranks=melt(ranks, id.vars = "param")

ggplot(lranks, aes(x = value) )+ geom_density()+
  #geom_histogram(binwidth = 50)+ 
  scale_x_continuous(limits=c(-1,401),expand = c(0,0))#+
  #ylim(c(0,400))#+
  #facet_wrap(.~ Var2, ncol = 3)
