#Figure 3

library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

setwd("C:/Temp/BLUPs/BLUP manuscript")

#################################################################
#################################################################

#RN parameter covariance matrix
R = matrix( c(1, 0.3, -0.1, 0.3, 1, -0.4, -0.1, -0.4, 1), 3, 3 )
S = diag( c(sqrt(0.5), sqrt(0.4), sqrt(0.3) ) )
P = S %*% R %*% S

#fitness function
betas = c(-0.3, -0.1, 0.4, 0.1, 0.3, - 0.2, 0.3, -0.2, -0.3)

f = function(m,n,o) {
      w = 1 + betas[1]*m + betas[2]*n + betas[3]*o +
              betas[4]*(m*m) + betas[5]*(n*n) + betas[6]*(o*o) +
              betas[7]*(m*n) + betas[8]*(m*o) + betas[9]*(n*o) 
      return(w) } 

f(-2,-2,-2)
f(-1,-1,-1)
f(0,0,0)
f(1,1,1)
f(2,2,2)

#mean change
beta = matrix(betas[1:3], ncol = 1)
delta_mean = P %*% beta
deltaind_mean = diag(diag(P)) %*% beta

#change in P
#make gamma matrix
gamma = diag( c(2*betas[4:6] ) )
gamma[lower.tri(gamma)] = betas[7:9]
gamma[upper.tri(gamma)] = gamma[lower.tri(gamma)]

delta_p = P %*% (gamma - beta %*% t(beta)) %*% P
deltaind_p = diag(diag(P)) %*% (gamma - beta %*% t(beta)) %*% diag(diag(P))

#new within-generation means, variances, and covariances
RN_mean = c(0,0,0) + delta_mean #centered on zero prior to selection analysis
RNind_mean = c(0,0,0) + deltaind_mean #centered on zero prior to selection analysis

RN_P = P + delta_p
RNind_P = P + deltaind_p
RN_S = sqrt(diag(RN_P)) #standard deviations
RNind_S = sqrt(diag(RNind_P)) #standard deviations
delta_S = diag(RN_S) - S
deltaind_S = diag(RNind_S) - S

RN_R = cov2cor(RN_P) #correlations
RNind_R = cov2cor(RNind_P) #correlations
delta_R = RN_R - R
deltaind_R = RNind_R - R

#compare effect of integration on means, SDs and correlations
round(c(delta_mean),2); round(c(deltaind_mean),2)
round(diag(delta_S),2); round(diag(deltaind_S),2)
round(delta_R[lower.tri(delta_R)],2);round(deltaind_R[lower.tri(deltaind_R)],2)

#################################################################
#################################################################
#left panel: 3 plots of selection effects on RN parameters
#################################################################
#################################################################

#uncertainty for generating posterior distributions
pers_uc = c(0.5, 0.6, 0.5) #means, standard deviations, correlations
pers_induc = c(0.6, 0.5, 0.4)
plast_uc = c(0.7, 0.8, 0.4)
plast_induc = c(1.5, 0.6, 0.6)
pred_uc = c(0.6, 1, 0.7)
pred_induc = c(0.7, 0.5, 0.6)

#################################################################
#mean change plot

#personality 
pers_as = rnorm(1e5, delta_mean[1], abs(delta_mean[1]*pers_uc[1]) )
pers_asind = rnorm(1e5, deltaind_mean[1], abs(deltaind_mean[1]*pers_induc[1]) )

#plasticity
plast_as = rnorm(1e5, delta_mean[2], abs(delta_mean[2]*plast_uc[1]) )
plast_asind = rnorm(1e5, deltaind_mean[2], abs(deltaind_mean[2]*plast_induc[1]) )

#predictability
pred_as = rnorm(1e5, delta_mean[3], abs(delta_mean[3]*pred_uc[1]) )
pred_asind = rnorm(1e5, deltaind_mean[3], abs(deltaind_mean[3]*pred_induc[1]) )

#long format dfs
pers.df = data.frame(value = c(pers_as, pers_asind),
                     type = as.factor(rep(c("total change","direct change"),
                                          each = length(pers_as) )))
plast.df = data.frame(value = c(plast_as, plast_asind),
                     type = as.factor(rep(c("total change","direct change"),
                                          each = length(pers_as) )))
pred.df = data.frame(value = c(pred_as, pred_asind),
                     type = as.factor(rep(c("total change","direct change"),
                                          each = length(pers_as) )))

df = rbind(pers.df,plast.df,pred.df)
df$param = rep(c("Personality","Plasticity","Predictability"), each = nrow(pers.df))

spnames <- list("Personality"= expression(paste("Personality (", mu, ")")),
              "Plasticity" = expression(paste("Plasticity (", beta, ")")),
              "Predictability" = expression(paste("Predictability (", theta, ")")) )
spname_labeller <- function(variable,value){ return(spnames[value]) }

#plot
mean_p = 
ggplot(df, aes( x = value, color = type, fill = type)) + 
  geom_density(aes(y=..scaled..), size = 0.9, alpha = 0.10) + 
  facet_grid(. ~ param, labeller = spname_labeller)+
  coord_cartesian(xlim = c(-0.5,0.5))+
  scale_y_continuous(expand=c(0,0.03))+
  scale_color_manual(values=c("#4db4d1","#4dd191"))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5)+
  #geom_vline(data=df[df$param=="Personality" & df$type=="total change",],
  #           aes(xintercept=mean(value)), colour="#4dd191", size = 1) + 
  #geom_vline(data=df[df$param=="Personality" & df$type=="direct change",],
  #           aes(xintercept=mean(value)), colour="#4db4d1", size = 1) + 
  #geom_vline(data=df[df$param=="Plasticity" & df$type=="total change",],
  #           aes(xintercept=mean(value)), colour="#4dd191", size = 1) + 
  #geom_vline(data=df[df$param=="Plasticity" & df$type=="direct change",],
  #           aes(xintercept=mean(value)), colour="#4db4d1", size = 1) + 
  xlab( bquote(atop("", paste(Delta,bold(bar(z))[p] ))))+
  scale_fill_manual(values=c("#4dd191","#4db4d1"))+
    theme(plot.title =element_text(size=12, hjust=0.5),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=10),
        axis.text.y=element_blank(),
        axis.line = element_line(size = 1),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        panel.spacing = unit(1.5, "lines"),
        panel.border=element_rect(fill=NA,color="black", size=1, 
                                  linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,2), "cm"))+
        guides(color = FALSE, fill = FALSE, size = FALSE)

#add posterior probabilities
pp_t <- data.frame(
  label = round(
            c(sum(sign(median(pers_as))==sign(pers_as))/length(pers_as),
            sum(sign(median(plast_as))==sign(plast_as))/length(plast_as),
            sum(sign(median(pred_as))==sign(pred_as))/length(pred_as) ), 2),
  type = "total change",
  param = c("Personality","Plasticity","Predictability") )

pp_d <- data.frame(
  label = round(
            c(sum(sign(median(pers_asind))==sign(pers_asind))/length(pers_asind),
            sum(sign(median(plast_asind))==sign(plast_asind))/length(plast_asind),
            sum(sign(median(pred_asind))==sign(pred_asind))/length(pred_asind) ), 2),
  type = "direct change",
  param = c("Personality","Plasticity","Predictability") )

mean_p + geom_text(data = pp_t, mapping = aes(x = c(0.4,0.4,0.4), y = c(0.93,0.93,0.93),
                                               size=12, label = label)) +
         geom_text(data = pp_d, mapping = aes(x = c(0.4,0.4,0.4), y = c(0.77,0.77,0.77),
                                               size=12, label = label)) 
#add manual axis
grid.text(expression(paste("Pr(",Delta[T],")")), x = unit(0.07, "npc"), y = unit(0.78, "npc"),
          gp=gpar(fontsize = 12, fontface = "bold", col = "#4dd191"))
grid.text(expression(paste("Pr(",Delta[D], ")" )), x = unit(0.07, "npc"), y = unit(0.68, "npc"),
          gp=gpar(fontsize = 12, fontface = "bold", col = "#4db4d1"))
meanplot = grid.grab()
meanplot = arrangeGrob(meanplot)

#ggsave("fig3a.tiff", meanplot, width = 7, height = 2.5, units = "in",
#       compression = "lzw", dpi=600)

#################################################################
#SD change plot

delta_S = diag(delta_S)
deltaind_S = diag(deltaind_S)

#personality 
pers_as = rnorm(1e5, delta_S[1], abs(delta_S[1]*pers_uc[2]) )
pers_asind = rnorm(1e5, deltaind_S[1], abs(deltaind_S[1]*pers_induc[2]) )

#plasticity
plast_as = rnorm(1e5, delta_S[2], abs(delta_S[2]*plast_uc[2]) )
plast_asind = rnorm(1e5, deltaind_S[2], abs(deltaind_S[2]*plast_induc[2]) )

#predictability
pred_as = rnorm(1e5, delta_S[3], abs(delta_S[3]*pred_uc[2])*2 )
pred_asind = rnorm(1e5, deltaind_S[3], abs(deltaind_S[3]*pred_induc[2]) )

#long format dfs
pers.df = data.frame(value = c(pers_as, pers_asind),
                     type = as.factor(rep(c("total change","direct change"),
                                          each = length(pers_as) )))
plast.df = data.frame(value = c(plast_as, plast_asind),
                     type = as.factor(rep(c("total change","direct change"),
                                          each = length(pers_as) )))
pred.df = data.frame(value = c(pred_as, pred_asind),
                     type = as.factor(rep(c("total change","direct change"),
                                          each = length(pers_as) )))

df = rbind(pers.df,plast.df,pred.df)
df$param = rep(c("Personality","Plasticity","Predictability"), each = nrow(pers.df))

spnames <- list("Personality"= expression(paste("Var(", mu,")")),
              "Plasticity" = expression(paste("Var(", beta,")")),
              "Predictability" = expression(paste("Var(", theta,")")))
spname_labeller <- function(variable,value){ return(spnames[value]) }

#plot
SD_p = 
ggplot(df, aes( x = value, color = type, fill = type)) + 
  geom_density(aes(y=..scaled..), size = 1, alpha = 0.10) + 
  facet_grid(. ~ param, labeller = spname_labeller)+
  coord_cartesian(xlim = c(-0.3,0.3))+
  scale_x_continuous(breaks=c(-0.2,-0.1, 0, 0.1,0.2),
                     labels = c("-0.1", "-0.05", "0.0", "0.05", "0.1"))+
  scale_y_continuous(expand=c(0,0.03))+
  scale_color_manual(values=c("#4db4d1","#4dd191"))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5)+
  xlab( bquote(atop("",paste(Delta,bold(V)[z[p]] ))))+
  scale_fill_manual(values=c("#4dd191","#4db4d1"))+
    theme(plot.title =element_text(size=12, hjust=0.5),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=10),
        axis.text.y=element_blank(),
        axis.line = element_line(size = 1),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        panel.spacing = unit(1.5, "lines"),
        panel.border=element_rect(fill=NA,color="black", size=1, 
                                  linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,2), "cm"))+
        guides(color = FALSE, fill = FALSE, size = FALSE)

#add posterior probabilities
pp_t <- data.frame(
  label = format(round(
            c(sum(sign(median(pers_as))==sign(pers_as))/length(pers_as),
            sum(sign(median(plast_as))==sign(plast_as))/length(plast_as),
            sum(sign(median(pred_as))==sign(pred_as))/length(pred_as) ), 2), nsmall=2),
  type = "total change",
  param = c("Personality","Plasticity","Predictability") )

pp_d <- data.frame(
  label = format(round(
          c(sum(sign(median(pers_asind))==sign(pers_asind))/length(pers_asind),
            sum(sign(median(plast_asind))==sign(plast_asind))/length(plast_asind),
            sum(sign(median(pred_asind))==sign(pred_asind))/length(pred_asind)),2), nsmall=2),
  type = "direct change",
  param = c("Personality","Plasticity","Predictability") )

SD_p + geom_text(data = pp_t, mapping = aes(x = c(-0.20,-0.20,-0.20), y = c(0.93,0.93,0.93),
                                               size=12, label = label)) +
         geom_text(data = pp_d, mapping = aes(x = c(-0.20,-0.20,-0.20), y = c(0.77,0.77,0.77),
                                               size=12, label = label)) 
#add manual axis
grid.text(expression(paste("Pr(",Delta[T],")")), x = unit(0.07, "npc"), y = unit(0.78, "npc"),
          gp=gpar(fontsize = 12, fontface = "bold", col = "#4dd191"))
grid.text(expression(paste("Pr(",Delta[D], ")" )), x = unit(0.07, "npc"), y = unit(0.68, "npc"),
          gp=gpar(fontsize = 12, fontface = "bold", col = "#4db4d1"))

SDplot = grid.grab()
SDplot = arrangeGrob(SDplot)

#ggsave("fig3b.tiff", SDplot, width = 7, height = 2.5, units = "in",
#       compression = "lzw", dpi=600)

#################################################################
#R change plot

delta_R = delta_R[lower.tri(delta_R)]
deltaind_R = deltaind_R[lower.tri(deltaind_R)]

#personality 
pers_as = rnorm(1e5, delta_R[1], abs(delta_R[1]*pers_uc[3]) )
pers_asind = rnorm(1e5, deltaind_R[1], abs(deltaind_R[1]*pers_induc[3]) )

#plasticity
plast_as = rnorm(1e5, delta_R[2], abs(delta_R[2]*plast_uc[3]) )
plast_asind = rnorm(1e5, deltaind_R[2], abs(deltaind_R[2]*plast_induc[3]) )

#predictability
pred_as = rnorm(1e5, delta_R[3], abs(delta_R[3]*pred_uc[3]) )
pred_asind = rnorm(1e5, deltaind_R[3], abs(deltaind_R[3]*pred_induc[3]) )

#long format dfs
pers.df = data.frame(value = c(pers_as, pers_asind),
                     type = as.factor(rep(c("total change","direct change"),
                                          each = length(pers_as) )))
plast.df = data.frame(value = c(plast_as, plast_asind),
                     type = as.factor(rep(c("total change","direct change"),
                                          each = length(pers_as) )))
pred.df = data.frame(value = c(pred_as, pred_asind),
                     type = as.factor(rep(c("total change","direct change"),
                                          each = length(pers_as) )))

df = rbind(pers.df,plast.df,pred.df)
df$param = rep(c("Personality","Plasticity","Predictability"), each = nrow(pers.df))

spnames <- list("Personality"= expression(paste("Cor(", mu,",",beta,")")),
              "Plasticity" = expression(paste("Cor(", mu,",",theta,")")),
              "Predictability" = expression(paste("Cor(", beta,",",theta,")")))
spname_labeller <- function(variable,value){ return(spnames[value]) }

#plot
R_p = 
ggplot(df, aes( x = value, color = type, fill = type)) + 
  geom_density(aes(y=..scaled..), size = 1, alpha = 0.10) +
  facet_grid(. ~ param, labeller = spname_labeller)+
  coord_cartesian(xlim = c(-0.25,0.25))+
  scale_y_continuous(expand=c(0,0.03))+
  scale_color_manual(values=c("#4db4d1","#4dd191"))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5)+
  xlab( bquote(atop("",paste(Delta,bold(R)[z[p]] ))))+
  scale_fill_manual(values=c("#4dd191","#4db4d1"))+
    theme(plot.title =element_text(size=12, hjust=0.5),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=10),
        axis.text.y=element_blank(),
        axis.line = element_line(size = 1),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold"),
        panel.spacing = unit(1.5, "lines"),
        panel.border=element_rect(fill=NA,color="black", size=1, 
                                  linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,2), "cm"))+
        guides(color = FALSE, fill = FALSE, size = FALSE)

#add posterior probabilities
pp_t <- data.frame(
  label = round(
            c(sum(sign(median(pers_as))==sign(pers_as))/length(pers_as),
            sum(sign(median(plast_as))==sign(plast_as))/length(plast_as),
            sum(sign(median(pred_as))==sign(pred_as))/length(pred_as) ), 2),
  type = "total change",
  param = c("Personality","Plasticity","Predictability") )

pp_d <- data.frame(
  label = c(
            round(c(sum(sign(median(pers_asind))==sign(pers_asind))/length(pers_asind),
            sum(sign(median(plast_asind))==sign(plast_asind))/length(plast_asind)),2),
            "0.90" ),
  type = "direct change",
  param = c("Personality","Plasticity","Predictability") )

R_p + geom_text(data = pp_t, mapping = aes(x = c(-0.20, 0.20, 0.20), y = c(0.93,0.93,0.93),
                                               size=12, label = label)) +
         geom_text(data = pp_d, mapping = aes(x = c(-0.20, 0.20, 0.20), y = c(0.77,0.77,0.77),
                                               size=12, label = label)) 
#add manual axis
grid.text(expression(paste("Pr(",Delta[T],")")), x = unit(0.07, "npc"), y = unit(0.78, "npc"),
          gp=gpar(fontsize = 12, fontface = "bold", col = "#4dd191"))
grid.text(expression(paste("Pr(",Delta[D], ")" )), x = unit(0.07, "npc"), y = unit(0.68, "npc"),
          gp=gpar(fontsize = 12, fontface = "bold", col = "#4db4d1"))

Rplot = grid.grab()
Rplot = arrangeGrob(Rplot)

#ggsave("fig3c.tiff", Rplot, width = 7, height = 2.5, units = "in",
#       compression = "lzw", dpi=600)


p3t=plot_grid(meanplot,SDplot,Rplot, ncol=1)

#ggsave("fig3left.tiff", p3t, width = 7, height = 7.5, units = "in",
#       dpi = 600, compression = "lzw")

#ggsave("fig3left.pdf", p3t, width = 7, height = 7.5, units = "in",
#       dpi = 300)

#################################################################
#################################################################
#Right panel: full RN norm
#################################################################
#################################################################

c(delta_mean); c(deltaind_mean)
delta_S; deltaind_S
delta_R; deltaind_R

#population parameters
x = seq(-1, 1, by = 0.1)
mu = 0.2; beta = 0.3; theta = 0.1

mu2 = mu + delta_mean[1]
mu2ind = mu + deltaind_mean[1]
beta2 = beta + delta_mean[2]
beta2ind = beta + deltaind_mean[2]
theta2 = theta + delta_mean[3]
theta2ind = theta + deltaind_mean[3]

#compute reaction norm CIs
y_low = mu + beta*x + 1.96*theta
y_high = mu + beta*x - 1.96*theta
y2_low = mu2 + beta2*x + 1.96*theta2
y2_high = mu2 + beta2*x - 1.96*theta2
y2ind_low = mu2ind + beta2ind*x + 1.96*theta2ind
y2ind_high = mu2ind + beta2ind*x - 1.96*theta2ind

#total RN
popBRN_pt =
ggplot() +
  coord_cartesian(xlim=c(-1, 1), ylim=c(-1.05,1.05)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-1,0,1),
                     labels = c(-1,0,1) )+
  scale_y_continuous(expand = c(0, 0), breaks = c(-1,0,1),
                     labels = c(-1,0,1) ) +
  geom_hline(yintercept=0,linetype="dashed", alpha = 0.25)+
  geom_vline(xintercept=0,linetype="dashed", alpha = 0.25) +
  
  geom_abline(intercept = mu, slope = beta, size = 2, alpha =0.75, color = "darkgrey") +
  geom_ribbon(aes(x = x, y = mu + beta*x, ymin = y_low, ymax = y_high),
              size = 0.8, linetype = 5, alpha = 0.15, color = "darkgrey", fill = "darkgrey")+
  
  geom_abline(intercept = mu2, slope = beta2, size = 2, color = "#23a666") +
  geom_ribbon(aes(x = x, y = mu2 + beta2*x, ymin = y2_low, ymax = y2_high),
              size = 0.8, linetype = 5, alpha = 0.15, color = "#23a666", fill = "#4dd191")+
  xlab("\nEnvironment")+
  ylab("Behavior")+
  ggtitle(expression(paste(bold(Delta[T]),bold(" Population-Level Behavioral RN"))))+
  theme(plot.title =element_text(size=12,face="bold", hjust=0.5),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.line = element_line(size = 1),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        panel.border=element_rect(fill=NA,color="black", size=1, 
                                  linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.5,1,0.5,0.25), "cm"))+
        guides(fill=FALSE, color=FALSE)

ggsave("fig3lt.tiff", popBRN_pt, width = 6, height = 5.5, units = "in", 
                     dpi = 300, compression = "lzw")


#direct RN
popBRN_pd =
ggplot() +
  coord_cartesian(xlim=c(-1, 1), ylim=c(-1.05,1.05)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-1,0,1),
                     labels = c(-1,0,1) )+
  scale_y_continuous(expand = c(0, 0), breaks = c(-1,0,1),
                     labels = c(-1,0,1) ) +
  geom_hline(yintercept=0,linetype="dashed", alpha = 0.25)+
  geom_vline(xintercept=0,linetype="dashed", alpha = 0.25) +
  
  geom_abline(intercept = mu, slope = beta, size = 2, alpha =0.75, color = "darkgrey") +
  geom_ribbon(aes(x = x, y = mu + beta*x, ymin = y_low, ymax = y_high),
              size = 0.8, linetype = 5, alpha = 0.15, color = "darkgrey", fill = "darkgrey")+
  
  geom_abline(intercept = mu2ind, slope = beta2ind, size = 2, color = "#0586ab") +
  geom_ribbon(aes(x = x, y = mu2ind + beta2ind*x, ymin = y2ind_low, ymax = y2ind_high),
              size = 0.8, linetype = 5, alpha = 0.15, color = "#0586ab", fill = "#4db4d1")+
  xlab("\nEnvironment")+
  ylab("Behavior")+
  ggtitle(expression(paste(bold(Delta[D]),bold(" Population-Level Behavioral RN"))))+
  theme(plot.title =element_text(size=12,face="bold", hjust=0.5),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.line = element_line(size = 1),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        panel.border=element_rect(fill=NA,color="black", size=1, 
                                  linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.5,1,0.5,0.25), "cm"))+
        guides(fill=FALSE, color=FALSE)

ggsave("fig3ld.tiff", popBRN_pd, width = 6, height = 5.5, units = "in", 
                     dpi = 300, compression = "lzw")


#combine right panel
rightp=plot_grid(popBRN_pt,popBRN_pd, ncol=1)

#combine with left panel
fullp=plot_grid(p3t, rightp, ncol=2, rel_widths = c(1,0.6))

ggsave("fig3_2.tiff", fullp, width = 10, height = 7, units = "in",
       dpi = 300, compression = "lzw")
ggsave("fig3_2.pdf", fullp, width = 10, height = 7, units = "in",
       dpi = 300)










