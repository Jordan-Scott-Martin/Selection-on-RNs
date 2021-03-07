#Figure 1

library(ggplot2)
library(grid)
library(pBrackets)

setwd("C:/Temp/BLUPs/BLUP manuscript")

#################################################################
#################################################################

#point estimate plot
x = seq(-1, 1, by = 0.1)
mu_j = -0.1
#mu_k = -0.2
beta_j = 0.37
#beta_k = -0.2
eta_j = 0.3
#eta_k = 0.3

y_j = mu_j + beta_j*x + rnorm(length(x), 0, eta_j)
y_jlow = mu_j + beta_j*x + 1.96*eta_j
y_jhigh = mu_j + beta_j*x - 1.96*eta_j

#y_k = mu_k + beta_k*x + rnorm(length(x), 0, eta_k)
#y_klow = mu_k + beta_k*x + 1.96*eta_k
#y_khigh = mu_k + beta_k*x - 1.96*eta_k

df = data.frame(y = y_j, x = x,
                low = y_jlow, high = y_jhigh,
                id = rep("j", each = length(x)))
saveRDS(df, "figa1df.RDS")
df = readRDS("figa1df.RDS")

ind_pointest =
ggplot(df, aes(x = x, y = y, color = id, fill = id) ) +
  coord_cartesian(xlim=c(-1, 1), ylim=c(-1.05,1.05)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-1,0,1),
                     labels = c(-1,0,1) )+
  scale_y_continuous(expand = c(0, 0), breaks = c(-1,0,1),
                     labels = c(-1,0,1) )+
  geom_hline(yintercept=0,linetype="dashed", alpha = 0.5)+
  geom_vline(xintercept=0,linetype="dashed", alpha = 0.5)+
  geom_ribbon(aes(ymin = low, ymax = high), size = 2, alpha = .15,
              color = NA, fill = "#7ed7f7")+
  geom_point(size = 4, alpha = 0.5, pch = 21, color = "black", fill = "#7ef7cd") +
  geom_abline(intercept = mu_j, slope = beta_j, size = 2, color = "#612da8") +
  geom_point(aes(x = 0, y = mu_j), pch = 21, size = 5, color = "black", fill="#eb007d") +
  scale_color_manual(values = c("#000000"))+
  scale_fill_manual(values = c("#87d0e0"))+
  xlab("\nEnvironment")+
  ylab("Behavior")+
  #ggtitle(expression(paste(bold("RN point estimates for individual "),italic("j"))))+
  theme(plot.title =element_text(size=12,face="bold", hjust=0.5),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12,angle=0,vjust=0.5),
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
        plot.margin = unit(c(0.5,1,0.5,0.5), "cm"))+
        guides(fill=FALSE, color=FALSE)

#add grobs

grobpers <- grobTree(textGrob(expression(paste(mu[j])), x=0.5,  y=0.39, hjust=0,
             gp=gpar(col="#c9026c", fontsize=25, fontface="bold")))

grobplast <- grobTree(textGrob(expression(paste(beta[j])), x=0.75,  y=0.62, hjust=0,
             gp=gpar(col="#7d00eb", fontsize=25, fontface="bold")))

grobpred <- grobTree(textGrob(expression(paste(theta[j])), x=0.10,  y=0.67, hjust=0,
             gp=gpar(col="#0098cf", fontsize=25, fontface="bold")))

#grobcurly <- bracketsGrob(0.1, 0.46, 0.1, 0.17, h=-0.05, lwd=4, col="#eb9500")

pa1=
ind_pointest + annotation_custom(grobpers) +  annotation_custom(grobplast) +
               annotation_custom(grobpred) #+ annotation_custom(grobcurly)

ggsave("figA1.tiff",pa1,width = 6, height = 5.5, units = "in", dpi = 300, compression = "lzw")

#################################################################
#################################################################

#uncertainty plots

pers_unc = rnorm(1e6, mu_j, 0.75)
plast_unc = rnorm(1e6, beta_j, 1)
pred_unc = rnorm(1e6, eta_j, 1.5)

figa2p1=
ggplot(data.frame(pers_unc), aes(x = pers_unc))+
  geom_density(color = "black", fill = "#c9026c", size = 1, alpha = 0.9)+
  geom_vline(xintercept = mu_j, size=1)+
  coord_cartesian(expand=FALSE, xlim=c(-3,3), ylim=c(0,1))+
  ggtitle(expression(paste("Personality", " (",mu[j], ")")))+
  #ggtitle(expression(paste(bold("Personality"), " (",bold(mu[j]), ")")))+
  theme(plot.title =element_text(size=12,face="bold", hjust=0.5),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=10),
        axis.text.y=element_blank(),
        axis.line = element_line(size = 1),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        panel.border=element_rect(fill=NA,color="black", size=1, 
                                  linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
        guides(fill=FALSE, color=FALSE)

figa2p2=
ggplot(data.frame(plast_unc), aes(x = plast_unc))+
  geom_density(color = "black", fill = "#7d00eb", size = 1, alpha = 0.9)+
  geom_vline(xintercept = beta_j, size=1)+
  coord_cartesian(expand=FALSE, xlim=c(-3,3), ylim=c(0,1))+
  ggtitle(expression(paste("Plasticity", " (",beta[j], ")")))+
  #ggtitle(expression(paste(bold("Plasticity"), " (",bold(beta[j]), ")")))+
  theme(plot.title =element_text(size=12,face="bold", hjust=0.5),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=10),
        axis.text.y=element_blank(),
        axis.line = element_line(size = 1),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        panel.border=element_rect(fill=NA,color="black", size=1, 
                                  linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
        guides(fill=FALSE, color=FALSE)

figa2p3=
ggplot(data.frame(pred_unc), aes(x = pred_unc))+
  geom_density(color = "black", fill = "#0098cf", size = 1, alpha = 0.9)+
  geom_vline(xintercept = eta_j, size=1)+
  coord_cartesian(expand=FALSE, xlim=c(-3,3), ylim=c(0,1))+
  ggtitle(expression(paste("Predictability", " (",theta[j], ")")))+
  #ggtitle(expression(paste(bold("Predictability"), " (",bold(theta[j]), ")")))+
  theme(plot.title =element_text(size=12,face="bold", hjust=0.5),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=10),
        axis.text.y=element_blank(),
        axis.line = element_line(size = 1),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        panel.border=element_rect(fill=NA,color="black", size=1, 
                                  linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
        guides(fill=FALSE, color=FALSE)

#################################################################

library(cowplot)

pa2=plot_grid(figa2p1,figa2p2,figa2p3, ncol=1)

title1 <- ggdraw() + 
  draw_label(
    expression(paste(bold("RN point estimates for individual "),italic("j"))),
    fontface = 'bold',  x = 0.58, y = 0.45, size = 12) +
    theme(plot.margin = margin(0, 0, 0, 0))

title2 <- ggdraw() + 
  draw_label(
    expression(paste(bold("RN uncertainty for individual "),italic("j"))),
    fontface = 'bold',  x = 0.5, y = 0.45, size = 12) +
    theme(plot.margin = margin(0, 0, 0, 0))

pa11=plot_grid(title1,pa1,ncol=1, rel_heights = c(0.1,1))
pa22=plot_grid(title2,pa2,ncol=1, rel_heights = c(0.1,1))

figa = plot_grid(pa11, pa22, rel_widths = c(1,0.4))
ggsave("fig1.tiff", figa, width = 8.5, height = 5.5, units = "in",
       dpi = 300, compression = "lzw")
ggsave("fig1.pdf", figa, width = 8.5, height = 5.5, units = "in",
       dpi = 300)
