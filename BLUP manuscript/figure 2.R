#Figure 2

library(ggplot2)
library(grid)
library(cowplot)
library(reshape)
library(RSA)


setwd("C:/Temp/BLUPs/BLUP manuscript")

#################################################################
#################################################################
#Top panel: 3 plots single trait effects  (stabilizing, disruptive, and balancing)
#################################################################
#################################################################

#stabilizing
z_pm = seq(-2,2,by=0.1)
w_st = 1 + 0*z_pm -0.3*(z_pm^2) 
plot(w_st ~ z_pm)

df_st = data.frame(z_pm, w_st)

#plot
st_p =
ggplot(df_st, aes(x = z_pm, y = w_st , color = w_st)) +
  geom_line(size = 2)+
  scale_colour_gradient( low="red", high="#f5b907" )+
  coord_cartesian(xlim=c(-2, 2), ylim = c(min(w_st) - 0.1, max(w_st) + 0.1)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-2,0,2),
                     labels = c(-2,0,2) )+
  scale_y_continuous(expand = c(0, 0), breaks = c(0,1),
                     labels = c(0,1) )+
  geom_hline(yintercept=1,linetype="dashed", alpha = 0.5)+
  geom_vline(xintercept=0,linetype="dashed", alpha = 0.5)+
  xlab(expression(paste(italic(z[pm]))))+
  ylab(expression(paste(italic(w))))+
  #ggtitle("Stabilizing selection\n")+
  theme(plot.title =element_text(size=12, hjust=0.5),
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
        plot.margin = unit(c(0.5,0.5,1,0.5), "cm"))+
        guides(fill=FALSE, color=FALSE)

###########################################################

#disruptive
z_pm = seq(-2,2,by=0.1)
w_dr = 1 + 0*z_pm + 0.3*(z_pm^2) 
plot(w_dr ~ z_pm)

df_dr = data.frame(z_pm, w_dr)

#plot
dr_p =
ggplot(df_dr, aes(x = z_pm, y = w_dr , color = w_dr)) +
  geom_line(size = 2)+
  scale_colour_gradient( low="red", high="#f5b907" )+
  coord_cartesian(xlim=c(-2, 2), ylim = c(min(w_dr) - 0.1, max(w_dr) + 0.1)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-2,0,2),
                     labels = c(-2,0,2) )+
  scale_y_continuous(expand = c(0, 0), breaks = c(1,2),
                     labels = c(1,2) )+
  geom_hline(yintercept=1,linetype="dashed", alpha = 0.5)+
  geom_vline(xintercept=0,linetype="dashed", alpha = 0.5)+
  xlab(expression(paste(italic(z[pm]))))+
  ylab(expression(paste(italic(w))))+
  #ggtitle("Disruptive selection\n")+
  theme(plot.title =element_text(size=12, hjust=0.5),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12,angle=0,vjust=0.5, color=NA),
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
        plot.margin = unit(c(0.5,0.5,1,0.5), "cm"))+
        guides(fill=FALSE, color=FALSE)


###########################################################

#balancing
z_pm = rep(seq(-2,2,by=0.1), each = 2)
state = rep(c(-1,1), length(seq(-2,2,by=0.1)) )
w_bl = 1 + 0*z_pm +0.3*(z_pm*state)

df_bl = data.frame(z_pm, w_bl,
                   state = rep(c("Context 1", "Context 2"), length(seq(-2,2,by=0.1)) ))

#plot
bl_p =
ggplot(df_bl, aes(x = z_pm, y = w_bl , color = w_bl, group = state)) +
  geom_line(size = 2)+
  scale_colour_gradient( low="red", high="#f5b907" )+
  coord_cartesian(xlim=c(-2, 2), ylim = c(min(w_bl) - 0.1, max(w_bl) + 0.1))+ 
  scale_x_continuous(expand = c(0, 0), breaks = c(-2,0,2),
                     labels = c(-2,0,2) )+
  scale_y_continuous(expand = c(0, 0), breaks = c(0,1,2),
                     labels = c(0,1,2) )+
  geom_hline(yintercept=1,linetype="dashed", alpha = 0.5)+
  geom_vline(xintercept=0,linetype="dashed", alpha = 0.5)+
  xlab(expression(paste(italic(z[pm]))))+
  ylab(expression(paste(italic(w))))+
  #ggtitle("Balancing selection\n")+
  theme(plot.title =element_text(size=12, hjust=0.5),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12,angle=0,vjust=0.5, color=NA),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.line = element_line(size = 1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12),
        panel.spacing = unit(1.5, "lines"),
        panel.border=element_rect(fill=NA,color="black", size=1, 
                                  linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.5,0.5,1,0.5), "cm"))+
        guides(fill=FALSE, color=FALSE)

#add text grobs for conditions

grob1 <- grobTree(textGrob("State 1", x= 0.12,  y = 0.25, hjust=0, rot = 27,
             gp=gpar(col="#131d87", fontsize=12, fontface="bold")))

grob2 <- grobTree(textGrob("State 2", x=0.7,  y = 0.41, hjust=0, rot = -26,
             gp=gpar(col="#131d87", fontsize=12, fontface="bold")))
bl_p2=
bl_p + annotation_custom(grob1) +  annotation_custom(grob2)


###########################################################
#combine plots
pt=plot_grid(st_p,dr_p,bl_p2, ncol=3)

#add left side title
title_st <- ggdraw() + 
  draw_label(
    expression(paste(bold("One trait"))),
    fontface = 'bold',  x = 0.5, y = 0.5, size = 12)

pt2 = plot_grid(title_st, pt, rel_heights=c(0.1,1), ncol = 1)

ggsave("fig2a.tiff",pt2, width = 10, height = 3, units = "in", dpi = 300, compression = "lzw")

#################################################################
#################################################################
#Bottom panel: 3 plots two trait effects  (bowl, saddle, dome)
#################################################################
#################################################################

dome=
  plotRSA(
    x =  -0.2,
    y =  -0.2,
   x2 = -0.2,
   y2 = -0.2,
   xy =  0.05,
   b0 = 1,
   type = "3d",
   ylab = expression(paste(italic(z[pm]))),
   xlab = expression(paste(italic(z[pn]))),
   zlab = expression(paste(italic(w))),
   label.rotation = list(z = 0),
   suppress.grid = TRUE,
   axes = NA,
   pal = colorRampPalette(c("red","red", "#ffce3d"))(50),
   legend = FALSE, coef = FALSE, param = FALSE,
   pad = -3.5, cex.tickLabel = 0.75
   )

bowl=
  plotRSA(
    x =  0.025,
    y = -0.025,
   x2 = 0.1,
   y2 = 0.1,
   xy = 0.1,
   b0 = 1,
   type = "3d",
   ylab = expression(paste(italic(z[pm]))),
   xlab = expression(paste(italic(z[pn]))),
   #zlab = expression(paste(italic(w))),
   zlab = " ",
   label.rotation = list(z = 0),
   suppress.grid = TRUE,
   axes = NA,
   pal = colorRampPalette(c("red","#ffce3d"))(50),
   legend = FALSE, coef = FALSE, param = FALSE,
   pad = -3.5, cex.tickLabel = 0.75
   )

saddle=
  plotRSA(
    x =  -0.2,
    y =  -0.1,
   x2 =  0.1,
   y2 =  -0.3,
   xy =  0.15,
   b0 = 1,
   type = "3d",
   ylab = expression(paste(italic(z[pm]))),
   xlab = expression(paste(italic(z[pn]))),
   #zlab = expression(paste(italic(w))),
   zlab = " ",
   label.rotation = list(z = 0),
   suppress.grid = TRUE,
   axes = NA,
   pal = colorRampPalette(c("red","#ffce3d"))(50),
   legend = FALSE, coef = FALSE, param = FALSE,
   pad = -3.5, cex.tickLabel = 0.75
   )

#################################################################
library(cowplot)

p1 = ggdraw(dome)
p2 = ggdraw(bowl)
p3 = ggdraw(saddle)

title1 <- ggdraw() + 
  draw_label("Dome", x = 0.5, y = 0.45, size = 12) +
  theme(plot.margin = margin(0, 0, 0, 0))

title2 <- ggdraw() + 
  draw_label("Bowl", x = 0.5, y = 0.45, size = 12) +
  theme(plot.margin = margin(0, 0, 0, 0))

title3 <- ggdraw() + 
  draw_label("Saddle", x = 0.5, y = 0.45, size = 12) +
    theme(plot.margin = margin(0, 0, 0, 0))

pcomb = plot_grid(p1, p2, p3, ncol = 3, scale = 0.9)
pcomb2 = plot_grid(title1,title2,title3, ncol = 3)
pcomb3 = plot_grid(pcomb2, pcomb, ncol = 1, rel_heights = c (0.1,1))
ggsave("fig2b.tiff", pcomb3, width = 11, height = 3.5, units = "in",
       dpi = 300)

#add left side title
title_mt <- ggdraw() + 
  draw_label(
    expression(paste(bold("Two traits"))),
    fontface = 'bold',  x = 0.5, y = 0.5, size = 12)

pb2 = plot_grid(title_mt, pcomb, ncol = 1, rel_heights = c(0.01,1), scale=1)


#add in top plots
pfig2 = plot_grid(pt2, pb2, ncol = 1 , rel_heights = c(0.7,1))
ggsave("fig2.tiff", pfig2, width = 11, height = 7, units = "in",
       dpi = 300, compression = "lzw")
ggsave("fig2.pdf", pfig2, width = 11, height = 7, units = "in",
       dpi = 300)

