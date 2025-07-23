############## 气态氮释放分析 ##############
library(readxl)
# setwd("D:/Workspace/book/test")
dat<-read.csv('D:/Workspace/book/Qingyuan/soil_data/Incubation/Gaseou N/kaihuang/gaseousN release_kai.csv',header=TRUE)
str(dat)

library(ggplot2)

ggplot(dat,aes(Temp,AC_N2,col=treatment))+
  geom_point(size=2)+
  facet_wrap(.~soillayer)

ggplot(dat,aes(Temp,N_N2O,col=treatment))+
  geom_point(size=2)+
  facet_wrap(plot~soillayer)

mythem<-theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1), 
              axis.text.x = element_text(colour = "black",size=12,angle = 0,hjust = .5,vjust =1),
              axis.ticks = element_line(linewidth = .5),
              axis.ticks.length = unit(2,"mm"),
              prism.ticks.length = unit(1,"mm"),
              axis.text.y = element_text(colour = "black",size = 12,hjust=.5), 
              axis.title.x =element_text(size=12), axis.title.y=element_text(colour = "black",size=12),
              legend.text = element_text(size=12), legend.title =element_text(size=12),
              legend.margin = margin(unit(c(0.5,1,0.5,1),'cm')),
              # legend.box.background = element_rect(fill ="white", size=0.6,colour = "black", linetype = "solid"),
              panel.background = element_rect(fill =NA, linewidth=0.6,colour = "black", linetype = "solid",
                                              inherit.blank=T),
              panel.grid=element_blank(),     #element_line(colour="gray",size=0.5),
              panel.grid.major = element_blank(),    #element_line(colour = "gray",size=0.6), 
              plot.background = element_rect(fill=NA,color=NA,linetype = "solid"),
              legend.key = element_rect(fill = NA,color=NA), 
              legend.background = element_rect(fill = NA), 
              plot.margin = unit(c(0.4,0.3,0.4,0.3), 'mm'),   #调整画图区域的间距，从上右下左调整
              strip.background = element_rect(fill = NA,color='black'), 
              strip.text = element_text(colour = "black",size = 12,hjust=.5),
              legend.position = "none")

################# N2O产生
####Oa+e
str(dat)
unique(dat$soillayer)
dat$soillayer<-ifelse(dat$soillayer=='Oa+e','O layer','0-10 cm')

dat<-within(dat,soillayer<-factor(dat$soillayer,levels = c('O layer','0-10 cm')))

names(dat)

{
  fit1_co <- summary(lm(log(P_N2O) ~ Temp,data=dat[dat$soillayer=='O layer' & dat$treatment=='control',]));fit1_co
  fit1_co$coefficients[[1]];fit1_co$coefficients[[2]]
  
  fit1_wo <- summary(lm(log(P_N2O) ~ Temp,data=dat[dat$soillayer=='O layer' & dat$treatment=='warmed',]));fit1_wo
  fit1_wo$coefficients[[1]];fit1_wo$coefficients[[2]]
  
  # 创建一个新数据框，用于绘制拟合曲线
  x_fit <- seq(5, 25, by = 0.01)
  y_fit.c <- exp(fit1_co$coefficients[[2]]*x_fit + fit1_co$coefficients[[1]])
  y_fit.w <- exp(fit1_wo$coefficients[[2]]*x_fit + fit1_wo$coefficients[[1]])
  
  fit_data <- data.frame(x = x_fit, control = y_fit.c, warmed = y_fit.w)
  
  library(ggpmisc);library(ggprism)
  p1_0<-ggplot(dat[dat$soillayer=='O layer',],
               aes(Temp,P_N2O,col=treatment))+
    geom_point(size=1,pch=1)+
    scale_color_manual(limits=c('control','warmed'),values = c('blue','red'),name=NULL)+
    geom_smooth(method = 'glm', 
                formula = y ~ x, 
                method.args = list(family = Gamma(link = "log")),  se=F,show.legend = T,linewidth=0.4)+
    facet_grid(.~soillayer,scales = 'free_y')+
    labs(x='Temperature (\u00B0C)',y=expression(paste('N'[2],'O production (ng N  g'^-1,' h'^-1~')')))+
    stat_poly_eq(aes(Temp,log(P_N2O),col=treatment,
                     label = paste(stat(rr.label),stat(p.value.label),sep="*\",\"~")),
                 formula =y~x,eq.x.rhs="x",coef.digits = 2,rr.digits = 2, 
                 eq.with.lhs = "italic(ln(y))~`=`~",hjust=0,  #给“y"换形式         
                 na.rm = FALSE,parse = TRUE,label.x.npc = 0.03,label.y.npc=c(0.90,0.80),size = 3)+
    scale_y_continuous(breaks=seq(0,6,2),limits=c(0,6),expand=c(0,0))+
    mythem+theme(legend.position = c(0.06,0.5));p1_0
  
  ####0-10cm
  fit1_co <- summary(lm(log(P_N2O) ~ Temp,data=dat[dat$soillayer=='0-10 cm' & dat$treatment=='control',]));fit1_co
  fit1_co$coefficients[[1]];fit1_co$coefficients[[2]]
  
  fit1_wo <- summary(lm(log(P_N2O) ~ Temp,data=dat[dat$soillayer=='0-10 cm' & dat$treatment=='warmed',]));fit1_wo
  fit1_wo$coefficients[[1]];fit1_wo$coefficients[[2]]
  
  # 创建一个新数据框，用于绘制拟合曲线
  x_fit <- seq(5, 25, by = 0.01)
  y_fit.c <- exp(fit1_co$coefficients[[2]]*x_fit + fit1_co$coefficients[[1]])
  y_fit.w <- exp(fit1_wo$coefficients[[2]]*x_fit + fit1_wo$coefficients[[1]])
  
  fit_data <- data.frame(x = x_fit, control = y_fit.c, warmed = y_fit.w)
  
  library(ggpmisc)
  p1_1<-ggplot(dat[dat$soillayer=='0-10 cm',],
               aes(Temp,P_N2O,col=treatment))+
    geom_point(size=1,pch=1)+
    scale_color_manual(limits=c('control','warmed'),values = c('blue','red'))+
    geom_smooth(method = 'glm', 
                formula = y ~ x, 
                method.args = list(family = Gamma(link = "log")),  se=F,show.legend = T,linewidth=0.4)+
    facet_grid(.~soillayer,scales = 'free_y')+
    labs(x='Temperature (\u00B0C)',y=expression(paste('N'[2],'O production (ng N  g'^-1,' h'^-1~')')))+
    stat_poly_eq(aes(Temp,log(P_N2O),col=treatment,
                     label = paste(stat(rr.label),stat(p.value.label),sep="*\",\"~")),
                 formula =y~x,eq.x.rhs="x",coef.digits = 3,rr.digits = 2, 
                 eq.with.lhs = "italic(ln(y))~`=`~",hjust=0,  #给“y"换形式         
                 na.rm = FALSE,parse = TRUE,label.x.npc = 0.03,label.y.npc=c(0.90,0.80),size = 3)+
    scale_y_continuous(breaks=seq(0,6,2),limits=c(0,6),expand=c(0,0))+
    mythem;p1_1
  
  
  ############# N2 产生
  ####Oa+e
  str(dat)
  fit1_co <- summary(lm(log(AC_N2) ~ Temp,data=dat[dat$soillayer=='O layer' &
                                                     dat$treatment=='control',]));fit1_co
  fit1_co$coefficients[[1]];fit1_co$coefficients[[2]]
  
  fit1_wo <- summary(lm(log(AC_N2) ~ Temp,data=dat[dat$soillayer=='O layer' & dat$treatment=='warmed',]));fit1_wo
  fit1_wo$coefficients[[1]];fit1_wo$coefficients[[2]]
  
  # 创建一个新数据框，用于绘制拟合曲线
  x_fit <- seq(5, 25, by = 0.01)
  y_fit.c <- exp(fit1_co$coefficients[[2]]*x_fit + fit1_co$coefficients[[1]])
  y_fit.w <- exp(fit1_wo$coefficients[[2]]*x_fit + fit1_wo$coefficients[[1]])
  
  fit_data <- data.frame(x = x_fit, control = y_fit.c, warmed = y_fit.w)
  
  library(ggpmisc);library(ggprism)
  p2_0<-ggplot(dat[dat$soillayer=='O layer',],
               aes(Temp,AC_N2,col=treatment))+
    geom_point(size=1,pch=1)+
    scale_color_manual(limits=c('control','warmed'),values = c('blue','red'))+
    geom_smooth(method = 'glm', 
                formula = y ~ x, 
                method.args = list(family = Gamma(link = "log")),  se=F,show.legend = T,linewidth=0.4)+
    facet_grid(.~soillayer,scales = 'free_y')+
    labs(x='Temperature (\u00B0C)',y=expression(paste('N'[2],' production (ng N  g'^-1,' h'^-1~')')))+
    stat_poly_eq(aes(Temp,log(AC_N2),col=treatment,
                     label = paste(stat(rr.label),stat(p.value.label),sep="*\",\"~")),
                 formula =y~x,eq.x.rhs="x",coef.digits = 2,rr.digits = 2, 
                 eq.with.lhs = "italic(ln(y))~`=`~",hjust=0,  #给“y"换形式         
                 na.rm = FALSE,parse = TRUE,label.x.npc = 0.03,label.y.npc=c(0.90,0.80),size = 3)+
    scale_y_continuous(breaks=seq(0,0.4,0.1),expand=c(0,0))+
    coord_cartesian(ylim=c(0,0.4),clip='on')+
    mythem;p2_0
  
  ####0-10cm
  fit1_co <- summary(lm(log(AC_N2) ~ Temp,data=dat[dat$soillayer=='0-10 cm' & dat$treatment=='control',]));fit1_co
  fit1_co$coefficients[[1]];fit1_co$coefficients[[2]]
  
  fit1_wo <- summary(lm(log(AC_N2) ~ Temp,data=dat[dat$soillayer=='0-10 cm' & dat$treatment=='warmed',]));fit1_wo
  fit1_wo$coefficients[[1]];fit1_wo$coefficients[[2]]
  
  # 创建一个新数据框，用于绘制拟合曲线
  x_fit <- seq(5, 25, by = 0.01)
  y_fit.c <- exp(fit1_co$coefficients[[2]]*x_fit + fit1_co$coefficients[[1]])
  y_fit.w <- exp(fit1_wo$coefficients[[2]]*x_fit + fit1_wo$coefficients[[1]])
  
  fit_data <- data.frame(x = x_fit, control = y_fit.c, warmed = y_fit.w)
  
  library(ggpmisc)
  str(dat)
  
  dat$plot_ab<-paste(dat$plot,dat$part,sep='')
  
  p2_1<-ggplot(dat[dat$soillayer=='0-10 cm',],
               aes(Temp,AC_N2,col=treatment))+
    geom_point(size=1,pch=1)+#geom_line(aes(group=plot_ab))+
    scale_color_manual(limits=c('control','warmed'),values = c('blue','red'))+
    geom_smooth(method = 'glm', 
                formula = y ~ x, 
                method.args = list(family = Gamma(link = "log")),  se=F,show.legend = T,linewidth=0.4)+
    facet_grid(.~soillayer,scales = 'free_y')+
    labs(x='Temperature (\u00B0C)',y=expression(paste('N'[2],' production (ng N  g'^-1,' h'^-1~')')))+
    stat_poly_eq(aes(Temp,log(AC_N2),col=treatment,
                     label = paste(stat(rr.label),stat(p.value.label),sep="*\",\"~")),
                 formula =y~x,eq.x.rhs="x",coef.digits = 2,rr.digits = 2, 
                 eq.with.lhs = "italic(ln(y))~`=`~",hjust=0,  #给“y"换形式         
                 na.rm = FALSE,parse = TRUE,label.x.npc = 0.03,label.y.npc=c(0.90,0.80),size = 3)+
    scale_y_continuous(breaks=seq(0,.4,.1),expand=c(0,0))+
    coord_cartesian(ylim=c(0,0.4),clip='on')+
    mythem;p2_1
  
  
  
  ############# N2O/N2 产生 ____T gradient
  ####Oa+e
  str(dat)
  dat$N2_N2O<-dat$AC_N2/dat$P_N2O
  p3_0<-ggplot(dat[dat$soillayer=='O layer',],
               aes(Temp,N2_N2O,col=treatment))+
    geom_point(size=1,pch=1)+
    scale_color_manual(limits=c('control','warmed'),values = c('blue','red'))+
    facet_grid(.~soillayer,scales = 'free_y')+
    # geom_smooth(method = 'lm',formula = y~x,se=F)+
    stat_poly_eq(aes(Temp,N2_N2O,col=treatment,
                     label = ifelse(stat(p.value.label)<0.05,paste(stat(rr.label),stat(p.value.label),sep="*\",\"~~~"),'ns')),
                 formula =y~x,eq.x.rhs="x",coef.digits = 2,rr.digits = 2,
                 eq.with.lhs = "italic(y)~`=`~",hjust=0,  #给“y"换形式
                 na.rm = FALSE,parse = TRUE,label.x.npc = 0.03,label.y.npc=c(0.90,0.80),size = 3)+
    labs(x='Temperature (\u00B0C)',y=expression(paste('N'[2],'/N'[2],'O')))+
    scale_y_continuous(breaks=seq(0,.4,.1),limits=c(0,.4),expand=c(0,0))+
    mythem;p3_0
  
  ####0-10cm
  p3_1<-ggplot(dat[dat$soillayer=='0-10 cm',],
               aes(Temp,N2_N2O,col=treatment))+
    geom_point(size=1,pch=1)+
    scale_color_manual(limits=c('control','warmed'),values = c('blue','red'))+
    facet_grid(.~soillayer,scales = 'free_y')+
    labs(x='Temperature (\u00B0C)',y=expression(paste('N'[2],'/N'[2],'O')))+
    # geom_smooth(method = 'lm',formula = y~x,se=F)+
    stat_poly_eq(aes(Temp,N2_N2O,col=treatment,
                     label = ifelse(stat(p.value.label)<0.05,paste(stat(rr.label),stat(p.value.label),sep="*\",\"~~~"),'ns')),
                 formula =y~x,eq.x.rhs="x",coef.digits = 2,rr.digits = 2,
                 eq.with.lhs = "italic(y)~`=`~",hjust=0,  #给“y"换形式
                 na.rm = FALSE,parse = TRUE,label.x.npc = 0.03,label.y.npc=c(0.90,0.80),size = 3)+
    scale_y_continuous(breaks=seq(0,.4,.1),limits=c(0,.4),expand=c(0,0))+
    mythem;p3_1
  
}


############# N2O/N2 产生______SM gradient
####Oa+e
str(dat)
dat$N2_N2O<-dat$AC_N2/dat$P_N2O

# 创建从蓝色到红色的颜色梯度函数
col_palette <- colorRampPalette(c("skyblue", "red"))

# 生成5个颜色
colors <- col_palette(5)

# 查看颜色
colors


p4_0<-ggplot(dat[as.factor(dat$soillayer)=='O layer',],
             aes(Dry_sm,N2_N2O,col=treatment,fill=treatment))+
  geom_point(size=1,pch=1)+  #aes(col=factor(Temp)),
  scale_color_manual(limits=c('control','warmed'),values = c('blue','red'))+
  facet_grid(.~soillayer,scales = 'free_y')+
  geom_smooth(method = 'lm',formula = y~x,se=F,linewidth=0.4)+
  stat_poly_eq(aes(label = paste(stat(eq.label),sep="*\",\"~")),
               formula =y~x,eq.x.rhs="x",coef.digits = 2,rr.digits = 2,
               eq.with.lhs = "italic(y)~`=`~",hjust=0,  #给“y"换形式
               na.rm = FALSE,parse = TRUE,label.x.npc = 0.07,label.y.npc=c(0.90,0.70),size = 3)+
  stat_poly_eq(aes(label = paste(stat(rr.label),stat(p.value.label),sep="*\",\"~")),
               formula =y~x,eq.x.rhs="x",coef.digits = 1,rr.digits = 2,
               eq.with.lhs = "italic(y)~`=`~",hjust=0,  #给“y"换形式
               na.rm = FALSE,parse = TRUE,label.x.npc = 0.07,label.y.npc=c(0.80,0.60),size = 3)+
  # geom_smooth(method = 'lm',formula = y~x,se=T,col='black')+
  # stat_poly_eq(aes(label = paste(stat(rr.label),stat(p.value.label),sep="*\",\"~")),
  #              formula =y~x+I(x^2),eq.x.rhs="x",coef.digits = 2,rr.digits = 2,
  #              eq.with.lhs = "italic(y)~`=`~",hjust=0,col='black',  #给“y"换形式
  #              na.rm = FALSE,parse = TRUE,label.x.npc = 0.07,label.y.npc=c(0.73),size = 3)+
  scale_y_continuous(breaks=seq(0,.4,.1),expand=c(0,0))+
  scale_x_continuous(breaks=seq(0.2,1.0,0.2),limits=c(0.2,1.0),expand=c(0.01,0.01))+
  coord_cartesian(ylim=c(0,0.4),clip='on')+
  labs(x=expression(paste('Soil moisture (g H'[2],'O g'^'-1',' soil)')),y=expression(paste('N'[2],'/N'[2],'O')))+
  mythem;p4_0

####0-10cm
p4_1<-ggplot(dat[as.factor(dat$soillayer)=='0-10 cm',],
             aes(Dry_sm,N2_N2O,col=treatment))+
  geom_point(size=1,pch=1)+  #aes(col=factor(Temp)),
  scale_color_manual(limits=c('control','warmed'),values = c('blue','red'))+
  facet_grid(.~soillayer,scales = 'free_y')+
  geom_smooth(method = 'lm',formula = y~x,se=F,linewidth=0.4)+
  stat_poly_eq(aes(label = paste(stat(eq.label),sep="*\",\"~")),
               formula =y~x,eq.x.rhs="x",coef.digits = 2,rr.digits = 2,
               eq.with.lhs = "italic(y)~`=`~",hjust=0,  #给“y"换形式
               na.rm = FALSE,parse = TRUE,label.x.npc = 0.07,label.y.npc=c(0.90,0.70),size = 3)+
  stat_poly_eq(aes(label = paste(stat(rr.label),stat(p.value.label),sep="*\",\"~")),
               formula =y~x,eq.x.rhs="x",coef.digits = 1,rr.digits = 2,
               eq.with.lhs = "italic(y)~`=`~",hjust=0,  #给“y"换形式
               na.rm = FALSE,parse = TRUE,label.x.npc = 0.07,label.y.npc=c(0.80,0.60),size = 3)+
  
  # geom_smooth(method = 'lm',formula = y~x,se=T,col='black')+
  # stat_poly_eq(aes(label = paste(stat(rr.label),stat(p.value.label),sep="*\",\"~")),
  #              formula =y~x+I(x^2),eq.x.rhs="x",coef.digits = 2,rr.digits = 2,
  #              eq.with.lhs = "italic(y)~`=`~",hjust=0,col='black',  #给“y"换形式
  #              na.rm = FALSE,parse = TRUE,label.x.npc = 0.07,label.y.npc=c(0.73),size = 3)+
  scale_x_continuous(breaks=seq(0.2,0.6,0.1),limits=c(0.2,0.6),expand=c(0.01,0.01))+
  scale_y_continuous(breaks=seq(0,.4,.1),expand=c(0,0))+
  coord_cartesian(ylim=c(0,0.4),clip='on')+
  labs(x=expression(paste('Soil moisture (g H'[2],'O g'^'-1',' soil)')),y=expression(paste('N'[2],'/N'[2],'O')))+
  mythem;p4_1

library(ggpubr)
figure<-ggarrange(p1_0+rremove("xlab")+theme(plot.margin = unit(c(0.1,0.13,0.1,0.03), 'cm')),
                  p2_0+rremove("xlab")+theme(plot.margin = unit(c(0.1,0.13,0.1,0.03), 'cm')),
                  p3_0+rremove("xlab")+theme(plot.margin = unit(c(0.1,0.13,0.1,0.03), 'cm')),
                  p4_0+rremove("xlab")+theme(plot.margin = unit(c(0.1,0.13,0.1,0.03), 'cm')),
                  p1_1+theme(plot.margin = unit(c(0.1,0.13,0.1,0.03), 'cm')),
                  p2_1+theme(plot.margin = unit(c(0.1,0.13,0.1,0.03), 'cm')),
                  p3_1+theme(plot.margin = unit(c(0.1,0.13,0.1,0.03), 'cm')),
                  p4_1+theme(plot.margin = unit(c(0.1,0.13,0.1,0.03), 'cm')),
                  labels = c("A", "B", "C","D","E","F","G","H"), 
                  font.label = list(size = 12, color = "black", face = "bold", family = NULL),
                  label.x = .032,label.y = 1.02,
                  ncol = 4, nrow = 2,align = "v",   ##"v"竖直对齐
                  widths = c(4,4,4,4), heights = c(3.2,3.45),
                  common.legend=FALSE);figure




setwd('D:/工作目录/202409/Manuscript_kai/Talk_20250503/Data and Code/Fig_S14')
ggsave("N2O+N2 production_new_20250721.pdf",figure, height = 5, width = 11)
