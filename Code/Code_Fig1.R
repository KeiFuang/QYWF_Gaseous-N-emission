library(ggplot2)
mythem<-theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1), 
              axis.text.x = element_text(colour = "black",size=10,angle = 0,hjust = .5,vjust =1),
              axis.ticks = element_line(linewidth = .25),
              axis.ticks.length = unit(2,"mm"),
              prism.ticks.length = unit(1,"mm"),
              axis.text.y = element_text(colour = "black",size = 10,hjust=.5), 
              axis.title.x =element_text(size=10), axis.title.y=element_text(colour = "black",size=10,hjust=.5),
              legend.text = element_text(size=10,hjust=0), legend.title =element_text(size=18),
              panel.background = element_rect(fill = NA,colour = NA ),   
              panel.border = element_rect(fill = NA,colour = "black",linewidth = 0.25,linetype = "solid" ),   
              panel.grid=element_blank(),    # element_line(colour="gray",size=30),
              panel.grid.major = element_blank(),    #element_line(colour = "gray",size=0.6), 
              plot.background = element_rect(fill=NA,color=NA,linetype = "solid"),
              legend.key = element_rect(fill = NA), 
              legend.background = element_rect(fill = NA), 
              plot.margin = unit(c(0.4,0.3,0.2,0.3),'cm'),   
              strip.background = element_rect(fill = NA,color='black'), 
              strip.text = element_text(colour = "black",size = 12,hjust=.5),
              legend.position = "none")

library(ggplot2);library(ggpubr);library(ggprism)

ggplot(dat[dat$group=='soil',],aes(depth,T_diff.m,group=year))+
  geom_hline(yintercept = 2,lty=2,linewidth=0.47,col='black')+
  geom_point(aes(fill=T_diff.m),position = position_jitter(width = 0.35),fill='red2',size=1.4,alpha=0.6,pch=21)+
  stat_mean(data=dat[dat$labs %in% c('sensor'),],aes(group =depth),pch=16,size=1.4,col='black')+
  geom_boxplot(data=dat[dat$labs %in% c('sensor'),],
               aes(group =depth),width=4,linewidth=0.28, outlier.shape = NULL,outlier.colour = NA,fill=NA,color='black',alpha=0.7)+ 
  coord_flip(clip="off")+mythem+
  scale_x_reverse(breaks = c(40,30,20,10,5,0),limits=c(43,-3),
                  labels=c("40","30","20","10","5",'O horizon'),expand=c(0,0))+
  labs(x='Soil depth (cm)',y=expression(paste('Temperature difference (\u00B0C)')))+
  # scale_fill_gradient2(low = "#0084D1",high = "#FF420E",mid = "yellow2",midpoint = 1.2)+  
  scale_y_continuous(limits=c(0,3),breaks=seq(0,3,1),minor_breaks=seq(0,3,0.2),expand=c(0,0))+
  guides(x="prism_offset_minor")->p1;p1

ggplot(dat[dat$labs %in% c('sensor','lolly')| dat$labs_layer=='sample_Organic',],
       aes(depth,diff_weight_moisture,group=year))+
  geom_hline(yintercept = 0,lty=2,linewidth=0.47,col='black')+
  geom_point(data=dat[dat$labs %in% c('sensor')| dat$labs_layer=='sample_Organic',],alpha=0.5,
             aes(shape=labs,y=100),size=1.4,outlier.shape = NULL,outlier.colour = NA,fill='white',color='black',alpha=0.6,notch=F)+ 
  geom_point(data=dat[dat$labs %in% c('sensor','lolly'),],aes(shape=labs),alpha=0.5,
             position = position_jitter(width = 0.35),fill='red2',size=1.4,pch=21)+
  geom_point(data=dat[dat$labs_layer=='sample_Organic',],alpha=0.5,
             position = position_jitter(width = 0.35),bg='red2',size=1.4,pch=21)+
  stat_mean(aes(group =layer),pch=16,size=1.4,col='black')+
  geom_boxplot(data=dat[dat$labs %in% c('sensor')| dat$labs_layer=='sample_Organic',],
               aes(group =layer),width=4, linewidth=0.28,outlier.shape = NULL,outlier.colour = NA,fill=NA,color='black',alpha=0.7,notch=F)+ 
  scale_shape_manual(limits=c('sensor','lolly'),labels=c('automatic','mannual'),
                     values=c(1,2),name='Methods')+
  coord_flip(clip="off")+mythem+theme(legend.position = 'none')+
  scale_x_reverse(breaks = c(40,30,20,10,5,0),limits=c(43,-3),
                  labels=c("40","30","20","10","5",'O horizon'),expand=c(0,0))+
  labs(x='Depth (cm)',y=expression(paste('Moisture difference (g H'[2],'O g'^-1,' soil)')))+
  scale_y_continuous(limits=c(-0.2,0.1),breaks=seq(-0.2,0.1,0.1),
                     minor_breaks=seq(-0.2,0.1,0.02),expand=c(0,0))+
  guides(x="prism_offset_minor")->p2;p2

figure<-ggarrange(p1,p2+rremove("y.text")+rremove("ylab"),
                  ncol = 2, nrow = 1,align = "h",   
                  labels = c("A", "B"),
                  label.x = c(0.89,0.88),label.y = 0.95,
                  font.label = list(size = 10, color ="black"),
                  widths = c(5.5,4.35), heights = c(6,6),
                  common.legend=F);figure


ggsave("D:/工作目录/202409/Manuscript_kai/Talk_20250503/Data and Code/Fig_1/Fig_1_warming effect on temp and moisture.pdf", 
       figure, width =6, height = 4, device=cairo_pdf) 

unique(dat$labs)

str(dat.m)

dat.m[dat.m$labs=='sample' & dat.m$soillayer=='Organic',c('diff_weight_moisture','diff_weight_moisture.se')]
dat.m[dat.m$labs=='sensor' & dat.m$soillayer=='5cm',c('diff_weight_moisture','diff_weight_moisture.se')]

