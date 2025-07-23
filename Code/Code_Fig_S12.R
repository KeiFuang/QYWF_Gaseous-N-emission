########## Ext.Data.Fig.9_Meta分析-增温N2O排放结果 ############
Sys.setlocale("LC_TIME","English")  
library(ggthemes);library(ggplot2)
mythem<-theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1), 
              axis.text.x = element_text(colour = "black",size=9,angle = 0,hjust = .5,vjust =1),
              axis.ticks = element_line(linewidth = .15),
              axis.ticks.length = unit(1,"mm"),
              prism.ticks.length = unit(0.5,"mm"),
              axis.text.y = element_text(colour = "black",size = 9,hjust=.5), 
              axis.title.x =element_text(size=9), axis.title.y=element_text(colour = "black",size=9),
              legend.text = element_text(size=9), legend.title =element_text(size=9),
              panel.background = element_rect(fill =NA, size=0.26,colour = "black", linetype = "solid",
                                              inherit.blank=T),
              panel.grid=element_blank(),     #element_line(colour="gray",size=0.5),
              panel.grid.major = element_blank(),    #element_line(colour = "gray",size=0.6), 
              plot.background = element_rect(fill=NA,color=NA,linetype = "solid"),
              legend.key = element_rect(fill = NA), 
              legend.background = element_rect(fill = NA), 
              plot.margin = unit(c(0.4,0.3,0.2,0.3),'cm'),   #调整画图区域的间距，从上右下左调整
              strip.background = element_rect(fill = NA,color='black'), 
              strip.text = element_text(colour = "black",size = 9,hjust=.5),
              legend.position = "none")

setwd('D:/Workspace/book/Qingyuan')
list.files('D:/Workspace/book/Qingyuan')
library(readxl)
mydat = read_excel("D:/Workspace/book/Qingyuan/meta_forest data.xlsx",sheet = 'workspace')

substrRight<- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))}

for(i in 1:length(mydat$Lat)){
  lat<-ifelse(substrRight(strsplit(mydat$Lat[i],'°')[[1]][2],1)=='N',1,-1)
  lon<-ifelse(substrRight(strsplit(mydat$Lon[i],'°')[[1]][2],1)=='E',1,-1)
  
  mydat$Lat[i]<-as.numeric(strsplit(mydat$Lat[i],'°')[[1]][1])*as.numeric(lat)
  mydat$Lon[i]<-as.numeric(strsplit(mydat$Lon[i],'°')[[1]][1])*as.numeric(lon)
  
  i=i+1
}


# site<-mysite[,c('Lat','Lon','SOC')]
mydat$Lat<-as.numeric(mydat$Lat)
mydat$Lon<-as.numeric(mydat$Lon)
str(mydat)

mydat$lnRR_v=mydat$N2O_w.sd^2/(mydat$Replications*mydat$N2O_w^2)+
  mydat$N2O_c.sd^2/(mydat$Replications*mydat$N2O_c^2)

mydat$WFPS=mydat$vwc_c/(1-mydat$BD/2.65)
library(metafor)

str(mydat)

mydat$pH<-as.numeric(mydat$pH)
library(car)
densityPlot(mydat$lnRR)
shapiro.test(mydat$lnRR) #P>0.05 符合正态

rma(lnRR,lnRR_v,data=mydat,method='REML')
rma(lnRR,lnRR_v,mods=~Country,data=mydat,method='REML') #Qm的P值<0.05,分类不妥




####发表偏倚
meta1<-rma(lnRR,lnRR_v,data=mydat,method='REML')

funnel(meta1, main="Funnel Plot",shade='white',back='white',
       ylim=c(0,2),xlab='Observed lnRR')
regtest(meta1)  #使用Egger线性回归法的结果也表明不存在发表偏倚（P=0.2811）


funnel(meta1, level=c(90, 95, 99), shade=c("white","gray", "darkgray"), refline=0, atransf=exp,at=log(c(.10, .25, .5, 1, 2, 4, 10)))
legend(1, 0.02, c("0.1 > p > 0.05", "0.05 > p >0.01", "< 0.01"), fill=c("white", "gray", "darkgray"))
# 轮廓增强漏斗图，该图从统计学99%、95%和90%三个水平层面进行检测。
# 从图5可看出，研究左右分布不均，主要分布在右侧，18项研究分布在无统计学意义区域中（白色区域），表明可能存在发表偏倚。


trf<-trimfill(meta1)
regtest(trf) 
funnel(trf, level=c(90, 95, 99), shade=c("white","gray", "darkgray"),
       refline=0, atransf=exp,at=log(c(.10, .25, .5, 1, 2, 4, 10)))
legend(1, 0.02, c("0.1 > p > 0.05", "0.05 > p >0.01", "< 0.01"), 
       fill=c("white", "gray", "darkgray"))

# 剪补后的轮廓增强漏斗图，由图可知，需要再增补一项研究（图中虚心圆圈）可以纠正漏斗图的不对称；
# 但增补的该项研究仍旧分布在无统计学意义区域（白色区域），
# 表明存在不具有统计学意义的研究未发表。因此，可以表明本Meta分析存在一定程度的发表偏倚。

pdf("funnel_plot_i.pdf", width = 12, height = 6.5 )
par(mfrow=c(1, 2))
funnel(meta1, main=NULL,shade='white',back='gray87',cex=1,
       ylim=c(0,2),xlim=c(-4,4),xlab='Observed lnRR',label = F)
polygon(x = c(-5, 4.3, 4.3, -5), 
        y = c(-0.08,-0.08,2.1,2.1),
        border = "black", lwd = 1)

funnel(trf, main=NULL,shade='white',back='gray87',cex=1,
       ylim=c(0,2),xlim=c(-4,4),xlab='Observed lnRR',label = F)
polygon(x = c(-5, 4.3, 4.3, -5), 
        y = c(-0.08,-0.08,2.1,2.1),
        border = "black", lwd = 1)

dev.off()

unique(mydat$Warming_methods)

#用metafor包画森林图#####
make_pct<-function(x)(exp(x)-1)*100#转化成百分比
library(dplyr)
# library(funModeling)
library(ggplot2)

mydat$dM=mydat$dM   #转换成WFPS
# mydat<-mydat[!mydat$Country=='Ireland',]  #对照排放为负，无法计算响应比

mydat1<-mydat[-which(is.na(mydat$lnRR)),c('ny','wy','Warming_methods','dT','WFPS','dM','lnRR','lnRR_v')]

unique(mydat1$ny)

str(mydat1)
mydat1$obs_yr<-ifelse(mydat1$ny==1,'ny_0-1',
                      ifelse(mydat1$ny==2,'ny_1-3','ny≥3'))

# mydat1$obs_yr<-factor(mydat1$obs_yr,
#                        levels = c("0-1","1-2","≥3"))


unique(mydat1$Warming_methods)
mydat1$Warming_methods<-factor(mydat1$Warming_methods,
                               levels = c("OTCs","geotherm","Infrared heaters","heating cables","Infrared heaters and heating cables"))


mydat1$dT_lab<-ifelse(mydat1$dT<=2,'T0-2',
                      ifelse(mydat1$dT<=4,'T2-4','T4-8'))

unique(mydat1$dT_lab)
mydat1$dT_lab<-factor(mydat1$dT_lab,
                      levels = c("T0-2","T2-4","T4-8"))
plot(mydat1$WFPS,mydat1$dM)

mydat1$WFPS_lab<-ifelse(mydat1$WFPS<=40,'W_20-40',
                        ifelse(mydat1$WFPS<=50,'W_40-50','W_50-60'))

plot(mydat1$dM)
str(mydat1)
#数据量较少，之划分为三类
plot(mydat1$dM)
mydat1$dM_lab<-ifelse(mydat1$dM<=0,'M_<0','M_>0')


plot(mydat1$wy)

mydat1$warm_yr<-ifelse(mydat1$wy<=2,'wy_0-2',
                       ifelse(mydat1$wy<=6,'wy_2-6','wy>6'))

mydat1$warm_yr<-factor(mydat1$warm_yr,
                       levels = c("wy_0-2","wy_2-6","wy>6"))


summarise<-dplyr::summarise  #由于函数所对应的包不知有一个，需手动指定


metalei1<-rma(lnRR,lnRR_v,mods=~dT_lab-1,data=mydat1,method='REML');metalei1 #Qm的P值<0.05,分类不妥
metalei2<-rma(lnRR,lnRR_v,mods=~WFPS_lab-1,data=mydat1,method='REML');metalei2 #Qm的P值<0.05,分类不妥
metalei3<-rma(lnRR,lnRR_v,mods=~dM_lab-1,data=mydat1,method='REML');metalei3 #Qm的P值<0.05,分类不妥
metalei4<-rma(lnRR,lnRR_v,mods=~warm_yr-1,data=mydat1,method='REML');metalei4 #Qm的P值<0.05,分类不妥
metalei5<-rma(lnRR,lnRR_v,mods=~obs_yr-1,data=mydat1,method='REML');metalei5 #Qm的P值<0.05,分类不妥
metalei6<-rma(lnRR,lnRR_v,mods=~Warming_methods-1,data=mydat1,method='REML');metalei6 #Qm的P值<0.05,分类不妥
metalei<-rma(lnRR,lnRR_v,data=mydat1,method='REML');metalei #计算总效应


str(mydat1)
metalei1.n<- mydat1 %>% group_by(dT_lab) %>% summarise(n=n())
metalei1.n<-na.omit(metalei1.n)
metalei1.df<-coef(summary(metalei1))%>%mutate(type="dT_lab",
                                              factor=levels(as.factor(mydat1$dT_lab)),
                                              size=metalei1.n$n)


metalei2.n<- mydat1 %>% group_by(WFPS_lab) %>% summarise(n=n())
metalei2.n<-na.omit(metalei2.n)
metalei2.df<-coef(summary(metalei2))%>%mutate(type="WFPS_lab",
                                              factor=levels(as.factor(mydat1$WFPS_lab)),
                                              size=metalei2.n$n)

metalei3.n<- mydat1 %>% group_by(dM_lab) %>% summarise(n=n())
metalei3.n<-na.omit(metalei3.n)
metalei3.df<-coef(summary(metalei3))%>%mutate(type="dM_lab",
                                              factor=levels(as.factor(mydat1$dM_lab)),
                                              size=metalei3.n$n)

metalei4.n<- mydat1 %>% group_by(warm_yr) %>% summarise(n=n())
metalei4.n<-na.omit(metalei4.n)
metalei4.df<-coef(summary(metalei4))%>%mutate(type="warm_yr",
                                              factor=levels(as.factor(mydat1$warm_yr)),
                                              size=metalei4.n$n)

metalei5.n<- mydat1 %>% group_by(obs_yr) %>% summarise(n=n())
metalei5.n<-na.omit(metalei5.n)
metalei5.df<-coef(summary(metalei5))%>%mutate(type="obs_yr",
                                              factor=levels(as.factor(mydat1$obs_yr)),
                                              size=metalei5.n$n)

metalei6.n<- mydat1 %>% group_by(Warming_methods) %>% summarise(n=n())
metalei6.n<-na.omit(metalei6.n)
metalei6.df<-coef(summary(metalei6))%>%mutate(type="Warming_methods",
                                              factor=levels(as.factor(mydat1$Warming_methods)),
                                              size=metalei6.n$n)

metalei.n<- mydat1 %>% group_by() %>% summarise(n=n())
metalei.n<-na.omit(metalei.n)
metalei.df<-coef(summary(metalei))%>%mutate(type="Warming_methods",
                                            factor='lnRR_N2O',
                                            size=metalei.n$n)

meta.df<-rbind(metalei1.df,metalei2.df,metalei3.df,metalei4.df,metalei5.df,metalei6.df,metalei.df)

unique(meta.df$factor)
meta.df$factor<-factor(meta.df$factor,
                       levels = c('M_>0','M_<0',
                                  "W_50-60","W_40-50","W_20-40",
                                  "wy>6","wy_2-6","wy_0-2",
                                  "ny≥3","ny_1-3","ny_0-1",
                                  "Infrared heaters and heating cables","heating cables","Infrared heaters","geotherm","OTCs",
                                  "T4-8","T2-4","T0-2","lnRR_N2O"))


unique(meta.df$type)
meta.df$type<-factor(meta.df$type,levels=c("dT_lab","WFPS_lab","dM_lab","warm_yr","obs_yr","Warming_methods"))

meta.df$estimateB<-make_pct(meta.df$estimate)#转化成百分数
meta.df$ci.ubB<-make_pct(meta.df$ci.ub) #提取变成百分数左边
meta.df$ci.lbB<-make_pct(meta.df$ci.lb)  #提取右边

fplot<-ggplot(data=meta.df[meta.df$type!='warm_yr',],aes(x=estimate,y=factor,col=factor(type)))+
  geom_point(size=3)+geom_errorbarh(aes(xmax=ci.ub,xmin=ci.lb),height=0.4)+
  geom_vline(aes(xintercept=0),linetype=2,size=0.45)+
  labs(x="Mean effect size",y="")+
  scale_x_continuous(limits = c(-1,3.5),breaks=seq(-1,3,1),expand = c(0,0))+
  geom_text(aes(label=paste('(',size,')')),x=-0.85,size=3)+
  theme_bw() +theme(legend.position="none" )+
  theme(panel.grid = element_blank())+
  geom_hline(aes(yintercept=2.5),linetype=1,size=0.1)+
  geom_hline(aes(yintercept=5.5),linetype=1,size=0.1)+
  geom_hline(aes(yintercept=8.5),linetype=1,size=0.1)+
  geom_hline(aes(yintercept=13.5),linetype=1,size=0.1)+
  geom_hline(aes(yintercept=16.5),linetype=1,size=0.1)+
  theme(legend.position="left");fplot

str(meta.df)
meta.df=meta.df[meta.df$type!='warm_yr',]
meta.df$factor<-factor(meta.df$factor,
                       levels = c('M_>0','M_<0',
                                  "W_50-60","W_40-50","W_20-40",
                                  "ny≥3","ny_1-3","ny_0-1",
                                  "Infrared heaters and heating cables","heating cables","Infrared heaters","geotherm","OTCs",
                                  "T4-8","T2-4","T0-2","lnRR_N2O"))
meta.df$numeric_vector<-as.numeric(meta.df$factor);meta.df$numeric_vector

library(ggprism)
fplot<-ggplot(data=meta.df,aes(x=numeric_vector,y=estimate))+
  annotate("rect",xmin=16.5,xmax=17.5,ymin=-3,ymax=-1,
           col="black",fill="white",alpha=.1)+
  annotate("rect",xmin=13.5,xmax=16.5,ymin=-3,ymax=-1,
           col="black",fill="red",alpha=.1)+
  annotate("rect",xmin=8.5,xmax=13.5,ymin=-3,ymax=-1,
           col="black",fill="gray",alpha=.1)+
  annotate("rect",xmin=5.5,xmax=8.5,ymin=-3,ymax=-1,
           col="black",fill="green",alpha=.1)+
  annotate("rect",xmin=2.5,xmax=5.5,ymin=-3,ymax=-1,
           col="black",fill="blue",alpha=.1)+
  annotate("rect",xmin=0.5,xmax=2.5,ymin=-3,ymax=-1,
           col="black",fill="yellow",alpha=.1)+
  geom_text(x=1, y=-2.5,label=expression(atop(paste(Delta, 'moisture'),paste('(%)'))),colour="black",inherit.aes = TRUE,
            check_overlap = T,angle=90,size=4)+
  geom_text(x=4, y=-2.5,label=expression(paste("moisture \n (%)")),colour="black",inherit.aes = TRUE,
            check_overlap = T,angle=90,size=4)+
  geom_text(x=7, y=-2.8,label="Duration",colour="black",inherit.aes = TRUE,
            check_overlap = T,angle=90,size=4)+
  geom_text(x=11, y=-2.5,label=expression(paste("Warming \n methods")),colour="black",inherit.aes = TRUE,
            check_overlap = T,angle=90,size=4)+
  geom_text(x=15, y=-2.8,label=expression(paste(Delta,'T')),colour="black",inherit.aes = TRUE,
            check_overlap = T,angle=90,size=4)+
  
  geom_point(size=1.6)+geom_errorbar(aes(ymax=ci.ub,ymin=ci.lb),width=0.24,size=0.5)+
  geom_hline(aes(yintercept=0),linetype=2,size=0.15)+
  labs(y="Mean effect size",x="")+
  scale_y_continuous(breaks=seq(-1,3,1),expand = c(0,0))+
  
  scale_x_discrete(limits=c(1:17),labels=c("≥0","<0",
                                           "50~60","40~50","25~40",
                                           "≥3 yr","1~3 yr","0~1 yr",
                                           "IR and heating cables","heating cables","IR","geotherm","OTCs",
                                           "4~8 \u00B0C","2~4 \u00B0C","0~2 \u00B0C","Overall"),expand=c(0,0))+
  coord_flip(clip = 'off',ylim = c(-1,3.5))+
  geom_text(aes(label=paste('(',size,')')),y=-0.75,size=3)+   #添加n
  geom_text(aes(y=estimate+0.35,               #添加p值
                label=ifelse(pval<0.05,'*','')),
            colour = "black", hjust=0,vjust=0.65,size=5)+
  geom_vline(aes(xintercept=2.5),linetype=1,size=0.1)+
  geom_vline(aes(xintercept=5.5),linetype=1,size=0.1)+
  geom_vline(aes(xintercept=8.5),linetype=1,size=0.1)+
  geom_vline(aes(xintercept=13.5),linetype=1,size=0.1)+
  geom_vline(aes(xintercept=16.5),linetype=1,size=0.1)+
  mythem+
  theme(axis.text.y = element_text(hjust = 1,vjust=0.5),
        axis.title.x = element_text(size=12));fplot

setwd('D:/工作目录/202409/Manuscript_kai/Talk_20241113/Data and Code/Extended Data Fig_9')
ggsave("Mean effect size.pdf",fplot, device=cairo_pdf,height = 7, width = 4.8)

#FE固定效应模型，mods=分类变量或连续变量
fixed1<-rma(lnRR,lnRR_v,mods=~MAT,data=mydat,method='FE');fixed1
fixed2<-rma(lnRR,lnRR_v,mods=~MAP,data=mydat,method='FE');fixed2
fixed3<-rma(lnRR,lnRR_v,mods=~N_dep,data=mydat,method='FE');fixed3
fixed4<-rma(lnRR,lnRR_v,mods=~pH,data=mydat,method='FE');fixed4
fixed5<-rma(lnRR,lnRR_v,mods=~Soil_temp,data=mydat,method='FE');fixed5
fixed6<-rma(lnRR,lnRR_v,mods=~WFPS,data=mydat,method='FE');fixed6

fixed7<-rma(lnRR,lnRR_v,mods=~dT,data=mydat,method='FE');fixed7
fixed8<-rma(lnRR,lnRR_v,mods=~dM,data=mydat,method='FE');fixed8
fixed9<-rma(lnRR,lnRR_v,mods=~SOC,data=mydat,method='FE');fixed9
fixed10<-rma(lnRR,lnRR_v,mods=~C_N,data=mydat,method='FE');fixed10
fixed11<-rma(lnRR,lnRR_v,mods=~BD,data=mydat, method='FE');fixed11
fixed12<-rma(lnRR,lnRR_v,mods=~lnRRsr,data=mydat,method='FE');fixed12


#REML最大似然法估算,混合随机效应模型，mods=分类变量或连续变量,
#Qt效应值的总体异质性，P<0.05说明异质性大；Qm已知因素引起的异质性，P<0.05即该因素对y有显著影响
fixed1<-rma(lnRR,lnRR_v,mods=~MAT,data=mydat,method='REML');fixed1
fixed2<-rma(lnRR,lnRR_v,mods=~MAP,data=mydat,method='REML');fixed2
fixed3<-rma(lnRR,lnRR_v,mods=~N_dep,data=mydat,method='REML');fixed3
fixed4<-rma(lnRR,lnRR_v,mods=~pH,data=mydat,method='REML');fixed4
fixed5<-rma(lnRR,lnRR_v,mods=~Soil_temp,data=mydat,method='REML');fixed5
fixed6<-rma(lnRR,lnRR_v,mods=~WFPS,data=mydat,method='REML');fixed6

fixed7<-rma(lnRR,lnRR_v,mods=~dT,data=mydat,method='REML');fixed7
fixed8<-rma(lnRR,lnRR_v,mods=~dM,data=mydat,method='REML');fixed8
fixed9<-rma(lnRR,lnRR_v,mods=~SOC,data=mydat,method='REML');fixed9
fixed10<-rma(lnRR,lnRR_v,mods=~C_N,data=mydat,method='REML');fixed10
fixed11<-rma(lnRR,lnRR_v,mods=~BD,data=mydat, method='REML');fixed11
fixed12<-rma(lnRR,lnRR_v,mods=~lnRRsr,data=mydat,method='REML');fixed12

fixed12$b   #相关系数

names(mydat)
mydat1<-mydat[-which(is.na(mydat$lnRR)),c('ny','wy','MAT','MAP','N_dep','pH','BD','Soil_temp','C_N','SOC','lnRRsr','RRn',
                                          'Warming_methods','dT','WFPS','dM','lnRR','lnRR_v')]



library(dplyr)
library(stringr)

str(mydat1)
meta.r2<-data.frame()
k=1
for(i in c('MAT','MAP','N_dep','pH','BD','Soil_temp','C_N','SOC','lnRRsr','dT','WFPS','dM')){
  i_vector=unlist(mydat1[,i])
  fixed_i<-rma(yi=mydat1$lnRR,vi=mydat1$lnRR_v,mods=~i_vector+0,data=mydat1,method='REML');fixed_i
  dat_i<-fixed_i$data   #提取fixed_i的数据框
  dat_i<-dat_i[!is.na(dat_i[,i]),]   #去除i列的NA行
  dat_i$yi=predict(fixed_i,unlist(dat_i[,i]))%>%as.data.frame()%>%pull(pred)   #预测
  
  meta.r2[k,1] = i
  
  fit_i=summary(lm(dat_i$lnRR~dat_i$yi))
  
  meta.r2[k,2] = round(fit_i$r.squared,3)   #R2
  meta.r2[k,3] = round(fit_i$coefficients[8],3)   #P-value
  
  k=k+1
}

str(fixed_i$data)

names(meta.r2)[1:3]=c('variables','R2','p_value')

#去除无用数据框
rm(i,fit_i,fixed_i,metalei1.df,metalei1.n,metalei2.df,metalei2.n,metalei3.df,metalei3.n,metalei4.df,metalei4.n,
   metalei5.df,metalei5.n,metalei6.df,metalei6.n,metalei7.df,metalei7.n,metalei8.df,metalei8.n,metalei9.df,metalei9.n,
   metalei10.df,metalei10.n,metalei11.df,metalei11.n,metalei12.df,metalei12.n,metalei.df,metalei.n,
   metalei,metalei1,metalei2,metalei3,metalei4,metalei5,metalei6,trf,fplot,dat_i,i_vector,k)

library(ggplot2);library(ggpmisc);library(ggprism)
str(mydat)

ggplot(mydat,aes(Lat,lnRR))+geom_point()+
  geom_smooth(method = 'lm',formula = y~x+I(x^2))+
  stat_poly_eq(formula =y~x+I(x^2),coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = paste(stat(eq.label),stat(rr.label),stat(p.value.label),
                                 sep ="*\",\"~~~")),#sep = "*\", \"*")),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 5)


mythem_sci<-theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1), 
                  axis.text.x = element_text(colour = "black",size=7,angle = 0,hjust = .5,vjust =1),
                  axis.ticks = element_line(linewidth = .15),
                  axis.ticks.length = unit(1,"mm"),
                  prism.ticks.length = unit(0.5,"mm"),
                  axis.text.y = element_text(colour = "black",size = 7), 
                  axis.title.x =element_text(size=7), axis.title.y=element_text(colour = "black",size=7),
                  legend.text = element_text(size=7), legend.title =element_text(size=7),
                  legend.margin = margin(unit(c(0.05,3,5,3),'mm')),
                  # legend.box.background = element_rect(fill ="white", size=0.6,colour = "black", linetype = "solid"),
                  panel.background = element_rect(fill =NA, linewidth=0.05,colour = "black", linetype = "solid",
                                                  inherit.blank=T),
                  panel.grid=element_blank(),     #element_line(colour="gray",size=0.5),
                  panel.grid.major = element_blank(),    #element_line(colour = "gray",size=0.6), 
                  plot.background = element_rect(fill=NA,color=NA,linetype = "solid"),
                  legend.key = element_rect(fill = NA), 
                  legend.background = element_rect(fill = NA), 
                  plot.margin = unit(c(0.1,0.1,0,0),'cm'),   #调整画图区域的间距，从上右下左调整
                  strip.background = element_rect(fill = NA,color='black',linewidth=0.05), 
                  strip.text = element_text(colour = "black",size = 7,hjust=.5),
                  legend.position = "none")

str(mydat)


########lnRR与环境因子的关系#########
pre1<-data.frame(x=seq(min(mydat$MAT),max(mydat$MAT),.5),
                 predict(fixed1,seq(min(mydat$MAT),max(mydat$MAT),.5)))

p1<-ggplot(mydat[mydat$lnRR<3,],aes(MAT,lnRR))+
  # geom_smooth(method = 'lm',formula = y~log(x),color='black',linewidth=0.3)+
  # geom_line(data=pre1,aes(x=x,y=pred),col='red')+
  # geom_ribbon(data=pre1,aes(x=x,y=pred,ymin=ci.lb,ymax=ci.ub),
  #             fill="red",alpha=0.2)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.20,
                                  paste(stat(eq.label),stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=1.5)+
  scale_x_continuous(limits=c(0,15),breaks = seq(0,15,5),expand = c(0,0))+
  scale_y_continuous(limits=c(-2,4),breaks = seq(-2,4,2),expand = c(0,0))+
  labs(x='MAT (\u00B0C)',y= expression(paste('lnRR of N'[2],'O emission')))+
  mythem;p1


pre2<-data.frame(x=seq(min(mydat$MAP),max(mydat$MAP),100),
                 predict(fixed2,seq(min(mydat$MAP),max(mydat$MAP),100)))

fixed2
p2<-ggplot(mydat,aes(MAP,lnRR))+
  geom_smooth(method = 'lm',formula = y~x,color='black',linewidth=0.3)+
  # geom_line(data=pre2,aes(x=x,y=pred),col='red')+
  # geom_ribbon(data=pre2,aes(x=x,y=pred,ymin=ci.lb,ymax=ci.ub),
  #             fill="red",alpha=0.2)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.40,
                                  paste(stat(eq.label),stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=-0.5,size = 2)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(aes(size=lnRR_v),pch=1)+geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=1.5)+
  scale_x_continuous(limits=c(0,1600),breaks = seq(0,1600,400),expand = c(0,0))+
  scale_y_continuous(limits=c(-1,2),breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x='MAP (mm)',y= expression(paste('lnRR of N'[2],'O emission')))+
  mythem;p2

pre3<-data.frame(x=seq(min(mydat$N_dep),max(mydat$N_dep),1),
                 predict(fixed3,seq(min(mydat$N_dep),max(mydat$N_dep),1)))

p3<-ggplot(mydat,aes(N_dep,lnRR))+
  # geom_smooth(method = 'lm',formula = y~log(x),color='black',linewidth=0.3)+
  # geom_line(data=pre3,aes(x=x,y=pred),col='red')+
  # geom_ribbon(data=pre3,aes(x=x,y=pred,ymin=ci.lb,ymax=ci.ub),
  #             fill="red",alpha=0.2)+
  stat_poly_eq(formula =y~log(x),coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.20,
                                  paste(stat(eq.label),stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=1.5)+
  # geom_vline(xintercept = 0,lty=2,col='red')+
  scale_x_continuous(limits=c(0,30),breaks = seq(0,30,10),expand = c(0,0))+
  scale_y_continuous(limits=c(-1,2),breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x=expression(paste('N deposition (kg N ha'^-1,'yr'^-1,')')),y= expression(paste('lnRR of N'[2],'O emission')))+mythem+
  coord_cartesian(clip="on",ylim=c(-2,4))+ guides(x="prism_offset_minor",y="prism_offset_minor");p3

pre4<-data.frame(x=seq(min(mydat$pH),max(mydat$pH),0.5),
                 predict(fixed4,seq(min(mydat$pH),max(mydat$pH),.5)))
p4<-ggplot(mydat,aes(pH,lnRR))+
  # geom_smooth(method = 'lm',formula = y~x,color='black',linewidth=0.3)+
  # geom_line(data=pre4,aes(x=x,y=pred),col='red')+
  # geom_ribbon(data=pre4,aes(x=x,y=pred,ymin=ci.lb,ymax=ci.ub),
  #             fill="red",alpha=0.2)+
  stat_poly_eq(formula =y~log(x),coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.20,
                                  paste(stat(eq.label),stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=1.5)+
  # geom_vline(xintercept = 0,lty=2,col='red')+
  scale_x_continuous(limits=c(1,10),breaks = seq(1,10,3),expand = c(0,0))+
  scale_y_continuous(limits=c(-1,2),breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x=expression(paste('pH')),y= expression(paste('lnRR of N'[2],'O emission')))+mythem+
  coord_cartesian(clip="on",ylim=c(-2,4))+ guides(x="prism_offset_minor",y="prism_offset_minor");p4
str(mydat)

pre5<-data.frame(x=seq(min(mydat$Soil_temp),max(mydat$Soil_temp),0.5),
                 predict(fixed5,seq(min(mydat$Soil_temp),max(mydat$Soil_temp),.5)))
p5<-ggplot(mydat,aes(Soil_temp,lnRR))+
  # geom_smooth(method = 'lm',formula = y~x,color='black',linewidth=0.3)+
  # geom_line(data=pre5,aes(x=x,y=pred),col='red')+
  # geom_ribbon(data=pre5,aes(x=x,y=pred,ymin=ci.lb,ymax=ci.ub),
  #             fill="red",alpha=0.2)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.50,
                                  paste(stat(eq.label),stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=1.5)+
  scale_x_continuous(limits=c(3,18),breaks = seq(3,18,5),expand = c(0,0))+
  scale_y_continuous(limits=c(-1,2),breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x='Soil temp (\u00B0C)',y= expression(paste('lnRR of N'[2],'O emission')))+
  mythem;p5

pre6<-data.frame(x=seq(min(na.omit(mydat$WFPS)),max(na.omit(mydat$WFPS)),1),
                 predict(fixed6,seq(min(na.omit(mydat$WFPS)),max(na.omit(mydat$WFPS)),1)))
str(mydat)
p6<-ggplot(mydat[mydat$lnRR<1,],aes(WFPS,lnRR))+
  geom_smooth(method = 'lm',formula = y~x+I(x^2),color='black',linewidth=0.3)+
  # geom_line(data=pre6,aes(x=x,y=pred),col='red')+
  # geom_ribbon(data=pre6,aes(x=x,y=pred,ymin=ci.lb,ymax=ci.ub),
  #             fill="red",alpha=0.2)+
  geom_point(data=mydat[mydat$lnRR>1,],pch=16,size=1.2)+
  stat_poly_eq(formula =y~x+I(x^2),coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.50,
                                  paste(stat(eq.label),stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=1.5)+
  scale_x_continuous(limits=c(0,80),breaks = seq(20,80,20),expand = c(0,0))+
  scale_y_continuous(limits=c(-1,2),breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x='WFPS (%)',y= expression(paste('lnRR of N'[2],'O emission')))+mythem+
  coord_cartesian(clip="on",ylim=c(-1,2))+ guides(x="prism_offset_minor",y="prism_offset_minor");p6

pre7<-data.frame(x=seq(min(na.omit(mydat$dT)),max(na.omit(mydat$dT)),.1),
                 predict(fixed7,seq(min(na.omit(mydat$dT)),max(na.omit(mydat$dT)),.1)))

p7<-ggplot(mydat,aes(dT,lnRR))+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.20,
                                  paste(stat(eq.label),stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=1.5)+
  scale_x_continuous(limits=c(0,9),breaks = seq(0,9,3),expand = c(0,0))+
  scale_y_continuous(limits=c(-1,2),breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x=expression(paste(Delta,' soil temp (\u00B0C)')),y= expression(paste('lnRR of N'[2],'O emission')))+mythem+
  coord_cartesian(clip="on",ylim=c(-2,4))+ guides(x="prism_offset_minor",y="prism_offset_minor");p7



mydat$dWFPS<-mydat$dM/(1-mydat$BD/2.65)
pre8<-data.frame(x=seq(min(na.omit(mydat$dWFPS)),max(na.omit(mydat$dWFPS)),1),
                 predict(fixed8,seq(min(na.omit(mydat$dWFPS)),max(na.omit(mydat$dWFPS)),1)))
str(mydat)
p8<-ggplot(mydat[mydat$lnRR<=3,],aes(dWFPS,lnRR))+
  geom_smooth(method = 'lm',formula = y~x,color='black',linewidth=0.3)+
  # geom_line(data=pre8,aes(x=x,y=pred),col='red')+
  # geom_ribbon(data=pre8,aes(x=x,y=pred,ymin=ci.lb,ymax=ci.ub),
  #             fill="red",alpha=0.2)+
  geom_point(data=mydat[mydat$lnRR>1,],pch=16,size=1.2)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.20,
                                  paste(stat(eq.label),stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=1.5)+
  scale_x_continuous(limits=c(-20,20),breaks = seq(-20,20,10),expand = c(0,0))+
  scale_y_continuous(limits=c(-1,2),breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x=expression(paste(Delta,' moisture (%)')),y= expression(paste('lnRR of N'[2],'O emission')))+mythem+
  coord_cartesian(clip="on",ylim=c(-1,2))+ guides(x="prism_offset_minor",y="prism_offset_minor");p8

pre9<-data.frame(x=seq(min(na.omit(mydat$SOC)),max(na.omit(mydat$SOC)),1),
                 predict(fixed9,seq(min(na.omit(mydat$SOC)),max(na.omit(mydat$SOC)),1)))

p9<-ggplot(mydat,aes(SOC,lnRR))+
  # geom_smooth(method = 'lm',formula = y~x,color='black',linewidth=0.3)+
  # geom_line(data=pre9,aes(x=x,y=pred),col='red')+
  # geom_ribbon(data=pre9,aes(x=x,y=pred,ymin=ci.lb,ymax=ci.ub),
  #             fill="red",alpha=0.2)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.20,
                                  paste(stat(eq.label),stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=1.5)+
  scale_x_continuous(limits=c(10,100),breaks = seq(10,100,30),expand = c(0,0))+
  scale_y_continuous(breaks = seq(-2,4,2),expand = c(0,0))+
  labs(x=expression(paste('SOC (g C kg'^-1,')')),y= expression(paste('lnRR of N'[2],'O emission')))+mythem+
  coord_cartesian(clip="on",ylim=c(-2,4))+ guides(x="prism_offset_minor",y="prism_offset_minor");p9

pre10<-data.frame(x=seq(min(na.omit(mydat$C_N)),max(na.omit(mydat$C_N)),1),
                  predict(fixed10,seq(min(na.omit(mydat$C_N)),max(na.omit(mydat$C_N)),1)))

p10<-ggplot(mydat,aes(C_N,lnRR))+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.20,
                                  paste(stat(eq.label),stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=1.5)+
  scale_x_continuous(limits=c(5,35),breaks = seq(5,35,10),expand = c(0,0))+
  scale_y_continuous(limits=c(-1,2),breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x=expression(paste('C/N')),y= expression(paste('lnRR of N'[2],'O emission')))+mythem+
  coord_cartesian(clip="on",ylim=c(-2,4))+ guides(x="prism_offset_minor",y="prism_offset_minor");p10


pre11<-data.frame(x=seq(min(na.omit(mydat$BD)),max(na.omit(mydat$BD)),.1),
                  predict(fixed11,seq(min(na.omit(mydat$BD)),max(na.omit(mydat$BD)),.1)))
p11<-ggplot(mydat,aes(BD,lnRR))+
  stat_poly_eq(formula =y~log(x),coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.20,
                                  paste(stat(eq.label),stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=1.5)+
  scale_x_continuous(limits=c(0.4,1.6),breaks = seq(0.4,1.6,0.4),expand = c(0,0))+
  scale_y_continuous(limits=c(-1,2),breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x=expression(paste('Bulk density')),y= expression(paste('lnRR of N'[2],'O emission')))+mythem+
  coord_cartesian(clip="on",ylim=c(-2,4))+ guides(x="prism_offset_minor",y="prism_offset_minor");p11

pre12<-data.frame(x=seq(min(na.omit(mydat$lnRRsr)),max(na.omit(mydat$lnRRsr)),.1),
                  predict(fixed12,seq(min(na.omit(mydat$lnRRsr)),max(na.omit(mydat$lnRRsr)),.1)))
p12<-ggplot(mydat,aes(lnRRsr,lnRR))+
  geom_smooth(method = 'lm',formula = y~x,color='black',linewidth=0.3)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.50,
                                  paste(stat(eq.label),stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=1.5)+
  scale_x_continuous(limits=c(-3,3),breaks = seq(-3,3,2),expand = c(0,0))+
  scale_y_continuous(limits=c(-1,2),breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x=expression(paste('lnRR of CO'[2],' emission')),y= expression(paste('lnRR of N'[2],'O emission')))+mythem+
  coord_cartesian(clip="on",ylim=c(-1,2))+ guides(x="prism_offset_minor",y="prism_offset_minor");p12

library(ggpubr)
figure<-ggarrange(p1+theme(plot.margin = unit(c(2,2,1,2),'mm')),p2+theme(plot.margin = unit(c(2,2,1,2),'mm')),
                  p3+theme(plot.margin = unit(c(2,2,1,2),'mm')),p4+theme(plot.margin = unit(c(2,5,1,2),'mm')),
                  
                  p5+theme(plot.margin = unit(c(0,2,1,2),'mm')),p6+theme(plot.margin = unit(c(0,2,1,2),'mm')),
                  p7+theme(plot.margin = unit(c(0,2,1,2),'mm')),p8+theme(plot.margin = unit(c(0,5,1,2),'mm')),
                  
                  p9+theme(plot.margin = unit(c(0,2,2,2),'mm')),p10+theme(plot.margin = unit(c(0,2,2,2),'mm')),
                  p11+theme(plot.margin = unit(c(0,2,2,2),'mm')),p12+theme(plot.margin = unit(c(0,5,2,2),'mm')),
                  labels =paste("(",letters[1:12],")",sep=""),
                  label.x = rep(c(0.005),12),label.y = 0.98,
                  ncol = 4, nrow = 3,align = "h",   ##"v"竖直对齐
                  font.label = list(size = 7, color ="black"),
                  widths = c(6,6,6,6.2), heights = c(5,5,5),legend = "bottom",
                  common.legend=T);figure

ggsave("N2O lnRR的影响因子1.pdf", figure, width =8, height =5.5,
       device=cairo_pdf) 

getwd()


########RRn与环境因子的关系#########
p1<-ggplot(mydat,aes(MAT,RRn))+
  # geom_smooth(method = 'lm',formula = y~x,color='black',linewidth=0.3)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.20,
                                  paste(stat(eq.label),stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0,size = 2.6)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+  geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=2)+
  scale_x_continuous(limits=c(0,15),breaks = seq(0,15,5),expand = c(0,0))+
  scale_y_continuous(limits=c(-1,2),breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x='MAT (\u00B0C)',y= expression(paste('RRn of N'[2],'O emission')))+
  guides(x="prism_offset_minor",y="prism_offset_minor")+
  mythem;p1

p2<-ggplot(mydat,aes(MAP,RRn))+
  geom_smooth(method = 'lm',formula = y~x,color='black',linewidth=0.3)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.20,
                                  paste(stat(eq.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2.6)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.50,
                                  paste(stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = 0.88,parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2.6)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+  geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=2)+
  scale_x_continuous(limits=c(0,1600),breaks = seq(0,1600,400),expand = c(0,0))+
  scale_y_continuous(limits=c(-1,2),breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x='MAP (mm)',y= expression(paste('RRn of N'[2],'O emission')))+
  guides(x="prism_offset_minor",y="prism_offset_minor")+
  mythem;p2

p3<-ggplot(mydat,aes(N_dep,RRn))+
  # geom_smooth(method = 'lm',formula = y~x,color='black',linewidth=0.3)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.20,
                                  paste(stat(eq.label),stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2.6)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+ geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=2)+
  # geom_vline(xintercept = 0,lty=2,col='red')+
  scale_x_continuous(limits=c(0,30),breaks = seq(0,30,10),expand = c(0,0))+
  scale_y_continuous(limits=c(-1,2),breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x=expression(paste('N deposition (kg N ha'^-1,'yr'^-1,')')),
       y= expression(paste('RRn of N'[2],'O emission')))+mythem+
  coord_cartesian(clip="on",ylim=c(-1,2))+ guides(x="prism_offset_minor",y="prism_offset_minor");p3


p4<-ggplot(mydat,aes(pH,RRn))+
  # geom_smooth(method = 'lm',formula = y~x,color='black',linewidth=0.3)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.10,
                                  paste(stat(eq.label),stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2.6)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+  geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=2)+
  # geom_vline(xintercept = 0,lty=2,col='red')+
  scale_x_continuous(limits=c(1,10),breaks = seq(1,10,3),expand = c(0,0))+
  scale_y_continuous(limits=c(-1,2),breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x=expression(paste('pH')),y= expression(paste('RRn of N'[2],'O emission')))+mythem+
  coord_cartesian(clip="on",ylim=c(-1,2))+ guides(x="prism_offset_minor",y="prism_offset_minor");p4

str(mydat)

p5<-ggplot(mydat,aes(Soil_temp,RRn))+
  # geom_smooth(method = 'lm',formula = y~x,color='black',linewidth=0.3)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.05,
                                  paste(stat(eq.label),stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2.6)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+ geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=2)+
  scale_x_continuous(limits=c(3,18),breaks = seq(3,18,5),expand = c(0,0))+
  scale_y_continuous(limits=c(-1,2),breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x='Soil temp (\u00B0C)',y= expression(paste('RRn of N'[2],'O emission')))+
  guides(x="prism_offset_minor",y="prism_offset_minor")+
  mythem;p5


p6<-ggplot(mydat,aes(WFPS,RRn))+
  # geom_smooth(method = 'lm',formula = y~x,color='black',linewidth=0.3)+
  stat_poly_eq(formula =y~x+I(x^2),coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.20,
                                  paste(stat(eq.label),stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2.6)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+  geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=2)+
  scale_x_continuous(limits=c(20,80),breaks = seq(20,80,20),expand = c(0,0))+
  scale_y_continuous(limits=c(-1,2),breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x='WFPS (%)',y= expression(paste('RRn of N'[2],'O emission')))+mythem+
  coord_cartesian(clip="on",ylim=c(-1,2))+ guides(x="prism_offset_minor",y="prism_offset_minor");p6

p7<-ggplot(mydat,aes(dT,RRn))+
  # geom_smooth(method = 'lm',formula = y~x,color='black',linewidth=0.3)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.20,
                                  paste(stat(eq.label),stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2.6)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+  geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=2)+
  scale_x_continuous(limits=c(0,9),breaks = seq(0,9,3),expand = c(0,0))+
  scale_y_continuous(limits=c(-1,2),breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x=expression(paste(Delta,' soil temp (\u00B0C)')),y= expression(paste('RRn of N'[2],'O emission')))+mythem+
  coord_cartesian(clip="on",ylim=c(-1,2))+ guides(x="prism_offset_minor",y="prism_offset_minor");p7


str(mydat)
p8<-ggplot(mydat,aes(dWFPS,RRn))+
  geom_smooth(method = 'lm',formula = y~x,color='black',linewidth=0.3)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.20,
                                  paste(stat(eq.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2.6)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.20,
                                  paste(stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = 0.88,parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2.6)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+  geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=2)+
  scale_x_continuous(limits=c(-20,20),breaks = seq(-20,20,10),expand = c(0,0))+
  scale_y_continuous(limits=c(-1,2),breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x=expression(paste(Delta,' WFPS (%)')),y= expression(paste('RRn of N'[2],'O emission')))+mythem+
  coord_cartesian(clip="on",ylim=c(-1,2))+ guides(x="prism_offset_minor",y="prism_offset_minor");p8


p9<-ggplot(mydat,aes(SOC,RRn))+
  # geom_smooth(method = 'lm',formula = y~x,color='black',linewidth=0.3)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.10,
                                  paste(stat(eq.label),stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2.6)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+  geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=2)+
  scale_x_continuous(limits=c(10,100),breaks = seq(10,100,30),expand = c(0,0))+
  scale_y_continuous(limits=c(-1,2),breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x=expression(paste('SOC (g C kg'^-1,')')),y= expression(paste('RRn of N'[2],'O emission')))+mythem+
  coord_cartesian(clip="on",ylim=c(-1,2))+ guides(x="prism_offset_minor",y="prism_offset_minor");p9


p10<-ggplot(mydat,aes(C_N,RRn))+
  # geom_smooth(method = 'lm',formula = y~x,color='black',linewidth=0.3)+
  stat_poly_eq(formula =y~poly(x,1),coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.20,
                                  paste(stat(eq.label),stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2.6)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+  geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=2)+
  scale_x_continuous(limits=c(5,35),breaks = seq(5,35,10),expand = c(0,0))+
  scale_y_continuous(limits=c(-1,2),breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x=expression(paste('C/N')),y= expression(paste('RRn of N'[2],'O emission')))+mythem+
  coord_cartesian(clip="on",ylim=c(-1,2))+ guides(x="prism_offset_minor",y="prism_offset_minor");p10


p11<-ggplot(mydat,aes(BD,RRn))+
  # geom_smooth(method = 'lm',formula = y~x,color='black',linewidth=0.3)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.20,
                                  paste(stat(eq.label),stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2.6)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+  geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=2)+
  scale_x_continuous(limits=c(0.4,1.6),breaks = seq(0.4,1.6,0.4),expand = c(0,0))+
  scale_y_continuous(limits=c(-1,2),breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x=expression(paste('Bulk density (g cm'^-3,')')),y= expression(paste('RRn of N'[2],'O emission')))+mythem+
  coord_cartesian(clip="on",ylim=c(-1,2))+ guides(x="prism_offset_minor",y="prism_offset_minor");p11

p12<-ggplot(mydat,aes(lnRRsr,RRn))+
  geom_smooth(method = 'lm',formula = y~x,color='black',linewidth=0.3)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.50,
                                  paste(stat(eq.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2.6)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.50,
                                  paste(stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = 0.88,parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2.6)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+  geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=2)+
  scale_x_continuous(limits=c(-1,3),breaks = seq(-1,3,1),expand = c(0,0))+
  scale_y_continuous(breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x=expression(paste('lnRR of CO'[2],' emission')),y= expression(paste('RRn of N'[2],'O emission')))+mythem+
  coord_cartesian(clip="on",ylim=c(-1,2))+ guides(x="prism_offset_minor",y="prism_offset_minor");p12

library(ggpubr)
figure<-ggarrange(p1+theme(plot.margin = unit(c(2,2,1,2),'mm')),p2+theme(plot.margin = unit(c(2,2,1,2),'mm')),
                  p3+theme(plot.margin = unit(c(2,2,1,2),'mm')),p4+theme(plot.margin = unit(c(2,5,1,2),'mm')),
                  
                  p5+theme(plot.margin = unit(c(0,2,1,2),'mm')),p6+theme(plot.margin = unit(c(0,2,1,2),'mm')),
                  p7+theme(plot.margin = unit(c(0,2,1,2),'mm')),p8+theme(plot.margin = unit(c(0,5,1,2),'mm')),
                  
                  p9+theme(plot.margin = unit(c(0,2,2,2),'mm')),p10+theme(plot.margin = unit(c(0,2,2,2),'mm')),
                  p11+theme(plot.margin = unit(c(0,2,2,2),'mm')),p12+theme(plot.margin = unit(c(0,5,2,2),'mm')),
                  labels =paste("(",letters[1:12],")",sep=""),
                  label.x = rep(c(0.005),12),label.y = 0.98,
                  ncol = 4, nrow = 3,align = "h",   ##"v"竖直对齐
                  font.label = list(size = 7, color ="black"),
                  widths = c(6,6,6,6.2), heights = c(5,5,5),legend = "bottom",
                  common.legend=T);figure

ggsave("N2O RRn的影响因子.pdf", figure, width =9, height =6,
       device=cairo_pdf) 


b1<-ggplot(mydat,aes(Soil_temp,dWFPS))+
  geom_smooth(method = 'lm',formula = y~x,color='black',linewidth=0.3)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.20,
                                  paste(stat(eq.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2.6)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.20,
                                  paste(stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = 0.88,parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2.6)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+  geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=2)+
  scale_y_continuous(breaks = seq(-20,20,10),expand = c(0,0))+
  scale_x_continuous(limits=c(3,18),breaks = seq(3,18,5),expand = c(0,0))+
  labs(y=expression(paste(Delta,' WFPS (%)')),x=expression(paste('Soil temp (°C)')))+mythem+
  coord_cartesian(clip="on",ylim=c(-20,20))+ guides(x="prism_offset_minor",y="prism_offset_minor");b1


b2<-ggplot(mydat,aes(dT,dWFPS))+
  # geom_smooth(method = 'lm',formula = y~x,color='black',linewidth=0.3)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.20,
                                  paste(stat(eq.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2.6)+
  # stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
  #              eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
  #              aes(label = ifelse(after_stat(p.value)<0.20,
  #                                 paste(stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
  #              label.x.npc = "left", label.y.npc = 0.88,parse = TRUE,
  #              vstep=0.2,hstep=0.1,size = 2.3)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+  geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=2)+
  scale_y_continuous(limits=c(-20,20),breaks = seq(-20,20,10),expand = c(0,0))+
  scale_x_continuous(limits=c(0,9),breaks = seq(0,9,3),expand = c(0,0))+
  labs(y=expression(paste(Delta,' WFPS (%)')),x=expression(paste(Delta,' Soil temp (°C)')))+mythem+
  coord_cartesian(clip="on",ylim=c(-20,20))+ guides(x="prism_offset_minor",y="prism_offset_minor");b2


b3<-ggplot(mydat,aes(MAP,dWFPS))+
  geom_smooth(method = 'lm',formula = y~x,color='black',linewidth=0.3)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.20,
                                  paste(stat(eq.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2.6)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.20,
                                  paste(stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = 0.88,parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2.6)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+  geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=2)+
  scale_y_continuous(limits=c(-20,20),breaks = seq(-20,20,10),expand = c(0,0))+
  scale_x_continuous(limits=c(0,1600),breaks = seq(0,1600,400),labels = seq(0,1600,400),expand = c(0,0))+
  labs(y=expression(paste(Delta,' WFPS (%)')),x=expression(paste('MAP (mm)')))+mythem+
  coord_cartesian(clip="on",ylim=c(-20,20))+ guides(x="prism_offset_minor",y="prism_offset_minor");b3

b4<-ggplot(mydat,aes(WFPS,dWFPS))+
  # geom_smooth(method = 'lm',formula = y~x,color='black',linewidth=0.3)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = ifelse(after_stat(p.value)<0.05,
                                  paste(stat(eq.label),sep ="*\",\"~~~"),'ns')),
               label.x.npc = "left", label.y.npc = "top",parse = TRUE,
               vstep=0.2,hstep=0.1,size = 2.6)+
  # stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
  #              eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
  #              aes(label = ifelse(after_stat(p.value)<0.05,
  #                                 paste(stat(rr.label),stat(p.value.label),sep ="*\",\"~~~"),'ns')),
  #              label.x.npc = "left", label.y.npc = 0.88,parse = TRUE,
  #              vstep=0.2,hstep=0.1,size = 2.3)+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=1.2)+  geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=2)+
  scale_y_continuous(limits=c(-20,20),breaks = seq(-20,20,10),expand = c(0,0))+
  scale_x_continuous(limits=c(20,80),breaks = seq(20,80,20),expand = c(0,0))+
  labs(y=expression(paste(Delta,' WFPS (%)')),x=expression(paste('WFPS (%)')))+mythem+
  coord_cartesian(clip="on",ylim=c(-20,20))+ guides(x="prism_offset_minor",y="prism_offset_minor");b4



###patchwork比ggarange更好用
library(patchwork)
p1+theme(plot.margin = unit(c(2,2,0,2),'mm'))+
  p2+rremove("y.text")+rremove("ylab")+theme(plot.margin = unit(c(2,2,0,2),'mm'))+
  p3+rremove("y.text")+rremove("ylab")+theme(plot.margin = unit(c(2,2,0,2),'mm'))+
  p4+rremove("y.text")+rremove("ylab")+theme(plot.margin = unit(c(2,5,0,2),'mm'))+
  
  p9+theme(plot.margin = unit(c(0.15,2,0,2),'mm'))+
  p10+rremove("y.text")+rremove("ylab")+theme(plot.margin = unit(c(0.15,2,0,2),'mm'))+
  p11+rremove("y.text")+rremove("ylab")+theme(plot.margin = unit(c(0.15,2,0,2),'mm'))+
  p12+rremove("y.text")+rremove("ylab")+theme(plot.margin = unit(c(0.15,5,0,2),'mm'))+
  
  p5+theme(plot.margin = unit(c(0.15,2,0,2),'mm'))+
  p7+rremove("y.text")+rremove("ylab")+theme(plot.margin = unit(c(0.15,2,0,2),'mm'))+
  p6+rremove("y.text")+rremove("ylab")+theme(plot.margin = unit(c(0.15,2,0,2),'mm'))+
  p8+rremove("y.text")+rremove("ylab")+theme(plot.margin = unit(c(0.15,5,0,2),'mm'))+
  
  b1+theme(plot.margin = unit(c(0.2,2,2,1.2),'mm'))+
  b2+rremove("y.text")+rremove("ylab")+theme(plot.margin = unit(c(0.2,2,2,1.2),'mm'))+
  b3+rremove("y.text")+rremove("ylab")+theme(plot.margin = unit(c(0.2,2,2,1.2),'mm'))+
  b4+rremove("y.text")+rremove("ylab")+theme(plot.margin = unit(c(0.2,5,2,1.2),'mm'))+
  
  plot_layout(ncol = 4,nrow=4)+
  plot_annotation(tag_levels = 'A')&    #前缀，连接符，后缀 tag_prefix = '(',tag_sep = '', tag_suffix = ')'
  theme(plot.tag.position = c(0.90, 0.90), 
        plot.tag = element_text(size = 9,vjust = 0,hjust=0,face="bold"))->figure;figure

ggsave("D:/工作目录/202409/Manuscript_kai/Talk_20250503/Data and Code/Fig. S12/Fig.S12_N2O RRn的影响因子SI.pdf", figure, width =7.5, height =7.5,
       device=cairo_pdf) 

############## RRn响应的地理分布 ##############
str(mydat)
unique(mydat$Site)

ggplot(mydat,aes(reorder(Site, Lon),RRn))+
  geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.3)+
  geom_point(pch=1,size=5)+  geom_point(data=mydat[mydat$Site=='Qingyuan',],pch=16,color='red',size=4.8)+
  scale_x_discrete(labels = c("Harvard Forest"="Harvard","Huntington Wildlife Forest"="Huntington", "Cloquet"="Cloquet",
                              "Ely"="Ely","Stillberg"="Stillberg","Dooary Forest"="Dooary","Achenkirch"="Achenkirch", 
                              "Maoxian"="Maoxian","Ningshan"="Ningshan","Wusutu"="Wusutu","Reykir"="Reykir","Qingyuan"="Qingyuan"),
                   expand = c(0.03,0.03))+
  scale_y_continuous(breaks = seq(-1,2,1),expand = c(0,0))+
  labs(x=NULL,y= expression(paste('RRn of N'[2],'O emission')))+mythem+
  theme(axis.title.y = element_text(size=14),axis.text.x = element_text(size=14),axis.text.y = element_text(size=14))+
  coord_cartesian(clip="on",ylim=c(-1,2))+ guides(y="prism_offset_minor")->p;p

ggsave("Supp.1_N2O RRn_b.pdf", p, width =12.2,height = 3.2,
       device=cairo_pdf) 
       