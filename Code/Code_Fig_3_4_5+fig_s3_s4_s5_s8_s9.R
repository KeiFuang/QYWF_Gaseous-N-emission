############## Combine data ############
############# mydata_h<-mydata_h[,-c(5,6,8:11,12,14,18,19)]
############# mydata_h<-mydata_h[as.Date(mydata_h$datehour)<='2024/1/1',]
# rm(list=ls());gc()

.libPaths(c("D:/Workspace/Rtrial/R library 2024",                        
            "C:/Program Files/R/R-4.1.3/library"))

pacman::p_load("ggplot2","ggpmisc","ggpubr","ggprism","doBy","car","readr","crayon","rstatix","plot3D","viridis",
               "lme4","MCMCglmm","tidyverse","directlabels","agricolae","lubridate","scales","ggpointdensity",
               "patchwork","reshape2","dplyr","animation")  

{Sys.setlocale("LC_TIME","English")  
  
  windowsFonts(Arial=windowsFont("Arial"),
               Times=windowsFont("Times New Roman"))
  
  mythem<-theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1), 
                axis.text.x = element_text(colour = "black",size=18,angle = 0,hjust = .5,vjust =1),
                axis.ticks = element_line(linewidth = .5),
                axis.ticks.length = unit(2,"mm"),
                prism.ticks.length = unit(1,"mm"),
                axis.text.y = element_text(colour = "black",size = 18,hjust=.5), 
                axis.title.x =element_text(size=18), axis.title.y=element_text(colour = "black",size=18,hjust=.5),
                legend.text = element_text(size=18,hjust=0), legend.title =element_text(size=18),
                # legend.margin = margin(unit(c(0.05,3,1,3),'mm')),
                # legend.box.background = element_rect(fill ="white", size=0.6,colour = "black", linetype = "solid"),
                panel.background = element_rect(fill =NA, linewidth=0.6,colour = "black", linetype = "solid",
                                                inherit.blank=T),
                panel.grid=element_blank(),     #element_line(colour="gray",size=0.5),
                panel.grid.major = element_blank(),    #element_line(colour = "gray",size=0.6), 
                plot.background = element_rect(fill=NA,color=NA,linetype = "solid"),
                legend.key = element_rect(fill = NA, colour = NA), 
                legend.background = element_rect(fill = NA), 
                plot.margin = unit(c(0.4,0.3,0.2,0.3),'cm'),   #调整画图区域的间距，从上右下左调整
                strip.background = element_rect(fill = NA,color='black'), 
                strip.text = element_text(colour = "black",size = 18,hjust=.5),
                legend.position = "none")
  
  mythem_sci<-theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1), 
                    axis.text.x = element_text(colour = "black",size=7,angle = 0,hjust = .5,vjust =1),
                    axis.ticks = element_line(linewidth = .15),
                    axis.ticks.length = unit(1,"mm"),
                    prism.ticks.length = unit(0.5,"mm"),
                    axis.text.y = element_text(colour = "black",size = 7), 
                    axis.title.x =element_text(size=7), axis.title.y=element_text(colour = "black",size=7,vjust =-1.5),
                    legend.text = element_text(size=7), legend.title =element_text(size=7),
                    # legend.margin = margin(unit(c(0.05,3,1,3),'mm')),
                    # legend.box.background = element_rect(fill ="white", size=0.6,colour = "black", linetype = "solid"),
                    panel.background = element_rect(fill =NA, linewidth=0.05,colour = "black", linetype = "solid",
                                                    inherit.blank=T),
                    panel.grid=element_blank(),     #element_line(colour="gray",size=0.5),
                    panel.grid.major = element_blank(),    #element_line(colour = "gray",size=0.6), 
                    plot.background = element_rect(fill=NA,color=NA,linetype = "solid"),
                    legend.key = element_rect(fill = NA, colour = NA), 
                    legend.background = element_rect(fill = NA), 
                    plot.margin = unit(c(0.1,0.1,0,0),'cm'),   #调整画图区域的间距，从上右下左调整
                    strip.background = element_rect(fill = NA,color='black',linewidth=0.05), 
                    strip.text = element_text(colour = "black",size = 7,hjust=.5),
                    legend.position = "none")
  
  ###计算数据均值，范围，se
  datFUN<-function(x){c(mean=round(mean(na.omit(x)),3),
                        range=paste(round(min(na.omit(x)),3),round(max(na.omit(x)),3),sep="~"),
                        n=length(na.omit(x)),
                        se=round(sd(na.omit(x))/sqrt(length(na.omit(x))),3),
                        sd=round(sd(na.omit(x)),3))
  }
  
  '%notin%' <- Negate('%in%')   #not in
}




##################### Data flow process #############
############# mydata_h -> mydata_c -> mydata_d -> F.pt -> F.treat -> warm_eff -> warm_eff.m
length(na.omit(mydata_h$system_id))   # 201481 records from 2018 to 2023

names(mydata_h)

mydata_c<-summaryBy(data=mydata_h,
                    air_temp + N2O_N + NO_N + N2_N + soil_temp + soil_moisture~date+system_id+plot,FUN=myfun)


#去除通量测量期间system_id=NA时的数据
ids<-unique(na.omit(mydata_c$system_id))

dates<-unique(mydata_c$date[mydata_c$system_id%in% ids])

str(mydata_c)

mydata_c$year<-format(mydata_c$date,"%Y")

mydata_c$plot<-NA
i=1
for(i in 1:length(mydata_c$date)){
  mydata_c$plot[i]<-
    if(mydata_c$year[i]!=2018){
      ifelse(mydata_c$system_id[i] %in% c('E1','E4','E7','E10','E13'),1,
             ifelse(mydata_c$system_id[i] %in% c('E2','E5','E8','E11','E14'),2,
                    ifelse(mydata_c$system_id[i] %in% c('E3','E6','E9','E12','E15'),3,
                           ifelse(mydata_c$system_id[i] %in% c('W1','W4','W7','W10','W13'),4,
                                  ifelse(mydata_c$system_id[i] %in% c('W2','W5','W8','W11','W14'),5,
                                         ifelse(mydata_c$system_id[i] %in% c('W3','W6','W9','W12','W15'),6,''))))))
    }else{
      ifelse(mydata_c$system_id[i] %in% c('E2','E3','E5','E6'),1,
             ifelse(mydata_c$system_id[i] %in% c('E12','E13','E15'),2,
                    ifelse(mydata_c$system_id[i] %in% c('E7','E9','E10'),3,
                           ifelse(mydata_c$system_id[i] %in% c('W1','W3','W5'),4,
                                  ifelse(mydata_c$system_id[i] %in% c('W12','W15','W16'),5,
                                         ifelse(mydata_c$system_id[i] %in% c('W6','W8','W9','W11'),6,''))))))
    }
  i=i+1
}

mydata_c$treatment<-ifelse(mydata_c$plot %in% c(1,4,5),'warmed','control')

names(mydata_h)
# 原始数据计算
mydata_d<-summaryBy(data=mydata_h,air_temp + N2O_N + NO_N + N2_N + soil_temp + soil_moisture~date+plot,FUN=myfun)


mydata_d<-mydata_d[,which(colnames(mydata_d) %in% c("date","plot","NO_N.m","N2O_N.m", "N2_N.m",
                                                    "soil_temp.m","soil_moisture.m","air_temp.m"))]



mydata_d$treatment<-ifelse(mydata_d$plot %in% c("1","4","5"),"warmed","control")
mydata_d$pair<-ifelse(mydata_d$plot %in% c("1","2"),"pair_1",
                      ifelse(mydata_d$plot %in% c("3","4"),"pair_2",
                             ifelse(mydata_d$plot %in% c("5","6"),"pair_3",NA)))

mydata_d$labs_warming<-
  ifelse(mydata_d$date>="2018/9/10"&mydata_d$date<="2018/11/15","warming",
         ifelse(mydata_d$date>="2019/6/8"&mydata_d$date<="2019/11/30","warming",
                ifelse(mydata_d$date>="2020/6/16"&mydata_d$date<="2020/8/26","warming",
                       ifelse(mydata_d$date>="2020/9/12"&mydata_d$date<="2020/12/6","warming",
                              ifelse(mydata_d$date>="2021/3/30"&mydata_d$date<="2021/12/7","warming",
                                     ifelse(mydata_d$date>="2022/3/26"&mydata_d$date<="2022/12/7","warming",
                                            ifelse(mydata_d$date>="2023/3/20"&mydata_d$date<="2023/12/7","warming",       
                                                   ifelse(mydata_d$date<="2018/9/9"&mydata_d$date>="2018/8/1","unwarming",
                                                          ifelse(mydata_d$date<="2019/6/7"&mydata_d$date>="2019/4/28","unwarming",
                                                                 ifelse(mydata_d$date<="2019/12/20"&mydata_d$date>="2019/12/1","unwarming",    
                                                                        ifelse(mydata_d$date<="2020/6/15"&mydata_d$date>="2020/3/29","unwarming",
                                                                               ifelse(mydata_d$date<="2020/12/20"&mydata_d$date>="2020/12/7","unwarming",  
                                                                                      ifelse(mydata_d$date<="2021/3/29"&mydata_d$date>="2021/3/8","unwarming",
                                                                                             ifelse(mydata_d$date<="2021/12/20"&mydata_d$date>="2021/12/8","unwarming",
                                                                                                    ifelse(mydata_d$date<="2022/3/26"&mydata_d$date>="2022/3/8","unwarming",  
                                                                                                           ifelse(mydata_d$date<="2023/3/20"&mydata_d$date>="2023/3/8","unwarming",
                                                                                                                  ifelse(mydata_d$date<="2023/12/20"&mydata_d$date>="2023/12/8","unwarming",NA
                                                                                                                  )))))))))))))))))

mydata_d<-mydata_d[which(mydata_d$date>="2018/8/1"),]
colnames(mydata_d)


unique(mydata_d$labs_warming)






###生长季土壤温度和湿度
mydata_h$month<-format(mydata_h$date,"%m")
unique(mydata_h$month)
mean(na.omit(mydata_h$soil_temp[mydata_h$month %in% c('05','06','07','08','09') & 
                                  mydata_h$treatment=='control' & mydata_h$year=='2019']))


###计算每个月每个样方的平均值
{str(mydata_d)
mydata_d$month<-format(mydata_d$date,"%m")
mydata_d$year<-format(mydata_d$date,"%Y")

colnames(mydata_d)
which(colnames(mydata_d) %in% c('NO_N.m','N2O_N.m','soil_temp.m','soil_moisture.m'))

mydata_d.m<-summaryBy(data=mydata_d,FUN=myfun,
                      mydata_d[,c(4,5,7,8)]~year+month+plot)  ###确认列数字


library(tidyverse);library(writexl)

mydata_d.m<-mydata_d.m[mydata_d.m$year %in% c(2021,2022,2023),]

mydata_d.m<-mydata_d.m[order(mydata_d.m$year,mydata_d.m$plot,mydata_d.m$month),]

###将不同的结果卸载一个xls文件
df = list(mydata_d.m,mydata_d) %>%
  set_names(c('mydata_d.m','mydata_d'))


getwd()
write_xlsx(df, "D:/工作目录/R plot/增温对NO+N2O的影响/mydata_monthly.xlsx")
}


###F.treat###
colnames(mydata_d)
which(colnames(mydata_d)%in% c("air_temp.m","N2O_N.m","NO_N.m", "soil_temp.m","soil_moisture.m","N2_N.m"))

names(mydata_d)[3:8]<-c("air_temp","N2O_N","NO_N","N2_N","soil_temp","soil_moisture")

mydata_d$WFPS<-mydata_d$soil_moisture/(1-0.7/2.65)

mydata_d$NXO=mydata_d$NO_N+mydata_d$N2O_N

mydata_d$gasN=mydata_d$NO_N+mydata_d$N2O_N+mydata_d$N2_N

mydata_d$NO_N2O=mydata_d$NO_N/mydata_d$N2O_N

mydata_d$year<-format(mydata_d$date,"%Y")

mydata_d$DOY<-as.numeric(format(mydata_d$date,"%j")) 

mydata_d$season<- ifelse(format(mydata_d$date,'%m') %in% c('03','04','05'),"spring",
                         ifelse(format(mydata_d$date,'%m') %in% c('06','07','08'),"summer",
                                ifelse(format(mydata_d$date,'%m') %in% c('09','10','11'),"fall","winter")))



mydata_d$warming<-
  ifelse(mydata_d$date>="2018/9/9" & mydata_d$date<="2018/11/15","Y",
         ifelse(mydata_d$date>="2019/6/9" & mydata_d$date<="2019/12/1","Y",
                ifelse(mydata_d$date>="2020/6/17" & mydata_d$date<="2020/8/20","Y",
                       ifelse(mydata_d$date>="2020/9/11" & mydata_d$date<="2020/12/7","Y",
                              ifelse(mydata_d$date>="2021/3/29"& mydata_d$date<="2021/12/7","Y",
                                     ifelse(mydata_d$date>="2022/3/26"& mydata_d$date<="2022/12/7","Y",
                                            ifelse(mydata_d$date>="2023/3/20"& mydata_d$date<="2023/12/8","Y",
                                                   ifelse(mydata_d$date>="2024/3/20"& mydata_d$date<="2024/12/8","Y",
                                                   "N"))))))))


str(mydata_d)
write.csv(mydata_d[,which(colnames(mydata_d) %in% 
                            c("date",'plot','treatment',
                              "NO_N","N2O_N",
                              "soil_temp","WFPS","air_temp"))],
          'D:/Workspace/book/Qingyuan/QY_Daily GHGs+NO flux_2019-2023.csv',row.names = F)


###新建文件夹(当前日期),存储图片
myname=format(Sys.time(),"%Y%m%d_%H%M");print(myname)    #get the local time
dir.create(paste("D:/工作目录/R plot/增温对NO+N2O的影响/",myname,sep=""))

setwd(paste("D:/工作目录/R plot/增温对NO+N2O的影响/",myname,sep=""))

getwd()

str(mydata_d)

mydata_d$N2_N2O.R=-1.1+0.13*mydata_d$WFPS


mydata_d$N2_N2N2O.R=1/(1+1/mydata_d$N2_N2O.R)


names(mydata_d)[c(3:8,14:17,21:22)]

F.treat<-summaryBy(data=mydata_d,mydata_d[c(3:8,14:17,21:22)]~date+treatment,FUN=myfun)

colnames(F.treat)
str(F.treat)

###设置date和pair范围
dates<-levels(as.factor(mydata_d$date))
colnames(mydata_d)

diff<-data.frame()   #warmed-control
warm_eff<-data.frame()   #(warmed-control)/control
k=1
max(dates)
for(i in dates){
  x=mydata_d[which(mydata_d$date==i &  mydata_d$treatment=="warmed"),]
  y=mydata_d[which(mydata_d$date==i & mydata_d$treatment=="control"),]
  if(length(x$date)==3 & length(y$date)==3){   #判断每天的向量为空集，空集不执行差值计算
    diff[k,1]=i
    diff[k,2]=mean(na.omit(x$air_temp))-mean(na.omit(y$air_temp))          #将空气温度进行差值运算
    diff[k,3]=mean(na.omit(x$N2O_N))-mean(na.omit(y$N2O_N))          #将N2O进行差值运算
    diff[k,4]=mean(na.omit(x$NO_N))-mean(na.omit(y$NO_N))          #将NO进行差值运算
    diff[k,5]=mean(na.omit(x$soil_temp))-mean(na.omit(y$soil_temp))          #将soil temp等进行差值运算
    diff[k,6]=mean(na.omit(x$soil_moisture))-mean(na.omit(y$soil_moisture))     #将soil moisture进行差值运算
    diff[k,7]="diff"                 # “差值”处理
    diff[k,8]=mean(na.omit(x$N2_N))-mean(na.omit(y$N2_N))               #N2_N差值运算
    diff[k,9]=mean(na.omit(x$NXO))-mean(na.omit(y$NXO))               #NO+N2O差值运算
    diff[k,10]=mean(na.omit(x$N2_N2O))-mean(na.omit(y$N2_N2O))               #N2_N2O比值差值运算
    diff[k,11]=mean(na.omit(x$WFPS))-mean(na.omit(y$WFPS))               #WFPS差值运算
    diff[k,12]=mean(na.omit(x$NO_N2O))-mean(na.omit(y$NO_N2O))               #N2_N2O比值差值运算
    diff[k,13]=mean(na.omit(x$N2_N2N2O))-mean(na.omit(y$N2_N2N2O))               #N2_N2N2O比值差值运算
    diff[k,14]=mean(na.omit(x$gasN))-mean(na.omit(y$gasN))               #N2_N2N2O比值差值运算
    
    diff[k,15:38]=NA
    
    ######计算warm_eff     
    warm_eff[k,1]=i
    warm_eff[k,2]=mean(na.omit(x$soil_temp))-mean(na.omit(y$soil_temp))    #Diff of waming soil temperature
    warm_eff[k,3]=mean(na.omit(x$air_temp))-mean(na.omit(y$air_temp))   #Diff of waming air temperature
    warm_eff[k,4]=mean(na.omit(x$soil_moisture))-mean(na.omit(y$soil_moisture))   #Diff of waming soil moisture
    warm_eff[k,5]=(mean(na.omit(x$soil_moisture))-mean(na.omit(y$soil_moisture)))/mean(na.omit(y$soil_moisture))   #Diff of waming soil moisture
    
    warm_eff[k,6]=(mean(na.omit(x$N2O_N))-mean(na.omit(y$N2O_N)))*100/mean(na.omit(y$N2O_N))         #N2O增温效果运算(w-c)/c*100%
    warm_eff[k,7]=(mean(na.omit(x$NO_N))-mean(na.omit(y$NO_N)))*100/mean(na.omit(y$NO_N))         #NO增温效果运算(w-c)/c*100%
    warm_eff[k,8]=(mean(na.omit(x$N2_N))-mean(na.omit(y$N2_N)))*100/mean(na.omit(y$N2_N))         #N2增温效果运算(w-c)/c*100%
    
    warm_eff[k,9]=mean(na.omit(y$soil_temp))    #background soil temp
    warm_eff[k,10]=mean(na.omit(y$soil_moisture))     #background soil moisture
    warm_eff[k,11]=mean(na.omit(x$N2O_N))/mean(na.omit(y$N2O_N))        #N2O通量比值w/c
    warm_eff[k,12]=mean(na.omit(x$NO_N))/mean(na.omit(y$NO_N))        #NO通量比值w/c
    warm_eff[k,13]=mean(na.omit(x$N2_N))/mean(na.omit(y$N2_N))        #N2通量比值w/c
    
    warm_eff[k,14]=mean(na.omit(x$WFPS))/mean(na.omit(y$WFPS))        #水分比值w/c
    warm_eff[k,15]=(warm_eff[k,6]+1)/warm_eff[k,14]    #每1%水分的w/c N2O
    warm_eff[k,16]=(warm_eff[k,7]+1)/warm_eff[k,14]    #每1%水分的w/c  NO
    warm_eff[k,17]=(warm_eff[k,8]+1)/warm_eff[k,14]    #每1%水分的w/c  N2
    
    warm_eff[k,18]=mean(na.omit(x$soil_temp))/mean(na.omit(y$soil_temp))                       #温度比值
    warm_eff[k,19]=(warm_eff[k,6]+1)/warm_eff[k,18]    #每1°下的w/c N2O
    warm_eff[k,20]=(warm_eff[k,7]+1)/warm_eff[k,18]    #每1°下的w/c NO
    warm_eff[k,21]=(warm_eff[k,8]+1)/warm_eff[k,18]    #每1°下的w/c N2
    warm_eff[k,22]=(mean(na.omit(x$NO_N2O))-mean(na.omit(y$NO_N2O)))*100/mean(na.omit(y$NO_N2O))  #NO/N2O 的增温效应
    warm_eff[k,23]=((mean(na.omit(x$N2O_N))+mean(na.omit(x$NO_N)))-(mean(na.omit(y$N2O_N))+mean(na.omit(y$NO_N))))*100/(mean(na.omit(y$N2O_N))+mean(na.omit(y$NO_N)))  #计算NO+N2O的增温效应
    warm_eff[k,24]=((mean(na.omit(x$N2O_N))+mean(na.omit(x$NO_N))+mean(na.omit(x$N2_N)))-(mean(na.omit(y$N2O_N))+mean(na.omit(y$NO_N))+mean(na.omit(x$N2_N))))*100/(mean(na.omit(y$N2O_N))+mean(na.omit(y$NO_N))+mean(na.omit(x$N2_N)))  #计算气态氮的增温效应
    
    warm_eff[k,25]=mean(na.omit(y$N2O_N))  #对照
    warm_eff[k,26]=mean(na.omit(y$NO_N))   #对照
    warm_eff[k,27]=mean(na.omit(y$N2_N))   #对照

    
    k=k+1  
  }
}

  
names(F.treat)

names(diff)[1:38]<-c("date","air_temp.m","N2O_N.m","NO_N.m","soil_temp.m",
                     "soil_moisture.m","treatment","N2_N.m","NXO.m","N2_N2O.R.m","WFPS.m","NO_N2O.m","N2_N2N2O.R.m", "gasN.m", 
                     "air_temp.n","air_temp.se","N2O_N.n","N2O_N.se",        
                     "NO_N.n", "NO_N.se", "soil_temp.n", "soil_temp.se",
                     "soil_moisture.n", "soil_moisture.se", "WFPS.n", "WFPS.se", "N2_N.n", "N2_N.se","NXO.n","NXO.se",
                     "NO_N2O.n","NO_N2O.se",  "N2_N2O.R.n", "N2_N2O.R.se",  "N2_N2N2O.R.n", "N2_N2N2O.R.se","gasN.n","gasN.se")  

names(warm_eff)[1:27]<-c("date","warm_st","warm_at","warm_sm","warm_sm.perc",
                         "eff_N2O","eff_NO","eff_N2",
                         "soil_temp","soil_moisture","ratio_N2O","ratio_NO","ratio_N2",
                         "wfps_ratio","ratio_N2O_w","ratio_NO_w","ratio_N2_w",
                         "temp_ratio","ratio_N2O_t","ratio_NO_t","ratio_N2_t","eff_NO_N2O",'eff_NO.N2O','eff_gasN',
                         'N2O_N','NO_N','N2_N')


names(diff);names(F.treat)

names(F.treat) %in% names(diff)

F.treat<-rbind(F.treat,diff);rm(diff)

str(warm_eff)
warm_eff$year<-format(as.Date(warm_eff$date),"%Y")

F.treat$year<-format(F.treat$date,'%Y')

myfun_mean<-function(x){m=mean(na.omit(x))}

warm_eff.m<-warm_eff

names(warm_eff.m)

str(warm_eff.m)
warm_eff.m$date<-as.Date(warm_eff.m$date)
warm_eff.m$year<-format(warm_eff.m$date,"%Y")
warm_eff.m$warming<-
  ifelse(warm_eff.m$date>="2018/9/9" & warm_eff.m$date<="2018/11/15","1",
         ifelse(warm_eff.m$date>="2019/6/9" & warm_eff.m$date<="2019/12/1","1",
                ifelse(warm_eff.m$date>="2020/6/17" & warm_eff.m$date<="2020/8/28","1",
                       ifelse(warm_eff.m$date>="2020/9/11" & warm_eff.m$date<="2020/12/7","1",
                              ifelse(warm_eff.m$date>="2021/3/30" & warm_eff.m$date<="2021/12/7","1",
                                     ifelse(warm_eff.m$date>="2022/3/26" & warm_eff.m$date<="2022/12/9","1",
                                            ifelse(warm_eff.m$date>="2023/3/20" & warm_eff.m$date<="2023/12/7","1",
                                                   ifelse(warm_eff.m$date>="2024/3/20" & warm_eff.m$date<="2024/12/7","1","0"))))))))

warm_eff.m$growing<-ifelse(format(warm_eff.m$date,'%m') %in% c('05','06','07','08','09'),'Y','N')

warm_eff.m<-within(warm_eff.m,year<-factor(warm_eff.m$year,levels = c('2024','2023',"2022","2021","2020","2019","2018")))


########### load QY rainfall Data
rainfall<-read.csv("D:/Workspace/book/Qingyuan/precipitation/QY air temp_precipitation_2024.csv",header=T,sep=",")
str(rainfall)

rainfall[,which(colnames(rainfall) %in% c("T_mean","Precipitation","VPD"))]<-
  lapply(rainfall[,which(colnames(rainfall) %in% c("T_mean","Precipitation","VPD"))], as.numeric)


rainfall$Date<-as.Date(rainfall$Date,tz="Asia/Taipei")
names(rainfall)[c(1,5)]<-c("date","precipitation");str(rainfall)

warm_eff.m<-merge(warm_eff.m,rainfall[,c(1,5)],by="date",all.x=T)
warm_eff.m$doy<-as.numeric(format(warm_eff.m$date,"%j"))
names(warm_eff.m)

#设置温度梯度分析
warm_eff.m$Tgradient<-ifelse(warm_eff.m$soil_temp<=10,'T1',
                             ifelse(warm_eff.m$soil_temp<=14,'T2',
                                    ifelse(warm_eff.m$soil_temp<=18,'T3','T4')))

warm_eff.m$month<-format(warm_eff.m$date,"%m")

warm_eff.m$season<-ifelse(warm_eff.m$month %in% c('05','06','07','08','09'),"growing",
                          ifelse(warm_eff.m$month %in% c('03','04'),"spring-thaw",
                                 ifelse(warm_eff.m$month %in% c('10','11'),"fall","winter")))


###### Fig 0 土壤0-10cm增温效果 ##########
unique(F.treat$treatment)
str(F.treat)
F.treat$DOY<-as.numeric(format(F.treat$date,"%j")) 

F.treat$season<-ifelse(format(F.treat$date,"%m") %in% c('05','06','07','08','09'),"growing",
                       ifelse(format(F.treat$date,"%m") %in% c('03','04'),"spring-thaw",
                              ifelse(format(F.treat$date,"%m") %in% c('10','11'),"fall","winter")))
F.treat$warming<-
  ifelse(F.treat$date>="2018/9/9" & F.treat$date<="2018/11/15","Y",
         ifelse(F.treat$date>="2019/6/9" & F.treat$date<="2019/12/1","Y",
                ifelse(F.treat$date>="2020/6/17" & F.treat$date<="2020/8/20","Y",
                       ifelse(F.treat$date>="2020/9/11" & F.treat$date<="2020/12/7","Y",
                              ifelse(F.treat$date>="2021/3/29"& F.treat$date<="2021/12/7","Y",
                                     ifelse(F.treat$date>="2022/3/26"& F.treat$date<="2022/12/7","Y",
                                            ifelse(F.treat$date>="2023/3/20"& F.treat$date<="2023/12/7","Y",
                                                   ifelse(F.treat$date>="2024/3/20"& F.treat$date<="2024/12/7","Y",
                                                   "N"))))))))


F.pt<-merge(F.treat,rainfall,by="date",all.x=T)   ###合并降水通量数据

library(ggprism)
pacman::p_load("ggplot2","ggpmisc","ggpubr","ggprism","doBy","car","readr","crayon","rstatix","plot3D",
              "lme4","MCMCglmm","tidyverse","directlabels","agricolae","lubridate","scales",
              "patchwork","reshape2","dplyr","animation")  #加载多个包

st<-ggplot(F.pt,aes(date,soil_temp.m,col=treatment))+
  annotate("rect",xmin=as.Date("2018-09-09"),
           xmax=as.Date("2018-11-16"),
           ymin=-10,ymax=30,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2019-06-08"),
           xmax=as.Date("2019-12-01"),
           ymin=-10,ymax=30,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2020-06-16"),
           xmax=as.Date("2020-8-26"),
           ymin=-10,ymax=30,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2020-09-12"),
           xmax=as.Date("2020-12-8"),
           ymin=-10,ymax=30,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2021-03-29"),
           xmax=as.Date("2021-12-7"),
           ymin=-10,ymax=30,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2022-03-26"),
           xmax=as.Date("2022-12-10"),
           ymin=-10,ymax=30,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2023-03-20"),
           xmax=as.Date("2023-12-08"),
           ymin=-10,ymax=30,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2024-03-20"),
           xmax=as.Date("2024-12-07"),
           ymin=-Inf,ymax=Inf,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2025-03-27"),
           xmax=as.Date(max(F.pt$date+1)),
           ymin=-Inf,ymax=Inf,alpha=.1,fill="red")+
  geom_hline(yintercept =2,linewidth =0.35,color="red",lty=2)+
  # geom_hline(yintercept =0,size =0.35,color="red",lty=2)+
  #geom_point(size=1.2)+
  geom_line(linewidth=0.45,show.legend = F)+
  geom_errorbar(aes(ymin=soil_temp.m-soil_temp.se, ymax=soil_temp.m+soil_temp.se),stat = "identity",
                linewidth=0.2,width=0,position=position_dodge(0),alpha=.2)+
  labs(y=expression("Soil temperature (\u00B0C)"),
       x=expression(paste("Date (yy/mm/dd)")))+
  scale_x_date(limits = c(as.Date("2018/8/1"),max(F.pt$date+1)),date_labels = "%Y/%m",
               date_minor_breaks="months",breaks="3 months",expand = c(0,0))+
  # scale_x_date(limits = c(as.Date("2018/8/1"),as.Date("2021/12/31")),date_labels = "%Y/%m",
  #              date_minor_breaks="months",breaks="3 months",expand = c(0,0))+
  scale_y_continuous(limits=c(-10,30),breaks=seq(-10,30,10),minor_breaks = seq(-10,30,2),expand=c(0,0))+
  scale_color_manual(values=c("red","blue","black"),limits=c("warmed","control","diff"),name="Treatment",
                     labels=c("Warmed","Control","Difference"))+
  geom_line(aes(y=35),linewidth=1.5,show.legend = T)+
  guides(y="prism_offset_minor",x="prism_offset_minor",
         color=guide_legend(override.aes=list(size=4),keyheight = 1,keywidth = 0.5,ncol = 1,byrow = FALSE))+
  theme(plot.subtitle = element_text(vjust = 1),plot.caption = element_text(vjust = 1), 
        plot.background   = element_rect(fill = NA,  colour = NA, linewidth = 1),
        panel.grid.major = element_blank(), panel.grid=element_blank(), 
        panel.background = element_rect(fill =NA,colour = "black", linewidth=0.6,linetype = "solid"),
        panel.border = element_rect(linewidth=0.6,fill=NA),
        axis.text.x = element_text(colour = "black",size=16,angle = 35,vjust=1,hjust=1), 
        axis.text.y = element_text(colour = "black",size=16), 
        axis.title = element_text(colour = "black",size=18),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black",size=16),
        legend.key = element_rect(fill = NA, color = NA), 
        legend.background = element_rect(fill = NA), 
        legend.position =c(1.05/10,9/10),
        axis.ticks = element_line(linewidth = 0.6),axis.ticks.length = unit(.2, 'cm'),
        plot.margin = unit(c(1,1,0,1),'lines'),
        strip.text.x = element_text(size=14, face="bold"),
        strip.background = element_rect(colour="black"));st


###only warming effect
p1<-ggplot(F.pt[F.pt$treatment!='diff',],aes(date,soil_temp.m,col=treatment))+
  annotate("rect",xmin=as.Date("2018-09-09"),
           xmax=as.Date("2018-11-16"),
           ymin=-10,ymax=30,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2019-06-08"),
           xmax=as.Date("2019-12-01"),
           ymin=-10,ymax=30,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2020-06-16"),
           xmax=as.Date("2020-8-26"),
           ymin=-10,ymax=30,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2020-09-12"),
           xmax=as.Date("2020-12-8"),
           ymin=-10,ymax=30,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2021-03-29"),
           xmax=as.Date("2021-12-7"),
           ymin=-10,ymax=30,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2022-03-26"),
           xmax=as.Date("2022-12-10"),
           ymin=-10,ymax=30,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2023-03-20"),
           xmax=as.Date("2023-12-08"),
           ymin=-10,ymax=30,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2024-03-20"),
           xmax=as.Date("2024-12-07"),
           ymin=-Inf,ymax=Inf,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2025-03-27"),
           xmax=as.Date(max(F.pt$date+1)),
           ymin=-Inf,ymax=Inf,alpha=.1,fill="red")+
  # geom_hline(yintercept =0,size =0.35,color="red",lty=2)+
  #geom_point(size=1.2)+
  geom_line(linewidth=0.45,show.legend = F)+
  geom_errorbar(aes(ymin=soil_temp.m-soil_temp.se, ymax=soil_temp.m+soil_temp.se),stat = "identity",
                linewidth=0.2,width=0,position=position_dodge(0),alpha=.2)+
  labs(y=expression("Soil temperature (\u00B0C)"),
       x=expression(paste("Date (yy/mm/dd)")))+
  scale_x_date(limits = c(as.Date("2018/7/1"),max(F.pt$date+1)), date_minor_breaks="months",
               breaks=c(as.Date('2018/7/1'),seq(as.Date('2019/1/1'),max(F.pt$date+1),'12 months')),
               date_label='%Y', expand = c(0,0))+
  scale_y_continuous(limits=c(-10,30),breaks=seq(-10,30,10),minor_breaks = seq(-10,30,2),expand=c(0,0))+
  scale_color_manual(values=c("red","blue","black"),limits=c("warmed","control","diff"),name="Treatment",
                     labels=c("Warmed","Control","Difference"))+
  geom_line(aes(y=35),linewidth=1.5,show.legend = T)+
  guides(y="prism_offset_minor",x="prism_offset_minor",
         color=guide_legend(override.aes=list(size=4),keyheight = 1,keywidth = 0.5,ncol = 1,byrow = FALSE))+
  theme(plot.subtitle = element_text(vjust = 1),plot.caption = element_text(vjust = 1), 
        plot.background   = element_rect(fill = NA,  colour = NA, linewidth = 1),
        panel.grid.major = element_blank(), panel.grid=element_blank(), 
        panel.background = element_rect(fill =NA,colour = "black", linewidth=0.6,linetype = "solid"),
        panel.border = element_rect(linewidth=0.6,fill=NA),
        axis.text.x = element_text(colour = "black",size=16,angle = 35,vjust=1,hjust=1), 
        axis.text.y = element_text(colour = "black",size=16), 
        axis.title = element_text(colour = "black",size=18),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black",size=16),
        legend.key = element_rect(fill = NA,color=NA), 
        legend.background = element_rect(fill = NA), 
        legend.position =c(1.00/10,9/10),
        axis.ticks = element_line(linewidth = 0.6),axis.ticks.length = unit(.2, 'cm'),
        plot.margin = unit(c(1,1,0,1),'lines'),
        strip.text.x = element_text(size=14, face="bold"),
        strip.background = element_rect(colour="black"));p1


d1<-ggplot(F.pt[F.pt$treatment=='diff',],aes(date,soil_temp.m,col=treatment))+
  annotate("rect",xmin=as.Date("2018-09-09"),
           xmax=as.Date("2018-11-16"),
           ymin=-Inf,ymax=Inf,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2019-06-08"),
           xmax=as.Date("2019-12-01"),
           ymin=-Inf,ymax=Inf,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2020-06-16"),
           xmax=as.Date("2020-8-26"),
           ymin=-Inf,ymax=Inf,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2020-09-12"),
           xmax=as.Date("2020-12-8"),
           ymin=-Inf,ymax=Inf,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2021-03-29"),
           xmax=as.Date("2021-12-7"),
           ymin=-Inf,ymax=Inf,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2022-03-26"),
           xmax=as.Date("2022-12-10"),
           ymin=-Inf,ymax=Inf,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2023-03-20"),
           xmax=as.Date("2023-12-08"),
           ymin=-Inf,ymax=Inf,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2024-03-20"),
           xmax=as.Date("2024-12-07"),
           ymin=-Inf,ymax=Inf,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2025-03-27"),
           xmax=as.Date(max(F.pt$date+1)),
           ymin=-Inf,ymax=Inf,alpha=.1,fill="red")+
  geom_hline(yintercept =2,linewidth =0.35,color="red",lty=2)+
  # geom_hline(yintercept =0,size =0.35,color="red",lty=2)+
  #geom_point(size=1.2)+
  geom_line(linewidth=0.45,show.legend = F)+
  geom_errorbar(aes(ymin=soil_temp.m-soil_temp.se, ymax=soil_temp.m+soil_temp.se),stat = "identity",
                linewidth=0.2,width=0,position=position_dodge(0),alpha=.2)+
  labs(y=expression("Warming effect (\u00B0C)"),
       x='Year')+
  scale_x_date(limits = c(as.Date("2018/7/1"),max(F.pt$date+1)), date_minor_breaks="months",
               breaks=c(as.Date('2018/7/1'),seq(as.Date('2019/1/1'),max(F.pt$date+1),'12 months')),
               date_label='%Y', expand = c(0,0))+
  scale_y_continuous(breaks=seq(-2,6,2),minor_breaks = seq(-2,6,0.5),expand=c(0,0))+
  scale_color_manual(values=c("black"),limits=c("diff"),name="Treatment",labels=c("Difference"))+
  geom_line(aes(y=35),linewidth=1.5,show.legend = T)+
  guides(y="prism_offset_minor",x="prism_offset_minor",
         color=guide_legend(override.aes=list(size=4),keyheight = 1,keywidth = 0.5,ncol = 1,byrow = FALSE))+
  theme(plot.subtitle = element_text(vjust = 1),plot.caption = element_text(vjust = 1), 
        plot.background   = element_rect(fill = NA,  colour = NA, linewidth = 0.2),
        panel.grid.major = element_blank(), panel.grid=element_blank(), 
        panel.background = element_rect(fill =NA,colour = "black", linewidth=0.26,linetype = "solid"),
        panel.border = element_rect(linewidth=0.6,fill=NA),
        axis.text.x = element_text(colour = "black",size=16,angle = 0,vjust=0,hjust=0.5), 
        axis.text.y = element_text(colour = "black",size=16), 
        axis.title = element_text(colour = "black",size=18),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black",size=16),
        legend.key = element_rect(fill = NA, color = NA), 
        legend.background = element_rect(fill = NA), 
        legend.position ='none',
        axis.ticks = element_line(linewidth = 0.6),axis.ticks.length = unit(.2, 'cm'),
        plot.margin = unit(c(0.8,1,1,1),'lines'),   #plot.margin = unit(rep(2,4),'lines'),
        strip.text.x = element_text(size=14, face="bold"),
        strip.background = element_rect(colour="black"))+
  coord_cartesian(clip="off",ylim=c(-2,6));d1

(p1+rremove('xlab')+rremove('x.text')+d1)+
  plot_layout(ncol=1,nrow=2,widths=c(5,5),heights=c(5,3))&
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        strip.text = element_text(size=14))->p;p
        # plot.margin = unit(c(4,2,1,3),"mm")


ggsave('Warming effect of temp.pdf',p, width =10, height =7, device=cairo_pdf)


{p1<-st+
  scale_y_continuous(breaks=seq(-10,30,10),minor_breaks = seq(-10,30,2),expand=c(0,0))+
  scale_x_date(limits = c(as.Date("2018/8/1"),max(F.pt$date+1)),date_labels = "%m",
               date_minor_breaks="months",breaks="3 months",expand = c(0,0),name='')+
  guides(y="prism_offset_minor",x="prism_offset_minor")+
  coord_cartesian(clip="off",ylim=c(-10,30))+
  theme(plot.subtitle = element_text(vjust = 1),plot.caption = element_text(vjust = 1), 
        plot.background   = element_rect(fill = NA,  colour = NA, size = 1),
        panel.grid.major = element_blank(), panel.grid=element_blank(), 
        panel.background = element_rect(fill =NA,colour = "black", size=0.6,linetype = "solid"),
        panel.border = element_rect(size=0.6,fill=NA),
        axis.text.x = element_text(colour = "black",size=16,angle = 0,vjust=0.5,hjust=0.5), 
        axis.text.y = element_text(colour = "black",size=16), 
        axis.title = element_text(colour = "black",size=18),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black",size=16),
        legend.key = element_rect(fill = NA,colour = NA), 
        legend.background = element_rect(fill = NA), 
        axis.ticks = element_line(linewidth = 0.6),
        axis.ticks.length = unit(.2, 'cm'),
        plot.margin = unit(c(0.8,1,2.5,1),'lines'),   #plot.margin = unit(rep(2,4),'lines'),
        strip.text.x = element_text(size=14, face="bold"),
        strip.background = element_rect(colour="black"))+
  annotate(geom="text",x=as.Date("2018/9/1"), y=-14.45,label="2018",colour="black",vjust=0,size=6)+
  annotate(geom="text",x=as.Date("2019/7/1"), y=-14.45,label="2019",colour="black",vjust=0,size=6)+
  annotate(geom="text",x=as.Date("2020/7/1"), y=-14.45,label="2020",colour="black",vjust=0,size=6)+
  annotate(geom="text",x=as.Date("2021/7/1"), y=-14.45,label="2021",colour="black",vjust=0,size=6)+
  annotate(geom="text",x=as.Date("2022/7/1"), y=-14.45,label="2022",colour="black",vjust=0,size=6)+
  annotate(geom="text",x=as.Date("2023/7/1"), y=-14.45,label="2023",colour="black",vjust=0,size=6)+
  annotate(geom="text",x=as.Date("2024/6/1"), y=-14.45,label="2024",colour="black",vjust=0,size=6)+
  ###加箭头
  geom_segment(aes(x=as.Date("2024/4/1"), xend=as.Date("2024/1/1"),y=-14, yend=-14),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2024/9/1"), xend=as.Date("2025/1/1"),y=-14, yend=-14),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2023/5/1"), xend=as.Date("2023/1/1"),y=-14, yend=-14),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2023/9/1"), xend=as.Date("2024/1/1"),y=-14, yend=-14),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2022/5/1"), xend=as.Date("2022/1/1"),y=-14, yend=-14),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2022/9/1"), xend=as.Date("2023/1/1"),y=-14, yend=-14),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2021/5/1"), xend=as.Date("2021/1/1"),y=-14, yend=-14),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2021/9/1"), xend=as.Date("2022/1/1"),y=-14, yend=-14),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2020/5/1"), xend=as.Date("2020/1/1"),y=-14, yend=-14),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2020/9/1"), xend=as.Date("2021/1/1"),y=-14, yend=-14),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2019/9/1"), xend=as.Date("2020/1/1"),y=-14, yend=-14), 
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2019/5/1"), xend=as.Date("2019/1/1"),y=-14, yend=-14),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2018/11/1"), xend=as.Date("2019/1/1"),y=-14, yend=-14), 
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2018/5/1"), xend=as.Date("2018/1/1"),y=-14, yend=-14),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  ###加年份分割线
  # geom_segment(x=as.Date("2018/1/1"), xend=as.Date("2018/1/1"),
  #              y=-150, yend=-210,colour="black",size=.5)+
  geom_segment(aes(x=as.Date("2019/1/1"), xend=as.Date("2019/1/1"),
                   y=-13, yend=-15),colour="black",size=.5)+
  geom_segment(aes(x=as.Date("2020/1/1"), xend=as.Date("2020/1/1"),
                   y=-13, yend=-15),colour="black",size=.5)+
  geom_segment(aes(x=as.Date("2021/1/1"), xend=as.Date("2021/1/1"),
                   y=-13, yend=-15),colour="black",size=.5)+
  geom_segment(aes(x=as.Date("2022/1/1"), xend=as.Date("2022/1/1"),
                   y=-13, yend=-15),colour="black",size=.5)+
  geom_segment(aes(x=as.Date("2023/1/1"), xend=as.Date("2023/1/1"),
                   y=-13, yend=-15),colour="black",size=.5)+
  geom_segment(aes(x=as.Date("2024/1/1"), xend=as.Date("2024/1/1"),
                   y=-13, yend=-15),colour="black",size=.5);p1

myname=paste("st",format(Sys.time(),"%Y%m%d"),sep="_");print(myname)    #获取当前时间

ggsave(paste(myname,"pdf",sep="."),p1, height = 6, width = 12,
       device=cairo_pdf)
}

str(F.pt)

sm<-ggplot(F.pt[F.pt$WFPS.m<=70 & F.pt$WFPS.m >= -10,],aes(date,WFPS.m,col=treatment))+
  annotate("rect",xmin=as.Date("2018-09-09"),
           xmax=as.Date("2018-11-16"),
           ymin=-10,ymax=70,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2019-06-08"),
           xmax=as.Date("2019-12-01"),
           ymin=-10,ymax=70,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2020-06-16"),
           xmax=as.Date("2020-8-26"),
           ymin=-10,ymax=70,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2020-09-12"),
           xmax=as.Date("2020-12-8"),
           ymin=-10,ymax=70,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2021-03-29"),
           xmax=as.Date("2021-12-7"),
           ymin=-10,ymax=70,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2022-03-26"),
           xmax=as.Date("2022-12-10"),
           ymin=-10,ymax=70,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2023-03-20"),
           xmax=as.Date("2023-12-08"),
           ymin=-10,ymax=70,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2024-03-20"),
           xmax=as.Date("2024-12-07"),
           ymin=-Inf,ymax=Inf,alpha=.1,fill="red")+
  annotate("rect",xmin=as.Date("2025-03-27"),
           xmax=as.Date(max(F.pt$date+1)),
           ymin=-Inf,ymax=Inf,alpha=.1,fill="red")+
  geom_hline(yintercept =0,size =0.25,color="red",lty=2)+
  geom_line(linewidth=0.45,show.legend = T)+
  geom_errorbar(aes(ymin=WFPS.m-WFPS.se, ymax=WFPS.m+WFPS.se),stat = "identity", 
                size=0.2,width=0,position=position_dodge(0),alpha=.2)+
  labs(y=expression("WFPS (0-10 cm, %)"),x='')+
  scale_x_date(limits = c(as.Date("2018/8/1"),max(F.pt$date+1)),date_labels = "%m",
               date_minor_breaks="months",breaks="3 months",expand = c(0,0))+
  # scale_x_date(limits = c(as.Date("2018/8/1"),as.Date("2021/12/31")),date_labels = "%Y/%m",
  #              date_minor_breaks="months",breaks="3 months",expand = c(0,0))+
  scale_y_continuous(breaks=seq(-10,70,20),minor_breaks = seq(-10,70,5),expand=c(0,0))+
  scale_color_manual(values=c("red","blue","black"),limits=c("warmed","control","diff"),name="Treatment",
                     labels=c("warmed","control","difference"))+
  guides(y="prism_offset_minor",x="prism_offset_minor")+
  theme(plot.subtitle = element_text(vjust = 1),plot.caption = element_text(vjust = 1), 
        plot.background   = element_rect(fill = NA,  colour = NA, size = 1),
        panel.grid.major = element_blank(), panel.grid=element_blank(), 
        panel.background = element_rect(fill =NA,colour = "black", size=0.6,linetype = "solid"),
        panel.border = element_rect(size=0.6,fill=NA),
        axis.text.x = element_text(colour = "black",size=16,angle = 0,vjust=0.5,hjust=0.5), 
        axis.text.y = element_text(colour = "black",size=16), 
        axis.title = element_text(colour = "black",size=18),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black",size=16),
        legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), 
        legend.position ="none",axis.ticks = element_line(linewidth = 0.6),
        axis.ticks.length = unit(.2, 'cm'),
        plot.margin = unit(c(0.8,1,2.5,1),'lines'),   #plot.margin = unit(rep(2,4),'lines'),
        strip.text.x = element_text(size=14, face="bold"),
        strip.background = element_rect(colour="black"))+
  coord_cartesian(clip="off",ylim=c(-10,70))+
  ###加标注文本
  annotate(geom="text",x=as.Date("2018/9/1"), y=-22,label="2018",colour="black",vjust=0,size=6)+
  annotate(geom="text",x=as.Date("2019/7/1"), y=-22,label="2019",colour="black",vjust=0,size=6)+
  annotate(geom="text",x=as.Date("2020/7/1"), y=-22,label="2020",colour="black",vjust=0,size=6)+
  annotate(geom="text",x=as.Date("2021/7/1"), y=-22,label="2021",colour="black",vjust=0,size=6)+
  annotate(geom="text",x=as.Date("2022/7/1"), y=-22,label="2022",colour="black",vjust=0,size=6)+
  annotate(geom="text",x=as.Date("2023/7/1"), y=-22,label="2023",colour="black",vjust=0,size=6)+
  annotate(geom="text",x=as.Date("2024/6/1"), y=-22,label="2024",colour="black",vjust=0,size=6)+
  ###加箭头
  geom_segment(aes(x=as.Date("2024/4/1"), xend=as.Date("2024/1/1"),y=-20, yend=-20),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2024/9/1"), xend=as.Date("2024/12/20"),y=-20, yend=-20),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2023/5/1"), xend=as.Date("2023/1/1"),y=-20, yend=-20),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2023/9/1"), xend=as.Date("2024/1/1"),y=-20, yend=-20),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2022/5/1"), xend=as.Date("2022/1/1"),y=-20, yend=-20),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2022/9/1"), xend=as.Date("2023/1/1"),y=-20, yend=-20),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2021/5/1"), xend=as.Date("2021/1/1"),y=-20, yend=-20),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2021/9/1"), xend=as.Date("2022/1/1"),y=-20, yend=-20),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2020/5/1"), xend=as.Date("2020/1/1"),y=-20, yend=-20),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2020/9/1"), xend=as.Date("2021/1/1"),y=-20, yend=-20),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2019/9/1"), xend=as.Date("2020/1/1"),y=-20, yend=-20), 
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2019/5/1"), xend=as.Date("2019/1/1"),y=-20, yend=-20),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2018/11/1"), xend=as.Date("2019/1/1"),y=-20, yend=-20), 
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2018/5/1"), xend=as.Date("2018/1/1"),y=-20, yend=-20),
               arrow=arrow(length = unit(0.03, "npc")),colour="black",size=.5,show.legend = F)+
  ###加年份分割线
  # geom_segment(x=as.Date("2018/1/1"), xend=as.Date("2018/1/1"),
  #              y=-150, yend=-210,colour="black",size=.5)+
  geom_segment(aes(x=as.Date("2019/1/1"), xend=as.Date("2019/1/1"),
                   y=-18, yend=-22),colour="black",size=.5)+
  geom_segment(aes(x=as.Date("2020/1/1"), xend=as.Date("2020/1/1"),
                   y=-18, yend=-22),colour="black",size=.5)+
  geom_segment(aes(x=as.Date("2021/1/1"), xend=as.Date("2021/1/1"),
                   y=-18, yend=-22),colour="black",size=.5)+
  geom_segment(aes(x=as.Date("2022/1/1"), xend=as.Date("2022/1/1"),
                   y=-18, yend=-22),colour="black",size=.5)+
  geom_segment(aes(x=as.Date("2023/1/1"), xend=as.Date("2023/1/1"),
                   y=-18, yend=-22),colour="black",size=.5)+
  geom_segment(aes(x=as.Date("2024/1/1"), xend=as.Date("2024/1/1"),
                   y=-18, yend=-22),colour="black",size=.5);sm


library(ggpubr)
figure<-ggarrange(st+rremove("x.text")+rremove("xlab"),sm,
                  labels = c("(a)", "(b)"),
                  label.x = 0.063,label.y = 0.97,
                  ncol = 1, nrow = 2,align = "v",   ##"v"竖直对齐
                  font.label = list(size = 20, color ="black"),
                  widths = c(6,6), heights = c(6,7.5),
                  common.legend=F);figure

myname=paste("st+sm",format(Sys.time(),"%Y%m%d"),sep="_");print(myname)    #获取当前时间

ggsave(paste(myname,"pdf",sep="."),figure, height = 9, width = 13.0,
       device=cairo_pdf)


rainfall$glab<-ifelse(rainfall$precipitation==0,'Y','N')

str(rainfall)
rainfall$T_mean<-as.numeric(rainfall$T_mean)

a<-ggplot(data=rainfall,aes(x=date,y=precipitation))+
  geom_line(data=rainfall,aes(x=as.Date('2017/1/1'),y=Inf,col=glab),size=1.3)+
  scale_color_manual(limits=c('Y','N'),values = c("red","blue"),labels=c("Air temp","Precipitation"))+
  geom_bar(position=position_dodge(),stat="identity",size=0.1,fill="blue",col="blue",alpha=.5)+
  geom_line(data=rainfall[which(rainfall$T_mean>= -25 & rainfall$T_mean<= 25),],
            aes(x=date,y=(T_mean+25)*2),size=0.8,color="red")+
  scale_x_date(limits = c(as.Date("2018/1/1"),as.Date("2024/12/27")),name=NULL,
               date_labels = rep(c('1','4','7','10'),6),
               date_minor_breaks="months",breaks="3 months",expand = c(0,0))+
  guides(x="prism_offset_minor",y="prism_offset_minor")+
  mythem+theme(plot.margin = unit(c(6,1,15,2),'mm'),legend.position = c(0.18,0.97),
               axis.text = element_text(colour = "black",size=16),
               axis.title = element_text(colour = "black",size=16),
               legend.text = element_text(size=16),
               legend.background = element_blank(),legend.title = element_blank())+  
  scale_y_continuous(breaks=seq(0,100,20),expand = c(0,0),
                     labels=seq(-25,25,10),
                     name=expression(paste('Air temperature (\u00B0C)')),
                     sec.axis=sec_axis(~./2-25,breaks=seq(-25,25,10),labels=seq(0,100,20),
                                       name=expression(paste('Precipitation (mm)'))))+
  coord_cartesian(clip="off",ylim=c(0,100))+   #在plot绘图区域外面标注的关键
  ###加标注文本
  annotate(geom="text",x=as.Date("2018/7/1"), y=0,label="2018",colour="black",vjust=3.2,size=7)+
  annotate(geom="text",x=as.Date("2019/7/1"), y=0,label="2019",colour="black",vjust=3.2,size=7)+
  annotate(geom="text",x=as.Date("2020/7/1"), y=0,label="2020",colour="black",vjust=3.2,size=7)+
  annotate(geom="text",x=as.Date("2021/7/1"), y=0,label="2021",colour="black",vjust=3.2,size=7)+
  annotate(geom="text",x=as.Date("2022/7/1"), y=0,label="2022",colour="black",vjust=3.2,size=7)+
  annotate(geom="text",x=as.Date("2023/7/1"), y=0,label="2023",colour="black",vjust=3.2,size=7)+
  annotate(geom="text",x=as.Date("2024/7/1"), y=0,label="2024",colour="black",vjust=3.2,size=7)+
  
  ###加箭头
  geom_segment(aes(x=as.Date("2024/5/1"), xend=as.Date("2024/1/1"),y=-8, yend=-8),
               arrow=arrow(length = unit(0.01, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2024/9/1"), xend=as.Date("2024/12/27"),y=-8, yend=-8),
               arrow=arrow(length = unit(0.01, "npc")),size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2023/5/1"), xend=as.Date("2023/1/1"),y=-8, yend=-8),
               arrow=arrow(length = unit(0.01, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2023/9/1"), xend=as.Date("2024/1/1"),y=-8, yend=-8),
               arrow=arrow(length = unit(0.01, "npc")),size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2022/5/1"), xend=as.Date("2022/1/1"),y=-8, yend=-8),
               arrow=arrow(length = unit(0.01, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2022/9/1"), xend=as.Date("2023/1/1"),y=-8, yend=-8),
               arrow=arrow(length = unit(0.01, "npc")),size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2021/5/1"), xend=as.Date("2021/1/1"),y=-8, yend=-8),
               arrow=arrow(length = unit(0.01, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2021/9/1"), xend=as.Date("2022/1/1"),y=-8, yend=-8),
               arrow=arrow(length = unit(0.01, "npc")),size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2020/5/1"), xend=as.Date("2020/1/1"),y=-8, yend=-8),
               arrow=arrow(length = unit(0.01, "npc")),colour="black",size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2020/9/1"), xend=as.Date("2021/1/1"),y=-8, yend=-8),
               arrow=arrow(length = unit(0.01, "npc")),size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2019/9/1"), xend=as.Date("2020/1/1"),y=-8, yend=-8), 
               arrow=arrow(length = unit(0.01, "npc")),size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2019/5/1"), xend=as.Date("2019/1/1"),y=-8, yend=-8),
               arrow=arrow(length = unit(0.01, "npc")),size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2018/9/1"), xend=as.Date("2019/1/1"),y=-8, yend=-8), 
               arrow=arrow(length = unit(0.01, "npc")),size=.5,show.legend = F)+
  geom_segment(aes(x=as.Date("2018/5/1"), xend=as.Date("2018/1/1"),y=-8, yend=-8),
               arrow=arrow(length = unit(0.01, "npc")),size=.5,show.legend = F)+
  ###加年份分割线
  geom_segment(x=as.Date("2018/1/1"), xend=as.Date("2018/1/1"),
               y=-6, yend=-10,size=.5)+
  geom_segment(x=as.Date("2019/1/1"), xend=as.Date("2019/1/1"),
               y=-6, yend=-10,size=.5)+
  geom_segment(x=as.Date("2020/1/1"), xend=as.Date("2020/1/1"),
               y=-6, yend=-10,size=.5)+
  geom_segment(x=as.Date("2021/1/1"), xend=as.Date("2021/1/1"),
               y=-6, yend=-10,size=.5)+
  geom_segment(x=as.Date("2022/1/1"), xend=as.Date("2022/1/1"),
               y=-6, yend=-10,size=.5)+
  geom_segment(x=as.Date("2023/1/1"), xend=as.Date("2023/1/1"),
               y=-6, yend=-10,size=.5)+
  geom_segment(x=as.Date("2024/1/1"), xend=as.Date("2024/1/1"),
               y=-6, yend=-10,size=.5)+
  geom_segment(x=as.Date("2024/12/27"), xend=as.Date("2024/12/27"),
               y=-6, yend=-10,size=.5);a

ggsave("QY temp-rainfall.pdf",a, height = 8, width = 12,
       device=cairo_pdf)


unique(F.pt$treatment)

getwd()

###### Fig_1_Seasonal NO and N2O fluxes under 6-year warming ##########
add<-F.pt[length(F.pt$soil_temp.m)+1,]
add$treatment="precipitation"
datai<-rbind(F.pt,add)   #制作降水图例


###log plot
F.pt$N2O_min=F.pt$N2O_N.m-F.pt$N2O_N.se
F.pt$N2O_max=F.pt$N2O_N.m+F.pt$N2O_N.se

dat<-F.pt[which(F.pt$date>="2018-7-8"),] #&F.pt$treatment!='diff'

str(dat)
unique(dat$labs_warming)
figure<-
  ggplot(dat[dat$year%in% c(2018,2019,2020,2021,2022,2023,2024) & dat$season=='growing',],
         aes(N2O_N.m,..density..,fill=treatment))+
  geom_histogram(col='black',binwidth = 2,alpha=0.5)+
  geom_vline(xintercept = 0,lty=2,col='black')+
  # geom_text(aes(label=as.character(round(..density..,2))),
  #           stat="bin",binwidth=0.5,boundary = 0,vjust=-0.5)+
  scale_x_continuous(limits = c(-10,40),expand=c(0.08,0.08),
                     name = expression(paste('N'[2],'O flux (μg N m'^-2,' h'^-1,')')))+
  scale_y_continuous(limits = c(0,0.17),breaks = seq(0,0.15,0.05),expand=c(0,0),name = 'Density')+
  scale_fill_manual(values=c("blue","red"),name="",
                    limits=c("control","warmed"),
                    labels=c("Control","Warmed"))+
  facet_grid(treatment~year)+mythem;figure

ggsave("频率分布-N2O flux.pdf",figure,height = 7, width = 16)

length(dat$N2O_N.m[dat$year%in% c(2023) & dat$season=='growing' & dat$N2O_N.m<0])/
  length(dat$N2O_N.m[dat$year%in% c(2023) & dat$season=='growing'])

length(dat$N2O_N.m[dat$year%in% c(2023) & dat$season=='growing' & dat$N2O_N.m<0.1])
#总体数据分布 2018: 3.9%,,2019: 7.8%,,2020: 19.4%,,2021: 16.4%,, 2022: 3.8%

length(dat$N2O_N.m[dat$year%in% c(2018) & dat$season=='growing' & dat$N2O_N.m<0.1 & dat$treatment=='warmed'])/
  length(dat$N2O_N.m[dat$year%in% c(2018) & dat$season=='growing' & dat$treatment=='warmed'])

length(dat$N2O_N.m[dat$year%in% c(2018) & dat$season=='growing' & dat$N2O_N.m<0.1 & dat$treatment=='warmed'])
#增温数据分布 2018: 1.3%（1天）,,2019: 10%（18天）,,2020: 26.2%（48天）,,2021: 25.7%（47天）,, 2022: 3.3%（6天）

length(dat$N2O_N.m[dat$year%in% c(2022) & dat$season=='growing' & dat$N2O_N.m<0.1 & dat$treatment=='control'])/
  length(dat$N2O_N.m[dat$year%in% c(2022) & dat$season=='growing' & dat$treatment=='control'])

length(dat$N2O_N.m[dat$year%in% c(2022) & dat$season=='growing' & dat$N2O_N.m<0.1 & dat$treatment=='control'])
#增温数据分布 2018: 6.6%（5天）,,2019: 5.6%（10d）,,2020: 12.6%（23d）,,2021: 7.1%（13d）,, 2022: 4.4%（8d）

dat$N2O_min[dat$N2O_min<0.1]=NA  #小于0.1不显示

unique(dat$labs_warming)
dat$N2O_N.m[dat$date=='2018/11/19']=NA #not completed observation

n_breaks=NA
for(i in 0:3){
  ni_breaks=seq(0.1,1,0.1)*10^i
  n_breaks=c(ni_breaks,n_breaks)
}


############### <Main_Fig 3> ###################
###### NO+N2O(log) 2018-2023 
F.pt$warming<-
  ifelse(F.pt$date>="2018/9/9" & F.pt$date<="2018/11/15","Y",
         ifelse(F.pt$date>="2019/6/9" & F.pt$date<="2019/12/1","Y",
                ifelse(F.pt$date>="2020/6/17" & F.pt$date<="2020/8/20","Y",
                       ifelse(F.pt$date>="2020/9/11" & F.pt$date<="2020/12/7","Y",
                              ifelse(F.pt$date>="2021/3/29"& F.pt$date<="2021/12/7","Y",
                                     ifelse(F.pt$date>="2022/3/26"& F.pt$date<="2022/12/7","Y",
                                            ifelse(F.pt$date>="2023/3/20"& F.pt$date<="2023/12/8","Y",
                                                   ifelse(F.pt$date>="2024/3/20"& F.pt$date<="2024/12/8","Y",
                                                        "N"))))))))


add<-F.pt[length(F.pt$soil_temp.m)+1,]
add$treatment="precipitation"
datai<-rbind(F.pt,add)   # make precipitation legend in Fig_3_<1>

{
  st<-ggplot(datai,aes(date,soil_temp.m,col=treatment))+
    annotate("rect",xmin=as.Date("2018-09-09"),
             xmax=as.Date("2018-11-16"),
             ymin=-10,ymax=30,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2019-06-08"),
             xmax=as.Date("2019-12-01"),
             ymin=-10,ymax=30,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2020-06-16"),
             xmax=as.Date("2020-08-26"),
             ymin=-10,ymax=30,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2020-09-13"),
             xmax=as.Date("2020-12-7"),
             ymin=-10,ymax=30,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2021-03-29"),
             xmax=as.Date("2021-12-7"),
             ymin=-10,ymax=30,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2022-03-26"),
             xmax=as.Date("2022-12-10"),
             ymin=-10,ymax=30,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2023-3-20"),
             xmax=as.Date("2023-12-08"),
             ymin=-10,ymax=30,alpha=.07,fill="red")+
    geom_hline(yintercept = 2,lty=2,linewidth=0.1,col="red")+
    geom_hline(yintercept = 0,lty=2,linewidth=0.1,col="red")+
    geom_line(aes(y=35),linewidth=0.7,show.legend = T)+   ####fake line
    geom_line(linewidth=0.15,show.legend = F)+
    scale_x_date(limits = c(as.Date("2018/8/1"),as.Date("2023/12/31")),
                 date_labels = "%Y/%m",position = "bottom",
                 date_minor_breaks="months",breaks="3 months",expand = c(0,0))+
    scale_y_continuous(limits=c(-10,30),breaks=seq(-10,30,10),minor_breaks=seq(-10,30,5),expand=c(0,0))+
    scale_color_manual(values=c("blue","red","black","grey"),limits=c("control","warmed","diff","precipitation"),name="",
                       labels=c("Control","Warmed","Difference","Precipitation"))+
    mythem_sci+
    guides(y="prism_offset_minor",x="prism_offset_minor",
           color=guide_legend(ncol = 4,keyheight = 3,keywidth = 0.35,size = 30,byrow = FALSE))+   #多行图例设置
    labs(x="Date",y=expression(atop(paste('Soil temp (\u00B0C)'))))+
    theme(plot.margin = unit(c(0.5,4,0.3,4),'lines'),
          axis.title = element_text(colour = "black",size=21),
          legend.position = c(0.20,0.97))+coord_cartesian(clip="on");st
  
  str(F.pt)
  F.pt$month<-as.numeric(format(F.pt$date,'%m'))
  
  ##土壤湿度
  sm<-ggplot(F.pt[which(F.pt$treatment!="diff" & F.pt$month %in% seq(3,11)),],aes(date,WFPS.m,col=treatment,group=interaction(year,treatment)))+
    annotate("rect",xmin=as.Date("2018-09-09"),
             xmax=as.Date("2018-11-16"),
             ymin=0,ymax=80,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2019-06-08"),
             xmax=as.Date("2019-12-01"),
             ymin=0,ymax=80,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2020-06-16"),
             xmax=as.Date("2020-08-26"),
             ymin=0,ymax=80,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2020-09-13"),
             xmax=as.Date("2020-12-7"),
             ymin=0,ymax=80,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2021-03-29"),
             xmax=as.Date("2021-12-7"),
             ymin=0,ymax=80,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2022-03-26"),
             xmax=as.Date("2022-12-10"),
             ymin=0,ymax=80,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2023-3-20"),
             xmax=as.Date("2023-12-08"),
             ymin=0,ymax=80,alpha=.07,fill="red")+
    geom_bar(aes(y=precipitation*2/3),position=position_dodge(),stat="identity",size=0.1,fill="grey",col="transparent",alpha=.45)+
    geom_line(linewidth=0.15,show.legend = F,na.rm = T)+
    geom_errorbar(aes(ymin=WFPS.m-WFPS.se, ymax=WFPS.m+WFPS.se),stat = "identity",
                  linewidth=0.05,width=0,position=position_dodge(0),alpha=.3,show.legend = F)+
    scale_colour_manual(values=c("blue","red","grey"),name="",
                        labels=c("control","warmed","precipitation"),
                        limits=c("control","warmed","diff"))+
    scale_x_date(limits = c(as.Date("2018/8/1"),as.Date("2023/12/31")),
                 date_labels = "%Y/%m",position = "bottom",
                 date_minor_breaks="months",breaks="3 months",name="Date",expand = c(0,0))+
    scale_y_continuous(limits=c(0,80),breaks=seq(0,80,20),expand=c(0,0),
                       labels = c(seq(0,60, by=20),""),
                       name = expression(atop(paste("WFPS"),paste("(0-10 cm, %)"))),
                       sec.axis=sec_axis(~.*3/2,breaks=seq(0,120,30),labels=seq(0,120,30),
                                         name=expression(paste('Precipitation (mm)'))))+
    guides(y="prism_offset_minor",x="prism_offset_minor")+mythem_sci+
    guides(y="prism_offset_minor",x="prism_offset_minor")+mythem_sci+
    theme(plot.margin = unit(c(0.5,4,0.3,4),'lines'),
          legend.position = "none")+coord_cartesian(clip="on");sm
  
  # sm<-ggplot(F.pt[which(F.pt$treatment!="diff" & F.pt$month %in% seq(3,11)),],aes(date,WFPS.m,col=treatment,group=year))+
  #   annotate("rect",xmin=as.Date("2018-09-09"),
  #            xmax=as.Date("2018-11-16"),
  #            ymin=0,ymax=80,alpha=.07,fill="red")+
  #   annotate("rect",xmin=as.Date("2019-06-08"),
  #            xmax=as.Date("2019-12-01"),
  #            ymin=0,ymax=80,alpha=.07,fill="red")+
  #   annotate("rect",xmin=as.Date("2020-06-16"),
  #            xmax=as.Date("2020-08-26"),
  #            ymin=0,ymax=80,alpha=.07,fill="red")+
  #   annotate("rect",xmin=as.Date("2020-09-13"),
  #            xmax=as.Date("2020-12-7"),
  #            ymin=0,ymax=80,alpha=.07,fill="red")+
  #   annotate("rect",xmin=as.Date("2021-03-29"),
  #            xmax=as.Date("2021-12-7"),
  #            ymin=0,ymax=80,alpha=.07,fill="red")+
  #   annotate("rect",xmin=as.Date("2022-03-26"),
  #            xmax=as.Date("2022-12-10"),
  #            ymin=0,ymax=80,alpha=.07,fill="red")+
  #   annotate("rect",xmin=as.Date("2023-3-20"),
  #            xmax=as.Date("2023-12-08"),
  #            ymin=0,ymax=80,alpha=.07,fill="red")+
  #   geom_bar(aes(y=precipitation*2/3),position=position_dodge(),stat="identity",size=0.1,fill="grey",col="transparent",alpha=.55)+
  #   geom_errorbar(aes(ymin=WFPS.m-WFPS.se, ymax=WFPS.m+WFPS.se),stat = "identity",
  #                 linewidth=0.015,width=0,position=position_dodge(0),alpha=.3,show.legend = F)+
  #   geom_line(linewidth=0.10,show.legend = F)+
  #   scale_colour_manual(values=c("blue","red","grey"),name="Treatment",
  #                       labels=c("control","warmed","difference"),
  #                       limits=c("control","warmed","diff"))+
  #   scale_x_date(limits = c(as.Date("2018/8/1"),as.Date("2023/12/31")),
  #                date_labels = "%Y/%m",position = "bottom",
  #                date_minor_breaks="months",breaks="3 months",name="Date",expand = c(0,0))+
  #   scale_y_continuous(limits=c(0,80),breaks=seq(0,80,20),expand=c(0,0),
  #                      labels = c(seq(0,60, by=20),""),
  #                      name = expression(atop(paste("WFPS"),paste("(0-10 cm, %)"))),
  #                      sec.axis=sec_axis(~.*3/2,breaks=seq(0,120,30),labels=seq(0,120,30),
  #                                        name=expression(paste('Precipitation (mm)'))))+
  #   
  #   # geom_hline(yintercept = 0,lty=2,size=0.4,col="grey")+
  #   guides(y="prism_offset_minor",x="prism_offset_minor")+mythem_sci+
  #   theme(plot.margin = unit(c(0.5,4,0.3,4),'lines'),
  #         legend.position = "Top")+coord_cartesian(clip="on");sm
  
  d0<-ggplot(F.pt[which(F.pt$treatment=="diff"  & F.pt$month %in% seq(3,11)),], aes(date,WFPS.m,col=treatment))+
    annotate("rect",xmin=as.Date("2018-09-09"),
             xmax=as.Date("2018-11-16"),
             ymin=-20,ymax=28,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2019-06-08"),
             xmax=as.Date("2019-12-01"),
             ymin=-20,ymax=28,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2020-06-16"),
             xmax=as.Date("2020-08-26"),
             ymin=-20,ymax=28,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2020-09-13"),
             xmax=as.Date("2020-12-7"),
             ymin=-20,ymax=28,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2021-03-29"),
             xmax=as.Date("2021-12-7"),
             ymin=-20,ymax=28,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2022-03-26"),
             xmax=as.Date("2022-12-10"),
             ymin=-20,ymax=28,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2023-03-20"),
             xmax=as.Date("2023-12-08"),
             ymin=-20,ymax=28,alpha=.07,fill="red")+
    scale_x_date(limits = c(as.Date("2018/8/1"),as.Date("2023/12/31")),
                 date_labels = "%Y/%m",position = "bottom",
                 date_minor_breaks="months",breaks="3 months",name="Date",expand = c(0,0))+
    scale_y_continuous(limits=c(-20,28),breaks=seq(-20,20,20),expand=c(0,0),
                       labels = c(seq(-20,20, by=20)),
                       name = expression(atop(paste("Diff of WFPS"),paste("(0-10 cm, %)"))))+
    geom_hline(yintercept = 0,col='red',lty='dashed',size=0.1,show.legend=T,na.rm=TRUE)+
    geom_point(size=0.1,stroke=0.25,col='gray1',alpha=0.7)+
    guides(y="prism_offset_minor",x="prism_offset_minor")+mythem_sci+
    theme(plot.margin = unit(c(0.5,4,0.3,4),'lines'),
          legend.position = "Top")+coord_cartesian(clip="on");d0
  
  p1<-ggplot(F.pt[which(F.pt$date>="2018/7/8"& F.pt$treatment!="diff"),],
             aes(date,NO_N.m,colour=treatment))+
    annotate("rect",xmin=as.Date("2018-09-09"),
             xmax=as.Date("2018-11-16"),
             ymin=0,ymax=75,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2019-06-08"),
             xmax=as.Date("2019-12-01"),
             ymin=0,ymax=75,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2020-06-16"),
             xmax=as.Date("2020-8-26"),
             ymin=0,ymax=75,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2020-09-13"),
             xmax=as.Date("2020-12-7"),
             ymin=0,ymax=75,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2021-03-29"),
             xmax=as.Date("2021-12-7"),
             ymin=0,ymax=75,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2022-03-26"),
             xmax=as.Date("2022-12-10"),
             ymin=0,ymax=75,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2023-03-20"),
             xmax=as.Date("2023-12-08"),
             ymin=0,ymax=75,alpha=.07,fill="red")+
    geom_point(size=0.1,stroke=0.25)+
    geom_errorbar(data=F.pt[which(F.pt$treatment!="diff"),],alpha=0.3,
                  aes(x=date,ymin=NO_N.m-NO_N.se,ymax=NO_N.m+NO_N.se),width=0,linewidth=.05,show.legend=F)+
    # geom_line(linewidth=0.15,show.legend=T,na.rm=TRUE)+
    scale_colour_manual(values=c("blue","red"),name=NULL,
                        limits=c("control","warmed"),
                        labels=c("Control","Warmed"))+
    scale_fill_manual(values=c("blue","red"),name=NULL,
                      limits=c("control","warmed"),
                      labels=c("Control","Warmed"))+
    scale_x_date(limits = c(as.Date("2018/8/1"),as.Date("2023/12/31")),date_labels = "%Y/%m",position = "bottom",
                 date_minor_breaks="months",breaks="3 months",name="Date",expand = c(0,0))+
    scale_y_continuous(limits=c(0,75),breaks=seq(0,60,20),expand=c(0,0),labels = c(seq(0,60, by=20)),
                       name = expression(atop(paste("NO flux"),
                                              paste("(μg "," N"," m"^-2," h"^-1,")"))))+
    guides(y="prism_offset_minor",x="prism_offset_minor")+mythem_sci+
    theme(plot.margin = unit(c(0.5,4,0.3,4),'lines'),
          legend.position = "Top")+
    coord_cartesian(clip="on");p1
  
  
  
  d1<-ggplot(F.pt[which(F.pt$date>="2018/7/8"& F.pt$treatment=="diff"),],
             aes(date,NO_N.m,colour=treatment))+
    annotate("rect",xmin=as.Date("2018-09-09"),
             xmax=as.Date("2018-11-16"),
             ymin=-30,ymax=40,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2019-06-08"),
             xmax=as.Date("2019-12-01"),
             ymin=-30,ymax=40,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2020-06-16"),
             xmax=as.Date("2020-8-26"),
             ymin=-30,ymax=40,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2020-09-13"),
             xmax=as.Date("2020-12-7"),
             ymin=-30,ymax=40,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2021-03-29"),
             xmax=as.Date("2021-12-7"),
             ymin=-30,ymax=40,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2022-03-26"),
             xmax=as.Date("2022-12-10"),
             ymin=-30,ymax=40,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2023-03-20"),
             xmax=as.Date("2023-12-08"),
             ymin=-30,ymax=40,alpha=.07,fill="red")+
    geom_hline(yintercept = 0,lty=2,col='red',linewidth=0.10,show.legend=T,na.rm=TRUE)+
    geom_point(size=0.1,stroke=0.25,col='gray1',alpha=0.7)+
    scale_x_date(limits = c(as.Date("2018/8/1"),as.Date("2023/12/31")),date_labels = "%Y/%m",position = "bottom",
                 date_minor_breaks="months",breaks="3 months",name="Date",expand = c(0,0))+
    scale_y_continuous(breaks=seq(-30,30,30),expand=c(0,0),
                       name = expression(atop(paste("Diff of NO flux"),
                                              paste("(μg "," N"," m"^-2," h"^-1,")"))))+
    guides(y="prism_offset_minor",x="prism_offset_minor")+mythem_sci+
    theme(plot.margin = unit(c(0.5,4,0.3,4),'lines'),
          legend.position = "Top")+
    coord_cartesian(clip="on",ylim=c(-30,40));d1
  
  p_log<-ggplot(dat[dat$treatment!='diff',],aes(date,N2O_N.m+10,colour=treatment))+
    annotate("rect",xmin=as.Date("2018-09-09"),
             xmax=as.Date("2018-11-16"),
             ymin=1,ymax=2500,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2019-06-08"),
             xmax=as.Date("2019-12-01"),
             ymin=1,ymax=2500,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2020-06-16"),
             xmax=as.Date("2020-8-26"),
             ymin=1,ymax=2500,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2020-09-13"),
             xmax=as.Date("2020-12-7"),
             ymin=1,ymax=2500,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2021-03-29"),
             xmax=as.Date("2021-12-7"),
             ymin=1,ymax=2500,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2022-03-26"),
             xmax=as.Date("2022-12-10"),
             ymin=1,ymax=2500,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2023-03-20"),
             xmax=as.Date("2023-12-08"),
             ymin=1,ymax=2500,alpha=.07,fill="red")+
    geom_hline(yintercept = 10,col='red',lty='dashed',linewidth=0.1,show.legend=T,na.rm=TRUE)+
    geom_point(size=0.1,stroke=0.25)+
    # geom_line(linewidth=0.15,show.legend=T,na.rm=TRUE)+
    geom_errorbar(aes(x=date,ymin=N2O_min+10,ymax=N2O_max+10),
                  width=0,linewidth=.05,show.legend=F,alpha=0.3)+###size=0.4原来
    scale_colour_manual(values=c("blue","red","grey"),name="",
                        limits=c("control","warmed","diff"),
                        labels=c("Control","Warmed","Warmed-contol"))+
    scale_y_log10(breaks=c(1,10,100,1000),guide = "prism_offset_minor",
                  minor_breaks = na.omit(n_breaks),
                  labels = c(expression(1),
                             expression(paste(10)),expression(paste(100)),
                             expression(paste(1000))),
                  name = expression(atop(paste("N"[2],"O flux"),paste("(+10 μg "," N"," m"^-2," h"^-1,")"))),expand = c(0,0))+
    scale_x_date(limits = c(as.Date("2018/8/1"),as.Date("2023/12/31")),date_labels = "%Y/%m",position = "bottom",
                 date_minor_breaks="months",breaks="3 months",name="Date",expand = c(0,0))+
    guides(y="prism_offset_minor",x="prism_offset_minor")+mythem_sci+
    theme(plot.margin = unit(c(0,4,3,4),'lines'),
          axis.ticks.length.x.top=unit(0.15,"cm"),
          prism.ticks.length.x.top=unit(0.12,"cm"),
          axis.text.x.top = element_blank(),
          legend.position = "Top")+
    coord_cartesian(clip="off",ylim = c(1,2500));p_log
  
  
  d2<-ggplot(dat[dat$treatment=='diff' & dat$N2O_N.m<50 & dat$N2O_N.m>-50,],aes(date,N2O_N.m,colour=treatment))+
    annotate("rect",xmin=as.Date("2018-09-09"),
             xmax=as.Date("2018-11-16"),
             ymin=-50,ymax=70,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2019-06-08"),
             xmax=as.Date("2019-12-01"),
             ymin=-50,ymax=70,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2020-06-16"),
             xmax=as.Date("2020-8-26"),
             ymin=-50,ymax=70,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2020-09-13"),
             xmax=as.Date("2020-12-7"),
             ymin=-50,ymax=70,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2021-03-29"),
             xmax=as.Date("2021-12-7"),
             ymin=-50,ymax=70,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2022-03-26"),
             xmax=as.Date("2022-12-10"),
             ymin=-50,ymax=70,alpha=.07,fill="red")+
    annotate("rect",xmin=as.Date("2023-03-20"),
             xmax=as.Date("2023-12-08"),
             ymin=-50,ymax=70,alpha=.07,fill="red")+
    geom_hline(yintercept = 0,col='red',lty='dashed',size=0.1,show.legend=T,na.rm=TRUE)+
    geom_point(size=0.1,stroke=0.25,col='gray1',alpha=0.7)+
    scale_y_continuous(breaks=seq(-50,50,50),minor_breaks = seq(-50,50,25),
                       name = expression(atop(paste("Diff of N"[2],"O flux"),paste("(μg "," N"," m"^-2," h"^-1,")"))),expand = c(0,0))+
    scale_x_date(limits = c(as.Date("2018/8/1"),as.Date("2023/12/31")),name=NULL,expand = c(0,0),
                 date_labels = c('8/1','11/1',rep(c('2/1','5/1','8/1','11/1'),4),c('2/1','5/1','8/1','11/1')),
                 breaks="3 months",date_minor_breaks="months",position = "bottom")+
    guides(y="prism_offset_minor",x="prism_offset_minor")+mythem_sci+
    theme(plot.margin = unit(c(0,4,3,4),'lines'),
          axis.ticks.length.x.top=unit(0.07,"cm"),
          prism.ticks.length.x.top=unit(0.04,"cm"),
          axis.text.x.top = element_blank(),
          legend.position = "Top")+
    coord_cartesian(clip="off",ylim = c(-50,70))+
    ###加标注文本
    annotate(geom="text",x=as.Date("2018/10/1"), y=-100,label="2018",colour="black",vjust=0,size=2.3)+
    annotate(geom="text",x=as.Date("2019/7/1"), y=-100,label="2019",colour="black",vjust=0,size=2.3)+
    annotate(geom="text",x=as.Date("2020/7/1"), y=-100,label="2020",colour="black",vjust=0,size=2.3)+
    annotate(geom="text",x=as.Date("2021/7/1"), y=-100,label="2021",colour="black",vjust=0,size=2.3)+
    annotate(geom="text",x=as.Date("2022/7/1"), y=-100,label="2022",colour="black",vjust=0,size=2.3)+
    annotate(geom="text",x=as.Date("2023/7/1"), y=-100,label="2023",colour="black",vjust=0,size=2.3)+
    ###加箭头
    geom_segment(aes(x=as.Date("2023/12/9"), xend=as.Date("2023/1/15"),y=-81, yend=-81),
                 arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
    geom_segment(aes(x=as.Date("2022/12/15"), xend=as.Date("2022/1/15"),y=-81, yend=-81),
                 arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
    geom_segment(aes(x=as.Date("2021/12/15"), xend=as.Date("2021/1/15"),y=-81, yend=-81),
                 arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
    geom_segment(aes(x=as.Date("2020/12/15"), xend=as.Date("2020/1/15"),y=-81, yend=-81),
                 arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
    geom_segment(aes(x=as.Date("2019/12/15"), xend=as.Date("2019/1/15"),y=-81, yend=-81),
                 arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
    geom_segment(aes(x=as.Date("2018/12/15"), xend=as.Date("2018/8/1"),y=-81, yend=-81),
                 arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F);d2
  
  figure<-ggarrange(st+rremove("x.text")+rremove("xlab")+
                      guides(color=guide_legend(override.aes=list(size=4),keyheight = 10,keywidth = 0.35,ncol = 4,byrow = FALSE))+
                      theme(axis.ticks.length.x.top=unit(-0.07,"cm"),
                            prism.ticks.length.x.top=unit(-0.04,"cm"),
                            axis.ticks.length.y=unit(0.07,"cm"),
                            prism.ticks.length.y=unit(0.04,"cm"),
                            axis.text.y.right = element_blank(),
                            axis.text.x = element_blank(),
                            axis.title.x = element_blank(),
                            # legend.text=element_text(colour = "black",size=16),
                            legend.position = c(0.38,1.00),
                            plot.margin = unit(c(0.2,0.2,0,0.2),'cm'),
                            plot.background = element_blank(), #拼合时图片互无遮挡
                            panel.background = element_rect(fill =NA, colour = "black", linetype = "solid")),
                    sm+rremove("x.text")+rremove("xlab")+
                      theme(axis.ticks.length.x.top=unit(-0.07,"cm"),
                            prism.ticks.length.x.top=unit(-0.04,"cm"),
                            axis.ticks.length.y=unit(0.07,"cm"),
                            prism.ticks.length.y=unit(0.04,"cm"),
                            axis.title.y.right = element_text(vjust=1,colour = "black"),
                            axis.text.y.right = element_text(hjust=0.5,colour = "black"),
                            axis.text.x = element_blank(),
                            axis.title.x = element_blank(),
                            plot.margin = unit(c(-0.1,0.2,0,0.2),'cm'),
                            plot.background = element_blank(), #拼合时图片互无遮挡
                            panel.background = element_rect(fill =NA, colour = "black", linetype = "solid")),
                    d0+rremove("x.text")+rremove("xlab")+
                      theme(axis.ticks.length.x.top=unit(-0.07,"cm"),
                            prism.ticks.length.x.top=unit(-0.04,"cm"),
                            axis.ticks.length.y=unit(0.07,"cm"),
                            prism.ticks.length.y=unit(0.04,"cm"),
                            axis.title.y.right = element_text(vjust=1,colour = "black"),
                            axis.text.y.right = element_text(hjust=0.5,colour = "black"),
                            axis.text.x = element_blank(),
                            axis.title.x = element_blank(),
                            plot.margin = unit(c(-0.1,0.2,0,0.2),'cm'),
                            plot.background = element_blank(), #拼合时图片互无遮挡
                            panel.background = element_rect(fill =NA, colour = "black", linetype = "solid")),
                    p1+rremove("x.text")+rremove("xlab")+
                      theme(axis.ticks.length.x.top=unit(-0.07,"cm"),
                            prism.ticks.length.x.top=unit(-0.04,"cm"),
                            axis.ticks.length.y=unit(0.07,"cm"),
                            prism.ticks.length.y=unit(0.04,"cm"),
                            axis.title.y.right = element_text(vjust=1,colour = "black"),
                            axis.text.y.right = element_text(hjust=0.5,colour = "black"),
                            axis.text.x = element_blank(),
                            axis.title.x = element_blank(),
                            plot.margin = unit(c(-0.1,0.2,0,0.2),'cm'),
                            plot.background = element_blank(), #拼合时图片互无遮挡
                            panel.background = element_rect(fill =NA, colour = "black", linetype = "solid")),
                    d1+rremove("x.text")+rremove("xlab")+
                      theme(axis.ticks.length.x.top=unit(-0.07,"cm"),
                            prism.ticks.length.x.top=unit(-0.04,"cm"),
                            axis.ticks.length.y=unit(0.07,"cm"),
                            prism.ticks.length.y=unit(0.04,"cm"),
                            axis.title.y.right = element_text(vjust=1,colour = "black"),
                            axis.text.y.right = element_text(hjust=0.5,colour = "black"),
                            axis.text.x = element_blank(),
                            axis.title.x = element_blank(),
                            plot.margin = unit(c(-0.1,0.2,0,0.2),'cm'),
                            plot.background = element_blank(), #拼合时图片互无遮挡
                            panel.background = element_rect(fill =NA, colour = "black", linetype = "solid")),
                    p_log+rremove("x.text")+rremove("xlab")+
                      theme(axis.ticks.length.x.top=unit(-0.07,"cm"),
                            prism.ticks.length.x.top=unit(-0.04,"cm"),
                            axis.ticks.length.y=unit(0.07,"cm"),
                            prism.ticks.length.y=unit(0.04,"cm"),
                            axis.title.y.right = element_text(vjust=1,colour = "black"),
                            axis.text.y.right = element_text(hjust=0.5,colour = "black"),
                            axis.text.x = element_blank(),
                            axis.title.x = element_blank(),
                            plot.margin = unit(c(-0.1,0.2,0,0.2),'cm'),
                            plot.background = element_blank(), #拼合时图片互无遮挡
                            panel.background = element_rect(fill =NA, colour = "black", linetype = "solid")),
                    d2+theme(axis.ticks.length.x.bottom=unit(0.07,"cm"),
                             prism.ticks.length.x.bottom=unit(0.04,"cm"),
                             axis.text.x.top = element_blank(),
                             axis.ticks.length.y=unit(0.07,"cm"),
                             prism.ticks.length.y=unit(0.04,"cm"),
                             legend.position="none", plot.margin = unit(c(-0.1,0.2,0.4,0.2),'cm'),
                             plot.background = element_blank(),   #拼合时图片互无遮挡
                             panel.background = element_rect(fill =NA, colour = "black", linetype = "solid")), 
                    labels = c("A", "B","C","D","E","F","G"), 
                    label.x = 0.123,label.y = c(0.90,1.05,1.05,1.05,1.05,1.05,1.05),                   
                    font.label = list(size = 7, face = "bold", color ="black"),
                    ncol = 1, nrow = 7,align = "v",   ##"v"竖直对齐
                    widths = 4, heights = c(3,3,2.2,3.5,2.2,3.5,3.6),
                    common.legend=F);figure
  
  ######10:14
  myname=paste("Fig_3_st+sw+no+n2o_log+diffs时间序列18-23",format(Sys.Date( ),"%Y-%m-%d"),sep="_");print(myname)    #获取当前时间
  ggsave(paste(myname,"pdf",sep="."),figure,height = 14, width = 14,units='cm')
}
  
###描述文字部分
unique(F.pt$season)
unique(F.pt$year)

datFUN(F.pt$NO_N.m[F.pt$treatment=='control' & F.pt$NO_N.m>0 & F.pt$year<=2023])   ###control NO mean flux
datFUN(F.pt$N2O_N.m[F.pt$treatment=='control'  & F.pt$year<=2023])   ###control N2O mean flux

datFUN(na.omit(F.pt$NO_N.m[F.pt$treatment=='control' & F.pt$year %in% c(2018:2023)]))
datFUN(na.omit(F.pt$N2O_N.m[F.pt$treatment=='control' & F.pt$year %in% c(2018:2023)]))

###1.1增温期间 对照处理生长季的释放速率
datFUN(na.omit(F.pt$NO_N.m[F.pt$treatment=='control' & F.pt$season=='growing' &
                              F.pt$warming=='Y' & F.pt$year %in% c(2018:2023)]))

datFUN(na.omit(F.pt$N2O_N.m[F.pt$treatment=='control' & F.pt$season=='growing' &
                              F.pt$warming=='Y' & F.pt$year %in% c(2018:2023)]))

###1.2增温期间 增温处理生长季的释放速率
datFUN(na.omit(F.pt$N2O_N.m[F.pt$treatment=='warmed' & F.pt$season=='growing' &
                              F.pt$warming=='Y' & F.pt$year %in% c(2018:2023)]))

###2.1增温期间 增温处理非生长季的释放速率
datFUN(na.omit(F.pt$N2O_N.m[F.pt$treatment=='control' & F.pt$seaso!='growing' &
                              F.pt$warming=='Y']))
###2.2增温期间 增温处理非生长季的释放速率
datFUN(na.omit(F.pt$N2O_N.m[F.pt$treatment=='warmed' & F.pt$season!='growing' &
                              F.pt$warming=='Y']))

mean(na.omit(F.pt$N2_N2O.R.m[F.pt$treatment=='control']))/(mean(na.omit(F.pt$N2_N2O.R.m[F.pt$treatment=='control']))+1)
range(na.omit(F.pt$N2_N2O.R.m[F.pt$treatment=='control']/(F.pt$N2_N2O.R.m[F.pt$treatment=='control']+1)))


mean(na.omit(F.pt$N2_N2O.R.m[F.pt$treatment=='warmed']))/(mean(na.omit(F.pt$N2_N2O.R.m[F.pt$treatment=='warmed']))+1)



t.test(mydata_h$N2O_N[mydata_h$year==2020 & mydata_h$treatment=='control'],
       mydata_h$N2O_N[mydata_h$year==2020 & mydata_h$treatment=='warmed'])


rm(x,y,st,sm,p1,p2,p1_log,p1diff_log,p_log,figure,datai,d2,d1,d0,add)

############# Fig.S3 ################

ggplot(F.pt[which(F.pt$treatment!="diff"),],aes(date,WFPS.m,col=treatment))+
  annotate("rect",xmin=as.Date("2018-09-09"),
           xmax=as.Date("2018-11-16"),
           ymin=0,ymax=80,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2019-06-08"),
           xmax=as.Date("2019-12-01"),
           ymin=0,ymax=80,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2020-06-16"),
           xmax=as.Date("2020-08-26"),
           ymin=0,ymax=80,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2020-09-13"),
           xmax=as.Date("2020-12-7"),
           ymin=0,ymax=80,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2021-03-29"),
           xmax=as.Date("2021-12-7"),
           ymin=0,ymax=80,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2022-03-26"),
           xmax=as.Date("2022-12-10"),
           ymin=0,ymax=80,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2023-3-20"),
           xmax=as.Date("2023-12-08"),
           ymin=0,ymax=80,alpha=.07,fill="red")+
  geom_bar(aes(y=precipitation*2/3),position=position_dodge(),stat="identity",size=0.1,fill="grey",col="transparent",alpha=.25)+
  geom_line(linewidth=0.15,show.legend = F)+
  geom_line(data = F.pt[which(F.pt$treatment!="diff" & F.pt$date<='2018/6/1'),],size=1,show.legend = T)+  #just for showing legend
  geom_errorbar(aes(ymin=WFPS.m-WFPS.se, ymax=WFPS.m+WFPS.se),stat = "identity",
                linewidth=0.05,width=0,position=position_dodge(0),alpha=.3,show.legend = F)+
  scale_colour_manual(values=c("blue","red","grey"),name="",
                      labels=c("control","warmed","precipitation"),
                      limits=c("control","warmed","diff"))+
  scale_x_date(limits = c(as.Date("2018/8/1"),as.Date("2023/12/31")),date_labels = "%Y/%m",position = "top",
               date_minor_breaks="months",breaks="3 months",name="Date",expand = c(0,0))+
  scale_y_continuous(limits=c(0,80),breaks=seq(0,80,20),expand=c(0,0),
                     labels = c(seq(0,60, by=20),"80"),
                     name = expression(atop(paste("WFPS"),paste("(0-10 cm, %)"))),
                     sec.axis=sec_axis(~.*3/2,breaks=seq(0,120,30),labels=seq(0,120,30),
                                       name=expression(paste('Precipitation (mm)'))))+
  guides(y="prism_offset_minor",x="prism_offset_minor")+mythem_sci+
  theme(plot.margin = unit(c(0.5,4,0.3,4),'lines'),legend.position=c(0.12,0.95))+coord_cartesian(clip="on")->p0;p0


str(F.treat)

###加和和比值之间的关系
str(F.pt)

nr_breaks=NA
for(i in seq(-2,2,1)){
  nir_breaks=seq(0.1,1,0.1)*10^i
  nr_breaks=na.omit(c(nir_breaks,nr_breaks))
}

print(nr_breaks)



names(F.pt)
p_R<-ggplot(F.pt[which(F.pt$date>="2018/7/8"& F.pt$treatment!="diff"),],
            aes(date,NO_N2O.m,colour=treatment))+
  annotate("rect",xmin=as.Date("2018-09-09"),
           xmax=as.Date("2018-11-16"),
           ymin=0.001,ymax=300,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2019-06-08"),
           xmax=as.Date("2019-12-01"),
           ymin=0.001,ymax=300,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2020-06-16"),
           xmax=as.Date("2020-8-26"),
           ymin=0.001,ymax=300,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2020-09-13"),
           xmax=as.Date("2020-12-7"),
           ymin=0.001,ymax=300,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2021-03-29"),
           xmax=as.Date("2021-12-7"),
           ymin=0.001,ymax=300,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2022-03-26"),
           xmax=as.Date("2022-12-10"),
           ymin=0.001,ymax=300,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2023-03-20"),
           xmax=as.Date("2023-12-08"),
           ymin=0.001,ymax=300,alpha=.07,fill="red")+
  geom_hline(yintercept = 1,lty=2, col='red', size=0.1,show.legend=T,na.rm=TRUE)+
  geom_point(size=0.2,stroke=0.25)+#geom_line(alpha=0.5,linewidth=0.1)+
  geom_errorbar(data=F.pt[which(F.pt$treatment!="diff"),],alpha=0.3,
                aes(x=date,ymin=NO_N2O.m,ymax=NO_N2O.m+NO_N2O.se),width=0,linewidth=.2,show.legend=F)+
  scale_colour_manual(values=c("blue","red"),name=NULL,
                      limits=c("control","warmed"),
                      labels=c("Control","Warmed"))+
  scale_fill_manual(values=c("blue","red"),name=NULL,
                    limits=c("control","warmed"),
                    labels=c("Control","Warmed"))+
  scale_x_date(limits = c(as.Date("2018/8/1"),as.Date("2023/12/31")),name=NULL,expand = c(0,0),
               date_labels = c('8/1','11/1',rep(c('2/1','5/1','8/1','11/1'),4),c('2/1','5/1','8/1')),breaks="3 months",
               date_minor_breaks="months",position = 'top')+
  scale_y_log10(breaks=c(0.001,0.01,0.1,1,10,100),guide = "prism_offset_minor",
                minor_breaks = na.omit(nr_breaks),
                labels = c(expression(0.001),expression(0.01),expression(paste(0.1)),expression(1),
                           expression(paste(10)),expression(paste(100))),expand=c(0,0),
                name = expression(paste('NO/N'[2],'O ratio')))+
  guides(y="prism_offset_minor",x="prism_offset_minor")+mythem_sci+
  coord_cartesian(clip="on",ylim = c(0.001,300));p_R

ggplot(F.treat[F.treat$treatment!='diff' & F.treat$soil_temp.m>=5,],aes(date,N2_N2O.R.m,col=treatment))+
  annotate("rect",xmin=as.Date("2018-09-09"),
           xmax=as.Date("2018-11-16"),
           ymin=0,ymax=11.5,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2019-06-08"),
           xmax=as.Date("2019-12-01"),
           ymin=0,ymax=11.5,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2020-06-16"),
           xmax=as.Date("2020-8-26"),
           ymin=0,ymax=11.5,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2020-09-13"),
           xmax=as.Date("2020-12-7"),
           ymin=0,ymax=11.5,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2021-03-29"),
           xmax=as.Date("2021-12-7"),
           ymin=0,ymax=11.5,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2022-03-26"),
           xmax=as.Date("2022-12-10"),
           ymin=0,ymax=11.5,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2023-03-20"),
           xmax=as.Date("2023-12-08"),
           ymin=0,ymax=11.5,alpha=.07,fill="red")+
  geom_point(size=0.2,stroke=0.25)+
  geom_errorbar(aes(x=date,ymin=N2_N2O.R.m-N2_N2O.R.se,ymax=N2_N2O.R.m+N2_N2O.R.se),
                width=0,linewidth=.05,show.legend=T,alpha=0.13)+###size=0.4原来
  labs(y=expression(paste('N'[2],'/N'[2],'O ratio')))+
  scale_colour_manual(values=c("blue","red","black"),name="",
                      limits=c("control","warmed","diff"),
                      labels=c("control","warmed","difference"))+
  scale_y_continuous(breaks=seq(0,10,2.5),guide = "prism_offset_minor",
                     minor_breaks = seq(0,10,1.25),expand=c(0,0))+
  scale_x_date(limits = c(as.Date("2018/8/1"),as.Date("2023/12/31")),position = "bottom",
               date_labels = c('8/1','11/1',rep(c('2/1','5/1','8/1','11/1'),4),c('2/1','5/1','8/1','11/1')),
               date_minor_breaks="months",breaks="3 months",name="Date",expand = c(0,0),
               sec.axis = dup_axis(trans = ~.,
                                   breaks=seq(as.Date("2018/8/1"),as.Date(max(F.pt$date+1)),'month'),name = NULL))+
  guides(y="prism_offset_minor",x="prism_offset_minor")+mythem_sci+
  theme(plot.margin = unit(c(2,4,0,4),'lines'),
        legend.position = "Top")+
  coord_cartesian(clip="off",ylim = c(0,11.5))+
  ###加标注文本
  annotate(geom="text",x=as.Date("2018/10/1"), y=-2.7,label="2018",colour="black",vjust=0,size=2.3)+
  annotate(geom="text",x=as.Date("2019/7/1"), y=-2.7,label="2019",colour="black",vjust=0,size=2.3)+
  annotate(geom="text",x=as.Date("2020/7/1"), y=-2.7,label="2020",colour="black",vjust=0,size=2.3)+
  annotate(geom="text",x=as.Date("2021/7/1"), y=-2.7,label="2021",colour="black",vjust=0,size=2.3)+
  annotate(geom="text",x=as.Date("2022/7/1"), y=-2.7,label="2022",colour="black",vjust=0,size=2.3)+
  annotate(geom="text",x=as.Date("2023/7/1"), y=-2.7,label="2023",colour="black",vjust=0,size=2.3)+
  ###加箭头
  geom_segment(aes(x=as.Date("2023/12/1"), xend=as.Date("2023/1/15"),y=-1.7, yend=-1.7),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(aes(x=as.Date("2022/12/15"), xend=as.Date("2022/1/15"),y=-1.7, yend=-1.7),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(aes(x=as.Date("2021/12/15"), xend=as.Date("2021/1/15"),y=-1.7, yend=-1.7),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(aes(x=as.Date("2020/12/15"), xend=as.Date("2020/1/15"),y=-1.7, yend=-1.7),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(aes(x=as.Date("2019/12/15"), xend=as.Date("2019/1/15"),y=-1.7, yend=-1.7),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(aes(x=as.Date("2018/12/15"), xend=as.Date("2018/8/1"),y=-1.7, yend=-1.7),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)->p2;p2   ###N2:N2O


ggplot(F.treat[F.treat$treatment!='diff' & F.treat$soil_temp.m>=5,],
       aes(date,N2_N2N2O.R.m,col=treatment))+
  annotate("rect",xmin=as.Date("2018-09-09"),
           xmax=as.Date("2018-11-16"),
           ymin=-Inf,ymax=Inf,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2019-06-08"),
           xmax=as.Date("2019-12-01"),
           ymin=-Inf,ymax=Inf,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2020-06-16"),
           xmax=as.Date("2020-8-26"),
           ymin=-Inf,ymax=Inf,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2020-09-13"),
           xmax=as.Date("2020-12-7"),
           ymin=-Inf,ymax=Inf,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2021-03-29"),
           xmax=as.Date("2021-12-7"),
           ymin=-Inf,ymax=Inf,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2022-03-26"),
           xmax=as.Date("2022-12-10"),
           ymin=-Inf,ymax=Inf,alpha=.07,fill="red")+
  annotate("rect",xmin=as.Date("2023-03-20"),
           xmax=as.Date("2023-12-08"),
           ymin=-Inf,ymax=Inf,alpha=.07,fill="red")+
  geom_point(size=0.2,stroke=0.25)+
  geom_errorbar(aes(x=date,ymin=ifelse(N2_N2N2O.R.m-N2_N2N2O.R.se<0, 0, N2_N2N2O.R.m-N2_N2N2O.R.se),
                    ymax=ifelse(N2_N2N2O.R.m+N2_N2N2O.R.se>3.35,3.35,N2_N2N2O.R.m+N2_N2N2O.R.se)),
                width=0,linewidth=.05,show.legend=T,alpha=0.13)+###size=0.4原来
  labs(y=expression(paste('N'[2],'/(N'[2],'+N'[2],'O) ratio')))+
  scale_colour_manual(values=c("blue","red","black"),name="",
                      limits=c("control","warmed","diff"),
                      labels=c("control","warmed","difference"))+
  scale_y_continuous(breaks=seq(0.4,1.2,0.2),guide = "prism_offset_minor",
                     minor_breaks = seq(0.4,1.2,0.1),expand=c(0,0))+
  scale_x_date(limits = c(as.Date("2018/8/1"),as.Date("2023/12/31")),position = "bottom",
               date_labels = c('8/1','11/1',rep(c('2/1','5/1','8/1','11/1'),4),c('2/1','5/1','8/1','11/1')),
               date_minor_breaks="months",breaks="3 months",name="Date",expand = c(0,0),
               sec.axis = dup_axis(trans = ~.,
                                   breaks=seq(as.Date("2018/8/1"),as.Date(max(F.pt$date+1)),'month'),name = NULL))+
  guides(y="prism_offset_minor",x="prism_offset_minor")+mythem_sci+
  theme(plot.margin = unit(c(2,4,0,4),'lines'),
        legend.position = "Top")+
  coord_cartesian(clip="off",ylim = c(0.4,1.1))+
  annotate(geom="text",x=as.Date("2018/10/1"), y=0.23,label="2018",colour="black",vjust=0,size=2.3)+
  annotate(geom="text",x=as.Date("2019/7/1"), y=0.23,label="2019",colour="black",vjust=-0,size=2.3)+
  annotate(geom="text",x=as.Date("2020/7/1"), y=0.23,label="2020",colour="black",vjust=-0,size=2.3)+
  annotate(geom="text",x=as.Date("2021/7/1"), y=0.23,label="2021",colour="black",vjust=-0,size=2.3)+
  annotate(geom="text",x=as.Date("2022/7/1"), y=0.23,label="2022",colour="black",vjust=-0,size=2.3)+
  annotate(geom="text",x=as.Date("2023/7/1"), y=0.23,label="2023",colour="black",vjust=-0,size=2.3)+
  ###加箭头
  geom_segment(data=F.treat[F.treat$date=='2022/6/2'&F.treat$treatment=='control',],
               aes(x=as.Date("2023/12/1"), xend=as.Date("2023/1/15"),y=.3, yend=.3),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(data=F.treat[F.treat$date=='2022/6/2'&F.treat$treatment=='control',],
               aes(x=as.Date("2022/12/15"), xend=as.Date("2022/1/15"),y=.3, yend=.3),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(data=F.treat[F.treat$date=='2022/6/2'&F.treat$treatment=='control',],
               aes(x=as.Date("2021/12/15"), xend=as.Date("2021/1/15"),y=.3, yend=.3),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(data=F.treat[F.treat$date=='2022/6/2'&F.treat$treatment=='control',],
               aes(x=as.Date("2020/12/15"), xend=as.Date("2020/1/15"),y=.3, yend=.3),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(data=F.treat[F.treat$date=='2022/6/2'&F.treat$treatment=='control',],
               aes(x=as.Date("2019/12/15"), xend=as.Date("2019/1/15"),y=.3, yend=.3),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)+
  geom_segment(data=F.treat[F.treat$date=='2022/6/2'&F.treat$treatment=='control',],
               aes(x=as.Date("2018/12/15"), xend=as.Date("2018/8/1"),y=.3, yend=.3),
               arrow=arrow(length = unit(0, "npc")),colour="black",size=.15,show.legend = F)->p3;p3 

figure<-ggarrange(p0+rremove("x.text")+rremove("xlab")+
                    guides(color=guide_legend(override.aes=list(size=3),keyheight = 5,keywidth = 0.15,ncol = 4,byrow = FALSE))+
                    theme(axis.ticks.length.x.top=unit(-0.07,"cm"),
                          prism.ticks.length.x.top=unit(-0.04,"cm"),
                          axis.ticks.length.y=unit(0.07,"cm"),
                          prism.ticks.length.y=unit(0.04,"cm"),
                          axis.title.y.right = element_text(vjust=1,colour = "black"),
                          axis.text.y.right = element_text(hjust=0.5,colour = "black"),
                          axis.text.x = element_blank(),
                          axis.title.x = element_blank(),
                          legend.position = c(0.28,0.98),
                          plot.margin = unit(c(0.2,0.2,0,0.2),'cm'),
                          plot.background = element_blank(), #拼合时图片互无遮挡
                          panel.background = element_rect(fill =NA, colour = "black", linetype = "solid")),
                  p_R+rremove("x.text")+rremove("xlab")+
                    theme(axis.ticks.length.x.top=unit(-0.07,"cm"),
                          prism.ticks.length.x.top=unit(-0.04,"cm"),
                          axis.ticks.length.y=unit(0.07,"cm"),
                          prism.ticks.length.y=unit(0.04,"cm"),
                          axis.text.x = element_blank(),
                          axis.title.x = element_blank(),
                          plot.margin = unit(c(0.07,0.2,0,0.2),'cm'),
                          plot.background = element_blank(), #拼合时图片互无遮挡
                          panel.background = element_rect(fill =NA, colour = "black", linetype = "solid")),
                  p2+theme(axis.ticks.length.x.top=unit(-0.07,"cm"),axis.ticks.length.x.bottom=unit(0.07,"cm"),
                             prism.ticks.length.x.top=unit(-0.07,"cm"),prism.ticks.length.x.bottom=unit(0.04,"cm"),
                             axis.text.x.top = element_blank(),
                             axis.title.x = element_blank(),
                             axis.ticks.length.y=unit(0.07,"cm"),
                             prism.ticks.length.y=unit(0.04,"cm"),
                             legend.position="none", plot.margin = unit(c(0,0.2,0.5,0.2),'cm'),
                             plot.background = element_blank(),   #拼合时图片互无遮挡
                             panel.background = element_rect(fill =NA, colour = "black", linetype = "solid")), 
                  labels = c("A", "B","C"), 
                  label.x = 0.13,label.y = c(0.92,0.98,0.98),                   
                  font.label = list(size = 7, face = "bold", color ="black"),
                  ncol = 1, nrow = 3,align = "v",   ##"v"竖直对齐
                  widths = 4, heights = c(3.2,3.5,4.6),
                  common.legend=F);figure

######10:14
myname=paste("sw+no_n2o_ratio+n2_n2o_ratio-23",format(Sys.Date( ),"%Y-%m-%d"),sep="_");print(myname)    #获取当前时间
ggsave(paste(myname,"pdf",sep="."),figure,height = 9, width = 14,units='cm')


figure<-ggarrange(p0+rremove("x.text")+rremove("xlab")+
                    guides(color=guide_legend(override.aes=list(size=3),keyheight = 5,keywidth = 0.35,ncol = 4,byrow = FALSE))+
                    theme(axis.ticks.length.x.top=unit(-0.07,"cm"),
                          prism.ticks.length.x.top=unit(-0.04,"cm"),
                          axis.ticks.length.y=unit(0.07,"cm"),
                          prism.ticks.length.y=unit(0.04,"cm"),
                          axis.title.y.right = element_text(vjust=1,colour = "black"),
                          axis.text.y.right = element_text(hjust=0.5,colour = "black"),
                          axis.text.x = element_blank(),
                          axis.title.x = element_blank(),
                          legend.position = c(0.28,0.98),
                          plot.margin = unit(c(0.2,0.2,0,0.2),'cm'),
                          plot.background = element_blank(), #拼合时图片互无遮挡
                          panel.background = element_rect(fill =NA, colour = "black", linetype = "solid")),
                  p_R+rremove("x.text")+rremove("xlab")+
                    theme(axis.ticks.length.x.top=unit(-0.07,"cm"),
                          prism.ticks.length.x.top=unit(-0.04,"cm"),
                          axis.ticks.length.y=unit(0.07,"cm"),
                          prism.ticks.length.y=unit(0.04,"cm"),
                          axis.title.y.right = element_text(vjust=1,colour = "black"),
                          axis.text.y.right = element_text(hjust=0.5,colour = "black"),
                          axis.text.x = element_blank(),
                          axis.title.x = element_blank(),
                          plot.margin = unit(c(0.07,0.2,0,0.2),'cm'),
                          plot.background = element_blank(), #拼合时图片互无遮挡
                          panel.background = element_rect(fill =NA, colour = "black", linetype = "solid")),
                  p3+theme(axis.ticks.length.x.top=unit(-0.07,"cm"),axis.ticks.length.x.bottom=unit(0.07,"cm"),
                           prism.ticks.length.x.top=unit(-0.07,"cm"),prism.ticks.length.x.bottom=unit(0.04,"cm"),
                           axis.text.x.top = element_blank(),
                           axis.title.x = element_blank(),
                           axis.ticks.length.y=unit(0.07,"cm"),
                           prism.ticks.length.y=unit(0.04,"cm"),
                           legend.position="none", plot.margin = unit(c(0,0.2,0.5,0.2),'cm'),
                           plot.background = element_blank(),   #拼合时图片互无遮挡
                           panel.background = element_rect(fill =NA, colour = "black", linetype = "solid")), 
                  labels = c("a", "b","c"), 
                  label.x = 0.13,label.y = c(0.92,0.98,0.98),                   
                  font.label = list(size = 7, face = "bold", color ="black"),
                  ncol = 1, nrow = 3,align = "none",   ##"v"竖直对齐
                  widths = 4, heights = c(3.2,3.5,4.6),
                  common.legend=F);figure

######10:14
myname=paste("sw+no_n2o_ratio+n2_n2n2o_ratio-23",format(Sys.Date( ),"%Y-%m-%d"),sep="_");print(myname)    #获取当前时间
ggsave(paste(myname,"pdf",sep="."),figure,height = 9, width = 14,units='cm')


##### the molar ratio of NO to N2O emissions was greater than 1:1 for 60% of the observations
length(F.treat$NO_N2O.m[F.treat$treatment=='control' & F.treat$NO_N2O.m>=1 & F.treat$year>=2019 & F.treat$year<=2023 ])/
  length(F.treat$NO_N2O.m[F.treat$treatment=='control' & F.treat$year>=2019 & F.treat$year<=2023])




###N2O/(n2o+n2)
datFUN(F.treat$N2_N2N2O.R.m[F.treat$treatment=='control' & F.treat$N2_N2N2O.R.m>0 & F.treat$season=='growing'])

datFUN(F.treat$N2_N2N2O.R.m[F.treat$treatment=='warmed' & F.treat$N2_N2N2O.R.m>0 & F.treat$season=='growing'])

t.test(F.treat$N2_N2O.R.m[F.treat$treatment=='control' & F.treat$N2_N2O.R.m>0 & F.treat$season=='growing' & F.treat$warming=='Y'],
       F.treat$N2_N2O.R.m[F.treat$treatment=='warmed' & F.treat$N2_N2O.R.m>0 & F.treat$season=='growing'& F.treat$warming=='Y'])

str(mydata_d)

################## <Main_Fig 4> ####################
###chamber level, cumulative flux  from 2019 to 2023, warming in growing seasons.

mydata_c$warming<-ifelse(mydata_c$date %in% mydata_d$date[mydata_d$warming=='Y'],'Y','N')

str(mydata_c)

cum_c<-data.frame()  ### calculate the cumulative emissions

for(i in c(paste('E',seq(1:15),sep=''),paste('W',seq(1:15),sep=''))){
  datai=mydata_c[mydata_c$system_id==i & mydata_c$date>='2019/1/1' & mydata_c$date<='2024/1/1',
                 c('date','warming','system_id','NO_N.m','N2O_N.m','N2_N.m','soil_temp.m')]
  datai=datai[!is.na(datai$date),] 
  datai<-datai[order(datai$date),]
  datai$NO_N.m[is.na(datai$NO_N.m)|datai$warming=='N']=0  #NA赋值给0
  datai$NO_N.cum=cumsum(datai$NO_N.m)*24*10^4*10^-9  #kg N ha-1
  
  datai$N2O_N.m[is.na(datai$N2O_N.m)|datai$warming=='N']=0  #NA赋值给0
  datai$N2O_N.cum=cumsum(datai$N2O_N.m)*24*10^4*10^-9  #kg N ha-1
  
  datai$N2_N.m[is.na(datai$N2_N.m)|datai$warming=='N'|datai$soil_temp.m<5]=0  #NA赋值给0
  datai$N2_N.cum=cumsum(datai$N2_N.m)*24*10^4*10^-9  #kg N ha-1
  
  cum_c<-rbind(cum_c,datai)
}

### define the plot for chambers
cum_c$plot<-
  ifelse(cum_c$system_id %in% c('E1','E4','E7','E10','E13'),1,
         ifelse(cum_c$system_id %in% c('E2','E5','E8','E11','E14'),2,
                ifelse(cum_c$system_id %in% c('E3','E6','E9','E12','E15'),3,
                       ifelse(cum_c$system_id %in% c('W1','W4','W7','W10','W13'),4,
                              ifelse(cum_c$system_id %in% c('W2','W5','W8','W11','W14'),5,
                                     ifelse(cum_c$system_id %in% c('W3','W6','W9','W12','W15'),6,''))))))

cum_c$treatment<-ifelse(cum_c$plot %in% c(1,4,5),'warmed','control')

ggplot(cum_c,aes(date,NO_N.cum,col=treatment))+
  geom_point(size=0.3)+
  scale_color_manual(limits=c('control','warmed'),values = c('blue','red'))

ggplot(cum_c,aes(date,N2O_N.cum,col=treatment))+
  geom_point(size=0.3)+
  scale_color_manual(limits=c('control','warmed'),values = c('blue','red'))

###计算在+2度的情况下，气态氮释放

cum_c1<-cum_c[which(cum_c$date == max(cum_c$date[cum_c$warming=='Y'])),]
names(cum_c1)

k=length(cum_c1$date)+1
for(i in c(paste('E',seq(1:15),sep=''),paste('W',seq(1:15),sep=''))){
  cum_c1[k,1]=cum_c$date[cum_c$system_id==i & cum_c$date == max(cum_c$date[cum_c$warming=='Y'])]
  cum_c1[k,2]=cum_c$warming[cum_c$system_id==i & cum_c$date == max(cum_c$date[cum_c$warming=='Y'])]
  cum_c1[k,3]=cum_c$system_id[cum_c$system_id==i & cum_c$date == max(cum_c$date[cum_c$warming=='Y'])]
  cum_c1[k,4]=cum_c$NO_N.m[cum_c$system_id==i & cum_c$date == max(cum_c$date[cum_c$warming=='Y'])]
  cum_c1[k,5]=cum_c$N2O_N.m[cum_c$system_id==i & cum_c$date == max(cum_c$date[cum_c$warming=='Y'])]
  cum_c1[k,6]=cum_c$N2O_N.m[cum_c$system_id==i & cum_c$date == max(cum_c$date[cum_c$warming=='Y'])]
  cum_c1[k,7]=cum_c$soil_temp.m[cum_c$system_id==i & cum_c$date == max(cum_c$date[cum_c$warming=='Y'])]+2
  cum_c1[k,8]=cum_c$NO_N.cum[cum_c$system_id==i & cum_c$date == max(cum_c$date[cum_c$warming=='Y'])]*(3.49)^(1/5)
  cum_c1[k,9]=cum_c$N2O_N.cum[cum_c$system_id==i & cum_c$date == max(cum_c$date[cum_c$warming=='Y'])]*(1.9)^(1/5)
  cum_c1[k,10]=cum_c$N2_N.cum[cum_c$system_id==i & cum_c$date == max(cum_c$date[cum_c$warming=='Y'])]*(2.4)^(1/5)
  cum_c1[k,11]=cum_c$plot[cum_c$system_id==i & cum_c$date == max(cum_c$date[cum_c$warming=='Y'])]
  cum_c1[k,12]='Q10_predicted'
  
  k=k+1
}

cum_c1$NON2O.cum=cum_c1$NO_N.cum+cum_c1$N2O_N.cum
cum_c1$gasN=cum_c1$NO_N.cum+cum_c1$N2O_N.cum+cum_c1$N2_N.cum

cum_c1<-cum_c1[-which(cum_c1$treatment=='Q10_predicted' & cum_c1$plot %in% c(1,4,5)),]

(mean(cum_c1$NO_N.cum[cum_c1$treatment=='warmed'])-mean(cum_c1$NO_N.cum[cum_c1$treatment=='control']))/
  mean(cum_c1$NO_N.cum[cum_c1$treatment=='control'])    # %cumulative decrease NO---19%

(mean(cum_c1$N2O_N.cum[cum_c1$treatment=='warmed'])-mean(cum_c1$N2O_N.cum[cum_c1$treatment=='control']))/
  mean(cum_c1$N2O_N.cum[cum_c1$treatment=='control'])  # %cumulative decrease N2O---15%

(mean(cum_c1$N2_N.cum[cum_c1$treatment=='warmed'])-mean(cum_c1$N2_N.cum[cum_c1$treatment=='control']))/
  mean(cum_c1$N2_N.cum[cum_c1$treatment=='control'])   #cumulative no change N2---


(mean(cum_c1$gasN[cum_c1$treatment=='warmed'])-mean(cum_c1$gasN[cum_c1$treatment=='control']))/
  mean(cum_c1$gasN[cum_c1$treatment=='control'])  # %cumulative decrease gaseous N---9%


(mean(cum_c1$NON2O.cum[cum_c1$treatment=='warmed'])-mean(cum_c1$NON2O.cum[cum_c1$treatment=='control']))/
  mean(cum_c1$NON2O.cum[cum_c1$treatment=='control'])   # %decrease NO+N2O---17%

my_comparisons <- list(c("control", "Q10_predicted"), c("Q10_predicted", "warmed"))

c1<-ggplot(cum_c1,aes(treatment,NO_N.cum,fill=treatment))+geom_point(pch=1)+
  geom_boxplot(aes(group =treatment),width=0.3,outlier.shape = 1,alpha=0.65)+ #geom_violin(fill=NA,width=0.25)+
  stat_mean(col='black',pch=16,size=2)+
  stat_compare_means(comparisons=my_comparisons,label.y = c(6, 6.3, 1.5),
                     method="t.test",
                     label="p.signif",tip.length = 0.01,
                     size=7)+
  stat_compare_means(comparisons=list(c('control','warmed')),
                     method="t.test",label.y = 1.2,label.y.npc = "bottom",
                     label="p.signif",tip.length = -0.01,vjust=1.5,
                     size=7)+
  scale_y_continuous(limits=c(0,8),breaks=seq(0,8,2),expand = c(0,0))+
  scale_fill_manual(limits=c('control','Q10_predicted','warmed'),values = c('blue','goldenrod','red'))+
  labs(x=NULL,y=expression(atop(paste('Cumulative NO'),paste('(kg N ha'^-1,')'))))+
  scale_x_discrete(limits=c('control','Q10_predicted','warmed'),labels=c('C','Q10','W'))+
  mythem;c1

c2<-ggplot(cum_c1,aes(treatment,N2O_N.cum,fill=treatment))+geom_point(pch=1)+
  geom_boxplot(aes(group =treatment),width=0.3,outlier.shape = 1,alpha=0.65)+ #geom_violin(fill=NA,width=0.25)+
  stat_mean(data=cum_c1,col='black',pch=16,size=2)+
  stat_compare_means(data=cum_c1[-which(cum_c1$treatment=='warmed'&cum_c1$N2O_N.cum>6),],
                     comparisons=my_comparisons,label.y = c(6.5, 7, 1.5),
                     method="t.test",
                     label="p.signif",tip.length = 0.01,
                     size=7)+
  stat_compare_means(data=cum_c1[-which(cum_c1$treatment=='warmed'&cum_c1$N2O_N.cum>6),],
                     comparisons=list(c('control','warmed')),
                     method="t.test",label.y = 1.2,label.y.npc = "bottom",
                     label="p.signif",tip.length = -0.01,vjust=1.5,
                     size=7)+
  scale_y_continuous(limits=c(0,8),breaks=seq(0,8,2),expand = c(0,0))+
  scale_fill_manual(limits=c('control','Q10_predicted','warmed'),values = c('blue','goldenrod','red'))+
  labs(x=NULL,y=expression(atop(paste('Cumulative N'[2],'O'),paste('(kg N ha'^-1,')'))))+
  scale_x_discrete(limits=c('control','Q10_predicted','warmed'),labels=c('C','Q10','W'))+
  mythem;c2

cum_c1$levels<-'Chamber level'
c3<-ggplot(cum_c1,aes(treatment,N2_N.cum,fill=treatment))+geom_point(pch=1)+facet_grid(levels~.)+
  geom_boxplot(aes(group =treatment),width=0.3,outlier.shape = 1,alpha=0.65)+ #geom_violin(fill=NA,width=0.25)+
  stat_mean(data=cum_c1[-which(cum_c1$treatment=='warmed'&cum_c1$N2_N.cum>15),],col='black',pch=16,size=2)+
  stat_compare_means(data=cum_c1[-which(cum_c1$treatment=='warmed'&cum_c1$N2_N.cum>12),],
                     comparisons=my_comparisons,label.y = c(12.5, 14, 1.5),
                     method="t.test",
                     label="p.signif",tip.length = 0.01,
                     size=7)+
  stat_compare_means(data=cum_c1[-which(cum_c1$treatment=='warmed'&cum_c1$N2_N.cum>12),],
                     comparisons=list(c('control','warmed')),
                     method="t.test",label.y = 0.8,label.y.npc = "bottom",
                     label="p.signif",tip.length = -0.01,vjust=1.5,
                     size=7)+
  scale_y_continuous(limits=c(0,18),breaks=seq(0,18,6),expand = c(0,0))+
  scale_fill_manual(limits=c('control','Q10_predicted','warmed'),values = c('blue','goldenrod','red'))+
  labs(x=NULL,y=expression(atop(paste('Cumulative N'[2]),paste('(kg N ha'^-1,')'))))+
  scale_x_discrete(limits=c('control','Q10_predicted','warmed'),labels=c('C','P','W'))+
  mythem;c3


###plot level, 不去除异常值
str(mydata_d)

cum_p<-data.frame()
for(i in unique(mydata_d$plot)){
  datai=mydata_d[mydata_d$plot==i & mydata_d$date>='2019/1/1',c('date','plot','NO_N','N2O_N','N2_N','warming','soil_temp')]
  datai=datai[!is.na(datai$date),] 
  datai<-datai[order(datai$date),]
  datai$NO_N[is.na(datai$NO_N)|datai$warming=='N']=0  #NA赋值给0
  datai$NO_N.cum=cumsum(datai$NO_N)*24*10^4*10^-9  #kg N ha-1
  
  datai$N2O_N[is.na(datai$N2O_N)|datai$warming=='N']=0  #NA赋值给0
  datai$N2O_N.cum=cumsum(datai$N2O_N)*24*10^4*10^-9  #kg N ha-1
  
  datai$NON2O.cum=datai$NO_N.cum+datai$N2O_N.cum
  
  datai$N2_N[is.na(datai$N2_N)|datai$warming=='N'|datai$soil_temp<5]=0  #NA赋值给0
  datai$N2_N.cum=cumsum(datai$N2_N)*24*10^4*10^-9  #kg N ha-1
  
  cum_p<-rbind(cum_p,datai)
}

cum_p$gasN=cum_p$NO_N.cum+cum_p$N2O_N.cum+cum_p$N2_N.cum

cum_p$treatment<-ifelse(cum_p$plot %in% c(1,4,5),'warmed','control')


######repeated measurement anov
# model=aov(Y ~ B * W + Error(Subject/W))，其中B是组间因子，W是组内因子，subject是实验对象的ID。
cum_p$NON2O<-cum_p$NO_N+cum_p$N2O_N

### NO
mod1<-aov(NO_N~treatment*date + Error(plot/date), cum_p[cum_p$warming=='Y',]);mod1
res1<-summary(mod1);res1

library(tidyr)
df<-do.call(rbind,res1$`Error: Within`)
df<-tibble::rownames_to_column(df,"data_code_market")
df<-tidyr::separate(df,data_code_market,into = c("factor"),sep="[.]")
df$factor
data.table::rbindlist(res1$`Error: Within`,use.names=TRUE, fill=TRUE, idcol="fenzu")

### N2O
mod2<-aov(N2O_N~treatment*date + Error(plot/date), cum_p[cum_p$warming=='Y',]);mod2
res2<-summary(mod2);res2


### N2
mod3<-aov(N2_N~treatment*date + Error(plot/date), cum_p[cum_p$warming=='Y',]);mod3
res3<-summary(mod3);res3


### gas N
mod3<-aov(gasN~treatment*date + Error(plot/date), cum_p[cum_p$warming=='Y',]);mod3
res3<-summary(mod3);res3


ggplot(cum_p[cum_p$warming=='Y',],aes(date,NO_N.cum,col=treatment))+
  geom_point(size=0.3)+
  scale_color_manual(limits=c('control','warmed'),values = c('blue','red'))

ggplot(cum_p[cum_p$warming=='Y',],aes(date,N2O_N.cum,col=treatment))+
  geom_point(size=0.3)+
  scale_color_manual(limits=c('control','warmed'),values = c('blue','red'))


ggplot(cum_p[cum_p$warming=='Y',],aes(date,N2_N.cum,col=treatment))+
  geom_point(size=0.3)+
  scale_color_manual(limits=c('control','warmed'),values = c('blue','red'))


cum_p1<-cum_p[which(cum_p$date == '2023/12/8'),]   ### cumulative emissions at this date

##增温期间NO和N2O累计释放
datFUN(cum_p1$NO_N.cum[cum_p1$treatment=='control'])
datFUN(cum_p1$NO_N.cum[cum_p1$treatment=='warmed'])

datFUN(cum_p1$N2O_N.cum[cum_p1$treatment=='control'])
datFUN(cum_p1$N2O_N.cum[cum_p1$treatment=='warmed'])

##累计释放减少百分比-NO
(mean(cum_p1$NO_N.cum[cum_p1$treatment=='warmed'])-mean(cum_p1$NO_N.cum[cum_p1$treatment=='control']))/
  mean(cum_p1$NO_N.cum[cum_p1$treatment=='control'])

##累计释放减少百分比-N2O
(mean(cum_p1$N2O_N.cum[cum_p1$treatment=='warmed'])-mean(cum_p1$N2O_N.cum[cum_p1$treatment=='control']))/
  mean(cum_p1$N2O_N.cum[cum_p1$treatment=='control'])

##累计释放减少百分比-N2
(mean(cum_p1$N2_N.cum[cum_p1$treatment=='warmed'])-mean(cum_p1$N2_N.cum[cum_p1$treatment=='control']))/
  mean(cum_p1$N2_N.cum[cum_p1$treatment=='control'])

##累计释放减少百分比-gasN
(mean(cum_p1$gasN[cum_p1$treatment=='warmed'])-mean(cum_p1$gasN[cum_p1$treatment=='control']))/
  mean(cum_p1$gasN[cum_p1$treatment=='control'])


###按年-逐日累积
mydata_c$year<-format(as.Date(mydata_c$date),"%Y")

cum_cy<-data.frame()

for(i in c(paste('E',seq(1:15),sep=''),paste('W',seq(1:15),sep=''))){
  for(j in 2019:2023){
    datai=mydata_c[mydata_c$system_id==i & mydata_c$year==j,
                   c('date','warming','system_id','plot','NO_N.m','N2O_N.m')]
    datai=datai[!is.na(datai$date),] 
    datai<-datai[order(datai$date),]
    # datai$NO_N.m[is.na(datai$NO_N.m)|datai$warming=='N']=0  #NA赋值给0
    datai$NO_N.m[is.na(datai$NO_N.m)]=0  #NA赋值给0
    datai$NO_N.cum=cumsum(datai$NO_N.m)*24*10^4*10^-9  #kg N ha-1
    
    # datai$N2O_N.m[is.na(datai$N2O_N.m)|datai$warming=='N']=0  #NA赋值给0
    datai$N2O_N.m[is.na(datai$N2O_N.m)]=0  #NA赋值给0
    datai$N2O_N.cum=cumsum(datai$N2O_N.m)*24*10^4*10^-9  #kg N ha-1
    cum_cy<-rbind(cum_cy,datai)
  }
}

cum_cy$NXO=cum_cy$NO_N.cum+cum_cy$N2O_N.cum

cum_cy$treatment<-ifelse(cum_cy$plot %in% c(1,4,5),'warmed','control')

ggplot(cum_cy,aes(date,NO_N.cum,col=treatment))+
  geom_point(size=0.3)+
  scale_color_manual(limits=c('control','warmed'),values = c('blue','red'))

ggplot(cum_cy,aes(date,N2O_N.cum,col=treatment))+
  geom_point(size=0.3)+
  scale_color_manual(limits=c('control','warmed'),values = c('blue','red'))

###按年累积排放
treatments<-c("control","warmed");treatments
# labs_warming<-unique(mydata_d$labs_warming);labs_warming
str(mydata_d)
str(rainfall)

mydata_d<-merge(mydata_d,rainfall[,c('date','T_mean','precipitation')],
                by=c('date'),all.x = T)

str(mydata_d)
# mydata_d$rainy<-ifelse(mydata_d$precipitation>=2,'Y','N')


mydata_d$season_period<-
  ifelse(format(mydata_d$date,"%m") %in% c('05','06','07','08','09','10'),'GS','NGS')

unique(mydata_d$warming)

emit_g<-data.frame()   #增温期间各plot有效排放尤其是气态氮
k=1
for(i in 1:6){
  for(t in c('all','growing','warming')){
    for(y in 2018:2023){
      if(t=='warming'){
        datai=mydata_d[which(mydata_d$plot==i & mydata_d$warming=='Y' & mydata_d$year==y),]
      }else if (t == 'all'){
        datai=mydata_d[which(mydata_d$plot==i & mydata_d$year==y &
                               mydata_d$warming %in% c('Y','N')),]
      }else if (t == 'growing'){
        datai=mydata_d[which(mydata_d$plot==i & mydata_d$year==y &
                               mydata_d$season_period=='GS'),]
      }
      if(length(datai$NO_N)>0){
        emit_g[k,1]=i   #plot num
        emit_g[k,2]=t    #labs_warming
        emit_g[k,3]=y   #year
        
        dataii<-datai[datai$soil_temp>=5 & datai$warming=='Y',]   #专为N2
        
        emit_g[k,4]=paste(round(min(na.omit(dataii$N2_N)),2),
                          round(max(na.omit(dataii$N2_N)),2),sep="<>")  #N2_field排放范围
        emit_g[k,5]=round(mean(na.omit(dataii$N2_N),na.rm=T),3)  #N2_field平均值
        emit_g[k,6]=length(dataii$date[complete.cases(dataii$N2_N)])  #N2_field实测days
        emit_g[k,7]=round(emit_g[k,5]*emit_g[k,6]*24*10^4/10^9,3)   #N2_field单位kg N ha-1 
        emit_g[k,8]=round(sd(na.omit(dataii$N2_N))/sqrt(emit_g[k,6])*24*10^4/10^9,3)
        
        
        emit_g[k,9]=paste(round(min(na.omit(datai$N2O_N)),2),
                          round(max(na.omit(datai$N2O_N)),2),sep="<>")  #N2O排放范围
        emit_g[k,10]=round(mean(na.omit(datai$N2O_N),na.rm=T),3)  #N2O平均值
        emit_g[k,11]=length(datai$date[complete.cases(datai$N2O_N)])  #N2O实测days
        emit_g[k,12]=round(emit_g[k,11]*emit_g[k,10]*24*10^4/10^9,3)   #单位kg N ha-1 
        emit_g[k,13]=round(sd(na.omit(datai$N2O_N),na.rm=T)/sqrt(emit_g[k,11])*24*10^4/10^9,3)
        
        emit_g[k,14]=paste(round(min(na.omit(datai$NO_N)),2),
                           round(max(na.omit(datai$NO_N)),2),sep="<>")   #NO排放范围
        emit_g[k,15]=round(mean(na.omit(datai$NO_N),na.rm=T),3)  #NO平均值
        emit_g[k,16]=length(datai$date[complete.cases(datai$NO_N)]) 
        emit_g[k,17]=round(emit_g[k,16]*emit_g[k,15]*24*10^4/10^9,3)  #单位kg N ha-1 
        emit_g[k,18]=round(sd(na.omit(datai$NO_N),na.rm=T)/sqrt(emit_g[k,16])*24*10^4/10^9,3)
        

        emit_g[k,19]=paste(round(min(na.omit(datai$soil_temp[complete.cases(datai$NO_N)])),2),
                           round(max(na.omit(datai$soil_temp[complete.cases(datai$NO_N)])),2),sep="<>")   
        emit_g[k,20]=round(mean(na.omit(datai$soil_temp[complete.cases(datai$NO_N)]),na.rm=T),3)  
        emit_g[k,21]=sd(na.omit(datai$soil_temp[complete.cases(datai$NO_N)]),na.rm=T)/emit_g[k,20]  
        
        emit_g[k,22]=paste(round(min(na.omit(datai$soil_moisture[complete.cases(datai$NO_N)])),2),
                           round(max(na.omit(datai$soil_moisture[complete.cases(datai$NO_N)])),2),sep="<>")   
        emit_g[k,23]=round(mean(na.omit(datai$soil_moisture[complete.cases(datai$NO_N)]),na.rm=T),3)  
        emit_g[k,24]=sd(na.omit(datai$soil_moisture[complete.cases(datai$NO_N)]),na.rm=T)/emit_g[k,23] 
        k=k+1
      }
    }
  }
}

######单位N2O_N,NO_N,(g N ha-1)

names(emit_g)[1:24]<-c('plot','labs_warming','year',
                       'N2_field_range','N2_field.m','days_N2_field','E_N2_field','N2_field.se',
                       'N2O_range','N2O.m','days_N2O','E_N2O','N2O.se',
                       'NO_range','NO.m','days_NO','E_NO','NO.se',
                       'temp_range','temp','temp.cv','VWC_range','VWC','VWC.cv')

emit_g$treatment<-ifelse(emit_g$plot %in% c("1","4","5"),"warmed","control")

unique(emit_g$labs_warming)  ###all the obs, growing season, warming period

emit_gt<-emit_g[emit_g$labs_warming=='all',]   #全部观测数据
emit_gg<-emit_g[emit_g$labs_warming=='growing',]   #生长季观测数据
emit_g<-emit_g[emit_g$labs_warming=='warming',]  #增温观测数据

str(emit_gt)   #'E_N2_field', 'E_N2O','E_NO'
emit_gt.yr<-summaryBy(data=emit_gt,emit_gt[,c('E_N2_field', 'E_N2O','E_NO')]~treatment+year,FUN=myfun)
emit_gg.yr<-summaryBy(data=emit_gg,emit_gg[,c('E_N2_field', 'E_N2O','E_NO')]~treatment+year,FUN=myfun)
emit_g.yr<-summaryBy(data=emit_g,emit_g[,c('E_N2_field', 'E_N2O','E_NO')]~treatment+year,FUN=myfun)

###将生长季排放和所有观测放一起all+growing, 计算冻融期间的排放


###计算3-4月，5-10月，11月冻融期间的排放
str(F.treat)
F.treat$season_cum<-
  ifelse(format(F.treat$date,'%m') %in% c('03','04'),'spring_freeze-thaw',
         ifelse(format(F.treat$date,'%m') %in% c('05','06','07','08','09','10'),'non_freeze-thaw',
                ifelse(format(F.treat$date,'%m') %in% c('11','12'),'winter_freeze-thaw', 'exclude'       
                )))

cum_t<-data.frame();k=1
for(i in 2019:2023){
  for(j in c('spring_freeze-thaw','non_freeze-thaw','winter_freeze-thaw')){
    for(t in c('control','warmed')){
      dat_ijt=F.treat[F.treat$treatment==t & F.treat$season_cum==j & F.treat$year==i & F.treat$warming=='Y',]
      cum_t[k,1]=i
      cum_t[k,2]=t
      cum_t[k,3]=j
      cum_t[k,4]=mean(na.omit(dat_ijt$NO_N.m))*length(na.omit(dat_ijt$NO_N.m))*24*10000/10^6
      cum_t[k,5]=mean(na.omit(dat_ijt$N2O_N.m))*length(na.omit(dat_ijt$N2O_N.m))*24*10000/10^6
      cum_t[k,6]=mean(na.omit(dat_ijt$soil_temp.m))
      cum_t[k,7]=mean(na.omit(dat_ijt$WFPS.m))
      k=k+1
    }
  }
}

names(cum_t)<-c('year','treatment','season_cum','NO_cum','N2O_cum','Temp','WFPS')


str(F.treat)
unique(F.treat$season)

# 'spring_freeze-thaw','non_freeze-thaw','winter_freeze-thaw'



cum_s<-data.frame();k=1
for(i in 2019:2023){
  for(j in c(c('growing season','dormant season'))){
    for(t in c('control','warmed')){
      if(j=='growing season'){
        dat_ijt<-F.treat[F.treat$treatment==t & F.treat$season_cum=='non_freeze-thaw' & F.treat$year==i,]
      }else{
        dat_ijt<-F.treat[F.treat$treatment==t & F.treat$season_cum!='non_freeze-thaw' & F.treat$year==i,]
      }
      
      cum_s[k,1]=i
      cum_s[k,2]=t
      cum_s[k,3]=j
      cum_s[k,4]=mean(na.omit(dat_ijt$NO_N.m))*length(na.omit(dat_ijt$NO_N.m))*24*10000/10^6
      cum_s[k,5]=mean(na.omit(dat_ijt$N2O_N.m))*length(na.omit(dat_ijt$N2O_N.m))*24*10000/10^6
      cum_s[k,6]=mean(na.omit(dat_ijt$soil_temp.m))
      cum_s[k,7]=mean(na.omit(dat_ijt$WFPS.m))
      k=k+1
    }
  }
}


names(cum_s)<-c('year','treatment','season_cum','NO_cum','N2O_cum','Temp','WFPS')



### n=15 chamber水平（增温期间）
treatments<-c("control","warmed");treatments
labs_warming<-unique(mydata_d$labs_warming);labs_warming
str(mydata_c)
unique(mydata_c$system_id)

unique(mydata_c$warming)  #分为Y+N

emit_c<-data.frame()   #增温期间各chamber有效排放尤其是气态氮
k=1

str(mydata_c)

mydata_c$season<-
  ifelse(mydata_c$date>='2018/5/1' & mydata_c$date<='2018/10/30','growing',
         ifelse(mydata_c$date>='2019/5/1' & mydata_c$date<='2019/10/30','growing',
                ifelse(mydata_c$date>='2020/5/1' & mydata_c$date<='2020/10/30','growing',
                       ifelse(mydata_c$date>='2021/5/1' & mydata_c$date<='2021/10/30','growing',
                              ifelse(mydata_c$date>='2022/5/1' & mydata_c$date<='2022/10/30','growing',
                                     ifelse(mydata_c$date>='2023/5/1' & mydata_c$date<='2023/10/30','growing',
                                            ifelse(mydata_c$date>='2024/5/1' & mydata_c$date<='2024/10/30','growing','dormant'
)))))))


for(t in c('growing','dormant')){
  for(y in 2018:2023){
    dat=mydata_c[which(mydata_c$season==t & mydata_c$year==y),]
    
    for(i in na.omit(unique(dat$system_id))){ 
      
      datai=dat[which(dat$system_id==i),]
      
      if(length(datai$NO_N.m)>0){
        emit_c[k,1]=i   #chamber num
        emit_c[k,2]=t    #season
        emit_c[k,3]=y   #year
        
        dataii<-datai[datai$soil_temp.m>=5 & datai$warming=='Y',]   #专为N2
        
        
        emit_c[k,4]=paste(round(min(na.omit(dataii$N2_N.m)),2),
                           round(max(na.omit(dataii$N2_N.m)),2),sep="<>")  #N2_field排放范围
        emit_c[k,5]=round(mean(na.omit(dataii$N2_N.m),na.rm=T),3)  #N2_field平均值
        emit_c[k,6]=length(dataii$date[complete.cases(dataii$N2_N.m)])  #N2_field实测days
        emit_c[k,7]=round(emit_c[k,5]*emit_c[k,6]*24*10^4/10^9,3)   #N2_field单位kg N ha-1 
        emit_c[k,8]=round(sd(na.omit(dataii$N2_N.m))/sqrt(emit_c[k,6])*24*10^4/10^9,3)
        
        
        emit_c[k,9]=paste(round(min(na.omit(datai$N2O_N.m)),2),
                           round(max(na.omit(datai$N2O_N.m)),2),sep="<>")  #N2O排放范围
        emit_c[k,10]=round(mean(na.omit(datai$N2O_N.m),na.rm=T),3)  #N2O平均值
        emit_c[k,11]=length(datai$date[complete.cases(datai$N2O_N.m)])  #N2O实测days
        emit_c[k,12]=round(emit_c[k,11]*emit_c[k,10]*24*10^4/10^9,3)   #单位kg N ha-1 
        emit_c[k,13]=round(sd(na.omit(datai$N2O_N.m),na.rm=T)/sqrt(emit_c[k,11])*24*10^4/10^9,3)
        
        emit_c[k,14]=paste(round(min(na.omit(datai$NO_N.m)),2),
                            round(max(na.omit(datai$NO_N.m)),2),sep="<>")   #NO排放范围
        emit_c[k,15]=round(mean(na.omit(datai$NO_N.m),na.rm=T),3)  #NO平均值
        emit_c[k,16]=length(datai$date[complete.cases(datai$NO_N.m)]) 
        emit_c[k,17]=round(emit_c[k,16]*emit_c[k,15]*24*10^4/10^9,3)  #单位kg N ha-1 
        emit_c[k,18]=round(sd(na.omit(datai$NO_N.m),na.rm=T)/sqrt(emit_c[k,16])*24*10^4/10^9,3)
        
        emit_c[k,19]=paste(round(min(na.omit(datai$soil_temp.m[complete.cases(datai$NO_N.m)])),2),
                            round(max(na.omit(datai$soil_temp.m[complete.cases(datai$NO_N.m)])),2),sep="<>")   
        emit_c[k,20]=round(mean(na.omit(datai$soil_temp.m[complete.cases(datai$NO_N.m)]),na.rm=T),3)  
        emit_c[k,21]=sd(na.omit(datai$soil_temp.m[complete.cases(datai$NO_N.m)]),na.rm=T)/emit_c[k,20]  
        
        emit_c[k,22]=paste(round(min(na.omit(datai$soil_moisture.m[complete.cases(datai$NO_N.m)])),2),
                            round(max(na.omit(datai$soil_moisture.m[complete.cases(datai$NO_N.m)])),2),sep="<>")   
        emit_c[k,23]=round(mean(na.omit(datai$soil_moisture.m[complete.cases(datai$NO_N.m)]),na.rm=T),3)  
        emit_c[k,24]=sd(na.omit(datai$soil_moisture.m[complete.cases(datai$NO_N.m)]),na.rm=T)/emit_c[k,23] 
        
        emit_c[k,25]=unique(datai$plot)
        emit_c[k,26]=unique(datai$treatment)
        
        k=k+1
        
      }
    }
  }
}

######单位N2O_N,NO_N,(g N ha-1)

names(emit_c)[1:26]<-c('chamber','season','year',
                        'N2_field_range','N2_field.m','days_N2_field','E_N2_field','N2_field.se',
                        'N2O_range','N2O.m','days_N2O','E_N2O','N2O.se',
                        'NO_range','NO.m','days_NO','E_NO','NO.se',
                        'temp_range','temp','temp.cv','VWC_range','VWC','VWC.cv','plot','treatment')

emit_c$treatment<-ifelse(emit_c$plot %in% c("1","4","5"),"warmed","control")

str(emit_c)

emit_c.treat<-summaryBy(data=emit_c,E_N2O+E_NO+E_N2_field~treatment+year+season,FUN=myfun)


library(tidyverse)
library(writexl)

###将不同的结果卸载一个xls文件
df = list(emit_gt,emit_gg,emit_g,emit_gt.yr,emit_gg.yr,emit_g.yr,cum_t, emit_c.treat,cum_s) %>%
  set_names(c('all obs','growing obs','warming obs','all_treat level',
              'growing_treat level','warming_treat level','warming_season','obs_season for chamber','obs_season for treatment'))

write_xlsx(df, "gasN emissions_plot_2018-2023.xlsx")   ###文章数据，

getwd()



###增温期间气态氮释放的响应
emit_g$labs_warming<-as.factor(emit_g$labs_warming)
  
  emit_g$gaseous_N<-emit_g$E_N2_field+emit_g$E_N2O+emit_g$E_NO
  
  emit_g<-within(emit_g,plot<-factor(plot,levels=c("1","4","5","2","3","6")))
  
  levels(emit_g$plot)
  
  str(emit_g)
  
  unique(emit_g$labs_warming)
  
  emit_g$E_N2_field.q10model<-emit_g$E_N2_field*2.6^(1/5)  # NCC N2q10 = 3.0 (QYM)
  emit_g$E_N2O.q10model<-emit_g$E_N2O*2.1^(1/5)  # NCC N2O q10 = 2.4 (QYM)
  emit_g$E_NO.q10model<-emit_g$E_NO*3.5^(1/5)  # NO q10 = 3.5 (QYM)
  
  dat_1<-emit_g[,c('year','plot','E_N2_field','E_N2O','E_NO')]
  dat_1$treatment<-ifelse(dat_1$plot %in% c(1,4,5),'warmed','control')
  
  dat_2<-emit_g[,c('year','plot','E_N2_field.q10model','E_N2O.q10model','E_NO.q10model')]
  
  names(dat_2)<-c('year','plot','E_N2_field','E_N2O','E_NO')
  dat_2$treatment<-'model'
  dat_2<-dat_2[dat_2$plot %in% c(2,3,6),]
  
  
  dat_fig1<-rbind(dat_1,dat_2)
  
  dat_fig1.t<-summaryBy(data=dat_fig1, E_N2_field+E_N2O+E_NO~year+treatment,FUN=myfun)
  
  

  ###多重比较control\warmed\model
  library(dplyr)
  result<- data.frame()
  for(i in 2018:2023){
    
    dfa<- dat_fig1 %>%
      filter(year == i)
    
    aov<- aov(E_NO~treatment, dfa)
    lsd<- LSD.test(aov, "treatment", p.adj = "none")
    
    lsdr1<- lsd$groups %>% as.data.frame() %>%
      mutate(treatment = rownames(lsd$groups),
             year= i) %>% dplyr::select(treatment, year,  groups)
    names(lsdr1)[3]<-'group_NO'
    
    
    
    aov<- aov(E_N2O~treatment, dfa)
    lsd<- LSD.test(aov, "treatment", p.adj = "none")
    
    lsdr2<- lsd$groups %>% as.data.frame() %>%
      mutate(treatment = rownames(lsd$groups),
             year= i) %>%
      dplyr::select(treatment, year,  groups)
    names(lsdr2)[3]<-'group_N2O'
    
    aov<- aov(E_N2_field~treatment, dfa)
    lsd<- LSD.test(aov, "treatment", p.adj = "none")
    
    lsdr3<- lsd$groups %>% as.data.frame() %>%
      mutate(treatment = rownames(lsd$groups),
             year= i) %>%
      dplyr::select(treatment, year,  groups)
    names(lsdr3)[3]<-'group_N2'
    
    
    lsdr<-merge(lsdr1,lsdr2,by=c('year','treatment'))
    lsdr<-merge(lsdr,lsdr3,by=c('year','treatment'))
    rm(lsdr1,lsdr2,lsdr3)
    
    result<- rbind(result, lsdr)
  }
  
  
  
  str(dat_fig1.t)
  str(result)
  dat_fig1.t<-merge(dat_fig1.t,result,by=c('year','treatment'),all=T)
  
  unique(dat_fig1.t$treatment)
  dat_fig1.t<-within(dat_fig1.t,treatment<-factor(treatment,levels = c("control","warmed","model")))
  
  summary(aov(E_NO~treatment*year+Error(plot/year), data=dat_fig1[dat_fig1$treatment!='model',]))
  summary(aov(E_N2O~treatment*year+Error(plot/year), data=dat_fig1[dat_fig1$treatment!='model',]))
  summary(aov(E_N2_field~treatment*year+Error(plot/year), data=dat_fig1[dat_fig1$treatment!='model',]))
  
  
  dat_fig1.t<-within(dat_fig1.t,treatment<-factor(treatment,levels = c('control','model','warmed')))
  
  ###将6年累积排放图并入
  write.csv(emit_g,'增温期间气态氮损失.csv',row.names = F)
  

### n=15 chamber水平（增温期间）
treatments<-c("control","warmed");treatments
labs_warming<-unique(mydata_d$labs_warming);labs_warming
str(mydata_c)
unique(mydata_c$system_id)

unique(mydata_c$warming)  #分为Y+N

emit_gc<-data.frame()   #增温期间各chamber有效排放尤其是气态氮
k=1

for(t in c('Y')){
  for(y in 2018:2023){
    dat=mydata_c[which(mydata_c$warming==t & mydata_c$year==y),]
    
    for(i in na.omit(unique(dat$system_id))){ 
      datai=dat[which(dat$system_id==i),]
      
      if(length(datai$NO_N.m)>0){
        emit_gc[k,1]=i   #chamber num
        emit_gc[k,2]=t    #labs_warming
        emit_gc[k,3]=y   #year
        
        dataii<-datai[datai$soil_temp.m>=5 & datai$warming=='Y',]   #专为N2
        
        
        emit_gc[k,4]=paste(round(min(na.omit(dataii$N2_N.m)),2),
                           round(max(na.omit(dataii$N2_N.m)),2),sep="<>")  #N2_field排放范围
        emit_gc[k,5]=round(mean(na.omit(dataii$N2_N.m),na.rm=T),3)  #N2_field平均值
        emit_gc[k,6]=length(dataii$date[complete.cases(dataii$N2_N.m)])  #N2_field实测days
        emit_gc[k,7]=round(emit_gc[k,5]*emit_gc[k,6]*24*10^4/10^9,3)   #N2_field单位kg N ha-1 
        emit_gc[k,8]=round(sd(na.omit(dataii$N2_N.m))/sqrt(emit_gc[k,6])*24*10^4/10^9,3)
        
        
        emit_gc[k,9]=paste(round(min(na.omit(datai$N2O_N.m)),2),
                           round(max(na.omit(datai$N2O_N.m)),2),sep="<>")  #N2O排放范围
        emit_gc[k,10]=round(mean(na.omit(datai$N2O_N.m),na.rm=T),3)  #N2O平均值
        emit_gc[k,11]=length(datai$date[complete.cases(datai$N2O_N.m)])  #N2O实测days
        emit_gc[k,12]=round(emit_gc[k,11]*emit_gc[k,10]*24*10^4/10^9,3)   #单位kg N ha-1 
        emit_gc[k,13]=round(sd(na.omit(datai$N2O_N.m),na.rm=T)/sqrt(emit_gc[k,11])*24*10^4/10^9,3)
        
        emit_gc[k,14]=paste(round(min(na.omit(datai$NO_N.m)),2),
                            round(max(na.omit(datai$NO_N.m)),2),sep="<>")   #NO排放范围
        emit_gc[k,15]=round(mean(na.omit(datai$NO_N.m),na.rm=T),3)  #NO平均值
        emit_gc[k,16]=length(datai$date[complete.cases(datai$NO_N.m)]) 
        emit_gc[k,17]=round(emit_gc[k,16]*emit_gc[k,15]*24*10^4/10^9,3)  #单位kg N ha-1 
        emit_gc[k,18]=round(sd(na.omit(datai$NO_N.m),na.rm=T)/sqrt(emit_gc[k,16])*24*10^4/10^9,3)
        
        emit_gc[k,19]=paste(round(min(na.omit(datai$soil_temp.m[complete.cases(datai$NO_N.m)])),2),
                            round(max(na.omit(datai$soil_temp.m[complete.cases(datai$NO_N.m)])),2),sep="<>")   
        emit_gc[k,20]=round(mean(na.omit(datai$soil_temp.m[complete.cases(datai$NO_N.m)]),na.rm=T),3)  
        emit_gc[k,21]=sd(na.omit(datai$soil_temp.m[complete.cases(datai$NO_N.m)]),na.rm=T)/emit_gc[k,20]  
        
        emit_gc[k,22]=paste(round(min(na.omit(datai$soil_moisture.m[complete.cases(datai$NO_N.m)])),2),
                            round(max(na.omit(datai$soil_moisture.m[complete.cases(datai$NO_N.m)])),2),sep="<>")   
        emit_gc[k,23]=round(mean(na.omit(datai$soil_moisture.m[complete.cases(datai$NO_N.m)]),na.rm=T),3)  
        emit_gc[k,24]=sd(na.omit(datai$soil_moisture.m[complete.cases(datai$NO_N.m)]),na.rm=T)/emit_gc[k,23] 
        
        emit_gc[k,25]=unique(datai$plot)
        emit_gc[k,26]=unique(datai$treatment)
        
        k=k+1
      }
    }
  }
}

######单位N2O_N,NO_N,(g N ha-1)

names(emit_gc)[1:26]<-c('chamber','labs_warming','year',
                        'N2_field_range','N2_field.m','days_N2_field','E_N2_field','N2_field.se',
                        'N2O_range','N2O.m','days_N2O','E_N2O','N2O.se',
                        'NO_range','NO.m','days_NO','E_NO','NO.se',
                        'temp_range','temp','temp.cv','VWC_range','VWC','VWC.cv','plot','treatment')

emit_gc$treatment<-ifelse(emit_gc$plot %in% c("1","4","5"),"warmed","control")

###增温期间
datFUN(emit_gc$E_N2_field[which(emit_gc$labs_warming=="Y" & emit_gc$treatment=="control"& emit_gc$year!='2018')])
datFUN(emit_gc$E_N2_field[which(emit_gc$labs_warming=="Y" & emit_gc$treatment=="warmed"& emit_gc$year!='2018')])
###促进了N2释放？？促进了冻融期间N2释放


##########日排放平均值和累积排放量的差异
levels(emit_gc$plot)

str(emit_gc)

unique(emit_gc$labs_warming)

emit_gc$N2_field.q10model<-emit_gc$N2_field.m*2.6^(1/5)  # N2q10 = 3.7 (QYM 30% WFPS)
emit_gc$N2O.q10model<-emit_gc$N2O.m *2.1^(1/5)  # N2O q10 = 2.4 (QYM 30% WFPS)
emit_gc$NO.q10model<-emit_gc$NO.m*3.49^(1/5)  # NO q10 = 3.6 (QYM 60% WFPS)

emit_gc$E_N2_field.q10model<-emit_gc$E_N2_field*2.6^(1/5)  # N2q10 = 3.7 (QYM 30% WFPS)
emit_gc$E_N2O.q10model<-emit_gc$E_N2O *2.1^(1/5)  # N2O q10 = 2.4 (QYM 30% WFPS)
emit_gc$E_NO.q10model<-emit_gc$E_NO*3.49^(1/5)  # NO q10 = 3.6 (QYM 60% WFPS)

emit_gc$NON2O.m<-emit_gc$N2O.m + emit_gc$NO.m

emit_gc$gas_N<-emit_gc$N2_field.m + emit_gc$N2O.m + emit_gc$NO.m

emit_gc$NON2O_N.q10model<-emit_gc$N2O.q10model+emit_gc$NO.q10model 

emit_gc$gaseous_N.q10model<-emit_gc$N2_field.q10model+ emit_gc$N2O.q10model+emit_gc$NO.q10model 


names(emit_gc)
dat_1<-emit_gc[,c('year','chamber','plot','N2_field.q10model','N2O.q10model','NO.q10model','gaseous_N.q10model','NON2O_N.q10model',
                  'E_N2_field.q10model','E_NO.q10model','E_N2O.q10model')]


dat_2<-emit_gc[,c('year','chamber','plot','N2_field.m','N2O.m','NO.m','gas_N','NON2O.m',
                  'E_N2_field','E_NO','E_N2O')]

dat_2$treatment<-ifelse(dat_2$plot %in% c(1,4,5),'warmed','control')

names(dat_1)<-c('year','chamber','plot','N2_field.m','N2O.m','NO.m','gas_N','NON2O.m','E_N2_field','E_NO','E_N2O')
dat_1<-dat_1[dat_1$plot %in% c(2,3,6),]
dat_1$treatment<-'model'


names(dat_1)
names(dat_2)
dat_fig0<-rbind(dat_1,dat_2)   ###

dat_fig0$chamber.treat<-paste(dat_fig0$chamber,dat_fig0$treatment,sep = '-')

dat_fig0$NO_stat.lab<-'Y'
dat_fig0$N2O_stat.lab<-'Y'
dat_fig0$NON2O_stat.lab<-'Y'
dat_fig0$N2_stat.lab<-'Y'
dat_fig0$gasN_stat.lab<-'Y'


###多重比较control\warmed\model
result0<- data.frame()
for(i in 2018:2023){
  
  dfa<- dat_fig0 %>%
    filter(year == i)
  
  x1c=dfa$chamber.treat[dfa$NO.m %in% boxplot(dfa$NO.m[dfa$treatment=='control'])$out]
  x2c=dfa$chamber.treat[dfa$N2O.m %in% boxplot(dfa$N2O.m[dfa$treatment=='control'])$out]
  x3c=dfa$chamber.treat[dfa$N2_field.m %in% boxplot(dfa$N2_field.m[dfa$treatment=='control'])$out]
  x4c=dfa$chamber.treat[dfa$gas_N %in% boxplot(dfa$gas_N[dfa$treatment=='control'])$out]
  x5c=dfa$chamber.treat[dfa$NON2O.m %in% boxplot(dfa$NON2O.m[dfa$treatment=='control'])$out]
  
  x1w=dfa$chamber.treat[dfa$NO.m %in% boxplot(dfa$NO.m[dfa$treatment=='warmed'])$out]
  x2w=dfa$chamber.treat[dfa$N2O.m %in% boxplot(dfa$N2O.m[dfa$treatment=='warmed'])$out]
  x3w=dfa$chamber.treat[dfa$N2_field.m %in% boxplot(dfa$N2_field.m[dfa$treatment=='warmed'])$out]
  x4w=dfa$chamber.treat[dfa$gas_N %in% boxplot(dfa$gas_N[dfa$treatment=='warmed'])$out]
  x5w=dfa$chamber.treat[dfa$NON2O.m %in% boxplot(dfa$NON2O.m[dfa$treatment=='warmed'])$out]
  
  x1m=dfa$chamber.treat[dfa$NO.m %in% boxplot(dfa$NO.m[dfa$treatment=='model'])$out]
  x2m=dfa$chamber.treat[dfa$N2O.m %in% boxplot(dfa$N2O.m[dfa$treatment=='model'])$out]
  x3m=dfa$chamber.treat[dfa$N2_field.m %in% boxplot(dfa$N2_field.m[dfa$treatment=='model'])$out]
  x4m=dfa$chamber.treat[dfa$gas_N %in% boxplot(dfa$gas_N[dfa$treatment=='model'])$out]
  x5m=dfa$chamber.treat[dfa$NON2O.m %in% boxplot(dfa$NON2O.m[dfa$treatment=='model'])$out]
  
  
  ###离群值进行标记
  dat_fig0[which(dat_fig0$year==i & dat_fig0$chamber.treat %in% unique(c(x1c,x1w,x1m))),'NO_stat.lab']='N'
  dat_fig0[which(dat_fig0$year==i & dat_fig0$chamber.treat %in% unique(c(x2c,x2w,x2m))),'N2O_stat.lab']='N'
  dat_fig0[which(dat_fig0$year==i & dat_fig0$chamber.treat %in% unique(c(x5c,x5w,x5m))),'NON2O_stat.lab']='N'
  dat_fig0[which(dat_fig0$year==i & dat_fig0$chamber.treat %in% unique(c(x3c,x3w,x3m))),'N2_stat.lab']='N'
  dat_fig0[which(dat_fig0$year==i & dat_fig0$chamber.treat %in% unique(c(x4c,x4w,x4m))),'gasN_stat.lab']='N'
  
  ###统计时去除离群值
  dfa1<-dfa[!dfa$chamber %in% unique(c(x1c,x1w,x1m)),]
  dfa2<-dfa[!dfa$chamber %in% unique(c(x2c,x2w,x2m)),]
  dfa5<-dfa[!dfa$chamber %in% unique(c(x5c,x5w,x5m)),]
  dfa3<-dfa[!dfa$chamber %in% unique(c(x3c,x3w,x3m)),]
  
  aov<- aov(NO.m~treatment, dfa1)
  lsd<- LSD.test(aov, "treatment", p.adj = "none")
  lsdr1<- lsd$groups %>% as.data.frame() %>%
    mutate(treatment = rownames(lsd$groups),
           year= i) %>%
    dplyr::select(treatment, year,  groups)
  names(lsdr1)[3]<-'group_NO'
  
  
  
  aov<- aov(N2O.m~treatment, dfa2)
  lsd<- LSD.test(aov, "treatment", p.adj = "none")
  lsdr2<- lsd$groups %>% as.data.frame() %>%
    mutate(treatment = rownames(lsd$groups),
           year= i) %>%
    dplyr::select(treatment, year,  groups)
  names(lsdr2)[3]<-'group_N2O'
  
  aov<- aov(N2_field.m~treatment, dfa3)
  lsd<- LSD.test(aov, "treatment", p.adj = "none")
  lsdr3<- lsd$groups %>% as.data.frame() %>%
    mutate(treatment = rownames(lsd$groups),
           year= i) %>%
    dplyr::select(treatment, year,  groups)
  names(lsdr3)[3]<-'group_N2'
  
  aov<- aov(gas_N~treatment, dfa3)
  lsd<- LSD.test(aov, "treatment", p.adj = "none")
  lsdr4<- lsd$groups %>% as.data.frame() %>%
    mutate(treatment = rownames(lsd$groups),
           year= i) %>%
    dplyr::select(treatment, year,  groups)
  names(lsdr4)[3]<-'group_gasN'
  
  aov<- aov(NON2O.m~treatment, dfa5)
  lsd<- LSD.test(aov, "treatment", p.adj = "none")
  lsdr5<- lsd$groups %>% as.data.frame() %>%
    mutate(treatment = rownames(lsd$groups),
           year= i) %>%
    dplyr::select(treatment, year,  groups)
  names(lsdr5)[3]<-'group_NON2O'
  
  lsdr<-merge(lsdr1,lsdr2,by=c('year','treatment'))
  lsdr<-merge(lsdr,lsdr3,by=c('year','treatment'))
  lsdr<-merge(lsdr,lsdr4,by=c('year','treatment'))
  lsdr<-merge(lsdr,lsdr5,by=c('year','treatment'))
  
  rm(lsdr1,lsdr2,lsdr3,lsdr4,lsdr5)
  
  result0<- rbind(result0, lsdr)
}


str(result0)
str(dat_fig0)
dat_ex<-dat_fig0  #将NA赋值给异常值,在求算处理水平均值
dat_ex$NO.m[dat_ex$NO_stat.lab=='N']=NA
dat_ex$N2O.m[dat_ex$N2O_stat.lab=='N']=NA
dat_ex$NON2O.m[dat_ex$NON2O_stat.lab=='N']=NA
dat_ex$N2_field.m[dat_ex$N2_stat.lab=='N']=NA
dat_ex$gas_N[dat_ex$gasN_stat.lab=='N']=NA

dat_fig0.t<-summaryBy(data=dat_ex, N2_field.m+N2O.m+NO.m+NON2O.m+gas_N~year+treatment,FUN=myfun)

str(dat_fig0.t)

dat_fig0.t<-merge(dat_fig0.t,result0,by=c('year','treatment'),all=T)

unique(dat_fig0.t$treatment)
dat_fig0.t<-within(dat_fig0.t,treatment<-factor(treatment,levels = c("control","warmed","model")))

summary(aov(NO.m~treatment*year+Error(year), data=dat_fig0[dat_fig0$treatment!='model',]))
summary(aov(N2O.m~treatment*year+Error(year), data=dat_fig0[dat_fig0$treatment!='model',]))
summary(aov(N2_field.m~treatment*year+Error(year), data=dat_fig0[dat_fig0$treatment!='model',]))

dat_fig0<-within(dat_fig0,treatment<-factor(treatment,levels = c("control","warmed","model")))
dat_fig0.t<-within(dat_fig0.t,treatment<-factor(treatment,levels = c("control","warmed","model")))

str(dat_fig0.t) 

str(dat_fig0)

dat<-dat_fig0[which(dat_fig0$treatment%in% c('control','warmed') & dat_fig0$year %in% seq(2018,2023)),
              c('gas_N','N2_field.m','NO.m','N2O.m','NON2O.m','year','treatment','chamber','plot',
                'NO_stat.lab','N2O_stat.lab','NON2O_stat.lab','N2_stat.lab','gasN_stat.lab')]

str(dat)
dat$id<-as.factor(seq(1,length(dat$year)))
dat$treatment<-as.factor(as.character(dat$treatment))
str(dat)
unique(dat$chamber)

###关系chamber带来的差异
summary(aov(scale(NO.m)~treatment*year+Error(chamber/year), 
            data=dat[dat$NO_stat.lab=='Y' & dat$year>=2019,]))

summary(aov(scale(N2O.m)~treatment*year+Error(chamber/year), 
            data=dat[dat$N2O_stat.lab=='Y' & dat$year>=2019,]))

summary(aov(scale(N2_field.m)~treatment*year+Error(chamber/year), 
            data=dat[dat$N2_stat.lab=='Y' & dat$year>=2019,]))

summary(aov(scale(gas_N)~treatment*year+Error(chamber/year), 
            data=dat[dat$gasN_stat.lab=='Y' & dat$year>=2019,]))

summary(aov(scale(NON2O.m)~treatment*year+Error(chamber/year), 
            data=dat[dat$NON2O_stat.lab=='Y'& dat$year>=2019,]))

########## 气态氮释放混合效应模型分析 #############
#####勿动该部分代码_20240530
str(dat_ex)

dat<-dat_ex[which(dat_ex$treatment !='model',),
            c('chamber','treatment','NO.m','N2O.m','N2_field.m','plot','year','NO_stat.lab','N2O_stat.lab','N2_stat.lab','NON2O_stat.lab')]


dat$ftreatment<-as.numeric(ifelse(dat$treatment=="control",0,2))

library(lme4);library(sjPlot)

### (1) plot level
str(emit_gg)
unique(emit_gg$treatment)
dat<-emit_gg[,c('treatment','NO.m','N2O.m','N2_field.m','plot','year')]

dat$E_NO<-dat$NO.m
dat$E_N2O<-dat$N2O.m
dat$E_N2_field<-dat$N2_field.m

dat$year=as.numeric(as.character(dat$year))
dat$ftreatment<-as.numeric(ifelse(dat$treatment=="control",0,2))

str(dat)
#NO  混合效应模型
fm1<-lmer(scale(E_NO)~ftreatment*scale(year)+(1|plot)+(1|year),
          data=dat[dat$year>=2018,]);fm1

summary(fm1)
tab_model(fm1)


#N2O
fm2<-lmer(scale(E_N2O)~ftreatment*scale(year)+(1|plot)+(1|year),
          data=dat[dat$year>=2018,]);fm2

summary(fm2)

tab_model(fm2)

#N2
fm3<-lmer(scale(E_N2_field)~ftreatment*scale(year)+(1|plot)+(1|year),
          data=dat[dat$year>=2019,]);fm3

summary(fm3)

tab_model(fm3)

#####(2) chamber level
str(dat_fig0)

dat<-dat_fig0[which(dat_fig0$treatment !='model',),
              c('chamber','treatment','NO.m','N2O.m','N2_field.m','NON2O.m','gas_N',
                'plot','year','NO_stat.lab','N2O_stat.lab','N2_stat.lab','NON2O_stat.lab')]

dat$year=as.numeric(as.character(dat$year))
dat$ftreatment<-as.numeric(ifelse(dat$treatment=="control",0,2))

str(dat)
########文章使用的分析方法
###NO  混合效应模型
fm1<-lmer(scale(NO.m)~ftreatment*scale(year)+(1|chamber)+(1|year),
          data=dat[dat$year>=2019 & dat$NO_stat.lab=='Y',]);fm1

summary(fm1)
tab_model(fm1)


###N2O remove the outlier
fm2<-lmer(scale(N2O.m)~ftreatment*scale(year)+(1|chamber)+(1|year),
          data=dat[dat$year>=2019 & dat$N2O_stat.lab=='Y',]);fm2

summary(fm2)

tab_model(fm2)

###N2 remove the outlier
fm3<-lmer(scale(N2_field.m)~ftreatment*scale(year)+(1|chamber)+(1|year),
          data=dat[dat$year>=2019 & dat$N2_stat.lab=='Y',]);fm3

summary(fm3)

tab_model(fm3)


###NO+N2O
fm3<-lmer(scale(NON2O.m)~ftreatment*scale(year)+(1|chamber)+(1|year),
          data=dat[dat$year>=2019 & dat$NON2O_stat.lab=='Y',]);fm3

summary(fm3)

tab_model(fm3)



###gasN
fm5<-lmer(scale(gas_N)~ftreatment*scale(year)+(1|chamber)+(1|year),
          data=dat[dat$year>=2019 & dat$NON2O_stat.lab=='Y',]);fm5

summary(fm5)

tab_model(fm5)


### NO
coefficients(summary(fm1))[ , "Estimate"]
presult<-car::Anova(fm1,type=2);presult
### N2O
coefficients(summary(fm2))[ , "Estimate"]
presult<-car::Anova(fm2,type=2);presult
### NO+N2O
coefficients(summary(fm3))[ , "Estimate"]
presult<-car::Anova(fm3,type=2);presult


############### <Main_Fig 4> ############

#######作图
unique(dat_fig0.t$treatment)

my_comparisons<-list(c('warmed','model'),c('control','warmed'))

dat_fig0$year<-as.factor(dat_fig0$year)

dat_fig0.t$year<-as.factor(dat_fig0.t$year)

{
  p1<-ggplot(dat_fig0.t,aes(year,NO.m.m))+
    geom_boxplot(data=dat_fig0,aes(year,NO.m,fill=treatment),width=0.5,alpha=0.7,size=0.2,outlier.size = 1,
                 position = position_dodge(0.8))+
    geom_point(data=dat_fig0[dat_fig0$NO_stat.lab=='Y',],aes(year,NO.m,group=treatment,shape=NO_stat.lab,alpha=NO_stat.lab),show.legend = F,
               position = position_jitterdodge(0.035, dodge.width = .8),size=1)+
    stat_summary(aes(group = treatment),fun = mean, geom = "point",
                 position =position_dodge(0.8), col = "black",size=1.2) +
    scale_shape_manual(limits=c('Y','N'),values=c(1,19))+
    scale_alpha_manual(limits=c('Y','N'),values=c(0.4,1))+
    geom_vline(xintercept = 1.5,lty=1,col='black',alpha=0.7,linewidth=0.45)+
    geom_text(data=dat_fig0.t[dat_fig0.t$treatment=='control' & dat_fig0.t$year=='2019',],col='black',
              aes(y=42,x='2020'),label=expression(paste('Treat: ', italic('P'),'= 0.011')),position=position_dodge(0.8),size=3.5,hjust=-.5)+
    geom_text(data=dat_fig0.t[dat_fig0.t$treatment=='control' & dat_fig0.t$year=='2019',],col='black',
              aes(y=38,x='2020'),label=expression(paste('Year: ', italic('P'),'< 0.001')),position=position_dodge(0.8),size=3.5,hjust=-.53)+
    geom_text(data=dat_fig0.t[dat_fig0.t$treatment=='control' & dat_fig0.t$year=='2019',],col='black',
              aes(y=34,x='2020'),label='Treat x Year:  ns',position=position_dodge(0.8),size=3.5,hjust=-.5)+
    scale_fill_manual(limits=c("control","warmed","model"),
                      values=c('royalblue1','red2','goldenrod'),name=NULL,
                      labels=c(expression(paste("Control")),
                               expression(paste("Warmed")),expression(paste("Q"[10]," expected"))))+
    scale_x_discrete(limits=c('2018','2019','2020','2021','2022','2023'),
                     labels=c(2018:2023),expand=c(0.05,0.05))+
    scale_y_continuous(limits=c(0,45),breaks=seq(0,45,15),minor_breaks = seq(0,45,5),expand=c(0,0))+
    guides(y="prism_offset_minor",fill = guide_legend(byrow = TRUE))+
    annotate("segment", y = 15, yend = 15, x = 1.46, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 30, yend = 30, x = 1.46, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 5, yend = 5, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 10, yend = 10, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 20, yend = 20, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 25, yend = 25, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 35, yend = 35, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 40, yend = 40, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    mythem+theme(legend.position=c(0.16,1.00), legend.justification=c(0,1),   #左上图例
                 # legend.key.size = unit(2, 'mm'),
                 legend.text = element_text(size=12),
                 legend.key.height = unit(0.6, 'cm'),
                 legend.key.width = unit(0.6, 'cm'))+
    labs(x='Year',y=expression(atop(paste('NO flux'),paste('(μg N m'^-2,' h'^-1,')'))));p1  
  
  
  stat.test <- dat_fig0  %>%
    filter(NO_stat.lab == 'Y')%>%
    group_by(year) %>%
    pairwise_t_test(NO.m ~ treatment,
                    p.adjust.method = "none") %>%    #两两T检验
    dplyr::select(year,group1,group2,n1,n2,p,p.signif) # Remove details
  stat.test
  
  stat.test <- stat.test %>% add_xy_position(x = "year")
  
  stat.test$g1_g2<-paste(stat.test$group1,stat.test$group2,sep='-')
  unique(stat.test$g1_g2)
  
  # t.test(dat_fig0$NO.m,)
  
  p1+ 
    stat_pvalue_manual(
      stat.test[stat.test$year!=2023&stat.test$p<0.050 & stat.test$g1_g2 %in% c("warmed-model"),],
      label = 'p.signif',
      y.position=25,label.size = 4.2,step.group.by='year',
      step.increase = 0.12,hide.ns = T, tip.length = 0.005)+
    stat_pvalue_manual(
      stat.test[stat.test$year!=2023&stat.test$p<0.056 & stat.test$g1_g2 %in% c("control-warmed"),],
      label = 'P = {round(p,3)}',
      y.position=20,label.size = 3.5,step.group.by='year',
      step.increase = 0.12,hide.ns = F, tip.length = 0.005)+
    stat_pvalue_manual(
      stat.test[stat.test$year==2023&stat.test$p<0.050 & stat.test$g1_g2%in% c("control-warmed","warmed-model"),], 
      label =  'p.signif',
      y.position=37,label.size = 4.2,prase=T,
      step.increase = 0.12,hide.ns = F, tip.length = 0.005)->p1;p1
  
  str(dat_fig0)
  dat_fig0$labs<-paste(dat_fig0$year,dat_fig0$treatment,sep='-')
  ggplot(dat_fig0,aes(year,NO.m,fill=treatment))+
    geom_point(data=dat_fig0,aes(year,NO.m,fill=treatment),show.legend = F,
               position = position_jitterdodge(0.035, dodge.width = .8),size=1,pch=1,alpha=0.3)+
    geom_vline(xintercept = 2018.5,lty=2,col='grey')+
    geom_boxplot(data=dat_fig0,aes(year,NO.m,fill=treatment,group=labs),position=position_dodge(0.8),outlier.colour ="black",width=0.65,
                 outlier.shape = 21,outlier.size = 0.7,alpha=.35)
  
  ggplot(dat_fig0,aes(year,N2O.m,fill=treatment))+
    geom_point(data=dat_fig0,aes(year,N2O.m,fill=treatment),show.legend = F,
               position = position_jitterdodge(0.035, dodge.width = .8),size=1,pch=1,alpha=0.3)+
    geom_vline(xintercept = 2018.5,lty=2,col='grey')+
    geom_boxplot(data=dat_fig0,aes(year,N2O.m,fill=treatment,group=labs),position=position_dodge(0.8),outlier.colour ="black",width=0.65,
                 outlier.shape = 21,outlier.size = 0.7,alpha=.35)
  
  ggplot(dat_fig0,aes(year,N2_field.m,fill=treatment))+
    geom_point(data=dat_fig0,aes(year,N2_field.m,fill=treatment),show.legend = F,
               position = position_jitterdodge(0.035, dodge.width = .8),size=1,pch=1,alpha=0.3)+
    geom_vline(xintercept = 2018.5,lty=2,col='grey')+
    geom_boxplot(data=dat_fig0,aes(year,N2_field.m,fill=treatment,group=labs),position=position_dodge(0.8),outlier.colour ="black",width=0.65,
                 outlier.shape = 21,outlier.size = 0.7,alpha=.35)
  
  
  p2<-ggplot(dat_fig0.t,aes(year,N2O.m.m))+
    geom_boxplot(data=dat_fig0,aes(x=year,y=N2O.m,fill=treatment),width=0.5,alpha=0.7,size=0.2,outlier.size = 1,
                 position = position_dodge(0.8))+
    geom_point(data=dat_fig0[dat_fig0$N2O_stat.lab=='Y',],
               aes(year,N2O.m,group=treatment,shape=N2O_stat.lab,alpha=N2O_stat.lab),show.legend = F,
               position = position_jitterdodge(0.035, dodge.width = .8),size=1)+
    stat_summary(aes(group = treatment),fun = mean, geom = "point",
                 position =position_dodge(0.8), col = "black",size=1.2) +
    scale_shape_manual(limits=c('Y','N'),values=c(1,19))+
    scale_alpha_manual(limits=c('Y','N'),values=c(0.4,1))+
    geom_vline(xintercept = 1.5,lty=1,col='black',alpha=0.7,linewidth=0.45)+
    geom_text(data=dat_fig0.t[dat_fig0.t$treatment=='control' & dat_fig0.t$year=='2019',],col='black',
              aes(y=56,x='2020'),label=expression(paste('Treat: ', italic('P'),'= 0.005')),position=position_dodge(0.8),size=3.5,hjust=-.5)+
    geom_text(data=dat_fig0.t[dat_fig0.t$treatment=='control' & dat_fig0.t$year=='2019',],col='black',
              aes(y=50,x='2020'),label=expression(paste('Year: ', italic('P'),'< 0.001')),position=position_dodge(0.8),size=3.5,hjust=-.53)+
    geom_text(data=dat_fig0.t[dat_fig0.t$treatment=='control' & dat_fig0.t$year=='2019',],col='black',
              aes(y=44,x='2020'),label='Treat x Year:  ns',position=position_dodge(0.8),size=3.5,hjust=-.5)+
    scale_fill_manual(limits=c("control","warmed","model"),
                      values=c('royalblue1','red2','goldenrod'),name=NULL,
                      labels=c(expression(paste("Control")),
                               expression(paste("Warmed")),expression(paste("Q"[10]," expected"))))+
    scale_x_discrete(limits=c('2018','2019','2020','2021','2022','2023'),
                     labels=c(2018:2023),expand=c(0.05,0.05))+
    scale_y_continuous(limits=c(0,60),breaks=seq(0,60,20),minor_breaks = seq(0,60,5),expand=c(0,0))+
    guides(y="prism_offset_minor",fill = guide_legend(byrow = TRUE))+
    annotate("segment", y = 20, yend = 20, x = 1.46, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 40, yend = 40, x = 1.46, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 5, yend = 5, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 10, yend = 10, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 15, yend = 15, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 25, yend = 25, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 30, yend = 30, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 35, yend = 35, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 45, yend = 45, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 50, yend = 50, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 55, yend = 55, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    mythem+theme(legend.position='none',   
                 # legend.key.size = unit(2, 'mm'),
                 legend.key.height = unit(0.7, 'cm'),
                 legend.key.width = unit(0.7, 'cm'))+
    labs(x='Year',y=expression(atop(paste('N'[2],'O flux'),paste('(μg N m'^-2~'h'^-1,')'))));p2 
  
  stat.test <- dat_fig0  %>%
    filter(N2O_stat.lab == 'Y')%>%
    group_by(year) %>%
    pairwise_t_test(N2O.m ~ treatment,
                    p.adjust.method = "none") %>%    #两两T检验
    dplyr::select(year,group1,group2,n1,n2,p,p.signif) # Remove details
  stat.test
  
  stat.test <- stat.test %>% add_xy_position(x = "year")
  
  stat.test$g1_g2<-paste(stat.test$group1,stat.test$group2,sep='-')
  
  p2+ 
    stat_pvalue_manual(
      stat.test[stat.test$year %notin% c(2022,2023) & stat.test$p<0.1 & stat.test$g1_g2 %in% c("control-warmed","warmed-model"),], 
      label = 'p.signif',
      y.position=26,label.size = 4.2,step.group.by='year',
      step.increase = 0.12,hide.ns = F, tip.length = 0.005)+
    stat_pvalue_manual(
      stat.test[stat.test$year %in% c(2022,2023) &stat.test$p<0.1 & stat.test$g1_g2 %in% c("control-warmed","warmed-model"),], 
      label = 'p.signif',
      y.position=50,label.size = 4.2,prase=T,
      step.increase = 0.12,hide.ns = F, tip.length = 0.005)->p2;p2
  
  p3<-ggplot(dat_fig0.t,aes(year,N2_field.m.m))+
    # geom_bar(aes(fill=treatment), position=position_dodge(0.8),stat="identity",size=0.17,col="black",alpha=.62,width=0.65)+
    # geom_errorbar(aes(ymin=N2_field.m.m-N2_field.m.se,ymax=N2_field.m.m+N2_field.m.se),
    #               position=position_dodge(0.8),stat="identity",width=0.15,size=0.5)+
    geom_boxplot(data=dat_fig0,aes(year,N2_field.m,fill=treatment),width=0.5,alpha=0.7,size=0.2,outlier.size = 1,
                 position = position_dodge(0.8))+
    geom_point(data=dat_fig0[dat_fig0$N2_stat.lab=='Y',],
               aes(year,N2_field.m,group=treatment,shape=N2_stat.lab,alpha=N2_stat.lab),show.legend = F,
               position = position_jitterdodge(0.035, dodge.width = .8),size=1)+
    stat_summary(aes(group = treatment),fun = mean, geom = "point",
                 position =position_dodge(0.8), col = "black",size=1.2) +
    scale_shape_manual(limits=c('Y','N'),values=c(1,19))+
    scale_alpha_manual(limits=c('Y','N'),values=c(0.4,1))+
    geom_vline(xintercept = 1.5,lty=1,col='black',alpha=0.7,linewidth=0.45)+
    # facet_wrap(.~year,nrow=1)+scale_y_continuous(limits=c(0,2.5),breaks=seq(0,2.4,0.6),expand=c(0,0))+
    # geom_text(aes(y=ifelse(N2_field.m.m>50,N2_field.m.m+18,N2_field.m.m+10),label=group_N2),
    #           position=position_dodge(0.8),size=5,hjust=0.5)+
    geom_text(data=dat_fig0.t[dat_fig0.t$treatment=='control' & dat_fig0.t$year=='2019',],col='black',
              aes(y=280,x='2020'),label='Treat:  ns',
              position=position_dodge(0.8),size=3.5,hjust=-1)+
    geom_text(data=dat_fig0.t[dat_fig0.t$treatment=='control' & dat_fig0.t$year=='2019',],col='black',
              aes(y=250,x='2020'),label=expression(paste('Year: ', italic('P'),'< 0.001')),
              position=position_dodge(0.8),size=3.5,hjust=-.65)+
    geom_text(data=dat_fig0.t[dat_fig0.t$treatment=='control' & dat_fig0.t$year=='2019',],col='black',
              aes(y=220,x='2020'),label=expression(paste('Treat x Year: ', italic('P'),'< 0.001')),
              position=position_dodge(0.8),size=3.5,hjust=-.43)+
    scale_fill_manual(limits=c("control","warmed","model"),
                      values=c('royalblue1','red2','goldenrod'),name=NULL,
                      labels=c(expression(paste("Control")),
                               expression(paste("Warmed")),expression(paste("Q"[10]," expected"))))+
    scale_x_discrete(limits=c('2018','2019','2020','2021','2022','2023'),
                     labels=c(2018:2023),expand=c(0.05,0.05))+
    scale_y_continuous(limits=c(0,300),breaks=seq(0,300,100),minor_breaks = seq(0,300,20),expand=c(0,0))+
    guides(y="prism_offset_minor",fill = guide_legend(byrow = TRUE))+
    annotate("segment", y = 100, yend = 100, x = 1.46, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 200, yend = 200, x = 1.46, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 20, yend = 20, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 40, yend = 40, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 60, yend = 60, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 80, yend = 80, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 120, yend = 120, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 140, yend = 140, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 160, yend = 160, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 180, yend = 180, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 220, yend = 220, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 240, yend = 240, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 260, yend = 260, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 280, yend = 280, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    mythem+theme(legend.position='none',
                 # legend.key.size = unit(2, 'mm'),
                 legend.key.height = unit(0.7, 'cm'),
                 legend.key.width = unit(0.7, 'cm'))+
    labs(x='Year',y=expression(atop(paste('N'[2],' flux'),paste('(μg N m'^-2~'h'^-1,')'))));p3
  
  
  stat.test <- dat_fig0  %>%
    filter(N2_stat.lab == 'Y')%>%
    group_by(year) %>%
    pairwise_t_test(N2_field.m ~ treatment,
                    p.adjust.method = "none") %>%    #两两T检验
    dplyr::select(year,group1,group2,n1,n2,p,p.signif) # Remove details
  stat.test
  
  stat.test <- stat.test %>% add_xy_position(x = "year")
  
  stat.test$g1_g2<-paste(stat.test$group1,stat.test$group2,sep='-')
  
  p3+ 
    stat_pvalue_manual(
      stat.test[stat.test$year == 2019 & stat.test$p<0.06 & stat.test$g1_g2 %in% c("control-warmed","warmed-model"),], 
      label = 'p.signif',
      y.position=135,label.size = 4.2,step.group.by='year',
      step.increase = 0.12,hide.ns = F, tip.length = 0.005)+
    stat_pvalue_manual(
      stat.test[stat.test$year!= 2019 &stat.test$p<0.06 & stat.test$g1_g2 %in% c("control-warmed","warmed-model"),],
      label = 'p.signif',
      y.position=70,label.size = 4.2,step.group.by='year',
      step.increase = 0.12,hide.ns = F, tip.length = 0.005)->p3;p3
  
  
  str(dat_fig0.t)
  unique(dat_fig0.t$year)
  str(dat_fig0)
  p4<-ggplot(dat_fig0.t,aes(year,gas_N.m))+
    geom_boxplot(data=dat_fig0,aes(year,gas_N,fill=treatment),width=0.5,alpha=0.7,size=0.2,outlier.size = 1,
                 position = position_dodge(0.8))+
    geom_point(data=dat_fig0[dat_fig0$gasN_stat.lab=='Y',],
               aes(year,gas_N,group=treatment,shape=gasN_stat.lab,alpha=gasN_stat.lab),show.legend = F,
               position = position_jitterdodge(0.035, dodge.width = .8),size=1)+
    stat_summary(aes(group = treatment),fun = mean, geom = "point",
                 position =position_dodge(0.8), col = "black",size=1.2) +
    scale_shape_manual(limits=c('Y','N'),values=c(1,19))+
    scale_alpha_manual(limits=c('Y','N'),values=c(0.4,1))+
    geom_vline(xintercept = 1.5,lty=1,col='black',alpha=0.7,linewidth=0.45)+
    geom_text(data=dat_fig0.t[dat_fig0.t$treatment=='control' & dat_fig0.t$year=='2019',],col='black',
              aes(y=280,x='2020'),label=expression(paste('Treat: ', italic('P'),'= 0.0104')),
              position=position_dodge(0.8),size=3.5,hjust=-0.57)+
    geom_text(data=dat_fig0.t[dat_fig0.t$treatment=='control' & dat_fig0.t$year=='2019',],col='black',
              aes(y=250,x='2020'),label=expression(paste('Year:  ',italic('P'),'< 0.001')),
              position=position_dodge(0.8),size=3.5,hjust=-0.65)+
    geom_text(data=dat_fig0.t[dat_fig0.t$treatment=='control' & dat_fig0.t$year=='2019',],col='black',
              aes(y=220,x='2020'),label=expression(paste('Treat x Year:  ',italic('P'),'= 0.003')),
              position=position_dodge(0.8),size=3.5,hjust=-0.45)+
    
    scale_fill_manual(limits=c("control","warmed","model"),
                      values=c('royalblue1','red2','goldenrod'),name=NULL,
                      labels=c(expression(paste("Control")),
                               expression(paste("Warmed")),expression(paste("Q"[10]," expected"))))+
    scale_x_discrete(limits=c('2018','2019','2020','2021','2022','2023'),
                     labels=c(2018:2023),expand=c(0.05,0.05))+
    scale_y_continuous(limits=c(0,300),breaks=seq(0,300,100),minor_breaks = seq(0,300,20),expand=c(0,0))+
    guides(y="prism_offset_minor",fill = guide_legend(byrow = TRUE))+
    annotate("segment", y = 100, yend = 100, x = 1.46, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 200, yend = 200, x = 1.46, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 20, yend = 20, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 40, yend = 40, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 60, yend = 60, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 80, yend = 80, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 120, yend = 120, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 140, yend = 140, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 160, yend = 160, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 180, yend = 180, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 220, yend = 220, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 240, yend = 240, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 260, yend = 260, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 280, yend = 280, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    mythem+theme(legend.position='none',
                 # legend.key.size = unit(2, 'mm'),
                 legend.key.height = unit(0.7, 'cm'),
                 legend.key.width = unit(0.7, 'cm'))+
    labs(x='Year',y=expression(atop(paste('Gaseous N flux'),paste('(μg N m'^-2~'h'^-1,')'))));p4
  
  
  stat.test <- dat_fig0  %>%
    filter(gasN_stat.lab == 'Y')%>%
    group_by(year) %>%
    pairwise_t_test(gas_N ~ treatment,
                    p.adjust.method = "none") %>%    #两两T检验
    dplyr::select(year,group1,group2,n1,n2,p,p.signif) # Remove details
  stat.test
  
  stat.test <- stat.test %>% add_xy_position(x = "year")
  
  stat.test$g1_g2<-paste(stat.test$group1,stat.test$group2,sep='-')
  
  p4+ 
    stat_pvalue_manual(
      stat.test[stat.test$year == 2019 & stat.test$p<0.06 & stat.test$g1_g2 %in% c("control-warmed","warmed-model"),], 
      label = 'p.signif',
      y.position=190,label.size = 4.2,step.group.by='year',
      step.increase = 0.12,hide.ns = F, tip.length = 0.005)+
    stat_pvalue_manual(
      stat.test[stat.test$year!= 2019 &stat.test$p<0.06 & stat.test$g1_g2 %in% c("control-warmed","warmed-model"),], 
      label = 'p.signif',
      y.position=175,label.size = 4.2,step.group.by='year',
      step.increase = 0.12,hide.ns = F, tip.length = 0.005)->p4;p4
  
  str(dat_fig0)
  p5<-ggplot(dat_fig0.t,aes(year,NON2O.m.m))+
    # geom_bar(aes(fill=treatment),position=position_dodge(0.8),stat="identity",size=0.17,col="black",alpha=.62,width=0.65)+
    # geom_errorbar(aes(ymin=NON2O.m.m-NON2O.m.se,ymax=NON2O.m.m+NON2O.m.se),
    #               position=position_dodge(0.8),stat="identity",width=0.15,size=0.5)+
    geom_boxplot(data=dat_fig0,aes(year,NON2O.m,fill=treatment),width=0.5,alpha=0.7,size=0.2,outlier.size = 1,
                 position = position_dodge(0.8))+
    geom_point(data=dat_fig0[dat_fig0$NON2O_stat.lab=='Y',],
               aes(year,NON2O.m,group=treatment,shape=NON2O_stat.lab,alpha=NON2O_stat.lab),show.legend = F,
               position = position_jitterdodge(0.035, dodge.width = .8),size=1)+
    stat_summary(aes(group = treatment),fun = mean, geom = "point",
                 position =position_dodge(0.8), col = "black",size=1.2) +
    scale_shape_manual(limits=c('Y','N'),values=c(1,19))+
    scale_alpha_manual(limits=c('Y','N'),values=c(0.4,1))+
    geom_vline(xintercept = 1.5,lty=1,col='black',alpha=0.7,linewidth=0.45)+
    # facet_wrap(.~year,nrow=1)+scale_y_continuous(limits=c(0,2.5),breaks=seq(0,2.4,0.6),expand=c(0,0))+
    # geom_text(aes(y=ifelse(N2O.m.m>15,N2O.m.m+0.25*N2O.m.m,N2O.m.m+5),label=group_N2O),position=position_dodge(0.8),size=5,hjust=0.5)+
    geom_text(data=dat_fig0.t[dat_fig0.t$treatment=='control' & dat_fig0.t$year=='2019',],col='black',
              aes(y=85,x='2020'),label=expression(paste('Treat: ', italic('P'),'< 0.001')),
              position=position_dodge(0.8),size=3.5,hjust=-.5)+
    geom_text(data=dat_fig0.t[dat_fig0.t$treatment=='control' & dat_fig0.t$year=='2019',],col='black',
              aes(y=77,x='2020'),label=expression(paste('Year: ', italic('P'),'= 0.001')),
              position=position_dodge(0.8),size=3.5,hjust=-.53)+
    geom_text(data=dat_fig0.t[dat_fig0.t$treatment=='control' & dat_fig0.t$year=='2019',],col='black',
              aes(y=69,x='2020'),label='Treat x Year:  ns',position=position_dodge(0.8),size=3.5,hjust=-.5)+
    scale_fill_manual(limits=c("control","warmed","model"),
                      values=c('royalblue1','red2','goldenrod'),name=NULL,
                      labels=c(expression(paste("Control")),
                               expression(paste("Warmed")),expression(paste("Q"[10]," expected"))))+
    scale_x_discrete(limits=c('2018','2019','2020','2021','2022','2023'),
                     labels=c(2018:2023),expand=c(0.05,0.05))+
    scale_y_continuous(limits=c(0,90),breaks=seq(0,90,30),minor_breaks = seq(0,90,10),expand=c(0,0))+
    guides(y="prism_offset_minor",fill = guide_legend(byrow = TRUE))+
    annotate("segment", y = 30, yend = 30, x = 1.46, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 60, yend = 60, x = 1.46, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 10, yend = 10, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 20, yend = 20, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 40, yend = 40, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 50, yend = 50, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 70, yend = 70, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    annotate("segment", y = 80, yend = 80, x = 1.475, xend = 1.5, colour = "black")+#创作视觉上的副本
    mythem+theme(legend.position='none',  
                 # legend.key.size = unit(2, 'mm'),
                 legend.key.height = unit(0.7, 'cm'),
                 legend.key.width = unit(0.7, 'cm'))+
    labs(x='Year',y=expression(atop(paste('NO + N'[2],'O flux'),paste('(μg N m'^-2~'h'^-1,')'))));p5 
  
  stat.test <- dat_fig0  %>%
    filter(NON2O_stat.lab == 'Y')%>%
    group_by(year) %>%
    pairwise_t_test(NON2O.m ~ treatment,
                    p.adjust.method = "none") %>%    #两两T检验
    dplyr::select(year,group1,group2,n1,n2,p,p.signif) # Remove details
  stat.test
  
  stat.test <- stat.test %>% add_xy_position(x = "year")
  
  stat.test$g1_g2<-paste(stat.test$group1,stat.test$group2,sep='-')
  
  p5+ 
    stat_pvalue_manual(
      stat.test[stat.test$year %notin% c(2022,2023) & stat.test$p<0.1 & stat.test$g1_g2 %in% c("control-warmed","warmed-model"),], 
      label = 'p.signif',
      y.position=40,label.size = 4.2,step.group.by='year',
      step.increase = 0.12,hide.ns = F, tip.length = 0.005)+
    stat_pvalue_manual(
      stat.test[stat.test$year %in% c(2022,2023) &stat.test$p<0.1 & stat.test$g1_g2 %in% c("control-warmed","warmed-model"),], 
      label = 'p.signif',
      y.position=68,label.size = 4.2,prase=T,
      step.increase = 0.08,hide.ns = F, tip.length = 0.005)->p5;p5
  
}

############## Normalization multi-year flux
str(dat_fig0)
unique(dat_fig0$chamber.treat)


datai.m <- dat_fig0[dat_fig0$year %in% seq(2019,2023),]

dat_fig.m<-data.frame()

k=1
for(i in unique(datai.m$chamber.treat)){
  dat_fig.m[k,1]=i
  dat_fig.m[k,2]=mean(datai.m$NO.m[datai.m$chamber.treat==i])
  dat_fig.m[k,3]=mean(datai.m$N2O.m[datai.m$chamber.treat==i])
  dat_fig.m[k,4]=mean(datai.m$NON2O.m[datai.m$chamber.treat==i])
  dat_fig.m[k,5]=mean(datai.m$N2_field.m[datai.m$chamber.treat==i])
  dat_fig.m[k,6]=mean(datai.m$gas_N[datai.m$chamber.treat==i])
  dat_fig.m[k,7]=unique(datai.m$chamber[datai.m$chamber.treat==i])
  dat_fig.m[k,8]=unique(datai.m$treatment[datai.m$chamber.treat==i])
  dat_fig.m[k,9]=unique(datai.m$plot[datai.m$chamber.treat==i])
  k=k+1
  
}

names(dat_fig.m)<-c('chamber.treat','NO.m','N2O.m','NON2O.m','N2.m','gasN.m','chamber','treatment','plot')

unique(dat_fig.m$treatment)

str(dat_fig.m)
###识别异常值
dat_fig.m$NO_stat.lab='Y'
dat_fig.m$N2O_stat.lab='Y'
dat_fig.m$NON2O_stat.lab='Y'
dat_fig.m$N2_stat.lab='Y'
dat_fig.m$gasN_stat.lab='Y'

for(i in unique(dat_fig.m$treatment)){
  dfa<- dat_fig.m %>%
    filter(treatment == i)
  x1c=dfa$chamber.treat[dfa$NO.m %in% boxplot(dfa$NO.m)$out]
  x2c=dfa$chamber.treat[dfa$N2O.m %in% boxplot(dfa$N2O.m)$out]
  x3c=dfa$chamber.treat[dfa$N2.m %in% boxplot(dfa$N2.m)$out]
  x4c=dfa$chamber.treat[dfa$gasN.m %in% boxplot(dfa$gasN.m)$out]
  x5c=dfa$chamber.treat[dfa$NON2O.m %in% boxplot(dfa$NON2O.m)$out]
  dat_fig.m$NO_stat.lab[dat_fig.m$chamber.treat %in% x1c]='N'
  dat_fig.m$N2O_stat.lab[dat_fig.m$chamber.treat %in% x2c]='N'
  dat_fig.m$N2_stat.lab[dat_fig.m$chamber.treat %in% x3c]='N'
  dat_fig.m$gasN_stat.lab[dat_fig.m$chamber.treat %in% x4c]='N'
  dat_fig.m$NON2O_stat.lab[dat_fig.m$chamber.treat %in% x5c]='N'
}


dat_fig.m.t<-summaryBy(data=dat_fig.m,FUN = myfun,
                       dat_fig.m[,c('NO.m','N2O.m','NON2O.m','N2.m','gasN.m')]~treatment)

str(dat_fig.m.t)
(dat_fig.m.t$NO.m.m[dat_fig.m.t$treatment=='warmed']-dat_fig.m.t$NO.m.m[dat_fig.m.t$treatment=='control'])/
  dat_fig.m.t$NO.m.m[dat_fig.m.t$treatment=='control']

(dat_fig.m.t$N2O.m.m[dat_fig.m.t$treatment=='warmed']-dat_fig.m.t$N2O.m.m[dat_fig.m.t$treatment=='control'])/
  dat_fig.m.t$N2O.m.m[dat_fig.m.t$treatment=='control']


(dat_fig.m.t$gasN.m.m[dat_fig.m.t$treatment=='warmed']-dat_fig.m.t$gasN.m.m[dat_fig.m.t$treatment=='control'])/
  dat_fig.m.t$gasN.m.m[dat_fig.m.t$treatment=='control']


my_comparisons<-list(c('warmed','model'),c('control','warmed'))

p1m<-ggplot(data=dat_fig.m,aes(treatment,NO.m,fill=treatment),)+#geom_point(pch=1)+
  geom_boxplot(width=0.4,alpha=0.7,size=0.2,outlier.size = 1,
               position = position_dodge(0.8))+
  geom_point(data=dat_fig.m[dat_fig.m$NO_stat.lab=='Y',],aes(treatment,NO.m,fill=treatment),
             position = position_jitterdodge(0.035, dodge.width = .8),size=1,pch=1,alpha=0.3)+
  stat_summary(data=dat_fig.m[dat_fig.m$NO_stat.lab=='Y',], fun = mean, geom = "point", fill = "black",size=1.2) +
  scale_y_continuous(breaks=seq(0,30,10),expand = c(0,0))+
  scale_fill_manual(limits=c('control','model','warmed'),                      
                    values=c('royalblue1','goldenrod','red2'),name=NULL)+
  stat_compare_means(data=dat_fig.m[dat_fig.m$NO_stat.lab=='Y',],aes(x=treatment,y=NO.m,
                                     label = paste0("p = ", after_stat(p.format))),
                     comparisons=my_comparisons,label.y = c(26.2, 21.7),
                     method="t.test",vjust=0.5,hide.ns = T,
                     label="p.signif",tip.length = 0.01,
                     size=4.2)+
  labs(x=NULL,y=expression(atop(paste('Normalized NO flux'),paste('(μg N m'^-2,' h'^-1,')'))))+
  scale_x_discrete(limits=c('control','warmed','model'),labels=NULL)+
  annotate("segment", y = -0.75, yend = 0, x = 2, xend = 2, colour = "black")+
  coord_cartesian(ylim = c(0,30),clip = 'off')+  guides(y="prism_offset_minor")+
  mythem+theme(axis.ticks.length.x = unit(0,'cm'),
               panel.background = element_rect(fill =NA, linewidth=0.4,colour = "black", linetype = "solid",inherit.blank=T),
               plot.margin = unit(c(0.25,0,0.15,0.5),'cm'));p1m


p2m<-ggplot(data=dat_fig.m,aes(treatment,N2O.m,fill=treatment),)+#geom_point(pch=1)+
  geom_boxplot(width=0.4,alpha=0.7,size=0.2,outlier.size = 1,
               position = position_dodge(0.8))+
  geom_point(data=dat_fig.m[dat_fig.m$NO_stat.lab=='Y',],aes(treatment,N2O.m,fill=treatment),
             position = position_jitterdodge(0.035, dodge.width = .8),size=1,pch=1,alpha=0.3)+
  stat_summary(data=dat_fig.m[dat_fig.m$N2O_stat.lab=='Y',], fun = mean, geom = "point", fill = "black",size=1.2) +
  scale_y_continuous(breaks=seq(0,30,10),expand = c(0,0))+
  scale_fill_manual(limits=c('control','model','warmed'),                      
                    values=c('royalblue1','goldenrod','red2'),name=NULL)+
  stat_compare_means(data=dat_fig.m[dat_fig.m$N2O_stat.lab=='Y',],
                     aes(x=treatment,y=N2O.m,label = paste0("p = ", after_stat(p.format))),
                     comparisons=my_comparisons,label.y = c(26.2,21.7),
                     method="t.test",vjust=0.5,hide.ns = T,
                     label="p.signif",tip.length = 0.01,
                     size=4.2)+
  labs(x=NULL,y=expression(atop(paste('Normalized N'[2],'O flux'),paste('(μg N m'^-2,' h'^-1,')'))))+
  scale_x_discrete(limits=c('control','warmed','model'),labels=NULL)+
  annotate("segment", y = -0.75, yend = 0, x = 2, xend = 2, colour = "black")+
  coord_cartesian(ylim = c(0,30),clip = 'off')+  guides(y="prism_offset_minor")+
  mythem+theme(axis.ticks.length.x = unit(0,'cm'),
               panel.background = element_rect(fill =NA, linewidth=0.4,colour = "black", linetype = "solid",inherit.blank=T),
               plot.margin = unit(c(0.25,0,0.15,0.5),'cm'));p2m

str(dat_fig.m)
p3m<-ggplot(data=dat_fig.m,aes(treatment,N2.m,fill=treatment),)+#geom_point(pch=1)+
  geom_boxplot(width=0.4,alpha=0.7,size=0.2,outlier.size = 1,
               position = position_dodge(0.8))+
  geom_point(data=dat_fig.m[dat_fig.m$N2_stat.lab=='Y',],aes(treatment,N2.m,fill=treatment),
             position = position_jitterdodge(0.035, dodge.width = .8),size=1,pch=1,alpha=0.3)+
  stat_summary(data=dat_fig.m[dat_fig.m$N2_stat.lab=='Y',], fun = mean, geom = "point", fill = "black",size=1.2) +
  scale_y_continuous(breaks=seq(0,90,30),expand = c(0,0))+
  scale_fill_manual(limits=c('control','model','warmed'),                      
                    values=c('royalblue1','goldenrod','red2'),name=NULL)+
  stat_compare_means(data=dat_fig.m[dat_fig.m$N2_stat.lab=='Y',],
                     aes(x=treatment,y=N2.m,label = paste0("p = ", after_stat(p.format))),
                     comparisons=my_comparisons,label.y = c(75,60),
                     method="t.test",vjust=0.5,hide.ns = T,
                     label="p.signif",tip.length = 0.01,
                     size=4.2)+
  labs(x=NULL,y=expression(atop(paste('Normalized  N'[2],' flux'),paste('(μg N m'^-2,' h'^-1,')'))))+
  scale_x_discrete(limits=c('control','warmed','model'),labels=NULL)+
  annotate("segment", y = -1.5, yend = 0, x = 2, xend = 2, colour = "black")+
  coord_cartesian(ylim = c(0,90),clip = 'off')+  guides(y="prism_offset_minor")+
  mythem+theme(axis.ticks.length.x = unit(0,'cm'),
               panel.background = element_rect(fill =NA, linewidth=0.4,colour = "black", linetype = "solid",inherit.blank=T),
               plot.margin = unit(c(0.25,0,0.65,0.5),'cm'));p3m

p4m<-ggplot(data=dat_fig.m,aes(treatment,gasN.m,fill=treatment),)+#geom_point(pch=1)+
  geom_boxplot(width=0.4,alpha=0.7,size=0.2,outlier.size = 1,
               position = position_dodge(0.8))+
  geom_point(data=dat_fig.m[dat_fig.m$gasN_stat.lab=='Y',],aes(treatment,gasN.m,fill=treatment),
             position = position_jitterdodge(0.035, dodge.width = .8),size=1,pch=1,alpha=0.3)+
  stat_summary(data=dat_fig.m[dat_fig.m$gasN_stat.lab=='Y',], fun = mean, geom = "point", fill = "black",size=1.2) +
  scale_y_continuous(breaks=seq(0,120,40),expand = c(0,0))+
  scale_fill_manual(limits=c('control','model','warmed'),                      
                    values=c('royalblue1','goldenrod','red2'),name=NULL)+
  stat_compare_means(data=dat_fig.m[dat_fig.m$gasN_stat.lab=='Y',],
                     aes(x=treatment,y=gasN.m,label = paste0("p = ", after_stat(p.format))),
                     comparisons=my_comparisons,label.y = c(98,84),
                     method="t.test",vjust=0.5,hide.ns = T,
                     label="p.signif",tip.length = 0.01,
                     size=4.2)+
  labs(x=NULL,y=expression(atop(paste('Normalized  gaseous N flux'),paste('(μg N m'^-2,' h'^-1,')'))))+
  scale_x_discrete(limits=c('control','warmed','model'),labels=NULL)+
  annotate("segment", y = -2, yend = 0, x = 2, xend = 2, colour = "black")+
  coord_cartesian(ylim = c(0,120),clip = 'off')+  guides(y="prism_offset_minor")+
  mythem+theme(axis.ticks.length.x = unit(0,'cm'),
               panel.background = element_rect(fill =NA, linewidth=0.4,colour = "black", linetype = "solid",inherit.blank=T),
               plot.margin = unit(c(0.25,0,0.65,0.5),'cm'))+
  annotate("text", label='2019-2023',y = -5.5,  x = 2, colour = "black",size=4.2);p4m

p5m<-ggplot(data=dat_fig.m,aes(treatment,NON2O.m,fill=treatment),)+#geom_point(pch=1)+
  geom_boxplot(width=0.4,alpha=0.7,size=0.2,outlier.size = 1,
               position = position_dodge(0.8))+
  geom_point(data=dat_fig.m[dat_fig.m$NO_stat.lab=='Y',],aes(treatment,NON2O.m,fill=treatment),
             position = position_jitterdodge(0.035, dodge.width = .8),size=1,pch=1,alpha=0.3)+
  stat_summary(data=dat_fig.m[dat_fig.m$NON2O_stat.lab=='Y',], fun = mean, geom = "point", fill = "black",size=1.2) +
  scale_y_continuous(breaks=seq(0,60,20),expand = c(0,0))+
  scale_fill_manual(limits=c('control','model','warmed'),                      
                    values=c('royalblue1','goldenrod','red2'),name=NULL)+
  stat_compare_means(data=dat_fig.m[dat_fig.m$NON2O_stat.lab=='Y',],
                     aes(x=treatment,y=NON2O.m,label = paste0("p = ", after_stat(p.format))),
                     comparisons=my_comparisons,label.y = c(53.2,45.7),
                     method="t.test",vjust=0.5,hide.ns = T,
                     label="p.signif",tip.length = 0.01,
                     size=4.2)+
  labs(x=NULL,y=expression(atop(paste('Normalized NO + N'[2],'O flux'),paste('(μg N m'^-2,' h'^-1,')'))))+
  scale_x_discrete(limits=c('control','warmed','model'),labels=NULL)+
  annotate("segment", y = -1.7, yend = 0, x = 2, xend = 2, colour = "black")+
  coord_cartesian(ylim = c(0,60),clip = 'off')+  guides(y="prism_offset_minor")+
  mythem+theme(axis.ticks.length.x = unit(0,'cm'),
               panel.background = element_rect(fill =NA, linewidth=0.4,colour = "black", linetype = "solid",inherit.blank=T),
               plot.margin = unit(c(0.25,0,0.65,0.5),'cm'))+
  annotate("text", label='2019-2023',y = -5,  x = 2, colour = "black",size=4.2);p5m

###累计排放柱状图  
{
  ###根据cum_c1最后一天的通量计算样方均值和se
  str(cum_c1)
  cum_c1$NO_stat.lab='Y'
  cum_c1$N2O_stat.lab='Y'
  cum_c1$N2_stat.lab='Y'
  cum_c1$gasN_stat.lab='Y'
  
  ###筛选异常值control\warmed\model
  for(i in unique(cum_c1$treatment)){
    
    dfa<- cum_c1 %>%
      filter(treatment == i)
    
    x1c=dfa$system_id[dfa$NO_N.m %in% boxplot(dfa$NO_N.m)$out]
    x2c=dfa$system_id[dfa$N2O_N.m %in% boxplot(dfa$N2O_N.m)$out]
    x3c=dfa$system_id[dfa$N2_N.m %in% boxplot(dfa$N2_N.m)$out]
    x4c=dfa$system_id[dfa$gasN %in% boxplot(dfa$gasN)$out]
    cum_c1$NO_stat.lab[cum_c1$system_id %in% x1c]='N'
    cum_c1$N2O_stat.lab[cum_c1$system_id %in% x2c]='N'
    cum_c1$N2_stat.lab[cum_c1$system_id %in% x3c]='N'
    cum_c1$gasN_stat.lab[cum_c1$system_id %in% x4c]='N'
  }
  
  dat_ex.t<-cum_c1  #将NA赋值给异常值,在求算处理水平均值
  dat_ex.t$NO_N.cum[dat_ex.t$NO_stat.lab=='N']=NA
  dat_ex.t$N2O_N.cum[dat_ex.t$N2O_stat.lab=='N']=NA
  dat_ex.t$N2_N.cum[dat_ex.t$N2_stat.lab=='N']=NA
  dat_ex.t$gasN[dat_ex.t$gasN_stat.lab=='N']=NA
  
  # cum_pL<-summaryBy(NO_N.cum+N2O_N.cum+N2_N.cum+gasN~plot+treatment,data=dat_ex.t,FUN=myfun)
  cum_pT<-summaryBy(NO_N.cum+N2O_N.cum+N2_N.cum+gasN~treatment,data=dat_ex.t,FUN=myfun)  #直接计算处理水平
  
  names(cum_pT)[c(2,4,5,7,8,10,11,13)]<-c(
    "NO_N.cum","NO_N.cum.se","N2O_N.cum","N2O_N.cum.se","N2_N.cum","N2_N.cum.se", 
    "gasN","gasN.se")
  
  ###多重比较control\warmed\model  20240323，，与T检验结果有出入。。。
  str(cum_pT)
  my_comparisons<-list(c('warmed','Q10_predicted'),c('control','warmed'))
  p1t<-ggplot(data=cum_pT,aes(treatment,NO_N.cum,fill=treatment),)+#geom_point(pch=1)+
    geom_bar(position=position_dodge(0.8),stat="identity",size=0.17,col="black",alpha=.62,width=0.65)+
    geom_errorbar(data=cum_pT,aes(treatment,NO_N.cum,
                                  ymin=NO_N.cum-NO_N.cum.se,
                                  ymax=NO_N.cum+NO_N.cum.se),
                  position=position_dodge(0.8),stat="identity",width=0.15,size=0.5)+
    geom_point(data=cum_c1,aes(treatment,NO_N.cum),
               position = position_jitterdodge(0.035, dodge.width = .8),size=1,pch=1,alpha=0.3)+
    scale_y_continuous(limits=c(0,8),breaks=seq(0,8,2),expand = c(0,0))+
    scale_fill_manual(limits=c('control','Q10_predicted','warmed'),                      
                      values=c('royalblue1','goldenrod','red2'),name=NULL)+
    stat_compare_means(data=cum_c1,aes(x=treatment,y=NO_N.cum,
                                       label = paste0("p = ", after_stat(p.format))),
                       comparisons=my_comparisons,label.y = c(6.1, 6.5, 7.2),
                       method="t.test",vjust=0.5,hide.ns = T,
                       label="p.signif",tip.length = 0.01,
                       size=8)+
    labs(x=NULL,y=expression(atop(paste('Cumulative NO'),paste('(kg N ha'^-1,')'))))+
    scale_x_discrete(limits=c('control','warmed','Q10_predicted'),labels=c('C','W','M'))+
    mythem;p1t
  
  
  p2t<-ggplot(cum_pT,aes(treatment,N2O_N.cum,fill=treatment))+#geom_point(pch=1)+
    geom_bar(position=position_dodge(0.8),stat="identity",size=0.17,col="black",alpha=.62,width=0.65)+
    geom_errorbar(aes(ymin=N2O_N.cum-N2O_N.cum.se,ymax=N2O_N.cum+N2O_N.cum.se),
                  position=position_dodge(0.8),stat="identity",width=0.15,size=0.5)+
    geom_point(data=cum_c1[cum_c1$N2O_N.cum<6,],aes(treatment,N2O_N.cum),              
               position = position_jitterdodge(0.035, dodge.width = .8),size=1,pch=1,alpha=0.3)+
    geom_point(data=cum_c1[cum_c1$N2O_N.cum>6 & cum_c1$treatment=='warmed',],aes(treatment,N2O_N.cum),
               position = position_jitterdodge(0.035, dodge.width = .8),size=1,pch=19,alpha=0.8)+
    stat_compare_means(data=cum_c1[-which(cum_c1$N2O_N.cum>6 & cum_c1$treatment=='warmed'),],
                       aes(x=treatment,y=N2O_N.cum,
                           label = paste0("p = ", after_stat(p.format))),
                       comparisons=my_comparisons,label.y = c(6.2, 6.9),
                       method="t.test",vjust=0.5,hide.ns = T,
                       label="p.signif",tip.length = 0.01,
                       size=8)+
    scale_y_continuous(limits=c(0,8),breaks=seq(0,8,2),expand = c(0,0))+
    scale_fill_manual(limits=c('control','Q10_predicted','warmed'),                      
                      values=c('royalblue1','goldenrod','red2'),name=NULL)+
    labs(x=NULL,y=expression(atop(paste('Cumulative N'[2],'O'),paste('(kg N ha'^-1,')'))))+
    scale_x_discrete(limits=c('control','warmed','Q10_predicted'),labels=c('C','W','M'))+
    mythem;p2t
  
  p3t<-ggplot(cum_pT,aes(treatment,N2_N.cum,fill=treatment))+#geom_point(pch=1)+
    geom_bar(position=position_dodge(0.8),stat="identity",size=0.17,col="black",alpha=.62,width=0.65)+
    geom_errorbar(aes(ymin=N2_N.cum-N2_N.cum.se,ymax=N2_N.cum+N2_N.cum.se),
                  position=position_dodge(0.8),stat="identity",width=0.15,size=0.5)+
    geom_point(data=cum_c1[cum_c1$N2_N.cum<16,],aes(treatment,N2_N.cum),
               position = position_jitterdodge(0.035, dodge.width = .8),size=1,pch=1,alpha=0.3)+
    geom_point(data=cum_c1[cum_c1$N2_N.cum>16,],aes(treatment,N2_N.cum),
               position = position_jitterdodge(0.035, dodge.width = .8),size=1,pch=19,alpha=0.8)+
    stat_compare_means(data=cum_c1[-which(cum_c1$N2_N.cum>16 & cum_c1$treatment=='warmed'),],
                       aes(x=treatment,y=N2_N.cum,
                           label = paste0("p = ", after_stat(p.format))),
                       comparisons=my_comparisons,label.y = c(13, 14, 15.5),
                       method="t.test",vjust=0.5,hide.ns = T,
                       label="p.signif",tip.length = 0.01,
                       size=8)+
    scale_y_continuous(limits=c(0,18),breaks=seq(0,18,6),expand = c(0,0))+
    scale_fill_manual(limits=c('control','Q10_predicted','warmed'),                      
                      values=c('royalblue1','goldenrod','red2'),name=NULL)+
    labs(x=NULL,y=expression(atop(paste('Cumulative N'[2],''),paste('(kg N ha'^-1,')'))))+
    scale_x_discrete(limits=c('control','warmed','Q10_predicted'),labels=c('C','W','M'))+
    mythem;p3t
  
  
  ###总气态氮释放
  str(cum_pT)
  str(cum_c1)

  
  pt<-ggplot(cum_pT,aes(treatment,gasN,fill=treatment))+
    geom_bar(position=position_dodge(0.8),stat="identity",size=0.17,col="black",alpha=.62,width=0.65)+
    geom_errorbar(aes(ymin=gasN-gasN.se,ymax=gasN+gasN.se),
                  position=position_dodge(0.8),stat="identity",width=0.15,size=0.5)+
    geom_point(data=cum_c1[cum_c1$gasN<=20,],aes(treatment,gasN),position = position_jitterdodge(0.035, dodge.width = .8),
               size=1,pch=1,alpha=0.3)+
    geom_point(data=cum_c1[cum_c1$gasN>20 & cum_c1$treatment=='warmed',],aes(treatment,gasN),position = position_jitterdodge(0.035, dodge.width = .8),
               size=1,pch=19,alpha=0.8)+
    stat_compare_means(data=cum_c1[-which(cum_c1$gasN>20 & cum_c1$treatment=='warmed'),],
                       aes(x=treatment,y=gasN,
                           label = paste0("p = ", after_stat(p.format))),
                       comparisons=my_comparisons,label.y = c(19.5, 21, 23.5),
                       method="t.test",vjust=0.5,hide.ns = T, na.rm = T,
                       label="p.signif",tip.length = 0.01,
                       size=8)+
    scale_y_continuous(limits=c(0,27),breaks=seq(0,27,9),expand = c(0,0))+
    scale_fill_manual(limits=c("control","Q10_predicted","warmed"),
                      values=c('royalblue1','goldenrod','red2'),name=NULL,
                      labels=c("Control","Predicted","Warmed"))+
    labs(x='Treatment',y=expression(atop(paste('Cumulative gaseous N'),paste('(kg N ha'^-1,')'))))+
    scale_x_discrete(limits=c('control','warmed','Q10_predicted'),labels=c('C','W','M'))+
    mythem;pt
  
  
  figure<-ggarrange(p1+rremove("xlab")+rremove("x.text")+
                      theme(plot.margin = unit(c(0.5,0,0.3,0.5),'cm'),
                            axis.text = element_text(size=12),axis.title = element_text(size=7)),
                    p1t+rremove("xlab")+rremove("x.text")+
                      theme(plot.margin = unit(c(0.5,0.5,0.3,0),'cm'),
                            axis.text = element_text(size=12),axis.title = element_text(size=7)),
                    p2+rremove("xlab")+rremove("x.text")+
                      theme(plot.margin = unit(c(0.25,0,0.3,0.5),'cm'),
                            axis.text = element_text(size=12),axis.title = element_text(size=7)),
                    p2t+rremove("xlab")+rremove("x.text")+
                      theme(plot.margin = unit(c(0.25,0.5,0.3,0),'cm'),
                            axis.text = element_text(size=12),axis.title = element_text(size=7)),
                    p3+rremove("xlab")+rremove("x.text")+
                      theme(plot.margin = unit(c(0.25,0,0.3,0.5),'cm'),
                            axis.text = element_text(size=12),axis.title = element_text(size=7)),
                    p3t+rremove("xlab")+rremove("x.text")+
                      theme(plot.margin = unit(c(0.25,0.5,0.3,0),'cm'),
                            axis.text = element_text(size=12),axis.title = element_text(size=7)),
                    
                    p4+theme(strip.text = element_blank(),strip.background = element_blank(),
                             axis.title = element_text(size=12),axis.text = element_text(size=7),
                             plot.margin = unit(c(0.25,0,0.15,0.5),'cm')),
                    pt+theme(strip.text = element_blank(),strip.background = element_blank(),
                             axis.title = element_text(size=12),axis.text = element_text(size=7),
                             plot.margin = unit(c(0.25,0.5,0.15,0),'cm')),
                    
                    labels = c("a", "b","c","d","e","f","g","h"),
                    label.x = c(0.035,0.14),label.y = c(0.98,0.98,0.98,0.98),
                    ncol = 2, nrow = 4,align = "v",   ##"v"竖直对齐
                    font.label = list(size = 20, color ="black"),
                    widths = c(16,4.3), heights = c(5,5,5,5.9),
                    common.legend=F);figure
  getwd()
  ggsave("Gaseous N response to warming.pdf", figure, width =16, height =12,
         device=cairo_pdf)
}


###2019-2023年平均释放速率
{
  ###根据每年的daily flux计算2019-2023年间均值和se（n = 5）
  str(dat_fig0)
  dat_fig0$E_NON2O<-dat_fig0$E_NO+dat_fig0$E_N2O
  dat_fig0$E_gasN<-dat_fig0$E_NO+dat_fig0$E_N2O+dat_fig0$E_N2_field
  
  dat_i=dat_fig0[dat_fig0$NO_stat.lab=='Y' & dat_fig0$year %in% seq(2019,2023),]
  dat_figt.1<-aggregate(dat_i[,c('NO.m','E_NO')], 
                        by=list(dat_i$chamber,dat_i$treatment,dat_i$plot),FUN = mean, na.rm=T)
  names(dat_figt.1)[1:5]<-c('chamber','treatment','plot','NO_N.m','E_NO')
  
  
  dat_i=dat_fig0[dat_fig0$N2O_stat.lab=='Y' & dat_fig0$year %in% seq(2019,2023),]
  dat_figt.2<-aggregate(dat_i[,c('N2O.m','E_N2O')], 
                        by=list(dat_i$chamber,dat_i$treatment,dat_i$plot),FUN = mean, na.rm=T)
  names(dat_figt.2)[1:5]<-c('chamber','treatment','plot','N2O_N.m','E_N2O')
  
  
  dat_i=dat_fig0[dat_fig0$NON2O_stat.lab=='Y' & dat_fig0$year %in% seq(2019,2023),]
  dat_figt.3<-aggregate(dat_i[,c('NON2O.m','E_NON2O','E_N2_field','E_gasN')], 
                        by=list(dat_i$chamber,dat_i$treatment,dat_i$plot),FUN = mean, na.rm=T)
  names(dat_figt.3)[1:7]<-c('chamber','treatment','plot','NON2O_N.m','E_NON2O','E_N2','E_gasN')
  
  dat_figt<-merge(dat_figt.1,dat_figt.2,by=c('chamber','treatment','plot'),all.x = T)
  dat_figt<-merge(dat_figt,dat_figt.3,by=c('chamber','treatment','plot'),all.x = T)
  rm(dat_figt.1,dat_figt.2,dat_figt.3)
  
  dat_figt$NO_stat.lab='Y'
  dat_figt$N2O_stat.lab='Y'
  dat_figt$NON2O_stat.lab='Y'
  
  ###筛选异常值control\warmed\model
  for(i in unique(dat_figt$treatment)){
    
    dfa<- dat_figt %>%
      filter(treatment == i)
    
    x1c=dfa$chamber[dfa$NO_N.m %in% boxplot(dfa$NO_N.m)$out]
    x2c=dfa$chamber[dfa$N2O_N.m %in% boxplot(dfa$N2O_N.m)$out]
    x5c=dfa$chamber[dfa$NON2O_N.m %in% boxplot(dfa$NON2O_N.m)$out]
    
    dat_figt$NO_stat.lab[dat_figt$chamber %in% x1c]='N'
    dat_figt$N2O_stat.lab[dat_figt$chamber %in% x2c]='N'
    dat_figt$NON2O_stat.lab[dat_figt$chamber %in% x5c]='N'
  }
  
  
  ###多重比较control\warmed\model  20240323，，与T检验结果有出入。。。
  str(dat_figt)
  unique(dat_figt$treatment)
  my_comparisons<-list(c('warmed','model'),c('control','warmed'))
  
  mean(dat_figt$E_NO[dat_figt$treatment=='warmed'])
  sd(dat_figt$E_NO[dat_figt$treatment=='warmed'])/sqrt(5)
  
  p1t<-ggplot(data=dat_figt,aes(treatment,E_NO,fill=treatment),)+#geom_point(pch=1)+
    geom_boxplot(data=dat_figt,width=0.4,alpha=0.7,size=0.2,outlier.size = 1,
                 position = position_dodge(0.8))+
    geom_point(data=dat_figt[dat_figt$NO_stat.lab=='Y',],show.legend = F,
               position = position_jitterdodge(0.035, dodge.width = .8),size=1,pch=1,alpha=0.4)+
    stat_summary(fun = mean, geom = "point", fill = "black",size=1.2) +
    scale_y_continuous(breaks=seq(0,1.5,.5),expand = c(0,0))+
    guides(y="prism_offset_minor")+
    scale_fill_manual(limits=c('control','warmed','model'),                      
                      values=c('royalblue1','red2','goldenrod'),name=NULL)+
    stat_compare_means(data=dat_figt[dat_figt$NO_stat.lab=='Y',],
                       aes(x=treatment,y=E_NO,label = paste0("p = ", after_stat(p.format))),
                       comparisons=my_comparisons,label.y = c(1.25,1.0),
                       method="t.test",vjust=0.5,hide.ns = T,
                       label="p.signif",tip.length = 0.01,size=4.2)+
    labs(x=NULL,y=expression(atop(paste('Cumulative NO flux'),paste('(kg N ha'^-1~'yr'^-1,')'))))+
    scale_x_discrete(limits=c('control','warmed','model'),labels=NULL)+
    mythem+theme(axis.ticks.length.x = unit(4,'mm'),               panel.background = element_rect(fill =NA, linewidth=0.4,colour = "black", linetype = "solid",inherit.blank=T),
                 axis.ticks.x = element_line(color='transparent'),
                 axis.text.x = element_blank())+
    coord_cartesian(ylim = c(0,1.5),clip = 'off')+
    annotate("segment", y = -0.04, yend = 0, x = 2, xend = 2, colour = "black");p1t
  
  mean(dat_figt$E_N2O[dat_figt$treatment=='warmed'])/mean(dat_figt$E_N2O[dat_figt$treatment=='control'])
  sd(dat_figt$E_N2O[dat_figt$treatment=='warmed'])/sqrt(5)
  
  p2t<-ggplot(data=dat_figt,aes(treatment,E_N2O,fill=treatment),)+#geom_point(pch=1)+
    geom_boxplot(data=dat_figt,width=0.4,alpha=0.7,size=0.2,outlier.size = 1,
                 position = position_dodge(0.8))+
    geom_point(data=dat_figt[dat_figt$N2O_stat.lab=='Y',],show.legend = F,
               position = position_jitterdodge(0.035, dodge.width = .8),size=1,pch=1,alpha=0.4)+
    stat_summary(fun = mean, geom = "point", fill = "black",size=1.2) +
    scale_y_continuous(breaks=seq(0,1.8,.6),expand = c(0,0))+
    guides(y="prism_offset_minor")+
    scale_fill_manual(limits=c('control','warmed','model'),                      
                      values=c('royalblue1','red2','goldenrod'),name=NULL)+
    stat_compare_means(data=dat_figt[dat_figt$N2O_stat.lab=='Y',],
                       aes(x=treatment,y=E_N2O,label = paste0("p = ", after_stat(p.format))),
                       comparisons=my_comparisons,label.y = c(1.6,1.3),
                       method="t.test",vjust=0.5,hide.ns = T,
                       label="p.signif",tip.length = 0.01,size=4.2)+
    labs(x=NULL,y=expression(atop(paste('Cumulative N'[2],'O flux'),paste('(kg N ha'^-1~'yr'^-1,')'))))+
    scale_x_discrete(limits=c('control','warmed','model'),labels=NULL)+
    mythem+theme(axis.ticks.length.x = unit(4,'mm'), panel.background = element_rect(fill =NA, linewidth=0.4,colour = "black", linetype = "solid",inherit.blank=T),
                 axis.ticks.x = element_line(color='transparent'),
                 axis.text.x = element_blank())+
    coord_cartesian(ylim = c(0,1.8),clip = 'off')+
    annotate("segment", y = -0.042, yend = 0, x = 2, xend = 2, colour = "black");p2t
  
  pt0<-ggplot(dat_figt,aes(treatment,E_NON2O,fill=treatment))+
    geom_boxplot(data=dat_figt,width=0.4,alpha=0.7,size=0.2,outlier.size = 1,
                 position = position_dodge(0.8))+
    geom_point(data=dat_figt,show.legend = F,
               position = position_jitterdodge(0.035, dodge.width = .8),size=1,pch=1,alpha=0.4)+
    stat_summary(fun = mean, geom = "point", fill = "black",size=1.2) +
    scale_y_continuous(breaks=seq(0,3,1),labels = c('0.0','1.0','2.0','3.0'),expand = c(0,0))+
    guides(y="prism_offset_minor")+
    scale_fill_manual(limits=c('control','warmed','model'),                      
                      values=c('royalblue1','red2','goldenrod'),name=NULL)+
    stat_compare_means(aes(x=treatment,y=E_NON2O,
                           label = paste0("p = ", after_stat(p.format))),
                       comparisons=my_comparisons,label.y = c(2.6, 2.15),
                       method="t.test",vjust=0.5,hide.ns = T, na.rm = T,
                       label="p.signif",tip.length = 0.01,size=4.2)+
    labs(x=NULL,y=expression(atop(paste('Cumulative NO + N'[2],'O flux'),paste('(kg N ha'^-1~'yr'^-1,')'))))+
    # scale_x_discrete(limits=c('control','warmed','model'),
    #                  labels=c(expression(paste("Control")),
    #                           expression(paste("Warmed")),expression(paste("Q"[10]," expected"))))+
    mythem+theme(axis.ticks.length.x = unit(4,'mm'), panel.background = element_rect(fill =NA, linewidth=0.4,colour = "black", linetype = "solid",inherit.blank=T),
                 axis.ticks.x = element_line(color='transparent'),
                 axis.text.x = element_blank())+
    coord_cartesian(ylim = c(0,3),clip = 'off')+
    annotate("segment", y = -0.10, yend = 0, x = 2, xend = 2, colour = "black")+
    annotate("text", label='2019-2023',y = -0.27,  x = 2, colour = "black",size=4.2);pt0
  
  
  figure<-ggarrange(p1+rremove("xlab")+rremove("x.text")+
                      theme(plot.margin = unit(c(0.25,0,0.13,0.5),'cm'),
                            panel.background = element_rect(fill =NA, linewidth=0.65,colour = "black", linetype = "solid",inherit.blank=T),
                            axis.text.y = element_text(size=12),axis.title.y = element_text(size=12)),
                    p1m+rremove("xlab")+rremove("x.text")+
                      theme(plot.margin = unit(c(0.25,0.05,0.36,0.25),'cm'),
                            axis.text.y = element_text(size=12),axis.title.y = element_text(size=12)),
                    p1t+rremove("xlab")+rremove("x.text")+
                      theme(plot.margin = unit(c(0.25,0.5,0,0.15),'cm'),
                            axis.text.y = element_text(size=12),axis.title.y = element_text(size=12)),
                    
                    p2+rremove("xlab")+rremove("x.text")+
                      theme(plot.margin = unit(c(0.25,0,0.13,0.5),'cm'),
                            panel.background = element_rect(fill =NA, linewidth=0.65,colour = "black", linetype = "solid",inherit.blank=T),
                            axis.text.y = element_text(size=12),axis.title.y = element_text(size=12)),
                    p2m+rremove("xlab")+rremove("x.text")+
                      theme(plot.margin = unit(c(0.25,0.05,0.36,0.25),'cm'),
                            axis.text.y = element_text(size=12),axis.title.y = element_text(size=12)),
                    p2t+rremove("xlab")+rremove("x.text")+
                      theme(plot.margin = unit(c(0.25,0.5,0,0.15),'cm'),
                            axis.text.y = element_text(size=12),axis.title.y = element_text(size=12)),
                    
                    p5+theme(strip.text = element_blank(),strip.background = element_blank(),
                             panel.background = element_rect(fill =NA, linewidth=0.65,colour = "black", linetype = "solid",inherit.blank=T),
                             axis.text.x = element_text(size=12),axis.title.x = element_text(size=12),
                             axis.text.y = element_text(size=12),axis.title.y = element_text(size=12),
                             plot.margin = unit(c(0.25,0,0.05,0.5),'cm')),
                    p5m+rremove("xlab")+
                      theme(plot.margin = unit(c(0.25,0.05,1.3,0.25),'cm'),
                            axis.text.y = element_text(size=12),axis.title.y = element_text(size=12)),
                    pt0+rremove("xlab")+
                      theme(strip.text = element_blank(),strip.background = element_blank(),
                            axis.text.y = element_text(size=12),axis.title.y = element_text(size=12),
                            plot.margin = unit(c(0.25,0.5,0.9,0.15),'cm')),
                    
                    labels = c("A", "B","C","D","E","F","G","H","I"),
                    label.x = c(0.16,0.40,0.38),label.y = c(0.95,0.95,0.95),
                    ncol = 3, nrow = 3,align = "none",   ##"v"竖直对齐
                    font.label = list(size = 12, color ="black"),
                    widths = c(10,3.5,3.8), heights = c(5,5,5.5),
                    common.legend=F);figure
  
  getwd()
  ggsave("Fig_4_NO+N2O response to warming1.pdf", figure, width =11, height =7.4,
         device=cairo_pdf)  #chamber尺度分析平均排放速率对增温的响应
} 


#### N2 emissions and total gaseous N emissions
str(dat_figt)

mean(dat_figt$E_N2[dat_figt$treatment=='control'])
sd(dat_figt$E_N2[dat_figt$treatment=='control'])/sqrt(15)

mean(dat_figt$E_gasN[dat_figt$treatment=='control'])
sd(dat_figt$E_gasN[dat_figt$treatment=='control'])/sqrt(15)

str(dat_figt)

{
  p3t<-ggplot(dat_figt,aes(treatment,E_N2,fill=treatment))+
    geom_boxplot(data=dat_figt,width=0.45,alpha=0.7,size=0.2,outlier.size = 1,
                 position = position_dodge(0.8))+
    geom_point(data=dat_figt,show.legend = F,
               position = position_jitterdodge(0.035, dodge.width = .8),size=1,pch=1,alpha=0.4)+
    stat_summary(fun = mean, geom = "point", fill = "black",size=1.2) +
    scale_y_continuous(breaks=seq(0,3,1),expand = c(0,0))+
    guides(y="prism_offset_minor")+
    scale_fill_manual(limits=c('control','warmed','model'),                      
                      values=c('royalblue1','red2','goldenrod'),name=NULL)+
    stat_compare_means(aes(x=treatment,y=E_N2,
                           label = paste0("p = ", after_stat(p.format))),
                       comparisons=my_comparisons,label.y = c(2.6, 2.15),
                       method="t.test",vjust=0.5,hide.ns = F, na.rm = T,
                       label="p.signif",tip.length = 0.01,size=4.2)+
    labs(x=NULL,y=expression(atop(paste('Cumulative N'[2],' flux'),paste('(kg N ha'^-1~'yr'^-1,')'))))+
    # scale_x_discrete(limits=c('control','warmed','model'),
    #                  labels=c(expression(paste("Control")),
    #                           expression(paste("Warmed")),expression(paste("Q"[10]," expected"))))+
    mythem+theme(axis.ticks.length.x = unit(4,'mm'),
                 panel.background = element_rect(fill =NA, linewidth=0.4,colour = "black", linetype = "solid",inherit.blank=T),
                 axis.ticks.x = element_line(color='transparent'),
                 axis.text.x = element_blank())+
    coord_cartesian(ylim = c(0,3),clip = 'off')+
    annotate("segment", y = -0.07, yend = 0, x = 2, xend = 2, colour = "black");p3t
  
  
  p4t<-ggplot(dat_figt,aes(treatment,E_gasN,fill=treatment))+
    geom_boxplot(data=dat_figt,width=0.45,alpha=0.7,size=0.2,outlier.size = 1,
                 position = position_dodge(0.8))+
    geom_point(data=dat_figt,show.legend = F,
               position = position_jitterdodge(0.035, dodge.width = .8),size=1,pch=1,alpha=0.4)+
    stat_summary(fun = mean, geom = "point", fill = "black",size=1.2) +
    scale_y_continuous(breaks=seq(0,6,2),expand = c(0,0))+
    guides(y="prism_offset_minor")+
    scale_fill_manual(limits=c('control','warmed','model'),                      
                      values=c('royalblue1','red2','goldenrod'),name=NULL)+
    stat_compare_means(aes(x=treatment,y=E_gasN,
                           label = paste0("p = ", after_stat(p.format))),
                       comparisons=my_comparisons,label.y = c(5, 4),
                       method="t.test",vjust=0.5,hide.ns = T, na.rm = T,
                       label="p.signif",tip.length = 0.01,size=4.2)+
    labs(x=NULL,y=expression(atop(paste('Cumulative gaseous N flux'),paste('(kg N ha'^-1~'yr'^-1,')'))))+
    # scale_x_discrete(limits=c('control','warmed','model'),
    #                  labels=c(expression(paste("Control")),
    #                           expression(paste("Warmed")),expression(paste("Q"[10]," expected"))))+
    mythem+theme(axis.ticks.length.x = unit(4,'mm'),
                 panel.background = element_rect(fill =NA, linewidth=0.4,colour = "black", linetype = "solid",inherit.blank=T),
                 axis.ticks.x = element_line(color='transparent'),
                 axis.text.x = element_blank())+
    coord_cartesian(ylim = c(0,6),clip = 'off')+
    annotate("segment", y = -0.15, yend = 0, x = 2, xend = 2, colour = "black")+
    annotate("text", label='2019-2023',y = -0.3,  x = 2, colour = "black",size=4.2);p4t
  
  
  figure<-ggarrange(p3+rremove("xlab")+rremove("x.text")+
                      theme(plot.margin = unit(c(0.5,0,0.3,0.5),'cm'),legend.position = c(0.28,0.83),
                            legend.text = element_text(size=12),
                            axis.text.y = element_text(size=12),axis.title.y = element_text(size=12)),
                    p3m+rremove("xlab")+rremove("x.text")+
                      theme(plot.margin = unit(c(0.5,0,0.5,0.5),'cm'),
                            axis.text.y = element_text(size=12),axis.title.y = element_text(size=12,vjust=3)),
                    p3t+rremove("xlab")+rremove("x.text")+
                      theme(plot.margin = unit(c(0.5,0.5,0.13,0.25),'cm'),
                            axis.text.y = element_text(size=12),axis.title.y = element_text(size=12)),
                    
                    p4+theme(strip.text = element_blank(),strip.background = element_blank(),
                             axis.title.x = element_text(size=12),axis.text.x = element_text(size=12),
                             axis.title.y = element_text(size=12),axis.text.y = element_text(size=12),
                             plot.margin = unit(c(0.25,0,1.55,0.5),'cm')),
                    p4m+rremove("xlab")+rremove("x.text")+
                      theme(plot.margin = unit(c(0.25,0,2.8,0.25),'cm'),
                            axis.text.y = element_text(size=12),axis.title.y = element_text(size=12)),
                    p4t+rremove("xlab")+
                      theme(strip.text = element_blank(),strip.background = element_blank(),
                            axis.text.y = element_text(size=12),axis.title.y = element_text(size=12),
                            plot.margin = unit(c(0.25,0.5,2.35,0.25),'cm')),
                    
                    labels = c("A", "B","C","D","E","F"),
                    label.x = c(0.17,0.45,0.39),label.y = c(0.94,0.94,0.94,0.97,0.97,0.97),
                    ncol = 3, nrow = 2,align = "none",   ##"v"竖直对齐
                    font.label = list(size = 12, color ="black"),
                    widths = c(10,3.5,3.5), heights = c(5,6.3),
                    common.legend=F);figure
  
  getwd()
  ggsave("Supplementary_Fig.5_gaseous N response to warming.pdf", figure, width =11, height =6,
         device=cairo_pdf)  #chamber尺度分析平均排放速率对增温的响应
}



###计算不同气态氮累积减少量
(cum_pT$NO_N.cum[cum_pT$treatment=='warmed']-cum_pT$NO_N.cum[cum_pT$treatment=='control'])/
  cum_pT$NO_N.cum[cum_pT$treatment=='control']

(cum_pT$N2O_N.cum[cum_pT$treatment=='warmed']-cum_pT$N2O_N.cum[cum_pT$treatment=='control'])/
  cum_pT$N2O_N.cum[cum_pT$treatment=='control']

(cum_pT$N2_N.cum[cum_pT$treatment=='warmed']-cum_pT$N2_N.cum[cum_pT$treatment=='control'])/
  cum_pT$N2_N.cum[cum_pT$treatment=='control']

(cum_pT$gasN[cum_pT$treatment=='warmed']-cum_pT$gasN[cum_pT$treatment=='control'])/
  cum_pT$gasN[cum_pT$treatment=='control']


################# Fig_S8_Seasonal contribution of gaseous N emission from 2019 to 2023 ###################

### output_emit_gc写出数据20241021 control
  mydata_c$season_period<-
    ifelse(format(mydata_c$date,"%m") %in% c('05','06','07','08','09','10'),'Plant growing season',
           'Plant dormant season')
  
  output_emit_gc<-data.frame()   #增温期间各chamber有效排放尤其是气态氮
  k=1
  
unique(mydata_c$season_period)
for(t in c("Plant dormant season","Plant growing season","All")){
  for(y in 2019:2023){
    
    if(t %in% c("Plant dormant season","Plant growing season")){
      dat=mydata_c[which(mydata_c$year==y & mydata_c$season_period==t),]
    }else{
      dat=mydata_c[which(mydata_c$year==y),]}
    
    for(i in na.omit(unique(dat$system_id))){ 
        datai=dat[which(dat$system_id==i),]
        dataii<-datai[datai$soil_temp.m>=5 & datai$warming=='Y',]   #专为N2
        
        if(length(datai$NO_N.m)>0){
          output_emit_gc[k,1]=i   #chamber num
          output_emit_gc[k,2]=t    #season
          output_emit_gc[k,3]=y   #year
          
          output_emit_gc[k,4]=round(mean(na.omit(datai$NO_N.m),na.rm=T),3)  #NO平均值
          output_emit_gc[k,5]=length(datai$date[complete.cases(datai$NO_N.m)]) 
          output_emit_gc[k,6]=round(output_emit_gc[k,4]*output_emit_gc[k,5]*24*10^4/10^9,3)  #单位kg N ha-1 

          output_emit_gc[k,7]=round(mean(na.omit(datai$N2O_N.m),na.rm=T),3)  #N2O平均值
          output_emit_gc[k,8]=length(datai$date[complete.cases(datai$N2O_N.m)])  #N2O实测days
          output_emit_gc[k,9]=round(output_emit_gc[k,7]*output_emit_gc[k,8]*24*10^4/10^9,3)   #单位kg N ha-1 

          # Plant dormant season不计算N2
          output_emit_gc[k,10]=ifelse(t=='Plant dormant season',NA,round(mean(na.omit(dataii$N2_N.m),na.rm=T),3))  #N2_field平均值
          output_emit_gc[k,11]=ifelse(t=='Plant dormant season',NA,length(dataii$date[complete.cases(dataii$N2_N.m)]))  #N2_field实测days
          output_emit_gc[k,12]=ifelse(t=='Plant dormant season',NA,round(output_emit_gc[k,10]*output_emit_gc[k,11]*24*10^4/10^9,3))   #N2_field单位kg N ha-1 

          output_emit_gc[k,13]=unique(datai$plot)
          output_emit_gc[k,14]=unique(datai$treatment)
          
          k=k+1
        }
      }
    }
}

rm(datai,dataii,dat)

######单位N2O_N,NO_N,(g N ha-1)  
 names(output_emit_gc)[1:14]<-c('chamber','season','year','NO.m','NO.n','E_NO','N2O.m','N2O.n','E_N2O','N2.m','N2.n','E_N2','plot','treatment')
 
 output_emit_gc$treatment<-ifelse(output_emit_gc$plot %in% c("1","4","5"),"warmed","control")
  
  
 output_emit_gc.E<-
   summaryBy(E_NO+E_N2O+E_N2~treatment+season+year, FUN = myfun,data=output_emit_gc)
 
 output_emit_gc.E$gasN<-output_emit_gc.E$E_NO.m+output_emit_gc.E$E_N2O.m+output_emit_gc.E$E_N2.m
 
 ###N2 释放量生长季
 datFUN(output_emit_gc.E$E_N2.m[output_emit_gc.E$treatment=='control' & output_emit_gc.E$season=='Plant growing season'])
 
 ### 气态氮 释放量生长季
 datFUN(output_emit_gc.E$gasN[output_emit_gc.E$treatment=='control' & output_emit_gc.E$season=='All'])
 
 ### NO 释放量全年
 datFUN(output_emit_gc.E$E_NO.m[output_emit_gc.E$treatment=='control' & output_emit_gc.E$season=='All'])
 
 ### N2O 释放量全年
 datFUN(output_emit_gc.E$E_N2O.m[output_emit_gc.E$treatment=='control' & output_emit_gc.E$season=='All'])
 
 library(tidyverse)
 library(writexl)
 
 ###将不同的结果卸载一个xls文件
 df = list(output_emit_gc,output_emit_gc.E) %>%
   set_names(c('Chamber','Treatment'))
 
 write_xlsx(df, "Output_gasN_20241027.xlsx")   ### Supplementary Data Table 1
 
 ###对照样地NO和N2O释放,文章数据
 str(F.treat)
 datFUN(F.treat$NO_N.m[F.treat$treatment=='control' & F.treat$year<=2023 & F.treat$NO_N.m>=0])
 datFUN(F.treat$N2O_N.m[F.treat$treatment=='control' & F.treat$year<=2023])
 
 
### plant growing and dormant periods emissions from 2021 to 2023
  
  str(mydata_c)
  
  ###(1)
  mydata_c$season_period<-
    ifelse(format(mydata_c$date,"%m") %in% c('05','06','07','08','09','10'),'Plant growing season',
           'Plant dormant season')
  
  emit_sf1<-data.frame()   #冻融期间排放（含处理）
  
  emit_sf_c<-data.frame()   #对照冻融期间排放 
  unique(mydata_c$warming)
  unique(mydata_c$season_period)
  
  k=1
  for(i in 2019:2023){
    for(j in unique(mydata_c$season_period)){
      dat=mydata_c[which(mydata_c$season_period==j &mydata_c$year==i & mydata_c$warming=='Y'),]
      dat.c=mydata_c[which(mydata_c$season_period==j &mydata_c$year==i),]
      
      for(t in na.omit(unique(dat$system_id))){ 
          datai=dat[which(dat$system_id==t),]
          datai.c=dat.c[which(dat.c$system_id==t),]
          
        if(length(datai$NO_N.m)>0){
          emit_sf1[k,1]=t   #chamber num
          emit_sf1[k,2]=j    #labs_warming
          emit_sf1[k,3]=i   #year
          
          emit_sf1[k,4]=paste(round(min(na.omit(datai$N2O_N.m)),2),
                            round(max(na.omit(datai$N2O_N.m)),2),sep="<>")  #N2O排放范围
          emit_sf1[k,5]=round(mean(na.omit(datai$N2O_N.m),na.rm=T),3)  #N2O平均值
          emit_sf1[k,6]=length(datai$date[complete.cases(datai$N2O_N.m)])  #N2O实测days
          emit_sf1[k,7]=round(emit_sf1[k,5]*emit_sf1[k,6]*24*10^4/10^9,3)   #单位kg N ha-1 
          emit_sf1[k,8]=round(sd(na.omit(datai$N2O_N.m),na.rm=T)/sqrt(emit_sf1[k,5])*24*10^4/10^9,3)
          
          emit_sf1[k,9]=paste(round(min(na.omit(datai$NO_N.m)),2),
                             round(max(na.omit(datai$NO_N.m)),2),sep="<>")   #NO排放范围
          emit_sf1[k,10]=round(mean(na.omit(datai$NO_N.m),na.rm=T),3)  #NO平均值
          emit_sf1[k,11]=length(datai$date[complete.cases(datai$NO_N.m)]) 
          emit_sf1[k,12]=round(emit_sf1[k,10]*emit_sf1[k,11]*24*10^4/10^9,3)  #单位kg N ha-1 
          emit_sf1[k,13]=round(sd(na.omit(datai$NO_N.m),na.rm=T)/sqrt(emit_sf1[k,10])*24*10^4/10^9,3)
          

          emit_sf1[k,14]=paste(round(min(na.omit(datai$soil_temp.m[complete.cases(datai$NO_N.m)])),2),
                             round(max(na.omit(datai$soil_temp.m[complete.cases(datai$NO_N.m)])),2),sep="<>")   
          emit_sf1[k,15]=round(mean(na.omit(datai$soil_temp.m[complete.cases(datai$NO_N.m)]),na.rm=T),3)  
          emit_sf1[k,16]=sd(na.omit(datai$soil_temp.m[complete.cases(datai$NO_N.m)]),na.rm=T)/emit_sf1[k,15]  
          
          emit_sf1[k,17]=paste(round(min(na.omit(datai$soil_moisture.m[complete.cases(datai$NO_N.m)])),2),
                             round(max(na.omit(datai$soil_moisture.m[complete.cases(datai$NO_N.m)])),2),sep="<>")   
          emit_sf1[k,18]=round(mean(na.omit(datai$soil_moisture.m[complete.cases(datai$NO_N.m)]),na.rm=T),3)  
          emit_sf1[k,19]=sd(na.omit(datai$soil_moisture.m[complete.cases(datai$NO_N.m)]),na.rm=T)/emit_sf1[k,18] 
          
          emit_sf1[k,20]=unique(datai$plot)
          emit_sf1[k,21]=unique(datai$treatment)
          
          
          emit_sf_c[k,1]=t   #chamber num
          emit_sf_c[k,2]=j    #labs_warming
          emit_sf_c[k,3]=i   #year
          emit_sf_c[k,4]=paste(round(min(na.omit(datai.c$N2O_N.m)),2),
                              round(max(na.omit(datai.c$N2O_N.m)),2),sep="<>")  #N2O排放范围
          emit_sf_c[k,5]=round(mean(na.omit(datai.c$N2O_N.m),na.rm=T),3)  #N2O平均值
          emit_sf_c[k,6]=length(datai.c$date[complete.cases(datai.c$N2O_N.m)])  #N2O实测days
          emit_sf_c[k,7]=round(emit_sf_c[k,5]*emit_sf_c[k,6]*24*10^4/10^9,3)   #单位kg N ha-1 
          emit_sf_c[k,8]=round(sd(na.omit(datai.c$N2O_N.m),na.rm=T)/sqrt(emit_sf_c[k,5])*24*10^4/10^9,3)
          
          emit_sf_c[k,9]=paste(round(min(na.omit(datai.c$NO_N.m)),2),
                              round(max(na.omit(datai.c$NO_N.m)),2),sep="<>")   #NO排放范围
          emit_sf_c[k,10]=round(mean(na.omit(datai.c$NO_N.m),na.rm=T),3)  #NO平均值
          emit_sf_c[k,11]=length(datai.c$date[complete.cases(datai.c$NO_N.m)]) 
          emit_sf_c[k,12]=round(emit_sf_c[k,10]*emit_sf_c[k,11]*24*10^4/10^9,3)  #单位kg N ha-1 
          emit_sf_c[k,13]=round(sd(na.omit(datai.c$NO_N.m),na.rm=T)/sqrt(emit_sf_c[k,10])*24*10^4/10^9,3)
          
          emit_sf_c[k,14]=unique(datai.c$plot)
          emit_sf_c[k,15]=unique(datai.c$treatment)
          
          k=k+1
        }
      }
    }
  }
  
  ######单位N2O_N,NO_N,(g N ha-1)
  
  names(emit_sf1)[1:21]<-c('chamber','season_period','year',
                         'N2O_range','N2O.m','days_N2O','E_N2O','N2O.se',
                         'NO_range','NO.m','days_NO','E_NO','NO.se',
                         'temp_range','temp','temp.cv','VWC_range','VWC','VWC.cv','plot','treatment')
  
  
  names(emit_sf_c)[1:15]<-c('chamber','season_period','year',
                           'N2O_range','N2O.m','days_N2O','E_N2O','N2O.se',
                           'NO_range','NO.m','days_NO','E_NO','NO.se',
                           'plot','treatment')
  
  emit_sf1$treatment<-ifelse(emit_sf1$plot %in% c("1","4","5"),"warmed","control")

  
 ###确定异常值for emit_sf1
  emit_sf1$NO_stat.lab='Y'
  emit_sf1$N2O_stat.lab='Y'

 ###筛选异常值control\warmed\model
  for(j in 2018:2023){
    for(i in unique(emit_sf1$treatment)){
      for(t in unique(emit_sf1$season_period)){
      dfa<- emit_sf1 %>%
        filter(treatment == i & year == j & season_period == t)

      if(length(dfa$E_NO)==15){
      x1c=dfa$chamber[dfa$E_NO %in% boxplot(dfa$E_NO)$out]
      x2c=dfa$chamber[dfa$E_N2O %in% boxplot(dfa$E_N2O)$out]

      emit_sf1$NO_stat.lab[emit_sf1$chamber %in% x1c & emit_sf1$year==j & emit_sf1$season_period==t]='N'
      emit_sf1$N2O_stat.lab[emit_sf1$chamber %in% x2c & emit_sf1$year==j & emit_sf1$season_period==t]='N'
      }
      }
    }
  }
  
  
 #########不同期间NO和N2O绘图
  str(emit_sf1)
  datai<-emit_sf1[emit_sf1$year %in% c(2019,2020,2021,2022,2023),]
  
  library(gghalves);library(ggprism)
  
 
  emit_sf1.range <- data.frame(
    treatment = rep("control", 4),
    year = as.factor(rep(2019,4)),
    E_NO = c(c(0,2),c(0,0.4)),
    E_N2O = c(c(0,2),c(0,4)),
    season_period = c(rep("Plant growing season",2), 
                      rep("Plant dormant season",2)))
  
  datai<-within(datai,season_period<-factor(season_period,
                                            levels = c('Plant growing season','Plant dormant season')))
  
  emit_sf1.range<-within(emit_sf1.range,season_period<-factor(season_period,
                                                              levels = c('Plant growing season','Plant dormant season')))
  
  datai<-within(datai,year<-factor(year, levels=c(2019,2020,2021,2022,2023)))
  
  p1<-ggplot(datai,aes(year, E_NO,col=treatment))+
    facet_wrap(.~season_period,scales = 'free_y')+
    geom_point(data=emit_sf1.range,col='transparent')+
    scale_y_continuous(n.breaks = 3,expand=c(0,0))+
    geom_boxplot(outlier.shape = NA,width = 0.55,alpha = 1,fill='white',linewidth = 0.2)+
    geom_point(position = position_jitterdodge(0.05, dodge.width = .55),
               pch=1,size=1,alpha=0.79,stroke=0.3,show.legend = T)+
    scale_color_manual(limits=c('control','warmed'),values=c('blue3','red4'))+
    ggpubr::stat_compare_means(data=datai[datai$NO_stat.lab=='Y',],
                               aes(group=treatment),show.legend = F,
                               label="p.signif",method = "t.test", #非参数检验
                               hide.ns = F,size=3.5,label.y.npc = 0.9,label.x = 1.5)+
    guides(y="prism_offset_minor")+labs(x='Year',y=expression(paste('NO emission (kg N ha'^-1,')')))+
    mythem_sci+theme(strip.background = element_rect(fill = alpha("grey", 0.3)));p1
  
  
 p2<-ggplot(datai,aes(year, E_N2O,col=treatment))+
    facet_wrap(.~season_period,scales = 'free_y')+
    scale_y_continuous(n.breaks = 3, expand=c(0,0.00001))+
    geom_point(data=emit_sf1.range,col='transparent')+
    geom_boxplot(outlier.shape = NA,width = 0.55,alpha = 1,fill='white',linewidth = 0.2)+
    geom_point(position = position_jitterdodge(0.05, dodge.width = .55),
               pch=1,size=1,alpha=0.79,stroke=0.3)+
    scale_color_manual(limits=c('control','warmed'),values=c('blue3','red4'))+
    ggpubr::stat_compare_means(data=datai[datai$N2O_stat.lab=='Y',],
                               aes(group=treatment),#show.legend = F,
                               tip.length = 0.03,bracket.size = 0.3,
                               label="p.signif",method = "t.test", #非参数检验
                               hide.ns = T,size=3.5,label.y.npc = 0.9,label.x = 1.5)+
   guides(y="prism_offset_minor")+labs(x='Year',y=expression(paste('N'[2],'O emission (kg N ha'^-1,')')))+
   mythem_sci+theme(strip.background = element_rect(fill = alpha("grey", 0.3)));p2
 
 (p1+rremove('x.text')+rremove('xlab')+
     theme(legend.position = c(0.07,0.65),legend.key.height=unit(0.3,"cm"),
           legend.key.width=unit(0.4,"cm"),legend.title = element_blank()))/
   (p2+theme(plot.margin = unit(c(0.5,0.1,0,0),'cm'),
             strip.background = element_blank(),strip.text = element_blank()))->p;p
 
 ###########绘制bar图
 str(datai)   #E_N2O和E_NO
 datai_vps<-summaryBy(E_NO+E_N2O~year+treatment+season_period, data=datai, FUN = myfun)
 

 datai_vps$R_NO=NA  #计算NO排放占比
 datai_vps$R_N2O=NA   #计算N2O排放占比
 
 i=1
 for(i in 1:length(datai_vps$year)){
   A_NO=sum(datai_vps$E_NO.m[datai_vps$treatment==datai_vps$treatment[i] & datai_vps$year==datai_vps$year[i]])
   A_N2O=sum(datai_vps$E_N2O.m[datai_vps$treatment==datai_vps$treatment[i] & datai_vps$year==datai_vps$year[i]])
   datai_vps$R_NO[i]=datai_vps$E_NO.m[i]*100/A_NO   #百分贡献
   datai_vps$R_N2O[i]=datai_vps$E_N2O.m[i]*100/A_N2O  #百分贡献
   
   
   i=i+1
 }
 
 str(datai_vps)
 
 ### NO的年排放及占比
 unique(datai_vps$year)
 
 datFUN(datai_vps$E_NO.m[datai_vps$season_period=='Plant growing season' & datai_vps$treatment =='control']+
   datai_vps$E_NO.m[datai_vps$season_period=='Plant dormant season'& datai_vps$treatment =='control'])
 
 datFUN(datai_vps$R_NO[datai_vps$season_period=='Plant growing season' & datai_vps$treatment =='control'])
 
 
 ### N2O的年排放及占比
 datFUN(datai_vps$R_N2O[datai_vps$season_period=='Plant growing season' & datai_vps$treatment =='control'])
 
 datFUN(datai_vps$E_N2O.m[datai_vps$season_period=='Plant growing season' & datai_vps$treatment =='control']+
          datai_vps$E_N2O.m[datai_vps$season_period=='Plant dormant season'& datai_vps$treatment =='control'])
 
 ### 增温之后
 datFUN(datai_vps$E_NO.m[datai_vps$season_period=='Plant growing season' & datai_vps$treatment =='warmed']+
          datai_vps$E_NO.m[datai_vps$season_period=='Plant dormant season'& datai_vps$treatment =='warmed'])
 
 datFUN(datai_vps$E_N2O.m[datai_vps$season_period=='Plant growing season' & datai_vps$treatment =='warmed']+
          datai_vps$E_N2O.m[datai_vps$season_period=='Plant dormant season'& datai_vps$treatment =='warmed'])
 
 
 datFUN(datai_vps$R_N2O[datai_vps$season_period=='Plant dormant season' & datai_vps$treatment =='control' &
                          datai_vps$year %in% c(2021,2022,2023)])
 
 
 datFUN((datai_vps$E_N2O.m[datai_vps$season_period=='Plant growing season' & datai_vps$treatment =='warmed']+
          datai_vps$E_N2O.m[datai_vps$season_period=='Plant dormant season'& datai_vps$treatment =='warmed'])/
          (datai_vps$E_N2O.m[datai_vps$season_period=='Plant growing season' & datai_vps$treatment =='control']+
             datai_vps$E_N2O.m[datai_vps$season_period=='Plant dormant season'& datai_vps$treatment =='control']))
 
 
 library(ggpattern)   ##直接从2023library复制过来
 
 unique(datai_vps$season_period)
 datai_vps<-within(datai_vps,
                   season_period<-factor(season_period,levels=c('Plant dormant season','Plant growing season')))
 
 ggplot(datai_vps,
        aes(treatment,R_NO,fill=treatment,group=season_period,pattern=season_period))+
   facet_wrap(.~year,strip.position='bottom',nrow=1)+
   scale_fill_manual(limits=c('control','warmed'),values=c('blue3','red3'),name=' Treatment')+
   scale_pattern_manual(limits=c('Plant growing season','Plant dormant season'),
                        values = c("stripe","none"),name=' Periods')+
   # geom_bar(stat="identity",position=position_stack(0),col='black',width=.55,alpha=0.005,linewidth=0.1)+
   scale_y_continuous(limits=c(0,100),breaks=seq(0,100,25),expand=c(0,0),
                      name=expression(atop(paste('Seasonal contribution'),
                                                paste('of NO emission (%)'))))+
   mythem_sci+
   theme(panel.spacing.x = unit(0,'mm'),strip.background = element_blank(),strip.text = element_blank(),
         axis.text.x = element_blank(),axis.ticks.length.x = unit(0,'mm'),
         axis.title.x = element_blank(),legend.position = 'top')+
   geom_col_pattern(aes(treatment,R_NO,fill=treatment,pattern=season_period),alpha=0.45,
                    position = position_stack(1),color='black',width=.55,linewidth=0.01,
                    pattern_spacing = 0.17,pattern_density = 0.5,pattern_size = 0.05, pattern_fill = "transparent")+
   guides(pattern = guide_legend(override.aes = list(fill = "white", pattern_fill = "black",pattern_density = 0.03,pattern_spacing = 0.03),
                                 keywidth = 1, keyheight = 1,nrow = 2),
          fill = guide_legend(override.aes = list(pattern = "none"),keywidth = 1, keyheight = 1,nrow = 2))->p3u;p3u
   
 
 
 ggplot(datai_vps,
        aes(treatment,R_N2O,fill=treatment,group=season_period,pattern=season_period))+
   geom_text(data=datai_vps[datai_vps$year=='2019' & datai_vps$treatment=='control' &
                              datai_vps$season_period=='Plant growing season',],
             aes(y= -5, label='2019'),hjust=-0.5,size=2.3)+
   geom_text(data=datai_vps[datai_vps$year=='2020' & datai_vps$treatment=='control' &
                              datai_vps$season_period=='Plant growing season',],
             aes(y= -5, label='2020'),hjust=-0.5,size=2.3)+
   geom_text(data=datai_vps[datai_vps$year=='2021' & datai_vps$treatment=='control' &
                              datai_vps$season_period=='Plant growing season',],
             aes(y= -5, label='2021'),hjust=-0.5,size=2.3)+
   geom_text(data=datai_vps[datai_vps$year=='2022' & datai_vps$treatment=='control' &
                              datai_vps$season_period=='Plant growing season',],
             aes(y= -5, label='2022'),hjust=-0.5,size=2.3)+
   geom_text(data=datai_vps[datai_vps$year=='2023' & datai_vps$treatment=='control' &
                              datai_vps$season_period=='Plant growing season',],
             aes(y= -5, label='2023'),hjust=-0.5,size=2.3)+
   facet_wrap(.~year,strip.position='bottom',nrow=1)+
   scale_fill_manual(limits=c('control','warmed'),values=c('blue3','red3'),name='Treatment')+
   scale_pattern_manual(limits=c('Plant growing season','Plant dormant season'),
                        values = c("stripe","none"),name=' Periods')+

   scale_y_continuous(breaks=seq(0,100,25),expand=c(0,0),
                      name=expression(atop(paste('Seasonal contribution'),
                      paste('of N'[2],'O emission (%)'))))+
   mythem_sci+coord_cartesian(clip="off",ylim=c(0,100))+
   theme(panel.spacing.x = unit(0,'mm'),strip.background = element_blank(),strip.text = element_blank(),
         axis.text.x = element_blank(),axis.ticks.length.x = unit(0,'mm'),
         plot.margin = unit(c(0.5,0.1,1,0),'cm'),
         axis.title.x = element_blank(),legend.position = 'none')+
   geom_col_pattern(aes(treatment,R_N2O,fill=treatment,pattern=season_period),alpha=0.45,
                    position = position_stack(1),color='black',width=.55,linewidth=0.01,
                    pattern_spacing = 0.17,pattern_density = 0.5,pattern_size = 0.05, pattern_fill = "transparent")+
   guides(pattern = guide_legend(override.aes = list(fill = "white")),
          fill = guide_legend(override.aes = list(pattern = "none")))->p3d;p3d
 
 
 
 
 p3u/p3d->p3;p3
 
 ((p1+rremove('x.text')+rremove('xlab')+
     theme(legend.position = c(0.09,0.65),legend.key.height=unit(0.3,"cm"),
           legend.key.width=unit(0.4,"cm"),legend.title = element_blank()))+p3u+
     plot_layout(ncol=2,nrow=1,widths=c(5,3.5),heights=c(5)))/
     (p2+theme(plot.margin = unit(c(0.5,0.1,1,0),'cm'),
               strip.background = element_blank(),strip.text = element_blank())+
        p3d+plot_layout(ncol=2,nrow=1,widths=c(5,3.5),heights=c(5,5)))->figure;figure
 
 ggsave("Fig_S8_NO+N2O emissions_box+bar 2019-2023.pdf", figure, width =7, height =5) 
 
 
 ### stat results of NO and N2O
 summary(aov(E_NO~treatment*factor(year)+Error(1/chamber),
             data=datai[datai$season_period=='Plant growing season' & datai$year %in% c(2019,2020,2021,2022,2023)&
                          datai$NO_stat.lab=='Y',]))
 
 summary(aov(E_N2O~treatment*factor(year)+Error(1/chamber),
             data=datai[datai$season_period=='Plant growing season' & datai$year %in% c(2019,2020,2021,2022,2023)&
                          datai$N2O_stat.lab=='Y',]))
 
 
 summary(aov(E_NO~treatment*factor(year)+Error(1/chamber),
             data=datai[datai$season_period=='Plant dormant season' & datai$year %in% c(2019,2020,2021,2022,2023)&
                          datai$NO_stat.lab=='Y',]))
 
 
 summary(aov(E_N2O~treatment*factor(year)+Error(1/chamber),
             data=datai[datai$season_period=='Plant dormant season' & datai$year %in% c(2019,2020,2021,2022,2023) &
                          datai$N2O_stat.lab=='Y',]))
 
  str(emit_sf_c)
  
  names(emit_sf_c)
  emit_sf_yr<-aggregate(emit_sf_c[,c("E_NO","E_N2O")], 
                        FUN=sum,by=list(emit_sf_c$treatment, emit_sf_c$year, 
                                        emit_sf_c$chamber))
  
  names(emit_sf_yr)[1:3]<-c("treatment",'year','chamber')
  emit_sf_yr<-summaryBy(data = emit_sf_yr, E_NO+E_N2O ~ year +treatment,
                       FUN = myfun)   ### extended table 1
  
  emit_sf_a<-summaryBy(data = emit_sf_c, E_NO+E_N2O ~ year +treatment+season_period,FUN = myfun)
  ###写出数据
  
  ###将不同的结果写出到一个xls文件
  #（1）增温期间日排放——呼吸室水平；（2）增温期间累计排放-呼吸室水平
  #（3）增温期间日排放-处理水平；（4）增温期间累计排放-处理水平
  #（5）增温期间生长季排放；（6）观测期间生长季排放
  df = list(dat_fig0,cum_c1,dat_fig0.t,cum_pT,emit_sf1,emit_sf_c) %>%
    set_names(c('warming_daily chamber level','warming_cumulative chamber level',
                'warming_daily treat level','warming_cumulative treat level',
                'warming_vegetation emit_2021-2023','Obs_vegetation emit_2021-2023'))
  
  library(writexl)
  write_xlsx(df, "gasN emissions_chamber_2018-2023.xlsx")

  
############## Fig_S9 typical NO and N2O emissions in 2023 #############
  p1<-ggplot(F.pt[F.pt$treatment!='diff',],aes(date,soil_temp.m,col=treatment))+
    annotate("rect",xmin=as.Date("2023/11/1"),
             xmax=as.Date("2023/12/9"),
             ymin=-Inf,ymax=Inf,alpha=.25,fill="grey")+
    annotate("rect",xmin=as.Date("2023/3/21"),
             xmax=as.Date("2023/11/1"),
             ymin=-Inf,ymax=Inf,alpha=.1,fill="green")+
    annotate("rect",xmin=as.Date("2023/3/21"),
             xmax=as.Date("2023/5/1"),
             ymin=-Inf,ymax=Inf,alpha=.1,fill="blue")+
    geom_hline(yintercept = 0,lty=2,size=0.5,colour='grey')+
    geom_vline(xintercept = as.Date('2023/5/1'),lty=2,size=0.15,colour='black')+
    geom_vline(xintercept = as.Date('2023/3/21'),lty=2,size=0.15,colour='black')+
    geom_vline(xintercept = as.Date('2023/11/1'),lty=2,size=0.15,colour='black')+
    geom_vline(xintercept = as.Date('2023/12/9'),lty=2,size=0.15,colour='black')+
    scale_y_continuous(breaks = seq(-5,25,10),minor_breaks = seq(-5,25,5),expand=c(0,0))+
    geom_point(size=0.7)+geom_line(linewidth=0.3)+
    scale_x_date(limits = c(as.Date("2023/3/1"),as.Date("2024/1/1")),#date_breaks = 'month',date_labels = '%m/%d',
                 breaks = c(as.Date("2023/3/1"),as.Date("2023/5/1"),as.Date("2023/7/1"),
                            as.Date("2023/9/1"),as.Date("2023/11/1"),as.Date("2024/1/1")),
                 labels = c('3/1','5/1','7/1','9/1','11/1','1/1'), #minor_breaks = "days",
                 expand=c(0.005,0.005))+
    scale_color_manual(limits=c("control","warmed","Precipitation"),
                       labels=c("Control","Warmed","Precipitation"),
                       values = c("blue","red",'gray'),name='Variables')+
    labs(x=NULL,y=expression('Soil temp (\u00B0C)'))+
    guides(y="prism_offset_minor",x="prism_offset_minor")+mythem_sci+
    coord_cartesian(clip="off",ylim=c(-5,25))+
    geom_text(data=F.pt[F.pt$date=='2023/4/15' & F.pt$treatment=='control',],
              aes(x=as.Date('2023/4/10'),y=-2.5,
                  label = 'Spring \n freeze-thaw'),col='black',size=4.2)+
    geom_text(data=F.pt[F.pt$date=='2023/4/15' & F.pt$treatment=='control',],
              aes(x=as.Date('2023/7/12'),y=-3,label = 'Plant growing season'),col='black',size=4.2)+
    geom_text(data=F.pt[F.pt$date=='2023/4/15' & F.pt$treatment=='control',],
              aes(x=as.Date('2023/11/20'),y=-3,label = 'Winter'),col='black',size=4.2)+
    theme(axis.text.x = element_text(size=20,angle = 0,hjust=0.5,vjust=0.5),
          legend.text = element_text(size=12),legend.title = element_text(size=12,face='bold'),
          legend.position = 'top',
          plot.margin = unit(c(0.5,0,0.25,1),'cm'));p1
  
  str(F.pt)
  w1<-ggplot(F.pt[F.pt$treatment!='diff',],aes(date,WFPS.m,col=treatment))+
    annotate("rect",xmin=as.Date("2023/11/1"),
             xmax=as.Date("2023/12/9"),
             ymin=-Inf,ymax=Inf,alpha=.25,fill="grey")+
    annotate("rect",xmin=as.Date("2023/3/21"),
             xmax=as.Date("2023/11/1"),
             ymin=-Inf,ymax=Inf,alpha=.1,fill="green")+
    annotate("rect",xmin=as.Date("2023/3/21"),
             xmax=as.Date("2023/5/1"),
             ymin=-Inf,ymax=Inf,alpha=.1,fill="blue")+
    geom_bar(data=F.pt[-which(F.pt$precipitation==0),],aes(x=date,y=precipitation*2/3),
             position=position_dodge(),stat="identity",size=0.4,fill="grey",col="grey",
             alpha=0.5,show.legend = T)+
    geom_hline(yintercept = 0,lty=1,size=1,colour='black')+
    geom_vline(xintercept = as.Date('2023/5/1'),lty=2,size=0.15,colour='black')+
    geom_vline(xintercept = as.Date('2023/3/21'),lty=2,size=0.15,colour='black')+
    geom_vline(xintercept = as.Date('2023/11/1'),lty=2,size=0.15,colour='black')+
    geom_vline(xintercept = as.Date('2023/12/9'),lty=2,size=0.15,colour='black')+
    geom_point(size=0.7)+geom_line(linewidth=0.3)+
    scale_x_date(limits = c(as.Date("2023/3/1"),as.Date("2024/1/1")),#date_breaks = 'month',date_labels = '%m/%d',
                 breaks = c(as.Date("2023/3/1"),as.Date("2023/5/1"),as.Date("2023/7/1"),
                            as.Date("2023/9/1"),as.Date("2023/11/1"),as.Date("2024/1/1")),
                 labels = c('3/1','5/1','7/1','9/1','11/1','1/1'), #minor_breaks = "days",
                 expand=c(0.005,0.005))+
    scale_y_continuous(limits=c(0,60),breaks=seq(0,60,20),expand=c(0,0),
                       labels = c(seq(0,60, by=20)),
                       name = expression(paste("WFPS (%)")),
                       sec.axis=sec_axis(~.*3/2,breaks=seq(0,90,30),
                                         name=expression(paste("Precipitation (mm)"))))+
    scale_color_manual(limits=c("control","warmed"),values = c("blue","red"))+
    labs(x="Date",y=expression('WFPS (%)'))+
    guides(y="prism_offset_minor",x="prism_offset_minor")+mythem_sci+
    theme(axis.text.x = element_text(size=20,angle = 0,hjust=0.5,vjust=0.5),
          plot.margin = unit(c(0.25,0,0.25,1),'cm'));w1
  
  str(F.pt)
  NO1<-ggplot(F.pt[F.pt$treatment!='diff',],aes(date,NO_N.m,col=treatment))+
    annotate("rect",xmin=as.Date("2023/11/1"),
             xmax=as.Date("2023/12/9"),
             ymin=-Inf,ymax=Inf,alpha=.25,fill="grey")+
    annotate("rect",xmin=as.Date("2023/3/21"),
             xmax=as.Date("2023/11/1"),
             ymin=-Inf,ymax=Inf,alpha=.1,fill="green")+
    annotate("rect",xmin=as.Date("2023/3/21"),
             xmax=as.Date("2023/5/1"),
             ymin=-Inf,ymax=Inf,alpha=.1,fill="blue")+
    geom_vline(xintercept = as.Date('2023/5/1'),lty=2,size=0.15,colour='black')+
    geom_vline(xintercept = as.Date('2023/3/21'),lty=2,size=0.15,colour='black')+
    geom_vline(xintercept = as.Date('2023/11/1'),lty=2,size=0.15,colour='black')+
    geom_vline(xintercept = as.Date('2023/12/9'),lty=2,size=0.15,colour='black')+
    geom_point(size=0.7)+geom_line(linewidth=0.3)+
    scale_x_date(limits = c(as.Date("2023/3/1"),as.Date("2024/1/1")),#date_breaks = 'month',date_labels = '%m/%d',
                 breaks = c(as.Date("2023/3/1"),as.Date("2023/5/1"),as.Date("2023/7/1"),
                            as.Date("2023/9/1"),as.Date("2023/11/1"),as.Date("2024/1/1")),
                 labels = c('3/1','5/1','7/1','9/1','11/1','1/1'), #minor_breaks = "days",
                 expand=c(0.005,0.005))+
    scale_color_manual(limits=c("control","warmed"),values = c("blue","red"),name=NULL)+
    scale_y_continuous(limits = c(0,90),breaks = seq(0,90,30),minor_breaks = seq(0,90,5),expand=c(0,0),
                       name = expression(atop(paste("NO flux"),
                                              paste("(μg "," N"," m"^-2," h"^-1,")"))))+   
    labs(x="Date",y=expression('WFPS (%)'))+
    guides(y="prism_offset_minor",x="prism_offset_minor")+mythem_sci+
    theme(axis.text.x = element_text(size=20,angle = 0,hjust=0.5,vjust=0.5),legend.position = 'none',
          plot.margin = unit(c(0.25,0.2,0.25,0.2),'cm'));NO1
  
  
  n_breaks_i=NA
  for(i in -1:3){
    ni_breaks_i=seq(0.1,1,0.1)*10^i
    n_breaks_i=na.omit(unique(c(ni_breaks_i,n_breaks_i)))}
  print(n_breaks_i)
  
  ###n2o
  N2O1<-ggplot(F.pt[F.pt$treatment!='diff',],aes(date,N2O_N.m+10,col=treatment))+
    annotate("rect",xmin=as.Date("2023/11/1"),
             xmax=as.Date("2023/12/9"),
             ymin=1,ymax=1000,alpha=.25,fill="grey")+
    annotate("rect",xmin=as.Date("2023/3/21"),
             xmax=as.Date("2023/11/1"),
             ymin=1,ymax=1000,alpha=.1,fill="green")+
    annotate("rect",xmin=as.Date("2023/3/21"),
             xmax=as.Date("2023/5/1"),
             ymin=1,ymax=1000,alpha=.1,fill="blue")+
    geom_hline(yintercept = 10,lty=2,size=0.5,colour='grey')+
    geom_vline(xintercept = as.Date('2023/5/1'),lty=2,size=0.15,colour='black')+
    geom_vline(xintercept = as.Date('2023/3/21'),lty=2,size=0.15,colour='black')+
    geom_vline(xintercept = as.Date('2023/11/1'),lty=2,size=0.15,colour='black')+
    geom_vline(xintercept = as.Date('2023/12/9'),lty=2,size=0.15,colour='black')+
    geom_point(size=0.7)+geom_line(linewidth=0.3)+
    scale_x_date(limits = c(as.Date("2023/3/1"),as.Date("2024/1/1")),#date_breaks = 'month',date_labels = '%m/%d',
                 breaks = c(as.Date("2023/3/1"),as.Date("2023/5/1"),as.Date("2023/7/1"),
                            as.Date("2023/9/1"),as.Date("2023/11/1"),as.Date("2024/1/1")),
                 labels = c('23/3/1','23/5/1','23/7/1','23/9/1','23/11/1','24/1/1'), #minor_breaks = "days",
                 expand=c(0.005,0.005))+
    scale_color_manual(limits=c("control","warmed","diff"),values = c("blue","red","black"))+
    scale_y_log10(breaks=c(1,10,100,1000),limits=c(1,1000),labels=c(1,10,100,1000),
                  minor_breaks = n_breaks_i,guide = "prism_offset_minor",expand = c(0,0),
                  name = expression(atop(paste("N"[2],"O flux"),
                                         paste("(+10 μg "," N"," m"^-2," h"^-1,")"))))+   
    labs(x="Date (yy/mm/dd)")+
    guides(y="prism_offset_minor",x="prism_offset_minor")+mythem_sci+
    theme(axis.text.x = element_text(size=12,angle = 0,hjust=0.5,vjust=0.5),
          axis.title.x = element_text(size=12,hjust=0.5,vjust=0.5),
          plot.margin = unit(c(0.25,0.2,0.5,0.2),'cm'));N2O1
  
  figure<-((p1+rremove('x.text')+rremove('xlab')+ theme(plot.margin = unit(c(0,0.4,0,0.4),'cm')))/
             (w1+rremove('x.text')+rremove('xlab')+ theme(plot.margin = unit(c(0.35,0.4,0,0.4),'cm')))/
             (NO1+rremove('x.text')+rremove('xlab')+ theme(plot.margin = unit(c(0.35,0.4,0,0.4),'cm')))/
             (N2O1+ theme(plot.margin = unit(c(0.35,0.4,0,0.4),'cm')))) +plot_layout(ncol = 1,nrow=4)+
    # plot_annotation(tag_levels = 'A', tag_prefix = '(',tag_sep = '', tag_suffix = ')')&    #前缀，连接符，后缀
    plot_annotation(tag_levels = 'A')&    #前缀，连接符，后缀
    theme(axis.text.y = element_text(size=12,hjust=0.5,vjust=0.5),
          axis.title.y = element_text(size=12,hjust=0.5,vjust=0.5),
          plot.tag.position = c(0.16, 0.86), 
          plot.tag = element_text(size = 12,vjust = 0,hjust=0,face="bold"));figure
  
  substr(Sys.time(),1,10)
  
  ggsave(paste(paste("Fig_4_typical NO and N2O fluxes in 2023",
                     substr(Sys.time(),1,10),sep=""),".pdf",sep=""),figure,width =8, height =9)
  
  
  rm(p1,p1t,p2,p2t,p3,p3t,p4,pt,x1c,x1w,x1m,x2c,x2w,x2m,x3c,x3w,x3m,x4c,x4w,x4m,x5c,x5w,x5m,
     p,p1,p2,p3,p3d,p3u,output_emit_gc,output_emit_gc.E,dat.c,datai_vps,datai.c,emit_sf1,emit_sf1.range,
     emit_sf_a,emit_sf_c,emit_sf_yr,A_NO,A_N2O,f1,f2,f3,f4,figure1,
     t,treatments,figure,fit_all,fm1,fm2,fm3,lsd,lsdr, a,aov,c1,c2,c3,df,mod1,mod2,mod3,
     my_comparisons,p,p2_log,p2diff_log,p5,pco2,presult,pt0,res1,res2,res3,w1,w2,wco2,sm1,stat.test,mean_wp,mean_wp.m,
     emit_g, emit_gc, emit_gg, emit_gg.yr, emit_g.yr, emit_gt, emit_gt.yr, emit_plot.ab,dat,cum_c,
     i,j,k,ids,dates,y,nr_breaks,ni_breaks,nir_breaks,myname,mylabs,labs_warming,n_breaks,dat_key,
     sm2,sm3,sm4,stst.test,dat_1,dat_2,cum_c1,cum_cy,cum_p,cum_p1,cum_pL,cum_pT,cum_t,dat_ex,dat_ex.t,dat_fig0,dat_fig0.t,
     dat_fig1,dat_fig1.t,dat_fig2,dat_fig2.t,dat_figt,dat_i,dat_ijt,datai,dataii,dfa,dfa1,dfa2,dfa3,dfa5,mydata_loop,rainfall,result,result0)
  
  
  
  mydata_q10.res
  
  
###### <Main_Fig 5> ##########
  mydata_q10.plot<-mydata_d[which(mydata_d$treatment!="diff" & mydata_d$year %in% seq(2019,2023)),]
  
  mydata_q10_all<-mydata_q10.plot;mydata_q10_all$year='All'
  
  unique(mydata_q10_all$year)
  mydata_q10.plot<-rbind(mydata_q10.plot,mydata_q10_all);rm(mydata_q10_all)
  unique(mydata_q10.plot$year)
  mydata_q10.plot$lab_season<-ifelse(format(mydata_q10.plot$date,"%m") %in% c('05','06','07','08','09','10'),
                                     'Plant growing season','Plant dormant season')
  
  str(mydata_q10.plot)
  
  mydata_q10.plot$lab_wfps<-round(mydata_q10.plot$WFPS)  #四舍五入
  mydata_q10.plot$lab_soiltemp<-round(mydata_q10.plot$soil_temp/0.5)*0.5  #四舍五入
  
  
  ### label the outlier from 2019-2023
  dates<-mydata_d$date[which(mydata_d$precipitation>=5)]
  mydata_q10.plot$labs_t<-ifelse(mydata_q10.plot$date>="2021/7/20" & mydata_q10.plot$date<="2021/8/1","N",
                                 ifelse(mydata_q10.plot$date>="2020/11/20" & mydata_q10.plot$date<="2020/12/20","N",
                                        ifelse(mydata_q10.plot$date %in% dates,"N","Y")))
  
  
  
  unique(mydata_q10.plot$year)
  
  #mydata_q10.plot$warming=="Y" &
  datai<-mydata_q10.plot[which(mydata_q10.plot$warming=="Y" & mydata_q10.plot$year %in% seq(2019,2023)),]
  
  dataj=datai;dataj$year='All'
  datai<-rbind(datai,dataj);remove(dataj)
  
  ###标记出95%的置信区间
  datai$statlab_wfps=NA
  datai$statlab_temp=NA
  
  i=1
  for(i in 1:length(datai$date)){
    x=datai$treatment[i]
    y=datai$lab_season[i]
    xmin=mean(na.omit(datai$WFPS[datai$year=='All' & datai$lab_season==y & datai$treatment==x]))-
      1.96*sd(na.omit(datai$WFPS[datai$year=='All' & datai$lab_season==y & datai$treatment==x]))
    xmax=mean(na.omit(datai$WFPS[datai$year=='All' & datai$lab_season==y & datai$treatment==x]))+
      1.96*sd(na.omit(datai$WFPS[datai$year=='All' & datai$lab_season==y & datai$treatment==x]))
    
    datai$statlab_wfps[i]=ifelse(datai$WFPS[i] >= xmin & datai$WFPS[i]<= xmax,'T','F')
    
    ymin=mean(na.omit(datai$soil_temp[datai$year=='All' & datai$lab_season==y & datai$treatment==x]))-
      1.96*sd(na.omit(datai$soil_temp[datai$year=='All' & datai$lab_season==y & datai$treatment==x]))
    ymax=mean(na.omit(datai$soil_temp[datai$year=='All' & datai$lab_season==y & datai$treatment==x]))+
      1.96*sd(na.omit(datai$soil_temp[datai$year=='All' & datai$lab_season==y & datai$treatment==x]))
    datai$statlab_temp[i]=ifelse(datai$soil_temp[i] >= ymin & datai$soil_temp[i]<= ymax,'T','F')
    
    
    i=i+1
  }
  
  
  ###wfps
  unique(datai$year)
  mydata_q10.res<-summaryBy(data=datai[datai$statlab_wfps=='T',],FUN=myfun,
                            datai[datai$statlab_wfps=='T',c("N2O_N","NO_N","soil_temp")]~lab_wfps+treatment+plot+year+lab_season)
  
  unique(mydata_d$year)
  
  mydata_q10.plot<-mydata_d[which(mydata_d$treatment!="diff" & mydata_d$year %in% seq(2019,2023)),]
  
  mydata_q10_all<-mydata_q10.plot;mydata_q10_all$year='All'
  
  unique(mydata_q10_all$year)
  mydata_q10.plot<-rbind(mydata_q10.plot,mydata_q10_all);rm(mydata_q10_all)
  unique(mydata_q10.plot$year)
  mydata_q10.plot$lab_season<-ifelse(format(mydata_q10.plot$date,"%m") %in% c('05','06','07','08','09','10'),
                                     'Plant growing season','Plant dormant season')
  
  str(mydata_q10.plot)
  
  mydata_q10.plot$lab_wfps<-round(mydata_q10.plot$WFPS)  #四舍五入
  mydata_q10.plot$lab_soiltemp<-round(mydata_q10.plot$soil_temp/0.5)*0.5  #四舍五入
  
  
  ### label the outlier from 2019-2023
  dates<-mydata_d$date[which(mydata_d$precipitation>=5)]
  mydata_q10.plot$labs_t<-ifelse(mydata_q10.plot$date>="2021/7/20" & mydata_q10.plot$date<="2021/8/1","N",
                                 ifelse(mydata_q10.plot$date>="2020/11/20" & mydata_q10.plot$date<="2020/12/20","N",
                                        ifelse(mydata_q10.plot$date %in% dates,"N","Y")))
  
  
  
  unique(mydata_q10.plot$year)
  
  datai<-mydata_q10.plot[which(mydata_q10.plot$warming=="Y" & mydata_q10.plot$year %in% seq(2019,2023)),]
  dataj<-datai;dataj$year='All'
  datai<-rbind(datai,dataj);rm(dataj)
  # datai$year='All'
  ###分处理和季节， 标记出95%的温湿度置信区间
  datai$statlab_wfps=NA
  datai$statlab_temp=NA
  
  i=1
  for(i in 1:length(datai$date)){
    x=datai$treatment[i]
    y=datai$lab_season[i]
    yr=datai$year[i]
    xmin=mean(na.omit(datai$WFPS[datai$year==yr & datai$lab_season==y & datai$treatment==x]))-
      1.96*sd(na.omit(datai$WFPS[datai$year==yr & datai$lab_season==y & datai$treatment==x]))
    xmax=mean(na.omit(datai$WFPS[datai$year==yr & datai$lab_season==y & datai$treatment==x]))+
      1.96*sd(na.omit(datai$WFPS[datai$year==yr & datai$lab_season==y & datai$treatment==x]))
    
    datai$statlab_wfps[i]=ifelse(datai$WFPS[i] >= xmin & datai$WFPS[i]<= xmax,'T','F')
    
    ymin=mean(na.omit(datai$soil_temp[datai$year==yr & datai$lab_season==y & datai$treatment==x]))-
      1.96*sd(na.omit(datai$soil_temp[datai$year==yr & datai$lab_season==y & datai$treatment==x]))
    ymax=mean(na.omit(datai$soil_temp[datai$year==yr & datai$lab_season==y & datai$treatment==x]))+
      1.96*sd(na.omit(datai$soil_temp[datai$year==yr & datai$lab_season==y & datai$treatment==x]))
    datai$statlab_temp[i]=ifelse(datai$soil_temp[i] >= ymin & datai$soil_temp[i]<= ymax,'T','F')
    
    i=i+1
  }
  
  
  ###wfps______mydata_q10.res+mydata_q10.res.all
  unique(datai$year)
  mydata_q10.w<-summaryBy(data=datai[datai$statlab_wfps=='T' & datai$year=='All',],FUN=myfun,
                          N2O_N+NO_N+soil_temp~lab_wfps+treatment+plot+year+lab_season)
  
  unique(mydata_q10.w$year)
  
  mydata_q10.wa<-summaryBy(data=datai[datai$statlab_wfps=='T' & datai$year=='All',],FUN=myfun,
                                N2O_N+NO_N+soil_temp~lab_wfps+treatment+plot+year)
  
  mydata_q10.wa<-na.omit(mydata_q10.wa)   #20240126
  
  
  ###soiltemp____mydata_q10.t
  mydata_q10.t<-
    summaryBy(data=datai[datai$year=='All',],FUN=myfun,
              N2O_N+NO_N+soil_temp~lab_soiltemp+lab_season+treatment+plot+year)
  
  
  mydata_q10.t<-na.omit(mydata_q10.t)
  
  
  # preparation of mydata_q10.res2

  unique(mydata_q10.plot$lab_season)
  unique(mydata_q10.t$year)
  {pt1<-
      ggplot(mydata_q10.t[mydata_q10.t$year=='All' & mydata_q10.t$lab_season=="Plant growing season"& mydata_q10.t$NO_N.n>=3,],
             aes(lab_soiltemp,NO_N.m,col=treatment))+
      geom_point(data=mydata_q10.plot[mydata_q10.plot$year=='All' & mydata_q10.plot$warming=='Y' &
                                        mydata_q10.plot$lab_season=="Plant growing season",],
                 aes(soil_temp,NO_N,colour = treatment),
                 size=1.2,alpha=.1,stroke=0)+
      geom_point(size=1.4,alpha=.5)+ geom_point(size=1.5,col="black",pch=1,stroke=0.1,alpha=0.4)+
      scale_y_continuous(breaks=seq(0,90,30),limits = c(0,90),expand=c(0,0))+
      geom_smooth(aes(colour = treatment),se=F,linewidth=.6,formula=y ~ exp(0.13*x),method = 'lm')+
      stat_poly_eq(formula =log(y)~x,eq.x.rhs="x",coef.digits = 2,rr.digits = 2, 
                   eq.with.lhs = "italic(ln(y))~`=`~",   #给“y"换形式         
                   aes(lab_soiltemp,NO_N.m,col=treatment,label = paste(stat(eq.label))),
                   parse = TRUE,label.x.npc = 'left', hjust = 0,
                   label.y.npc = c(0.98,0.82),size = 3.5)+
      stat_poly_eq(formula =log(y)~x,eq.x.rhs="x",coef.digits = 2,rr.digits = 2, 
                   eq.with.lhs = "italic(ln(y))~`=`~",   #给“y"换形式         
                   aes(lab_soiltemp,NO_N.m,col=treatment,label =paste(stat(rr.label),stat(p.value.label),sep ="*\",\"~~")),
                   parse = TRUE,label.x.npc = 'left', hjust = 0,
                   label.y.npc = c(0.90,0.76),size = 3.5)+
      scale_x_continuous(breaks=seq(0,24,8),limits = c(0,24),expand=c(0,0),
                         guide = "prism_offset",minor_breaks = seq(0,24,1))+
      scale_colour_manual(limits=c("control","warmed"),values=c("blue","red"),
                          labels=c("control","warmed"),name=NULL)+
      scale_fill_manual(limits=c("control","warmed"),values=c("blue","red"),
                        labels=c("control","warmed"),name=NULL)+
      labs(x=expression("Soil temperature (\u00B0C)"),
           y=expression(paste("NO flux (μg"," N"," m"^-2," h"^-1,")")))+
      guides(x="prism_offset_minor",y="prism_offset_minor")+ #mythem+
      theme(plot.subtitle = element_text(vjust = 1),plot.caption = element_text(vjust = 1), 
            panel.grid.major = element_blank(), panel.grid=element_blank(), 
            panel.background = element_blank(),
            panel.border = element_rect(size=0.6,fill=NA),
            axis.text.x = element_text(colour = "black",size=12), 
            axis.text.y = element_text(colour = "black",size=12), 
            axis.title = element_text(colour = "black",size=12),
            axis.ticks = element_line(linewidth = .5),
            axis.ticks.length = unit(2,"mm"),
            prism.ticks.length = unit(1,"mm"),
            legend.title = element_text(colour = "black",size=12),
            legend.text = element_text(colour = "black",size=12),
            legend.key = element_rect(fill = NA,color=NA), 
            legend.background = element_rect(fill = NA), 
            legend.position =c(0.26,0.6),
            plot.background = element_blank(),
            plot.margin = unit(c(1,1,0,0),'lines'),
            panel.spacing.x = unit(5, 'mm'),
            strip.text = element_text(size=12, face="bold"),
            strip.background = element_rect(colour="black",fill=NA));pt1
    
    
    
  n_breaks_i=NA
  for(i in -1:3){
    ni_breaks_i=seq(0.1,1,0.1)*10^i
    n_breaks_i=na.omit(unique(c(ni_breaks_i,n_breaks_i)))}
  print(n_breaks_i)
  
    pt2<-
      ggplot(mydata_q10.t[mydata_q10.t$year=='All' & mydata_q10.t$lab_season=="Plant growing season" & mydata_q10.t$N2O_N.n>=3,],
             aes(lab_soiltemp,N2O_N.m,col=treatment))+
      # geom_vline(xintercept = 0,col='black',lty=2,size=0.35)+
      geom_vline(xintercept = 15,col='black',lty=2,size=0.35)+
      geom_point(data=mydata_q10.plot[mydata_q10.plot$year=='All' & mydata_q10.plot$warming=='Y' &
                                        mydata_q10.plot$lab_season=="Plant growing season",],
                 aes(soil_temp,N2O_N,colour = treatment),
                 size=1.2,alpha=.1,stroke=0)+
      geom_point(size=1.4,alpha=.5)+ geom_point(size=1.5,col="black",pch=1,stroke=0.1,alpha=0.4)+
      scale_y_log10(breaks=c(0.01,0.1,1,10,100,1000),limits=c(0.01,1000),labels=c(0.01,0.1,1,10,100,1000),
                    minor_breaks = n_breaks_i,guide = "prism_offset_minor",expand = c(0,0))+
      geom_smooth(data=mydata_q10.t[mydata_q10.t$year=='All' & mydata_q10.t$lab_season=="Plant growing season" &
                                      mydata_q10.t$lab_soiltemp<= 15,],
                  aes(colour = treatment),se=F,linewidth=.6,formula=y ~ x,method = 'lm')+
      geom_smooth(data=mydata_q10.t[mydata_q10.t$year=='All' & mydata_q10.t$lab_season=="Plant growing season" &
                                      mydata_q10.t$lab_soiltemp>= 15,],
                  aes(colour = treatment),se=F,linewidth=.6,formula=y ~ x,method = 'lm')+
      stat_poly_eq(data=mydata_q10.t[mydata_q10.t$year=='All' & mydata_q10.t$lab_season=="Plant growing season" &
                                       mydata_q10.t$lab_soiltemp<= 15& mydata_q10.t$N2O_N.n>=3,],
                   formula =y~x,eq.x.rhs="x",coef.digits = 2,rr.digits = 2, 
                   eq.with.lhs = "italic(ln(y))~`=`~",   #给“y"换形式         
                   aes(col=treatment,label = paste(stat(eq.label))),
                   parse = TRUE,label.x.npc = 'left', hjust = 0,
                   label.y.npc = c(0.98,0.82),size = 3.5)+
      stat_poly_eq(data=mydata_q10.t[mydata_q10.t$year=='All' & mydata_q10.t$lab_season=="Plant growing season" & 
                                       mydata_q10.t$lab_soiltemp<= 15& mydata_q10.t$N2O_N.n>=3,],
                   formula =y~x,eq.x.rhs="x",coef.digits = 2,rr.digits = 2, 
                   eq.with.lhs = "italic(ln(y))~`=`~",   #给“y"换形式         
                   aes(col=treatment,label =paste(stat(rr.label),stat(p.value.label),sep ="*\",\"~~")),
                   parse = TRUE,label.x.npc = 'left', hjust = 0,
                   label.y.npc = c(0.90,0.76),size = 3.5)+
      stat_poly_eq(data=mydata_q10.t[mydata_q10.t$year=='All' & mydata_q10.t$lab_season=="Plant growing season" & 
                                       mydata_q10.t$lab_soiltemp>= 15& mydata_q10.t$N2O_N.n>=3,],
                   formula =y~x,eq.x.rhs="x",coef.digits = 2,rr.digits = 2, 
                   eq.with.lhs = "italic(ln(y))~`=`~",   #给“y"换形式         
                   aes(col=treatment,label = paste(stat(eq.label))),
                   parse = TRUE,label.x.npc = 0.35, hjust = 0,
                   label.y.npc = c(0.28,0.12),size = 3.5)+
      stat_poly_eq(data=mydata_q10.t[mydata_q10.t$year=='All' & mydata_q10.t$lab_season=="Plant growing season" &
                                       mydata_q10.t$lab_soiltemp>= 15 & mydata_q10.t$N2O_N.n>=3,],
                   formula =y~x,eq.x.rhs="x",coef.digits = 2,rr.digits = 2, 
                   eq.with.lhs = "italic(ln(y))~`=`~",   #给“y"换形式         
                   aes(col=treatment,label =paste(stat(rr.label),stat(p.value.label),sep ="*\",\"~~")),
                   parse = TRUE,label.x.npc = 0.35, hjust = 0,
                   label.y.npc = c(0.20,0.04),size = 3.5)+
      scale_x_continuous(breaks=seq(0,24,8),limits = c(0,24),expand=c(0,0),
                         guide = "prism_offset",minor_breaks = seq(0,24,1))+
      scale_colour_manual(limits=c("control","warmed"),values=c("blue","red"),
                          labels=c("control","warmed"),name=NULL)+
      scale_fill_manual(limits=c("control","warmed"),values=c("blue","red"),
                        labels=c("control","warmed"),name=NULL)+
      labs(x=expression("Soil temperature (\u00B0C)"),
           y=expression(paste("N"[2],"O flux (μg"," N"," m"^-2," h"^-1,")")))+
      guides(x="prism_offset_minor",y="prism_offset_minor")+ #mythem+
      theme(plot.subtitle = element_text(vjust = 1),plot.caption = element_text(vjust = 1), 
            panel.grid.major = element_blank(), panel.grid=element_blank(), 
            panel.background = element_blank(),
            panel.border = element_rect(size=0.6,fill=NA),
            axis.text.x = element_text(colour = "black",size=12), 
            axis.text.y = element_text(colour = "black",size=12), 
            axis.title = element_text(colour = "black",size=12),
            axis.ticks = element_line(linewidth = .5),
            axis.ticks.length = unit(2,"mm"),
            prism.ticks.length = unit(1,"mm"),
            legend.title = element_text(colour = "black",size=12),
            legend.text = element_text(colour = "black",size=12),
            legend.key = element_rect(fill = NA), 
            legend.background = element_rect(fill = NA), 
            legend.position ='none',
            plot.background = element_blank(),
            plot.margin = unit(c(1,1,0,0),'lines'),
            panel.spacing.x = unit(5, 'mm'),
            strip.text = element_text(size=12, face="bold"),
            strip.background = element_rect(colour="black",fill=NA));pt2
    
    ####非生长季温度与NO和N2O
    pt1n<-
      ggplot(mydata_q10.t[mydata_q10.t$year=='All' & mydata_q10.t$lab_season=='Plant dormant season'& mydata_q10.t$NO_N.n>=3,],
             aes(lab_soiltemp,NO_N.m,col=treatment))+
      geom_vline(xintercept = 0,col='black',lty=2,size=0.35)+
      geom_point(data=mydata_q10.plot[mydata_q10.plot$year=='All' & mydata_q10.plot$warming=='Y' &
                                        mydata_q10.plot$lab_season=='Plant dormant season',],
                 aes(soil_temp,NO_N,colour = treatment),
                 size=1.2,alpha=.1,stroke=0)+
      geom_point(size=1.4,alpha=.5)+ geom_point(size=1.5,col="black",pch=1,stroke=0.1,alpha=0.4)+
      scale_y_continuous(breaks=seq(0,90,30),limits = c(0,90),expand=c(0,0))+
      geom_smooth(data=mydata_q10.t[mydata_q10.t$year=='All' & mydata_q10.t$lab_season=='Plant dormant season' & 
                                      mydata_q10.t$lab_soiltemp>= 0 & mydata_q10.t$NO_N.n>=3,],
                  aes(colour = treatment),se=F,linewidth=.6,formula=y ~ exp(0.10*x),method = 'lm')+
      stat_poly_eq(data=mydata_q10.t[mydata_q10.t$year=='All' & mydata_q10.t$lab_season=='Plant dormant season',],
                   formula =log(y)~x,eq.x.rhs="x",coef.digits = 2,rr.digits = 2, 
                   eq.with.lhs = "italic(ln(y))~`=`~",   #给“y"换形式         
                   aes(lab_soiltemp,NO_N.m,col=treatment,label = paste(stat(eq.label))),
                   parse = TRUE,label.x.npc = 0.34, hjust = 0,
                   label.y.npc = c(0.98,0.82),size = 3.5)+
      stat_poly_eq(data=mydata_q10.t[mydata_q10.t$year=='All' & mydata_q10.t$lab_season=='Plant dormant season',],
                   formula =log(y)~x,eq.x.rhs="x",coef.digits = 2,rr.digits = 2, 
                   eq.with.lhs = "italic(ln(y))~`=`~",   #给“y"换形式         
                   aes(lab_soiltemp,NO_N.m,col=treatment,label =paste(stat(rr.label),stat(p.value.label),sep ="*\",\"~~")),
                   parse = TRUE,label.x.npc = 0.34, hjust = 0,
                   label.y.npc = c(0.90,0.76),size = 3.5)+
      scale_x_continuous(breaks=seq(-5,10,5),limits = c(-5,10),expand=c(0,0),
                         guide = "prism_offset",minor_breaks = seq(-5,10,1))+
      scale_colour_manual(limits=c("control","warmed"),values=c("blue","red"),
                          labels=c("control","warmed"),name=NULL)+
      scale_fill_manual(limits=c("control","warmed"),values=c("blue","red"),
                        labels=c("control","warmed"),name=NULL)+
      labs(x=expression("Soil temperature (\u00B0C)"),
           y=expression(paste("NO flux (μg"," N"," m"^-2," h"^-1,")")))+
      guides(x="prism_offset_minor",y="prism_offset_minor")+ #mythem+
      theme(plot.subtitle = element_text(vjust = 1),plot.caption = element_text(vjust = 1), 
            panel.grid.major = element_blank(), panel.grid=element_blank(), 
            panel.background = element_blank(),
            panel.border = element_rect(size=0.6,fill=NA),
            axis.text.x = element_text(colour = "black",size=12), 
            axis.text.y = element_text(colour = "black",size=12), 
            axis.title = element_text(colour = "black",size=12),
            axis.ticks = element_line(linewidth = .5),
            axis.ticks.length = unit(2,"mm"),
            prism.ticks.length = unit(1,"mm"),
            legend.title = element_text(colour = "black",size=12),
            legend.text = element_text(colour = "black",size=12),
            legend.key = element_rect(fill = NA), 
            legend.background = element_rect(fill = NA), 
            legend.position ='none',
            plot.background = element_blank(),
            plot.margin = unit(c(1,1,0,0),'lines'),
            panel.spacing.x = unit(5, 'mm'),
            strip.text = element_text(size=12, face="bold"),
            strip.background = element_rect(colour="black",fill=NA));pt1n
    
    pt2n<-
      ggplot(mydata_q10.t[mydata_q10.t$year=='All' & mydata_q10.t$lab_season=='Plant dormant season' & mydata_q10.t$N2O_N.n>=3,],
             aes(lab_soiltemp,N2O_N.m,col=treatment))+
      # geom_line(data = dat_break,color = 'blue',fill = 'blue',lty=1,size=1.5)+
      geom_vline(xintercept = 0,col='black',lty=2,size=0.35)+
      geom_vline(xintercept = 15,col='black',lty=2,size=0.35)+
      
      geom_point(data=mydata_q10.plot[mydata_q10.plot$year=='All' & mydata_q10.plot$warming=='Y' &
                                        mydata_q10.plot$lab_season=='Plant dormant season',],
                 aes(soil_temp,N2O_N,colour = treatment),
                 size=1.2,alpha=.1,stroke=0)+
      geom_point(size=1.4,alpha=.5)+ geom_point(size=1.5,col="black",pch=1,stroke=0.1,alpha=0.4)+
      scale_y_log10(breaks=c(0.01,0.1,1,10,100,1000),limits=c(0.01,1000),labels=c(0.01,0.1,1,10,100,1000),
                    minor_breaks = n_breaks_i,guide = "prism_offset_minor",expand = c(0,0))+
      geom_smooth(data=mydata_q10.t[mydata_q10.t$year=='All' & mydata_q10.t$lab_season=='Plant dormant season' &
                                      mydata_q10.t$lab_soiltemp>= 0 & mydata_q10.t$N2O_N.n>=3,],
                  aes(colour = treatment),
                  se=F,linewidth=.6,formula=y ~ x+I(x^2),method = 'lm')+
      stat_poly_eq(data=mydata_q10.t[mydata_q10.t$year=='All' & mydata_q10.t$lab_season=='Plant dormant season'&
                                       mydata_q10.t$lab_soiltemp>= 0,],
                   formula =y ~ x+I(x^2),eq.x.rhs='x',coef.digits = 2,rr.digits = 2, 
                   eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
                   aes(lab_soiltemp,N2O_N.m,col=treatment,label = paste(stat(eq.label))),
                   parse = TRUE,label.x.npc = 0.34, hjust = 0,
                   label.y.npc = c(0.22,0.10),size = 3.5)+
      stat_poly_eq(data=mydata_q10.t[mydata_q10.t$year=='All' & mydata_q10.t$lab_season=='Plant dormant season'&
                                       mydata_q10.t$lab_soiltemp>= 0,],
                   formula =y ~ x+I(x^2),eq.x.rhs="x",coef.digits = 2,rr.digits = 2, 
                   eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
                   aes(lab_soiltemp,N2O_N.m,col=treatment,label =paste(stat(rr.label),stat(p.value.label),sep ="*\",\"~~")),
                   parse = TRUE,label.x.npc = 0.34, hjust = 0,
                   label.y.npc = c(0.16,0.04),size = 3.5)+
      scale_x_continuous(breaks=seq(-5,10,5),limits = c(-5,10),expand=c(0,0),
                         guide = "prism_offset",minor_breaks = seq(-5,10,1))+
      scale_colour_manual(limits=c("control","warmed"),values=c("blue","red"),
                          labels=c("control","warmed"),name=NULL)+
      scale_fill_manual(limits=c("control","warmed"),values=c("blue","red"),
                        labels=c("control","warmed"),name=NULL)+
      guides(x="prism_offset_minor",y="prism_offset_minor")+
      labs(x=expression("Soil temperature (\u00B0C)"),
           y=expression(paste("N"[2],"O flux (μg"," N"," m"^-2," h"^-1,")")))+ # mythem+
      theme(plot.subtitle = element_text(vjust = 1),plot.caption = element_text(vjust = 1), 
            panel.grid.major = element_blank(), panel.grid=element_blank(), 
            panel.background = element_blank(),
            panel.border = element_rect(size=0.6,fill=NA),
            axis.text.x = element_text(colour = "black",size=12), 
            axis.text.y = element_text(colour = "black",size=12), 
            axis.title = element_text(colour = "black",size=12),
            axis.ticks = element_line(linewidth = .5),
            axis.ticks.length = unit(2,"mm"),
            prism.ticks.length = unit(1,"mm"),
            legend.title = element_text(colour = "black",size=12),
            legend.text = element_text(colour = "black",size=12),
            legend.key = element_rect(fill = NA), 
            legend.background = element_rect(fill = NA), 
            legend.position ='none',
            plot.background = element_blank(),
            plot.margin = unit(c(1,1,0,0),'lines'),
            panel.spacing.x = unit(5, 'mm'),
            strip.text = element_text(size=12, face="bold"),
            strip.background = element_rect(colour="black",fill=NA));pt2n
   }
    
  ### 生长季NXO与湿度的关系 
  {pw1<-ggplot(mydata_q10.w[mydata_q10.w$year=='All' & mydata_q10.w$lab_season=='Plant growing season' & mydata_q10.w$NO_N.n>=3,],
               aes(lab_wfps,NO_N.m,col=treatment))+
      geom_point(data=datai[datai$statlab_wfps=='T' & datai$year=='All' & datai$labs_t=='Y' &
                              datai$lab_season=='Plant growing season',],
                 aes(WFPS,NO_N,colour = treatment),
                 size=1.2,alpha=.1,stroke=0)+
      geom_point(size=1.4,alpha=.5)+ geom_point(size=1.5,col="black",pch=1,stroke=0.1,alpha=0.4)+
      geom_smooth(aes(lab_wfps,NO_N.m),method="lm",se=F,linewidth=.6,
                  formula=y~x+I(x^2),show.legend=T)+
      scale_x_continuous(breaks=seq(15,55,10),limits = c(15,55),minor_breaks = seq(15,55,5),expand=c(0,0))+
      scale_y_continuous(breaks=seq(0,90,30),limits = c(0,90),expand=c(0,0))+
      scale_colour_manual(limits=c("control","warmed"),values=c("blue","red"),
                          labels=c("Control","Warmed"),name=NULL)+
      scale_fill_manual(limits=c("control","warmed"),values=c("blue","red"),
                        labels=c("Control","Warmed"),name=NULL)+
      stat_poly_eq(formula =y~x+I(x^2),coef.digits = 2,eq.x.rhs="x", vstep=0.11,
                   eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
                   aes(label = paste(after_stat(eq.label),sep ="*\",\"~~")),
                   label.x.npc = "left", label.y.npc = c(0.98,0.82),show.legend=T,
                   parse = TRUE,size = 3.5)+
      stat_poly_eq(formula =y~x+I(x^2),coef.digits = 2,eq.x.rhs="x", vstep=0.11,
                   eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
                   aes(label = paste(..rr.label..,..p.value.label..,sep ="*\",\"~~")),
                   label.x.npc = "left", label.y.npc = c(0.90,0.76),show.legend=T,
                   parse = T,size = 3.5)+
      labs(x=expression("WFPS (%)"),
           y=expression(paste("NO flux (μg"," N"," m"^-2," h"^-1,")")))+
      guides(x="prism_offset_minor",y="prism_offset_minor")+ #mythem+
      theme(plot.subtitle = element_text(vjust = 1),plot.caption = element_text(vjust = 1), 
            panel.grid.major = element_blank(), panel.grid=element_blank(), 
            panel.background = element_blank(),
            panel.border = element_rect(size=0.6,fill=NA),
            axis.text.x = element_text(colour = "black",size=12), 
            axis.text.y = element_text(colour = "black",size=12), 
            axis.title = element_text(colour = "black",size=12),
            axis.ticks = element_line(linewidth = .5),
            axis.ticks.length = unit(2,"mm"),
            prism.ticks.length = unit(1,"mm"),
            legend.title = element_text(colour = "black",size=12),
            legend.text = element_text(colour = "black",size=12),
            legend.key = element_rect(fill = NA), 
            legend.background = element_rect(fill = NA), 
            legend.position ="none",
            plot.background = element_blank(),
            plot.margin = unit(c(2,1,2,1),'lines'),
            panel.spacing.x = unit(5, 'mm'),
            strip.text = element_text(size=12, face="bold"),
            strip.background = element_rect(colour="black",fill=NA));pw1
    
    pw2<-
      ggplot(mydata_q10.w[mydata_q10.w$year=='All' & mydata_q10.w$lab_season=='Plant growing season' & mydata_q10.w$N2O_N.n>=3,],
             aes(lab_wfps,N2O_N.m,col=treatment))+
      geom_point(data=datai[datai$statlab_wfps=='T' & datai$year=='All' & datai$labs_t=='Y' &
                              datai$lab_season=='Plant growing season',],
                 aes(WFPS,N2O_N,colour = treatment),
                 size=1.2,alpha=.1,stroke=0)+
      geom_point(size=1.4,alpha=.5)+ geom_point(size=1.5,col="black",pch=1,stroke=0.1,alpha=0.4)+
      geom_smooth(linewidth=.6,formula=y~x,method = 'lm',se=F)+
      scale_x_continuous(breaks=seq(15,55,10),limits = c(15,55),minor_breaks = seq(15,55,5),expand=c(0,0))+
      # scale_y_continuous(breaks=seq(-10,110,30),limits = c(-10,110),expand=c(0,0))+
      scale_y_log10(breaks=c(0.01,0.1,1,10,100,1000),limits=c(0.01,1000),labels=c(0.01,0.1,1,10,100,1000),
                    minor_breaks = n_breaks_i,guide = "prism_offset_minor",expand = c(0,0))+
      scale_colour_manual(limits=c("control","warmed"),values=c("blue","red"),
                          labels=c("Control","Warmed"),name=NULL)+
      scale_fill_manual(limits=c("control","warmed"),values=c("blue","red"),
                        labels=c("Control","Warmed"),name=NULL)+
      stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", vstep=0.11,
                   eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
                   aes(label = paste(after_stat(eq.label),sep ="*\",\"~~")),
                   label.x.npc = "left", label.y.npc = c(0.98,0.82),show.legend=T,
                   parse = TRUE,size = 3.5)+
      stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", vstep=0.11,
                   eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
                   aes(label = paste(..rr.label..,..p.value.label..,sep ="*\",\"~~")),
                   label.x.npc = "left", label.y.npc = c(0.90,0.76),show.legend=T,
                   parse = T,size = 3.5)+
      labs(x=expression("WFPS (%)"),
           y=expression(paste("N"[2],"O flux (μg"," N"," m"^-2," h"^-1,")")))+
      guides(x="prism_offset_minor",y="prism_offset_minor")+ #mythem+
      theme(plot.subtitle = element_text(vjust = 1),plot.caption = element_text(vjust = 1), 
            panel.grid.major = element_blank(), panel.grid=element_blank(), 
            panel.background = element_blank(),
            panel.border = element_rect(size=0.6,fill=NA),
            axis.text.x = element_text(colour = "black",size=12), 
            axis.text.y = element_text(colour = "black",size=12), 
            axis.title = element_text(colour = "black",size=12),
            axis.ticks = element_line(linewidth = .5),
            axis.ticks.length = unit(2,"mm"),
            prism.ticks.length = unit(1,"mm"),
            legend.title = element_text(colour = "black",size=12),
            legend.text = element_text(colour = "black",size=12),
            legend.key = element_rect(fill = NA), 
            legend.background = element_rect(fill = NA), 
            legend.position ="none",
            plot.background = element_blank(),
            plot.margin = unit(c(0,1,1,1),'lines'),
            panel.spacing.x = unit(5, 'mm'),
            strip.text = element_text(size=12, face="bold"),
            strip.background = element_rect(colour="black",fill=NA));pw2
    
    #非生长季NXO与湿度的关系 
    pw1n<-
      ggplot(mydata_q10.w[mydata_q10.w$year=='All' & mydata_q10.w$lab_season=='Plant dormant season' & mydata_q10.w$NO_N.n>=3,],
             aes(lab_wfps,NO_N.m,col=treatment))+
      geom_point(data=datai[datai$statlab_wfps=='T' & datai$year=='All' & datai$labs_t=='Y' &
                              datai$lab_season=='Plant dormant season',],
                 aes(WFPS,NO_N,colour = treatment),
                 size=1.2,alpha=.1,stroke=0)+
      # geom_label(aes(label=season))+
      geom_point(size=1.4,alpha=.5)+ geom_point(size=1.5,col="black",pch=1,stroke=0.1,alpha=0.4)+
      geom_smooth(aes(lab_wfps,NO_N.m),method="lm",se=F,linewidth=.8,
                  formula=y ~ x+I(x^2),show.legend=T)+
      scale_x_continuous(breaks=seq(15,55,10),limits = c(15,55),minor_breaks = seq(15,55,5),expand=c(0,0))+
      scale_y_continuous(breaks=seq(0,90,30),limits = c(0,90),expand=c(0,0))+
      scale_colour_manual(limits=c("control","warmed"),values=c("blue","red"),
                          labels=c("Control","Warmed"),name=NULL)+
      scale_fill_manual(limits=c("control","warmed"),values=c("blue","red"),
                        labels=c("Control","Warmed"),name=NULL)+
      stat_poly_eq(formula =y~x+I(x^2),coef.digits = 2,eq.x.rhs="x", vstep=0.11,
                   eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
                   aes(label = paste(after_stat(eq.label),sep ="*\",\"~~")),
                   label.x.npc = "left", label.y.npc = c(0.98,0.83),show.legend=T,
                   parse = TRUE,size = 3.5)+
      stat_poly_eq(formula =y~x+I(x^2),coef.digits = 2,eq.x.rhs="x", vstep=0.11,
                   eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
                   aes(label = paste(..rr.label..,..p.value.label..,sep ="*\",\"~~")),
                   label.x.npc = "left", label.y.npc = c(0.90,0.76),show.legend=T,
                   parse = T,size = 3.5)+
      labs(x=expression("WFPS (%)"),
           y=expression(paste("NO flux (μg"," N"," m"^-2," h"^-1,")")))+
      guides(x="prism_offset_minor",y="prism_offset_minor")+ #mythem+
      theme(plot.subtitle = element_text(vjust = 1),plot.caption = element_text(vjust = 1), 
            panel.grid.major = element_blank(), panel.grid=element_blank(), 
            panel.background = element_blank(),
            panel.border = element_rect(size=0.6,fill=NA),
            axis.text.x = element_text(colour = "black",size=12), 
            axis.text.y = element_text(colour = "black",size=12), 
            axis.title = element_text(colour = "black",size=12),
            axis.ticks = element_line(linewidth = .5),
            axis.ticks.length = unit(2,"mm"),
            prism.ticks.length = unit(1,"mm"),
            legend.title = element_text(colour = "black",size=12),
            legend.text = element_text(colour = "black",size=12),
            legend.key = element_rect(fill = NA), 
            legend.background = element_rect(fill = NA), 
            legend.position ="none",
            plot.background = element_blank(),
            plot.margin = unit(c(2,1,2,1),'lines'),
            panel.spacing.x = unit(5, 'mm'),
            strip.text = element_text(size=12, face="bold"),
            strip.background = element_rect(colour="black",fill=NA));pw1n
    
    pw2n<-
      ggplot(mydata_q10.w[mydata_q10.w$year=='All' & mydata_q10.w$lab_season=='Plant dormant season' & mydata_q10.w$N2O_N.n>=3,],
             aes(lab_wfps,N2O_N.m,col=treatment))+
      geom_point(data=datai[datai$statlab_wfps=='T' & datai$year=='All' & datai$labs_t=='Y' &
                              datai$lab_season=='Plant dormant season',],
                 aes(WFPS,N2O_N,colour = treatment),size=1.2,alpha=.1,stroke=0)+
      geom_point(size=1,alpha=.5)+ geom_point(size=1.2,col="black",pch=1,stroke=0.2,alpha=0.4)+
      geom_smooth(linewidth=.8,formula=y~x,method = 'lm',se=F)+
      scale_x_continuous(breaks=seq(15,55,10),limits = c(15,55),minor_breaks = seq(15,55,5),expand=c(0,0))+
      scale_y_log10(breaks=c(0.1,1,10,100,1000),limits=c(0.1,1000),labels=c(0.1,1,10,100,1000),
                    minor_breaks = n_breaks_i,guide = "prism_offset_minor",expand = c(0,0))+
      scale_colour_manual(limits=c("control","warmed"),values=c("blue","red"),
                          labels=c("Control","Warmed"),name=NULL)+
      scale_fill_manual(limits=c("control","warmed"),values=c("blue","red"),
                        labels=c("Control","Warmed"),name=NULL)+
      stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", vstep=0.11,
                   eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
                   aes(label = paste(after_stat(eq.label),sep ="*\",\"~~")),hjust=0,
                   label.x.npc = "left", label.y.npc = c(0.22,0.10),show.legend=T,
                   parse = TRUE,size = 3.5)+
      stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", vstep=0.11,
                   eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
                   aes(label = paste(..rr.label..,..p.value.label..,sep ="*\",\"~~")),hjust=0,
                   label.x.npc = "left", label.y.npc = c(0.16,0.04),show.legend=T,
                   parse = T,size = 3.5)+
      labs(x=expression("WFPS (%)"),
           y=expression(paste("N"[2],"O flux (μg"," N"," m"^-2,"h"^-1,")")))+
      guides(x="prism_offset_minor",y="prism_offset_minor")+ #mythem+
      theme(plot.subtitle = element_text(vjust = 1),plot.caption = element_text(vjust = 1), 
            panel.grid.major = element_blank(), panel.grid=element_blank(), 
            panel.background = element_blank(),
            panel.border = element_rect(size=0.6,fill=NA),
            axis.text.x = element_text(colour = "black",size=12), 
            axis.text.y = element_text(colour = "black",size=12), 
            axis.title = element_text(colour = "black",size=12),
            axis.ticks = element_line(linewidth = .5),
            axis.ticks.length = unit(2,"mm"),
            prism.ticks.length = unit(1,"mm"),
            legend.title = element_text(colour = "black",size=12),
            legend.text = element_text(colour = "black",size=12),
            legend.key = element_rect(fill = NA), 
            legend.background = element_rect(fill = NA), 
            legend.position ="none",
            plot.background = element_blank(),
            plot.margin = unit(c(0,1,1,1),'lines'),
            panel.spacing.x = unit(5, 'mm'),
            strip.text = element_text(size=12, face="bold"),
            strip.background = element_rect(colour="black",fill=NA));pw2n 
  }
  
  figure1<-ggarrange(pt1+theme(plot.margin = unit(c(0.5,0.5,0.2,0.5),'lines'),
                               axis.text.x = element_blank(),axis.title.x = element_blank())+
                       facet_wrap('Plant growing season'~.),
                     pt2+theme(plot.margin = unit(c(0.5,0.5,0.5,0.5),'lines')),
                     labels = c("A", ""),
                     label.x = 0.01,label.y = c(1.0,0.99),
                     ncol = 1, nrow = 2,align = 'v',   ##"v"竖直对齐
                     font.label = list(size = 12, color ="black"),
                     widths = c(6), heights = c(4.5,5),
                     common.legend=F);figure1
  
  figure2<-ggarrange(pt1n+rremove("ylab")+rremove("y.text")+
                       theme(plot.margin = unit(c(0.5,0.5,0.2,0.5),'lines'),
                             axis.text.x = element_blank(),axis.title.x = element_blank())+
                       facet_wrap('Plant dormant season'~.),
                     pt2n+rremove("ylab")+rremove("y.text")+
                       theme(plot.margin = unit(c(0.5,0.5,0.5,0.5),'lines')),
                     labels = c(" ", " "),
                     label.x = 0.87,label.y = c(0.86,0.95),
                     ncol = 1, nrow = 2,align = 'v',   ##"v"竖直对齐
                     font.label = list(size = 12, color ="black"),
                     widths = c(6), heights = c(4.5,5),
                     common.legend=F);figure2
  
  figure3<-ggarrange(pw1+rremove("ylab")+rremove("y.text")+
                       theme(plot.margin = unit(c(0.5,0.5,0.2,2),'lines'),
                             axis.text.x = element_blank(),axis.title.x = element_blank())+
                       facet_wrap('Plant growing season'~.),
                     pw2+rremove("ylab")+rremove("y.text")+
                       theme(plot.margin = unit(c(0.5,0.5,0.5,2),'lines')),
                     labels = c("B",""),
                     label.x = 0,label.y = c(1.0,0.85),
                     ncol = 1, nrow = 2,align = 'v',   ##"v"竖直对齐
                     font.label = list(size = 12),
                     widths = c(6), heights = c(4.5,5),
                     common.legend=F);figure3
  
  figure4<-ggarrange(pw1n+rremove("ylab")+rremove("y.text")+
                       theme(plot.margin = unit(c(0.5,3,0.2,0.5),'lines'),
                             axis.text.x = element_blank(),axis.title.x = element_blank())+
                       facet_wrap('Plant dormant season'~.),
                     pw2n+rremove("ylab")+rremove("y.text")+
                       theme(plot.margin = unit(c(0.5,3,0.5,0.5),'lines')),
                     labels = c(" "," "),
                     label.x = 0.87,label.y = c(0.86,0.95),
                     ncol = 1, nrow = 2,align = 'v',   ##"v"竖直对齐
                     font.label = list(size = 12, color ="black"),
                     widths = c(6), heights = c(4.5,5),
                     common.legend=F);figure4
  
  figure<-ggarrange(figure1,figure2,figure3,figure4,
                    ncol = 4, nrow = 1,align = 'none',   ##"v"竖直对齐
                    widths = c(6,4.6,5.1,5.4), heights = c(5),
                    common.legend=F);figure
  
  
  
  ### lnRR of NO + N2O
  dat<-warm_eff.m[which(warm_eff.m$warming=="1" & 
                          warm_eff.m$year %in% c('2023','2022','2021','2020','2019')),]
  
  dat$warm_WFPS<-dat$warm_sm/(1-0.7/2.65)
  
  dat$season<-ifelse(dat$month %in% c('03','04','05'),'spring',
                     ifelse(dat$month %in% c('06','07','08'),'summer',
                            ifelse(dat$month %in% c('09','10','11'),'fall','winter')))
  
  
  dat<-within(dat,season<-factor(season,levels=c("spring","summer","fall","winter")))
  
  ###
  # dat<-dat[dat$warm_WFPS>= -10 & dat$warm_WFPS<=5,]
  
  dat<-dat[dat$season!='winter',]   #冬季数据较少
  
  dat$WFPS<-dat$soil_moisture/(1-0.7/2.65)
  
  ### ln RR
  {a0<-ggplot(dat[dat$month %in% c('06','07','08'),],
              aes(WFPS,log(eff_NO/100+1)))+
      geom_hline(yintercept =0,size =0.4,color="red",lty=2)+
      # geom_vline(xintercept =0,size =0.4,color="red",lty=2)+
      scale_x_continuous(limits=c(15,55),breaks=seq(15,55,10),expand=c(0,0))+
      scale_y_continuous(limits=c(-1,1),breaks=seq(-1,1,1),expand=c(0,0))+
      guides(y="prism_offset_minor",x="prism_offset_minor")+
      # geom_pointdensity(aes(WFPS,log(eff_NO/100+1)),size=1.5)+
      # scale_color_viridis() +
      geom_point(size=1.5,col='green4',alpha=0.36)+
      geom_point(size=1.5,pch=1,alpha=0.68)+ facet_grid(.~"NO", labeller=label_parsed)+
      geom_smooth(method="lm",col='transparent',formula = y~x,se=T,alpha=0.23,show.legend = T)+
      geom_smooth(method="lm",formula = y~x,se=F,show.legend = T,col='green3')+   
      stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x",col='green4',
                   eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式
                   aes(label = ifelse(..p.value..<0.05,
                                      paste(..rr.label..,..p.value.label..,sep ="*\",\"~~"),'ns')),     #..eq.label..,
                   parse = TRUE,label.x.npc = 0.02, label.y.npc = c(0.933),size = 3.5)+  #vstep=0.08,
      mythem+
      theme(legend.title = element_blank(),
            legend.position ="none",
            legend.text = element_text(colour = "black",size=6,vjust=0.5),
            legend.background = element_rect(colour = "black",size=0.1),
            panel.border = element_blank(),
            panel.background = element_rect(fill =NA, linewidth=0.5,colour = "black", linetype = "solid",inherit.blank=T),
            plot.background = element_blank(),
            axis.title.y = element_text(colour = "black",size=12,vjust=2),
            axis.title.x = element_text(colour = "black",size=12,vjust=0.5),
            axis.text.x = element_text(colour = "black",size=12,vjust=0.5),
            axis.text.y = element_text(colour = "black",size=12,vjust=0.5),
            plot.margin = unit(c(0.5,0.25,0.2,1.8),'lines'))+
      coord_cartesian(clip="off")+
      labs(x=expression("WFPS (%)"),
           y=expression(paste("lnRR of NO flux")));a0
    
    b0<-ggplot(dat[dat$month %in% c('06','07','08') & dat$N2O_N>0 ,],
               aes(WFPS,log(eff_N2O/100+1)))+
      geom_hline(yintercept =0,size =0.4,color="red",lty=2)+
      # geom_vline(xintercept =0,size =0.4,color="red",lty=2)+
      scale_x_continuous(limits=c(15,55),breaks=seq(15,55,10),expand=c(0,0))+
      scale_y_continuous(limits=c(-4,4),breaks=seq(-4,4,4),expand=c(0,0))+
      guides(y="prism_offset_minor",x="prism_offset_minor")+
      # geom_pointdensity(size=1.5)+
      # scale_color_viridis() +
      geom_point(size=1.5,col='green4',alpha=0.36)+
      geom_point(size=1.5,pch=1,alpha=0.68)+ facet_grid(.~"N[2]*O", labeller=label_parsed)+
      geom_smooth(method="lm",col='transparent',formula = y~x,se=T,alpha=0.23,show.legend = T)+
      geom_smooth(method="lm",formula = y~x,se=F,show.legend = T,col='green3')+   
      stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x",col='green4',
                   eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式
                   aes(label = ifelse(..p.value..<0.05,
                                      paste(..rr.label..,..p.value.label..,sep ="*\",\"~~"),'ns')),     #..eq.label..,
                   parse = TRUE,label.x.npc = 0.02, label.y.npc = c(0.933),size = 3.5)+  #vstep=0.08,
      mythem+
      theme(legend.title = element_blank(),
            legend.position ="none",
            legend.text = element_text(colour = "black",size=6,vjust=0.5),
            legend.background = element_rect(colour = "black",size=0.1),
            panel.border = element_blank(),
            plot.background = element_blank(),panel.background = element_rect(fill =NA, linewidth=0.5,colour = "black", linetype = "solid",inherit.blank=T),
            axis.title.y = element_text(colour = "black",size=12,vjust=-1),
            axis.title.x = element_text(colour = "black",size=12,vjust=0.5),
            axis.text.x = element_text(colour = "black",size=12,vjust=0.5),
            axis.text.y = element_text(colour = "black",size=12,vjust=0.5),
            plot.margin = unit(c(0.5,0.7,0.2,0.15),'lines'))+
      coord_cartesian(clip="off")+
      labs(x=expression("WFPS (%)"),
           y=expression(paste("lnRR of N"[2],"O flux")));b0
    
    a1<-ggplot(dat[dat$month %in% c('06','07','08'),],
               aes(warm_sm.perc*100,log(eff_NO/100+1)))+
      geom_hline(yintercept =0,size =0.4,color="red",lty=2)+
      # geom_vline(xintercept =0,size =0.4,color="red",lty=2)+
      scale_x_continuous(limits=c(-24,8),breaks=seq(-24,8,8),expand=c(0,0))+
      scale_y_continuous(limits=c(-1,1),breaks=seq(-1,1,1),expand=c(0,0))+
      guides(y="prism_offset_minor",x="prism_offset_minor")+
      geom_point(size=1.5,col='green4',alpha=0.36)+
      geom_point(size=1.5,pch=1,alpha=0.68)+ facet_grid(.~"NO", labeller=label_parsed)+
      geom_smooth(method="lm",col='transparent',formula = y~x,se=T,alpha=0.23,show.legend = T)+
      geom_smooth(method="lm",formula = y~x,se=F,show.legend = T,col='green3')+   
      stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x",col='green4',
                   eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式
                   aes(label = ifelse(..p.value..<0.05,
                                      paste(..rr.label..,..p.value.label..,sep ="*\",\"~~"),'ns')),     #..eq.label..,
                   parse = TRUE,label.x.npc = 0.02, label.y.npc = c(0.933),size = 3.5)+  #vstep=0.08,
      mythem+
      theme(legend.title = element_blank(),
            legend.position ="none",
            legend.text = element_text(colour = "black",size=6,vjust=0.5),
            legend.background = element_rect(colour = "black",size=0.1),
            panel.border = element_blank(), panel.background = element_rect(fill =NA, linewidth=0.5,colour = "black", linetype = "solid",inherit.blank=T),
            plot.background = element_blank(),
            axis.title.y = element_text(colour = "black",size=12,vjust=0),
            axis.title.x = element_text(colour = "black",size=12,vjust=0.5),
            axis.text.x = element_text(colour = "black",size=12,vjust=0.5),
            axis.text.y = element_text(colour = "black",size=12,vjust=0.5),
            plot.margin = unit(c(0.5,0.2,0.2,0.8),'lines'))+
      coord_cartesian(clip="off")+
      labs(x=expression(paste(Delta,"WFPS (%)")),
           y=expression(paste("lnRR of NO flux")));a1
    
    b1<-ggplot(dat[dat$month %in% c('06','07','08') & dat$N2O_N>0,],
               aes(warm_sm.perc*100,log(eff_N2O/100+1)))+
      geom_hline(yintercept =0,size =0.4,color="red",lty=2)+
      # geom_vline(xintercept =0,size =0.4,color="red",lty=2)+
      scale_x_continuous(limits=c(-24,8),breaks=seq(-24,8,8),expand=c(0,0))+
      scale_y_continuous(limits=c(-4,4),breaks=seq(-4,4,4),expand=c(0,0))+
      guides(y="prism_offset_minor",x="prism_offset_minor")+
      geom_point(size=1.5,col='green4',alpha=0.36)+
      geom_point(size=1.5,pch=1,alpha=0.68)+ facet_grid(.~"N[2]*O", labeller=label_parsed)+
      geom_smooth(method="lm",col='transparent',formula = y~x,se=T,alpha=0.23,show.legend = T)+
      geom_smooth(method="lm",formula = y~x,se=F,show.legend = T,col='green3')+   
      stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x",col='green4',
                   eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式
                   aes(label = ifelse(..p.value..<0.05,
                                      paste(..rr.label..,..p.value.label..,sep ="*\",\"~~"),'ns')),     #..eq.label..,
                   parse = TRUE,label.x.npc = 0.02, label.y.npc = c(0.933),size = 3.5)+  #vstep=0.08,
      mythem+
      theme(legend.title = element_blank(),
            legend.position ="none",
            legend.text = element_text(colour = "black",size=6,vjust=0.5),
            legend.background = element_rect(colour = "black",size=0.1),
            panel.border = element_blank(), panel.background = element_rect(fill =NA, linewidth=0.5,colour = "black", linetype = "solid",inherit.blank=T),
            plot.background = element_blank(),
            axis.title.y = element_text(colour = "black",size=12,vjust=0),
            axis.title.x = element_text(colour = "black",size=12,vjust=0.5),
            axis.text.x = element_text(colour = "black",size=12,vjust=0.5),
            axis.text.y = element_text(colour = "black",size=12,vjust=0.5),
            plot.margin = unit(c(0.5,2.8,0.2,0.45),'lines'))+
      coord_cartesian(clip="off")+
      labs(x=expression(paste(Delta,"WFPS (%)")),
           y=expression(paste("lnRR of N"[2],"O flux")));b1
  }

  p_new<-ggarrange(a0+theme(strip.text = element_blank(),strip.background = element_blank()),
                   a1+rremove('ylab')+rremove('y.text')+
                     theme(strip.text = element_blank(),strip.background = element_blank()),
                   b0+theme(strip.text = element_blank(),strip.background = element_blank()),
                   b1+rremove('ylab')+rremove('y.text')+
                     theme(strip.text = element_blank(),strip.background = element_blank()),
                   labels = c("C","","D",""),
                   label.x = 0.02,label.y = c(1.02),
                   ncol = 4, nrow = 1,align = "none",   ##"v"竖直对齐
                   font.label = list(size = 12, color ="black"),
                   widths = c(1.0,0.79,0.91,0.93), heights = c(1,1),legend = 'bottom',
                   common.legend=F);p_new
  
  p_new_main<-ggarrange(figure,p_new,
                        ncol = 1, nrow = 2,align = "none",   ##"v"竖直对齐
                        widths = c(1,0.6), heights = c(0.8,0.4),legend = 'bottom',
                        common.legend=F);p_new_main
  
  ggsave("Fig 5_NO+N2O与水分的关系_main1.pdf", p_new_main, width =12, height =9,
         device=cairo_pdf) 
  
  
############## Fig_S9. The relationship of NO and N2O fluxes with soil temperature #######################
  ###用All绘制温度和水分对NO和N2O排放的影响
  unique(mydata_q10.res$year)
  
  mydata_q10.res.all<-summaryBy(data=datai[datai$statlab_wfps=='T',],FUN=myfun,
                                datai[datai$statlab_wfps=='T',c("N2O_N","NO_N","soil_temp")]~lab_wfps+treatment+plot+year)
  
  # mydata_q10.res<-na.omit(mydata_q10.res)   #20240126
  
  
  ###soiltemp
  #[datai$statlab_temp=='T',]，datai$statlab_temp=='T'
  str(datai)
  
  mydata_q10.res1<-
    summaryBy(data=datai[datai$statlab_temp=='T',],FUN=myfun,
              N2O_N+NO_N+WFPS~lab_soiltemp+treatment+plot+year)
  
  
  mydata_q10.res1<-na.omit(mydata_q10.res1)
  
  
  
  
  ##计算NO温度的指数相关性
  exp.no<-data.frame()
  exp.no.m<-data.frame()
  
  exp.n2o.t<-data.frame()
  exp.n2o.tm<-data.frame()
  rm(exp)
  
  str(mydata_q10.res1)
  k=1
  for(i in c('All')){
    for(j in c("control","warmed")){
        dai=mydata_q10.res1[which(mydata_q10.res1$treatment==j),]
      #NO温度敏感性计算,用0.5度间隔拟合 
      myexp<-summary(lm(log(NO_N.m)~lab_soiltemp,data=dai[dai$NO_N.m>0,]));myexp
      exp.no[k,1]<-i
      exp.no[k,2]<-j
      exp.no[k,3]<-myexp$coefficients[[1]]   #y=ax+b中的b
      exp.no[k,4]<-myexp$coefficients[[2]]   #y=ax+b中的a
      exp.no[k,5]<-myexp$r.squared           #R2
      exp.no[k,6]<-myexp$coefficients[[8]]   #P
      x=seq(min(dai$lab_soiltemp),max(dai$lab_soiltemp),0.5)
      exp.no.mi<-data.frame(soil_temp.m=x,
                            NO_N.m=exp(myexp$coefficients[[2]]*x+myexp$coefficients[[1]]),
                            year=rep(i,length(x)),
                            treatment=rep(j,length(x)),
                            a=myexp$coefficients[[2]],
                            b=myexp$coefficients[[1]])
      exp.no.m<-rbind(exp.no.m,exp.no.mi)
      
      #N2O温度敏感性计算    
      myexp<-summary(lm(log(N2O_N.m+10)~lab_soiltemp,
                        data=dai[which(dai$lab_soiltemp>7.5),]));myexp
      exp.n2o.t[k,1]<-i
      exp.n2o.t[k,2]<-j
      exp.n2o.t[k,3]<-myexp$coefficients[[1]]   #y=ax+b中的b
      exp.n2o.t[k,4]<-myexp$coefficients[[2]]   #y=ax+b中的a
      exp.n2o.t[k,5]<-myexp$r.squared           #R2
      exp.n2o.t[k,6]<-myexp$coefficients[[8]]   #P
      k=k+1
      
      x=seq(min(dai$lab_soiltemp),max(dai$lab_soiltemp),0.5)
      exp.n2o.tmi<-data.frame(soil_temp.m=x,
                              N2O_N.m=exp(myexp$coefficients[[2]]*x+myexp$coefficients[[1]])-10,
                              year=rep(i,length(x)),
                              treatment=rep(j,length(x)),
                              a=myexp$coefficients[[2]],
                              b=myexp$coefficients[[1]])
      exp.n2o.tm<-rbind(exp.n2o.tm,exp.n2o.tmi)
      
      rm(exp.no.mi);rm(exp.n2o.tmi)
    }
  }
  
  names(exp.no)[1:6]<-c("year","treatment","b","a","R2","P_value")
  names(exp.n2o.t)[1:6]<-c("year","treatment","b","a","R2","P_value")
  
  exp.no$q_10<-exp(10*exp.no$a)
  
  exp.no$eq_txt<-paste(round(exp(exp.no$b),2),'e^',
                       paste(round(exp.no$a,3),'*x'),sep="~")
  
  exp.n2o.t$q_10<-exp(10*exp.n2o.t$a)
  
  exp.n2o.t$eq_txt<-paste(round(exp(exp.n2o.t$b),2),'e^',
                          paste(round(exp.n2o.t$a,3),'*x'),sep="~")
  
  
  str(mydata_q10.res1)
  
  str(mydata_q10.res)
  
  unique(mydata_q10.plot$year)
  
  unique(mydata_q10.plot$lab_season)
  
  pt1<-ggplot(mydata_q10.res1[mydata_q10.res1$year=='All',],aes(lab_soiltemp,NO_N.m,col=treatment))+
    # facet_grid(.~year,scale="free_y")+
    geom_point(data=mydata_q10.plot[mydata_q10.plot$year %in% c(2019:2023)  & mydata_q10.plot$warming=='Y',],
               aes(soil_temp,NO_N,colour = treatment),
               size=1.2,alpha=.05,stroke=0,show.legend = F)+
    geom_point(size=1.4,alpha=.5)+ geom_point(size=1.5,col="black",pch=1,stroke=0.1,alpha=0.4)+
    geom_smooth(method = 'glm', 
                formula = y ~ x, linewidth=0.6,
                method.args = list(family = Gamma(link = "log")), alpha = 0.6,se=F)+     
    scale_x_continuous(breaks=seq(-5,25,10),limits = c(-6,26),expand=c(0,0),
                       guide = "prism_offset",minor_breaks = seq(-5,25,5))+
    guides(x="prism_offset_minor",y="prism_offset_minor")+
    scale_y_continuous(breaks=seq(0,80,20),limits = c(0,80),expand = c(0,0))+
    scale_colour_manual(values=c("blue","red"),name=NULL,limits=c("control", "warmed"),
                        labels=c("control","warmed"))+
    geom_text(data=exp.no[which(exp.no$treatment=="control" & exp.no$year=='All'),],
              aes(x=rep(-5,length(exp.no$year[exp.no$treatment=="control" & exp.no$year=='All'])),
                  y=rep(75,length(exp.no$year[exp.no$treatment=="control" & exp.no$year=='All'])),
                  label =paste('~italic(y)',eq_txt,sep="==")), #paste(txt,round(q_10,2),sep="==")),
              size=4.2, vjust=0.4, hjust=0, inherit.aes = FALSE,
              color = "blue", parse = TRUE)+
    geom_text(data=exp.no[which(exp.no$treatment=="warmed" & exp.no$year=='All'),],
              aes(x=rep(-5,length(exp.no$year[exp.no$treatment=="control" & exp.no$year=='All'])),
                  y=rep(65,length(exp.no$year[exp.no$treatment=="control" & exp.no$year=='All'])),
                  label =paste('~italic(y)',eq_txt,sep="==")), #paste(txt,round(q_10,2),sep="==")),
              size=4.2, vjust=0.4, hjust=0, inherit.aes = FALSE,
              color = "red", parse = TRUE)+
    geom_text(data=exp.no[which(exp.no$treatment=="control" & exp.no$year=='All'),],
              aes(x=rep(4.5,length(exp.no$year[exp.no$treatment=="control" & exp.no$year=='All'])),
                  y=rep(75,length(exp.no$year[exp.no$treatment=="control" & exp.no$year=='All'])),
                  label =paste('~~~italic(R^2)',round(R2,2),sep="==")), #paste(txt,round(q_10,2),sep="==")),
              size=4.2, vjust=0.4, hjust=0, inherit.aes = FALSE,
              color = "blue", parse = TRUE)+
    geom_text(data=exp.no[which(exp.no$treatment=="warmed" & exp.no$year=='All'),],
              aes(x=rep(4.5,length(exp.no$year[exp.no$treatment=="control" & exp.no$year=='All'])),
                  y=rep(65,length(exp.no$year[exp.no$treatment=="control" & exp.no$year=='All'])),
                  label =paste('~~~italic(R^2)',round(R2,2),sep="==")), #paste(txt,round(q_10,2),sep="==")),
              size=4.2, vjust=0.4, hjust=0, inherit.aes = FALSE,
              color = "red", parse = TRUE)+
    geom_text(data=exp.no[which(exp.no$treatment=="control" & exp.no$year=='All'),],
              aes(x=rep(-5,length(exp.no$year[exp.no$treatment=="control" & exp.no$year=='All'])),
                  y=rep(70,length(exp.no$year[exp.no$treatment=="control" & exp.no$year=='All'])),
                  label =paste('~italic(Q[10])',round(q_10,2),sep="==")), #paste(txt,round(q_10,2),sep="==")),
              size=4.2, vjust=0.4, hjust=0, inherit.aes = FALSE,
              color = "blue", parse = TRUE)+
    geom_text(data=exp.no[which(exp.no$treatment=="warmed" & exp.no$year=='All'),],
              aes(x=rep(-5,length(exp.no$year[exp.no$treatment=="control" & exp.no$year=='All'])),
                  y=rep(60,length(exp.no$year[exp.no$treatment=="control" & exp.no$year=='All'])),
                  label =paste('~italic(Q[10])',round(q_10,2),sep="==")), #paste(txt,round(q_10,2),sep="==")),
              size=4.2, vjust=0.4, hjust=0, inherit.aes = FALSE,
              color = "red", parse = TRUE)+
    labs(x=expression("Soil temperature (\u00B0C)"),
         y=expression(paste("NO flux (μg"," N"," m"^-2," h"^-1,")")))+
    theme(plot.subtitle = element_text(vjust = 1),plot.caption = element_text(vjust = 1), 
          panel.grid.major = element_blank(), panel.grid=element_blank(), 
          panel.background = element_blank(),
          panel.border = element_rect(size=0.6,fill=NA),
          axis.text.x = element_text(colour = "black",size=12), 
          axis.text.y = element_text(colour = "black",size=12), 
          axis.title = element_text(colour = "black",size=12),
          axis.ticks.length = unit(1.7, 'mm'),
          axis.ticks = element_line(linewidth = 0.35),
          legend.title = element_text(colour = "black",size=12),
          legend.text = element_text(colour = "black",size=12),
          legend.key = element_rect(fill = NA,color=NA), 
          legend.background = element_rect(fill = NA), 
          legend.position =c(0.20,0.55),
          plot.background = element_rect(linetype = "solid",colour=NA),
          plot.margin = unit(c(2,1,0,1),'lines'),
          panel.spacing.x = unit(5, 'mm'),
          strip.text = element_text(size=12, face="bold"),
          strip.background = element_rect(colour="black",fill=NA));pt1
 
  
  library(segmented)
  #fit simple linear regression model
  str(mydata_q10.res1)
  fit <- lm(log(N2O_N.m,10) ~ lab_soiltemp, data=mydata_q10.res1[mydata_q10.res1$year=='All' & mydata_q10.res1$lab_soiltemp>= -1,])

  segmented.fit <- segmented(fit, seg.Z = ~lab_soiltemp, npsi=2)
  
  #view summary of segmented model
  summary(segmented.fit)   
  
  pt2<-
    ggplot(mydata_q10.res1[mydata_q10.res1$year=='All',],
           aes(lab_soiltemp,N2O_N.m,colour = treatment))+
    geom_vline(xintercept = 5,col='black',lty=2,size=0.35)+
    geom_vline(xintercept = 15,col='black',lty=2,size=0.35)+
    
    geom_point(data=mydata_q10.plot[mydata_q10.plot$year %in% c(2018:2023) & mydata_q10.plot$warming=='Y',],
               aes(soil_temp,N2O_N,colour = treatment),
               size=1.2,alpha=.1,stroke=0)+
    geom_point(size=1.4,alpha=.5)+ geom_point(size=1.5,col="black",pch=1,stroke=0.1,alpha=0.4)+
    scale_y_log10(breaks=c(0.01,0.1,1,10,100,1000),limits=c(0.01,1000),labels=c(0.01,0.1,1,10,100,1000),
                  minor_breaks = n_breaks_i,guide = "prism_offset_minor",expand = c(0,0))+
    stat_poly_eq(data=mydata_q10.res1[mydata_q10.res1$lab_soiltemp<=5  & mydata_q10.res1$year=='All',],
                 formula =y~x,eq.x.rhs="x",coef.digits = 3,rr.digits = 2, 
                 eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
                 aes(lab_soiltemp,N2O_N.m,colour = treatment,label = paste(stat(eq.label))),
                 parse = TRUE,label.x.npc = 0.02, hjust = 0,
                 label.y.npc = c(0.22,0.10),size = 4.2)+
    stat_poly_eq(data=mydata_q10.res1[mydata_q10.res1$lab_soiltemp<=5  & mydata_q10.res1$year=='All',],
                 formula =y~x,eq.x.rhs="x",coef.digits = 3,rr.digits = 2, 
                 eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
                 aes(lab_soiltemp,N2O_N.m,colour = treatment,label =paste(stat(rr.label),stat(p.value.label),sep ="*\",\"~~")),
                 parse = TRUE,label.x.npc = 0.02, hjust = 0,
                 label.y.npc = c(0.16,0.04),size = 4.2)+
    
    geom_smooth(data=mydata_q10.res1[mydata_q10.res1$lab_soiltemp<=5  & mydata_q10.res1$lab_soiltemp>=-1 &
                                       mydata_q10.res1$year=='All',],aes(colour = treatment),
                se=F,linewidth=.6,formula=y~x,method = 'lm')+
    geom_smooth(data=mydata_q10.res1[mydata_q10.res1$lab_soiltemp>=5 & mydata_q10.res1$lab_soiltemp<=15 &
                                       mydata_q10.res1$year=='All',],aes(colour = treatment),
                se=F,linewidth=.6,formula=y~x,method = 'lm')+
    geom_smooth(data=mydata_q10.res1[mydata_q10.res1$lab_soiltemp>=15 & mydata_q10.res1$year=='All',],aes(colour = treatment),
                se=F,linewidth=.6,formula=y~x,method = 'lm')+
    stat_poly_eq(data=mydata_q10.res1[mydata_q10.res1$lab_soiltemp>=5 & mydata_q10.res1$lab_soiltemp<=15 & mydata_q10.res1$year=='All',],
                 formula =y~x,eq.x.rhs="x",coef.digits = 3,rr.digits = 2, 
                 eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
                 aes(lab_soiltemp,N2O_N.m,col=treatment,label = paste(stat(eq.label))),
                 parse = TRUE,label.x.npc = 0.35, hjust = 0,
                 label.y.npc = c(0.97,0.85),size = 4.2)+
    stat_poly_eq(data=mydata_q10.res1[mydata_q10.res1$lab_soiltemp>=5 & mydata_q10.res1$lab_soiltemp<=15 & mydata_q10.res1$year=='All',],
                 formula =y~x,eq.x.rhs="x",coef.digits = 3,rr.digits = 2, 
                 eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
                 aes(lab_soiltemp,N2O_N.m,col=treatment,label =paste(stat(rr.label),stat(p.value.label),sep ="*\",\"~~")),
                 parse = TRUE,label.x.npc = 0.35, hjust = 0,
                 label.y.npc = c(0.91,0.79),size = 4.2)+
    stat_poly_eq(data=mydata_q10.res1[mydata_q10.res1$lab_soiltemp>=15 & mydata_q10.res1$year=='All',],
                 formula =y~x,eq.x.rhs="x",coef.digits = 3,rr.digits = 2, 
                 eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
                 aes(lab_soiltemp,N2O_N.m,col=treatment,label = paste(stat(eq.label))),
                 parse = TRUE,label.x.npc = 0.67, hjust = 0,
                 label.y.npc = c(0.22,0.10),size = 4.2)+
    stat_poly_eq(data=mydata_q10.res1[mydata_q10.res1$lab_soiltemp>=15 & mydata_q10.res1$year=='All',],
                 formula =y~x,eq.x.rhs="x",coef.digits = 3,rr.digits = 2, 
                 eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
                 aes(lab_soiltemp,N2O_N.m,col=treatment,label =paste(stat(rr.label),stat(p.value.label),sep ="*\",\"~~")),
                 parse = TRUE,label.x.npc = 0.67, hjust = 0,
                 label.y.npc = c(0.16,0.04),size = 4.2)+
    scale_x_continuous(breaks=seq(-5,25,10),limits = c(-6,26),expand=c(0,0),
                       guide = "prism_offset",minor_breaks = seq(-5,25,5))+
    scale_colour_manual(limits=c("control","warmed"),values=c("blue","red"),
                        labels=c("control","warmed"),name=NULL)+
    scale_fill_manual(limits=c("control","warmed"),values=c("blue","red"),
                      labels=c("control","warmed"),name=NULL)+
    labs(x=expression("Soil temperature (\u00B0C)"),
         y=expression(paste("N"[2],"O flux (μg"," N"," m"^-2," h"^-1,")")))+
    guides(x="prism_offset_minor",y="prism_offset_minor")+
    theme(plot.subtitle = element_text(vjust = 1),plot.caption = element_text(vjust = 1), 
          panel.grid.major = element_blank(), panel.grid=element_blank(), 
          panel.background = element_blank(),
          panel.border = element_rect(size=0.6,fill=NA),
          axis.text.x = element_text(colour = "black",size=12), 
          axis.text.y = element_text(colour = "black",size=12), 
          axis.title = element_text(colour = "black",size=12),
          axis.ticks.length = unit(1.7, 'mm'),
          axis.ticks = element_line(linewidth = 0.35),
          legend.title = element_text(colour = "black",size=12),
          legend.text = element_text(colour = "black",size=12),
          legend.key = element_rect(fill = NA), 
          legend.background = element_rect(fill = NA), 
          legend.position ='none',
          plot.background = element_blank(),
          plot.margin = unit(c(1,1,0,0),'lines'),
          panel.spacing.x = unit(5, 'mm'),
          strip.text = element_text(size=12, face="bold"),
          strip.background = element_rect(colour="black",fill=NA));pt2
  
  figure<-ggarrange(pt1+theme(plot.margin = unit(c(1,0.5,0.1,0.1),'lines')),
                    pt2+theme(plot.margin = unit(c(1,0.5,0.1,0.1),'lines')),
                    labels = c("A", "B"),
                    label.x = 0.83,label.y = c(0.93,0.93),
                    ncol = 2, nrow = 1,align = 'v',   ##"v"竖直对齐
                    font.label = list(size = 12, color ="black"),
                    widths = c(6,6), heights = c(5,5),
                    common.legend=F);figure
  
  ggsave("Fig_S9_NO+N2O_T_all.pdf",figure, device=cairo_pdf,height = 5.2, width = 11.5)
  
 
