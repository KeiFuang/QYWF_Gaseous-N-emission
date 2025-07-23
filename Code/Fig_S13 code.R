{pacman::p_load("ggplot2","ggpmisc","ggpubr","ggprism","doBy","car","readr","crayon","rstatix",
                "lme4","MCMCglmm","tidyverse","sjPlot","directlabels","agricolae","lubridate",
                "patchwork","reshape2","dplyr","animation")  #加载多个包
  Sys.setlocale("LC_TIME","English")  
  
  windowsFonts(Arial=windowsFont("Arial"),
               Times=windowsFont("Times New Roman"))
  
  mythem<-theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1), 
                axis.text.x = element_text(colour = "black",size=21,angle = 0,hjust = .5,vjust =1),
                axis.ticks = element_line(linewidth = .5),
                axis.ticks.length = unit(2,"mm"),
                prism.ticks.length = unit(1,"mm"),
                axis.text.y = element_text(colour = "black",size = 21,hjust=.5), 
                axis.title.x =element_text(size=21), axis.title.y=element_text(colour = "black",size=21),
                legend.text = element_text(size=21), legend.title =element_text(size=21),
                legend.margin = margin(unit(c(0.05,3,5,3),'mm')),
                # legend.box.background = element_rect(fill ="white", size=0.6,colour = "black", linetype = "solid"),
                panel.background = element_rect(fill =NA, linewidth=0.6,colour = "black", linetype = "solid",
                                                inherit.blank=T),
                panel.grid=element_blank(),     #element_line(colour="gray",size=0.5),
                panel.grid.major = element_blank(),    #element_line(colour = "gray",size=0.6), 
                plot.background = element_rect(fill=NA,color=NA,linetype = "solid"),
                legend.key = element_rect(fill = NA), 
                legend.background = element_rect(fill = NA), 
                plot.margin = unit(c(5,5,2,2), 'mm'),   #调整画图区域的间距，从上右下左调整
                strip.background = element_rect(fill = NA,color='black'), 
                strip.text = element_text(colour = "black",size = 20,hjust=.5),
                legend.position = "none")
  
  ###计算数据均值，范围，se
  datFUN<-function(x){c(mean=round(mean(na.omit(x)),3),
                        range=paste(round(min(na.omit(x)),3),round(max(na.omit(x)),3),sep="~"),
                        n=length(na.omit(x)),
                        se=round(sd(na.omit(x))/sqrt(length(na.omit(x))),3),
                        sd=round(sd(na.omit(x)),3))
  }
}

library(plotly);library(grid)
list.files('D:/Workspace/book/Qingyuan/soil_data/Incubation/Gaseou N')
xi<-readxl::read_excel('D:/Workspace/book/Qingyuan/soil_data/Incubation/Gaseou N/xidan_forests.xlsx',sheet='kai_recal',
                       col_names=TRUE,col_types=NULL)

str(xi)
ggplot(xi,aes(WFPS,N2_N2O))+
  geom_smooth(method='lm',formula = y~poly(x,1),se=T,col='black',alpha=0.3,linewidth=0.6)+
  geom_point(size=2,pch=1,stroke=1)+
  # scale_y_continuous(breaks=seq(0,30,10),expand=c(0.05,0.05))+
  scale_x_continuous(limits=c(20,140),breaks=seq(20,140,40),expand=c(0.05,0.05))+
  labs(x='WFPS (%)',y=expression(paste('N'[2],'/N'[2],'O ratio')))+
  mythem+theme(panel.grid.major = element_line(linewidth = 0.1,colour = 'gray'),
               panel.grid.minor = element_line(linewidth = 0.1,colour = 'gray'),
               axis.title = element_text(size=9),
               axis.text =  element_text(size=9),
               plot.margin = unit(c(1,1,1,1),'cm'))+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = paste(stat(eq.label),sep ="*\",\"~~~")),#sep = "*\", \"*")),
               label.x.npc = "left", label.y.npc = 0.95,parse = TRUE,
               vstep=0.2,hstep=0.1,size = 4.2)+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = paste(stat(rr.label),stat(p.value.label),sep ="*\",\"~~~")),#sep = "*\", \"*")),
               label.x.npc = "left", label.y.npc = 0.85,parse = TRUE,
               vstep=0.2,hstep=0.1,size = 4.2)+
  coord_cartesian(clip = 'on',ylim = c(0,20))->p1;p1


setwd('D:/工作目录/202409/Manuscript_kai/Talk_20250503/Data and Code/Supplementary')

ggsave("N2_N2O ratio with moisture gradient.pdf", p1, width =9, height =7) 

xi$N2_N2O_mod<- -1.1+0.13*xi$WFPS

xi$res_N2_N2O<- xi$N2_N2O-xi$N2_N2O_mod

xi$res_N2<- xi$N2_N2O_mod*xi$N2O-xi$N2


ggplot(xi,aes(WFPS,res_N2_N2O))+
  geom_smooth(method='lm',formula = y~poly(x,1),se=T,col='black',alpha=0.3)+
  geom_point(size=2,pch=1,stroke=1)


ggplot(xi,aes(res_N2_N2O))+
  geom_smooth(method='lm',formula = y~poly(x,1),se=T,col='black',alpha=0.3)+
  geom_point(size=2,pch=1,stroke=1)


library(ggplot2)

ggplot(xi, aes(x = res_N2_N2O)) +
  geom_histogram(
    aes(y = ..density..),  # 转换为密度
    bins = 5,            # 调整分箱数
    fill = "lightblue",
    color = "black",
    alpha = 0.7
  ) +
  geom_density(           # 添加密度曲线
    color = "red",
    linewidth = 1
  ) +
  labs(
    x = "res_N2_N2O",
    y = "Density",
    title = "Distribution of res_N2_N2O"
  ) +
  theme_minimal()


# 首先计算必要的参数
# 计算均值和标准差
mean_val <- mean(xi$res_N2_N2O, na.rm = TRUE)
sd_val <- sd(xi$res_N2_N2O, na.rm = TRUE)
n <- length(na.omit(xi$res_N2_N2O))

# 计算95%置信区间
conf_int <- qnorm(c(0.025, 0.975), mean = mean_val, sd = sd_val)

# 创建临时图获取bin宽度
temp_plot <- ggplot(xi, aes(x = res_N2_N2O)) + geom_histogram(bins = 5)
bin_width <- diff(ggplot_build(temp_plot)$data[[1]]$x)[1]

# 正式绘图
p2 <- ggplot(xi, aes(x = res_N2_N2O)) +
  stat_function(
    fun = function(x) {
      dnorm(x, mean = mean_val, sd = sd_val) * n * bin_width
    },
    color = "black",
    linewidth = 0.6
  ) +
  geom_vline(xintercept = 0,lty=2,col='red')+
  # 添加95%置信区间阴影
  # geom_area(
  #   stat = "function",
  #   fun = function(x) {
  #     dnorm(x, mean = mean_val, sd = sd_val) * n * bin_width
  #   },
  #   fill = "gray80",
  #   alpha = 0.5,
  #   xlim = conf_int
  # ) +
  # 添加置信区间边界线
  # geom_vline(
  #   xintercept = conf_int,
  #   linetype = "dashed",
  #   color = "red",
  #   linewidth = 0.5
  # ) +
  labs(
    x = expression(atop(
      paste('N'[2], '/N'[2], 'O ratio residual'),
      paste('(Observed-Simulated)')
    )),
    y = "Count"
  ) +
  scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) +
  # scale_x_continuous(limits = c(-5, 5), expand = c(0, 0)) +
  mythem+
  theme(
    panel.grid.major = element_line(linewidth = 0.1, colour = 'gray'),
    panel.grid.minor = element_line(linewidth = 0.1, colour = 'gray'),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 9),
    plot.margin = unit(c(1, 1, 1, 1), 'cm'),
    panel.background = element_rect(fill =NA, linewidth=0.6,colour = "black", linetype = "solid",
                                    inherit.blank=T)
  );p2

# # 显示图形
# print(p2)
# 
# # 可选：在图形上添加置信区间数值标签
# p2 + annotate(
#   "text",
#   x = mean_val,
#   y = 3.5,
#   label = paste0("95% CI: [", round(conf_int[1], 2), ", ", round(conf_int[2], 2), "]"),
#   size = 3
# )

ggplot(xi,aes(WFPS,res_N2_N2O))+
  geom_smooth(method='lm',formula = y~poly(x,1),se=T,col='black',alpha=0.3,linewidth=0.6)+
  geom_point(size=2,pch=1,stroke=1)+
  scale_x_continuous(limits=c(20,140),breaks=seq(20,140,40),expand=c(0.05,0.05))+
  labs(x='WFPS (%)',y=expression(atop(paste('N'[2],'/N'[2],'O ratio residual'),
                                      paste('(Observed-Simulated)'))))+
  mythem+theme(panel.grid.major = element_line(linewidth = 0.1,colour = 'gray'),
               panel.grid.minor = element_line(linewidth = 0.1,colour = 'gray'),
               axis.title = element_text(size=9),
               axis.text =  element_text(size=9),
               plot.margin = unit(c(1,1,1,1),'cm'))+
  stat_poly_eq(formula =y~x,coef.digits = 2,eq.x.rhs="x", 
               eq.with.lhs = "italic(y)~`=`~",   #给“y"换形式         
               aes(label = paste(stat(rr.label),stat(p.value.label),sep ="*\",\"~~~")),#sep = "*\", \"*")),
               label.x.npc = "left", label.y.npc = 0.95,parse = TRUE,
               vstep=0.2,hstep=0.1,size = 4.2)+
  coord_cartesian(clip = 'on',ylim = c(-5,5))->p3;p3

(p1+ labs(tag = "A") +theme(axis.title.x = element_text(size=12),
                            axis.title.y = element_text(size=12),
                            axis.text.x = element_text(size=12),
                            axis.text.y = element_text(size=12),
                            plot.tag = element_text(face = "bold",size=12),
                            plot.margin = unit(c(2,0,1,1.5), 'mm'))+
    p2+ labs(tag = "B") +theme(axis.title.x = element_text(size=12),
                               axis.title.y = element_text(size=12),
                               axis.text.x = element_text(size=12),
                               axis.text.y = element_text(size=12),
                               plot.tag = element_text(face = "bold",size=12),
                               plot.margin = unit(c(2,0,1,1.5), 'mm'))+
    p3 + labs(tag = "C") +theme(axis.title.x = element_text(size=12),
                                axis.title.y = element_text(size=12),
                                axis.text.x = element_text(size=12),
                                axis.text.y = element_text(size=12),
                                plot.tag = element_text(face = "bold",size=12),
                                plot.margin = unit(c(2,0,1,1.5), 'mm')))->p;p

setwd('D:/工作目录/202409/Manuscript_kai/Talk_20250503/Data and Code/Fig_S13')

ggsave("N2_N2O ratio with moisture gradient_new.pdf", p, width =10, height =3.5) 
