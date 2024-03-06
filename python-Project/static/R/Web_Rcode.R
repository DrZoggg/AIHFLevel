
###################################################'@=================START!!!!!!!!
library(DynNom)
library(tidyverse)
library(dplyr)
library(ggbreak)
library(ggforce)
library(xgboost)
library(mlr3)
library(mlr3viz)
library(mlr3learners)
library(mlr3verse)
library(mlr3tuning)
library(data.table)
library(mlr3proba)
library(rms)
library(magrittr)
library(mlr3extralearners)
library(survival)
library(scales)
library(ggplot2)
library(fmsb)
library(ggradar)


library(getopt)


load(file = 'Web_preparedata.rda')

##############接收python传来的数据#############

command=matrix(c( 
  'age','A', 1, 'numeric',
  'arrhythmia','H', 1, 'numeric', 
  'CAD','C', 1, 'numeric', 
  'CKD','K', 1, 'numeric', 
  'Cr','R', 1, 'numeric',
  'EF','E', 1, 'numeric', 
  'GFR','G', 1, 'numeric', 
  'LY_Per','L', 1, 'numeric', 
  'MCHC','M', 1, 'numeric', 
  'SV','S', 1, 'numeric', 
  'TN','T', 1, 'numeric', 
  'TBIL','B', 1, 'numeric',
  'Timepoint','I', 1, 'numeric'),byrow=T,ncol=4)
args=getopt(command)
#print(args)


Age <- args$age
Arrhythmia <- args$arrhythmia
CAD <- args$CAD
CKD_stage <- args$CKD
Cr <- args$Cr
EF <- args$EF
GFR  <- args$GFR
LY_Per  <- args$LY_Per
MCHC  <- args$MCHC
SV <- args$SV
TBIL <- args$TBIL
TN <- args$TN
Timepoint <- args$Timepoint



####################################'@1.输入指标数值
#Age = 79
#Arrhythmia = 0
#CAD = 1
#CKD_stage = 1
#Cr = 55.3
#EF = 45
#GFR = 94
#LY_Per = 8.7
#MCHC = 332.9
#SV = 41
#TBIL = 24
#TN = 0.02
#Timepoint = 600


#####################################'@2.根据模型预测患者-----(风险等级)
df_test <- data.frame(
  ACM = 300,
  ACM_Censor = 1,
  Age = Age,
  Arrhythmia = Arrhythmia,
  CAD = CAD,
  CKD_stage = CKD_stage,
  Cr = Cr,
  EF = EF,
  GFR = GFR,
  LY_Per = LY_Per,
  MCHC = MCHC,
  SV = SV,
  TBIL = TBIL,
  TN = TN
)

dat.mat <- df_test%>%as.matrix()
dtrain <- list(data=dat.mat[,3:ncol(dat.mat),drop=F],
               label=dat.mat[,'ACM']*(-(-1)^(as.numeric(dat.mat[,'ACM_Censor']))))
Dtrain <- xgb.DMatrix(dtrain$data,label=dtrain$label)
RS2 = predict(fit2$model,newdata = Dtrain,type='risk')
RS2 = cbind(df_test[,1:2],RS=RS2)
RS2 = RS2%>%mutate(
  risk = case_when(
    RS2$RS > 1.548 ~ 'High risk',
    RS2$RS > 0.435 & RS2$RS <= 1.548 ~ 'Intermidiate risk',
    RS2$RS <= 0.435 ~ 'Low risk'
  )
)
Risk = RS2$risk



#####################################'@3.计算(未来生存率)(全因死亡风险)(置信区间se)
m.pred <- survest(coxpbc, newdata = RS2,times = seq(0,939,1))
surv <- m.pred$surv
time <- m.pred$time
se.pred <- m.pred$std.err


#####################################'@4.可视化(各指标对全因死亡的贡献度)

std1 <- function(x){
  return ((x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T)))
}
source('ZG_shap.score.rank.R')
source('ZG_shap.prep.R')
source('ZG_plot.shap.summary.R')
dat.mat <- df_test%>%as.matrix()
dtrain <- list(data=dat.mat[,3:ncol(dat.mat),drop=F],
               label=dat.mat[,'ACM']*(-(-1)^(as.numeric(dat.mat[,'ACM_Censor']))))
Dtrain <- xgb.DMatrix(dtrain$data,label=dtrain$label)
shap_result = shap.score.rank(xgb_model = fit2$model,
                              X_train = Dtrain,
                              shap_approx = F)
shap_long_hd = shap.prep(shap  = shap_result, X_train = df_test , top_n = 12)

impdf <- shap_long_hd%>%dplyr::select(variable,value,mean_value)
impdf[,3] <- std1(impdf[,'mean_value',drop=F])
impdf_radar <- impdf%>%dplyr::select(variable,mean_value)%>%
  pivot_wider(
    names_from = "variable",
    values_from = "mean_value")
#impdf_radar <- impdf_radar%>%cbind(mod='A',.)
impdf_radar <- rbind(Min = 0,impdf_radar)%>%rbind(Max =1,.)


##########可视化1
pdf('imp_radar.pdf',
    width = 5,height = 4.5)
par(mar = c(.2, .2, .2, .2))
radarchart(impdf_radar,axistype=1,axislabcol = '#2C8EAE' ,
           vlcex= .9,palcex= 2,pcol = '#00AFBB',pfcol = scales::alpha("#00AFBB", 0.5),
           plwd = 4.5 ,plty = 1,pty= 30 ,cglcol = "grey70",cglty = 1,
           cglwd = .5)
dev.off()



##########可视化2
x_bound <- max(abs(shap_long_hd$value))
imp_bar <- ggplot(data = shap_long_hd)+coord_flip() +
  labs(title  = "Feature Importance", x = "",y = '',color = "Feature Value")+
  geom_col(aes(x = variable, y = value, fill = mean_value),width=.7,color='black') +
  # geom_text(data = unique(shap_long_hd[, c("variable", "value"), with = F]),
  #           aes(x = variable, y=value, label = sprintf("%.4f", value)),
  #           size = 3, alpha = 0.7,hjust = -0.2,fontface = "bold") +
  geom_text(data = unique(shap_long_hd[, c("variable", "value"), with = F]),
            aes( x = variable,  label = sprintf("%.4f", value)),y=0.03,
            size = 3, alpha = 0.7,hjust = -0.2,fontface = "bold") +
  theme_bw(base_rect_size = 1.4) +
  guides(fill = guide_legend(title.position = 'right'))+
  scale_fill_gradient(low="white", high="#5AB4C2") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12,colour = 'grey10',face = 'bold',vjust= 2),
        plot.title = element_text(size = 11.5,colour = 'darkred',face = 'bold',hjust= .5),
        axis.text = element_text(size=9,face = 'bold'),
        axis.ticks = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line( color = 'black',size =.6 ),
        panel.grid.major  = element_line(color = "grey95", size = .5, linetype = "solid"),
        panel.grid.minor = element_line(color = "white", size = .3, linetype = "solid"),
        panel.background = element_rect(fill='white'),
        legend.key.height = unit(26,'mm'),
        legend.key.width  = unit(1.3,'mm'),
        legend.position="none") +
  geom_hline(yintercept = 0,color = scales::alpha('grey70',1),size=.5,linetype = "solid") +
  scale_y_continuous(limits = c(-x_bound, x_bound),expand = c(0,0.1)) +
  scale_x_discrete(limits = rev(levels(shap_long_hd$variable)),expand = c(0.07,0.1))


ggsave(imp_bar,filename = c('imp_bar.pdf'),
       width = 5.4,height = 5.8)


#####################################'@5.可视化(全因死亡风险KM曲线)
KMdf<-surv%>%t()%>%as.data.frame();KMdf$time<-time;rownames(KMdf)<-time;colnames(KMdf)<-c('Surv','time')
Est_surv <- KMdf[which(KMdf$time==Timepoint),'Surv']


KMplot <- ggplot(KMdf,aes(x=time,y=Surv)) +
  geom_line(size=1.4,col = scales::alpha('#2C8EAE',1) )+
  geom_vline(xintercept = Timepoint,linetype = "dashed",size = .85,
             col = scales::alpha('#D1D2D4',1))+
  annotate('segment',x=Timepoint+180,xend = Timepoint,
           y=Est_surv+0.015,yend = Est_surv+0.0009,
           size=1.5,arrow=arrow(),
           alpha=.85,color='#ff8c3e')+
  annotate('text',x = 310,y = range(KMdf$Surv)[1]-0.02,
           label=paste0(Timepoint,'-day \nAll-cause mortality = ',
                        percent(round(1-Est_surv,5),0.001)),
           fontface='bold.italic',colour="grey40",size= 4.5)+
  theme_bw(base_rect_size = 1.8)+
  scale_color_manual(values = c('#E36258','#9FD5D6'))+
  scale_y_continuous(limits = c(0.15,1))+
  labs(y='Survival Probability',x='Follow-up Time (Day)')+
  ggtitle('HF&RD-CPS system')+
  scale_y_continuous(limits = c(range(KMdf$Surv)[1]-0.05,1),
                     labels = scales::percent_format(accuracy = 1))+
  theme(legend.position = c('none'),
        legend.background = element_blank(),
        legend.key=element_rect(fill='transparent', colour='transparent'),
        axis.text = element_text(size = 9,colour = 'grey15',face = 'bold'),
        axis.title.x = element_text(size = 11,colour = 'grey15',face = 'bold',vjust= .01),
        axis.title.y = element_text(size = 11,colour = 'grey15',face = 'bold',vjust= 2),
        plot.title = element_text(size=12.5,hjust=0.5,colour = 'darkred',face = 'bold'),
        panel.grid.major  = element_line(color = "grey95", size = .5, linetype = "solid"),
        panel.grid.minor = element_line(color = "white", size = .3, linetype = "solid"),
        #panel.background = element_rect(fill='#F3F6F6'),
        panel.background = element_rect(fill='white'),
        panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid"))

ggsave(KMplot,filename = 'KMplot.pdf',
       width = 5.3,height = 4.4)



#####################################'@6.可视化(未来生存率的95%可信区间)
se.pred <- m.pred$std.err
Est_surv <- KMdf[which(KMdf$time==Timepoint),'Surv']
Est_surv_6 <- KMdf[which(KMdf$time==180),'Surv']
Est_surv_12 <- KMdf[which(KMdf$time==365),'Surv']
Est_surv_24 <- KMdf[which(KMdf$time==730),'Surv']
Est_surv_30 <- KMdf[which(KMdf$time==900),'Surv']

SurvPro<-surv%>%t()%>%as.data.frame();SurvPro$time<-time;rownames(SurvPro)<-time;SurvPro$se<-se.pred[1:length(se.pred)]
colnames(SurvPro)<-c('Surv','time','se')
SurvPro <- SurvPro%>% filter(time==c(180)|
                               time==c(365)|
                               time==c(730)|
                               time==c(900)|
                               time==c(Timepoint)
)
SurvPro$time <- factor(as.character(SurvPro$time),levels = c(180,365,730,900,Timepoint))

Survplot <- ggplot(SurvPro,aes(y=Surv,x=time,fill=time))+
  geom_col(position=position_dodge(0.8),width=.5,size = .8,
           color = 'black')+
  geom_errorbar(aes(x=time,ymin=Surv-se,ymax=Surv+se),width=0.3,
                size=.6,color = 'black')+
  theme_bw(base_rect_size = 1.5)+
  labs(x = 'Time Point (month)', y = 'Survival Probability')+
  ggtitle('95% Credible Interval')+
  theme(axis.title.x = element_text(size = 12,colour = 'grey15',face = 'bold',vjust= .01),
        axis.title.y = element_text(size = 12,colour = 'grey15',face = 'bold',vjust= 2),
        axis.text.x = element_text(size=10,angle = 0,hjust = 1,face = 'bold'),
        axis.text.y = element_text(size=10,face = 'bold'),
        legend.position = 'none',
        panel.grid.major  = element_line(color = "grey95", size = .5, linetype = "solid"),
        panel.grid.minor = element_line(color = "white", size = .3, linetype = "solid"),
        #panel.background = element_rect(fill='#F3F6F6'),
        panel.background = element_rect(fill='white'),
        plot.title = element_text(size=14,hjust=0.5,colour = 'darkred',face = 'bold'),
        strip.text = element_text(size=12),
        axis.ticks.x = element_blank())+
  scale_fill_manual(values = c('#B0997F','#93a8ac','#ffc857','#61a5c2','#119da4','#FF6666'))+
  scale_y_continuous(expand = c(0,0.09))+
  facet_zoom(ylim = c(sort(SurvPro$Surv)[1]-0.03, .92),zoom.size=2)


ggsave(Survplot,filename = 'Survplot.pdf',
       width = 6.2,height = 4.6)




print(Risk)
print( round(SurvPro['180','Surv'],4)  )
print( round(SurvPro['365','Surv'],4)  )
print( round(SurvPro['730','Surv'],4)  )
print( round(SurvPro['900','Surv'],4)  )


####################################'@end!!!!!!!!


