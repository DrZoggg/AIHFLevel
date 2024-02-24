



if (T) {
  options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  options(BioC_mirror = "http://mirrors.ustc.edu.cn/bioc/") 
  library(e1071)
  library(parallel)
  library(preprocessCore)
  library(sva)
  library(ggplot2)
  library(survminer)
  library(survival)
  library(rms)
  library(randomForest)
  library(pROC)
  library(glmnet)
  library(pheatmap)
  library(timeROC)
  library(vioplot)
  library(corrplot)
  library(forestplot)
  library(survivalROC)
  library(beeswarm)
  library(Biostrings)
  library(pheatmap)
  library(stringr)
  library(dplyr)
  library(RColorBrewer)
  library(tibble)
  library(cowplot)
  library(ggcorrplot)
  library(tidyverse)
  library(plyr)
  library(circlize)
  library(SummarizedExperiment)
  library(ggsci)
  library(ComplexHeatmap)
  library(pastecs)
  library(fs)
  library(pastecs)
  library(scales)
  library(ggkm)
  library(survcomp)
  library(pbapply)
  library(survival)
  library(pec)
  library(caret)
  library(leaps)
  library(randomForestSRC)
  library(mlr3)
  library(mlr3viz)
  library(mlr3learners)
  library(mlr3verse)
  library(mlr3tuning)
  library(data.table)
  library(mlr3proba)
  library("magrittr")
  library(mlr3extralearners)
  library(iml)
  library("DALEX")
  library(DALEXtra)
  library(dcurves)
  library(MASS)
  library(survival)
  library(risksetROC)
  library(ggDCA)
  library(beepr)
  library(BRRR)
  options(stringsAsFactors = FALSE)
}



#############################################################################################
##################################'@Web
#############################################################################################
library(xgboost)
rm(list=ls())
load( 'val_dd_list.rda');rm(fit,val_dd_list,task_test,task_train,task_meta )
load('HF_metadata.rda');df<-dataset3 %>% dplyr::select(-c('Name','ID'));rm(dataset3)
df$ID <-  c(paste0('G', c(1:712)));row.names(df) <- df$ID;df <- df %>% dplyr::select(-c('ID'))

df <- df %>% mutate(
  CKD_stage = case_when(
    CKD_stage %in% 'I' ~ 1,CKD_stage %in% 'II' ~ 2,
    CKD_stage %in% 'IIIa'| CKD_stage %in% 'IIIb' ~ 3,
    CKD_stage %in% 'IV' ~ 4,CKD_stage %in% 'V' ~ 5))
df <- df[,c("ACM","ACM_Censor",rid)]
df <- apply(df, 2,as.numeric)%>%as.data.frame()
expr.surv <- df
expr.surv$ACM <- expr.surv$ACM*30

lrn_xgboost <- lrn("surv.xgboost",
                   max_depth = 3,eta = 0.01,gamma = 0,subsample = 0.75,colsample_bytree = 0.8,
                   min_child_weight = 3,lambda=1,alpha=0,nrounds = 500);measures <- c("surv.cindex")
resampling_HD <- rsmp("holdout",ratio = 0.7);set.seed(1771)
rr_HD <- resample(as_task_surv(expr.surv,time="ACM",event="ACM_Censor"),lrn_xgboost,resampling_HD,store_models=T)
index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
fit2 <- rr_HD$learners[[1]]

rr_HD$resampling$train_set(1)
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task$select(rid);task
task_train <- task$filter(rr_HD$resampling$train_set(1));task_train
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task$select(rid);task
task_test <- task$filter(rr_HD$resampling$test_set(1));task_test
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task$select(rid);task
task_meta <- task;task_meta

fit2 <- rr_HD$learners[[1]]
fit2$predict(task_test)$score(msrs(measures))
fit2$predict(task_train)$score(msrs(measures))
fit2$predict(task_meta)$score(msrs(measures))

val_dd_list2 <- list()
val_dd_list2[[1]]=task_train$data();val_dd_list2[[2]]=task_test$data();val_dd_list2[[3]]=task_meta$data()
names(val_dd_list2) <- c('Train_1','Test_2','Meta_3')

rs2 <- lapply(val_dd_list2,function(x){
  
  dat.mat <- x%>%as.matrix()
  dtrain <- list(data=dat.mat[,3:ncol(dat.mat),drop=F],
                 label=dat.mat[,'ACM']*(-(-1)^(as.numeric(dat.mat[,'ACM_Censor']))))
  Dtrain <- xgb.DMatrix(dtrain$data,label=dtrain$label)
  RS = predict(fit2$model,newdata = Dtrain,type='risk')
  cbind(x[,1:2],RS=RS)
})

C_index <- sapply(rs2,function(x){
  concordance.index(scale(x$RS),
                    surv.time = x$ACM,
                    surv.event = x$ACM_Censor,
                    method = "noether")$c.index
});C_index



pbc <- rs2$Meta_3
dd <- datadist(pbc);options(datadist="dd");options(na.action="na.delete")
coxpbc <- cph(formula = Surv(ACM,ACM_Censor) ~RS,data=pbc,x=T,y=T,surv = T,na.action=na.delete)

pbc <- rs2$Meta_3
dd <- datadist(pbc);options(datadist="dd");options(na.action="na.delete")
coxpbc <- cph(formula = Surv(ACM,ACM_Censor) ~RS,data=pbc,x=T,y=T,surv = T,na.action=na.delete)

save(coxpbc,file = 'Web_coxpbc.rda')
save(fit2,file = 'Web_fit2.rda')
save(coxpbc,fit2,expr.surv,file = 'Web_preparedata.rda')
rm(list=ls())

###################################################'@=================START!!!!!!!!
setwd('Web')
dev.off()
rm(list=ls())
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
load(file = 'Web_preparedata.rda')

72
0
1
4
101.3
51.0
23.490
28.3
328.00
42
2.79
0.35

####################################'@1.IUPUT
Age = 72
Arrhythmia = 0
CAD = 1
CKD_stage = 4
Cr = 101.3
EF = 51.0
GFR = 23.490
LY_Per = 28.3
MCHC = 328.00
SV = 42
TBIL = 2.79
TN = 0.35
Timepoint = 600

filename = 70555.5


#####################################'@2.Prognostic-stratification
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



#####################################'@3.Calculation-of-long-term-survival-probabilities
m.pred <- survest(coxpbc, newdata = RS2,times = seq(0,939,1))
surv <- m.pred$surv
time <- m.pred$time
se.pred <- m.pred$std.err


#####################################'@4.Contribution

std1 <- function(x) {
  return((x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T)))
}

shap.score.rank <- function(xgb_model = xgb_mod, shap_approx = TRUE,
                            X_train = mydata$train_mm) {
  require(xgboost)
  require(data.table)
  shap_contrib <- predict(xgb_model, X_train,
    predcontrib = TRUE, approxcontrib = shap_approx
  )
  shap_contrib <- as.data.table(shap_contrib)
  shap_contrib[, BIAS := NULL]
  cat("make SHAP score by decreasing order\n\n")
  mean_shap_score <- colMeans(abs(shap_contrib))[order(colMeans(abs(shap_contrib)), decreasing = T)]
  return(list(
    shap_score = shap_contrib,
    mean_shap_score = (mean_shap_score)
  ))
}

shap.prep <- function(shap = shap_result, X_train = mydata$train_mm, top_n) {
  require(ggforce)
  std1 <- function(x) {
    return((x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T)))
  }
  # descending order
  if (missing(top_n)) top_n <- dim(X_train)[2] # by default, use all features
  if (!top_n %in% c(1:dim(X_train)[2])) stop("supply correct top_n")
  require(data.table)
  shap_score_sub <- as.data.table(shap$shap_score)
  shap_score_sub <- shap_score_sub[, names(shap$mean_shap_score)[1:top_n], with = F]
  shap_score_long <- melt.data.table(shap_score_sub, measure.vars = colnames(shap_score_sub))

  # feature values: the values in the original dataset
  fv_sub <- as.data.table(X_train)[, names(shap$mean_shap_score)[1:top_n], with = F]
  # standardize feature values
  fv_sub_long <- melt.data.table(fv_sub, measure.vars = colnames(fv_sub))
  fv_sub_long[, stdfvalue := std1(value), by = "variable"]
  # SHAP value: value
  # raw feature value: rfvalue;
  # standarized: stdfvalue
  names(fv_sub_long) <- c("variable", "rfvalue", "stdfvalue")
  shap_long2 <- cbind(shap_score_long, fv_sub_long[, c("rfvalue", "stdfvalue")])
  shap_long2[, mean_value := mean(abs(value)), by = variable]
  setkey(shap_long2, variable)
  return(shap_long2)
}


plot.shap.summary <- function(data_long) {
  x_bound <- max(abs(data_long$value))
  require("ggforce") # for `geom_sina`
  plot1 <- ggplot(data = data_long) +
    coord_flip() +
    labs(
      title = "SHapley Additive exPlanations", x = "",
      y = "SHAP value (Impact on Model Output)",
      color = "Feature Value"
    ) +
    # sina plot:
    geom_sina(aes(x = variable, y = value, color = stdfvalue), size = 1.7) +
    # print the mean absolute value:
    geom_text(
      data = unique(data_long[, c("variable", "mean_value"), with = F]),
      aes(x = variable, y = -Inf, label = sprintf("%.4f", mean_value)),
      size = 3, alpha = 0.7,
      hjust = -0.2,
      fontface = "bold"
    ) + # bold
    # theme_classic(base_rect_size = 1.5) +
    theme_bw(base_rect_size = 1.65) +
    guides(fill = guide_legend(title.position = "right")) +
    # # add a "SHAP" bar notation
    # annotate("text", x = -Inf, y = -Inf, vjust = -0.2, hjust = 0, size = 3,
    #          label = expression(group("|", bar(SHAP), "|"))) +
    # scale_color_gradient(low="#FFCC33", high="#6600CC",
    #                      breaks=c(0,1), labels=c("Low","High")) +
    scale_color_gradient2(
      low = "#5AB4C2", mid = "white", high = "#E36258",
      labels = c("Low", "", "High"), breaks = c(0, 0.5, 1),
      midpoint = 0.5 # breaks=c(0,1), labels=c("Low",'0',"High")
    ) +
    theme(
      axis.title.x = element_text(size = 10, colour = "black", face = "bold", vjust = .01),
      axis.title.y = element_text(size = 12, colour = "grey15", face = "bold", vjust = 2),
      plot.title = element_text(size = 12, colour = "darkred", face = "bold", hjust = .6),
      axis.text = element_text(size = 9, face = "bold"),
      axis.ticks = element_blank(),
      axis.line.y = element_blank(),
      # axis.line.x = element_line( color = 'black',size = .8 ),

      panel.grid.major = element_line(color = "white", size = .8, linetype = "solid"),
      panel.grid.minor = element_line(color = "white", size = .5, linetype = "solid"),

      # panel.background = element_rect(fill='white'),
      panel.background = element_rect(fill = scales::alpha("#F0F8FF", .7)),
      legend.key.height = unit(27, "mm"),
      legend.key.width = unit(1.3, "mm"),
      legend.position = "right"
    ) +
    geom_hline(yintercept = 0, color = scales::alpha("grey80", 1), size = .65, linetype = "dashed") + # the vertical line
    # geom_hline(yintercept = 0,color = scales::alpha('white',1),size= 1.0,linetype = "solid") + # the vertical line

    scale_y_continuous(limits = c(-x_bound, x_bound)) +
    # reverse the order of features
    scale_x_discrete(limits = rev(levels(data_long$variable)))

  return(plot1)
}

dat.mat <- df_test %>% as.matrix()
dtrain <- list(
  data = dat.mat[, 3:ncol(dat.mat), drop = F],
  label = dat.mat[, "ACM"] * (-(-1)^(as.numeric(dat.mat[, "ACM_Censor"])))
)
Dtrain <- xgb.DMatrix(dtrain$data, label = dtrain$label)
shap_result <- shap.score.rank(
  xgb_model = fit2$model,
  X_train = Dtrain,
  shap_approx = F
)
shap_long_hd <- shap.prep(shap = shap_result, X_train = df_test, top_n = 12)

impdf <- shap_long_hd %>% dplyr::select(variable, value, mean_value)
impdf[, 3] <- std1(impdf[, "mean_value", drop = F])
impdf_radar <- impdf %>%
  dplyr::select(variable, mean_value) %>%
  pivot_wider(
    names_from = "variable",
    values_from = "mean_value"
  )
impdf_radar <- rbind(Min = 0, impdf_radar) %>% rbind(Max = 1, .)


#########
png(paste0(filename, "_imp_radar.png"), width = 4.3, height = 4, res = 800, units = "cm", family = "sans")
par(mar = c(.2, .2, .2, .2))
radarchart(impdf_radar,
  axistype = 1, axislabcol = "#2C8EAE", calcex = .37,
  vlcex = .32, palcex = 2.5, pcol = "#00AFBB", pfcol = scales::alpha("#00AFBB", 0.5),
  plwd = 1.2, plty = 1, pty = 30, cglcol = "grey90", cglty = 1,
  cglwd = .4
)
dev.off()


pdf(paste0(filename, "_imp_radar.pdf"), width = 2.3, height = 3)
par(mar = c(.2, .2, .2, .2))
radarchart(impdf_radar,
  axistype = 1, axislabcol = "#2C8EAE", calcex = .37,
  vlcex = .32, palcex = 2.5, pcol = "#00AFBB", pfcol = scales::alpha("#00AFBB", 0.5),
  plwd = 1.2, plty = 1, pty = 30, cglcol = "grey90", cglty = 1,
  cglwd = .4
)
dev.off()

#########
x_bound <- max(abs(shap_long_hd$value))
par(mar = c(0, 0, 0, 0))
imp_bar <- ggplot(data = shap_long_hd) +
  coord_flip() +
  labs(title = "Feature Importance", y = "", color = "Feature Value") +
  geom_col(aes(x = variable, y = value, fill = mean_value), width = .7, color = "black", size = .25) +
  geom_text(
    data = unique(shap_long_hd[, c("variable", "value"), with = F]),
    aes(x = variable, label = sprintf("%.4f", value)), y = 0.03,
    size = .9, alpha = 0.7, hjust = -0.2, fontface = "bold"
  ) +
  theme_bw(base_rect_size = .6) +
  guides(fill = guide_legend(title.position = "right")) +
  scale_fill_gradient(low = "white", high = "#5AB4C2") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(size = 3.3, colour = "darkred", face = "bold", hjust = .5, vjust = -3),
    axis.text.x = element_text(size = 2.4, face = "bold", vjust = 5),
    axis.text.y = element_text(size = 2.4, face = "bold", hjust = 1),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(size = .2),
    axis.ticks.length.x = unit(0.04, "cm"),
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
    panel.grid.major = element_line(color = "grey95", size = .25, linetype = "solid"),
    panel.grid.minor = element_line(color = "white", size = .85, linetype = "solid"),
    panel.background = element_rect(fill = "white"),
    legend.key.height = unit(26, "mm"),
    legend.key.width = unit(1.3, "mm"),
    legend.position = "none",
    plot.margin = ggplot2::margin(
      t = 0.0, # 顶部边缘距离
      r = 0.1, # 右边边缘距离
      b = 0.0, # 底部边缘距离
      l = 0.1, # 左边边缘距离
      unit = "cm"
    )
  ) +
  geom_hline(yintercept = 0, color = scales::alpha("grey70", 1), size = .25, linetype = "solid") +
  scale_y_continuous(limits = c(-x_bound, x_bound), expand = c(0, 0.1)) +
  scale_x_discrete(limits = rev(levels(shap_long_hd$variable)), expand = c(0.07, 0.1))


# png( paste0(filename, "_imp_bar.png" ),width = 10,height = 10.6,res = 800,units = "cm")
# imp_bar
# dev.off()
png(paste0(filename, "_imp_bar.png"), width = 3.8, height = 4, res = 700, units = "cm")
imp_bar
dev.off()

pdf(paste0(filename, "_imp_bar.pdf"), width = 3.1, height = 4)
imp_bar
dev.off()

#####################################' @5.long-term-KM
KMdf <- surv %>%
  t() %>%
  as.data.frame()
KMdf$time <- time
rownames(KMdf) <- time
colnames(KMdf) <- c("Surv", "time")
Est_surv <- KMdf[which(KMdf$time == Timepoint), "Surv"]

KMplot <- ggplot(KMdf, aes(x = time, y = Surv)) +
  geom_line(size = .7, col = scales::alpha("#2C8EAE", 1)) +
  geom_vline(
    xintercept = Timepoint, linetype = "dashed", size = .4,
    col = scales::alpha("#D1D2D4", 1)
  ) +
  annotate("segment",
    x = Timepoint + 180, xend = Timepoint,
    y = Est_surv + 0.015, yend = Est_surv + 0.0009,
    size = .6, arrow = arrow(length = unit(0.1, "inches")),
    alpha = .85, color = "#ff8c3e"
  ) +
  annotate("text",
    x = 310, y = range(KMdf$Surv)[1] - 0.02,
    label = paste0(
      Timepoint, "-day \nAll-cause mortality = ",
      percent(round(1 - Est_surv, 5), 0.001)
    ),
    fontface = "bold.italic", colour = "grey40", size = 1.7
  ) +
  theme_bw(base_rect_size = .7) +
  scale_color_manual(values = c("#E36258", "#9FD5D6")) +
  scale_y_continuous(limits = c(0.15, 1)) +
  labs(y = "Survival Probability", x = "Follow-up Time (Day)") +
  ggtitle("HF&RD-CPS system") +
  scale_y_continuous(
    limits = c(range(KMdf$Surv)[1] - 0.05, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme(
    legend.position = c("none"),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    # axis.ticks = element_blank(),
    axis.ticks = element_line(size = .3),
    axis.ticks.length = unit(0.06, "cm"),
    plot.margin = ggplot2::margin(
      t = 0.0,
      r = 0.1,
      b = 0.0,
      l = 0.0,
      unit = "cm"
    ),
    axis.text.x = element_text(size = 3.4, colour = "grey15", face = "bold", vjust = 2),
    axis.text.y = element_text(size = 3.4, colour = "grey15", face = "bold", hjust = 1.5),
    axis.title.x = element_text(size = 4.2, colour = "grey15", face = "bold", vjust = 4),
    axis.title.y = element_text(size = 4.2, colour = "grey15", face = "bold", vjust = -4),
    plot.title = element_text(size = 4.2, hjust = 0.5, colour = "darkred", face = "bold", vjust = -3),
    panel.grid.major = element_line(color = "grey95", size = .25, linetype = "solid"),
    panel.grid.minor = element_line(color = "grey95", size = .15, linetype = "solid"),
    panel.background = element_rect(fill = "white")
  )

png(paste0(filename, "_KMplot.png"), width = 5.25, height = 4, res = 700, units = "cm")
KMplot
dev.off()


pdf(paste0(filename, "_KMplot.pdf"), width = 3.95, height = 4)
KMplot
dev.off()


####################################' @6.可视化(未来生存率的95%可信区间)
se.pred <- m.pred$std.err
Est_surv <- KMdf[which(KMdf$time == Timepoint), "Surv"]
Est_surv_6 <- KMdf[which(KMdf$time == 180), "Surv"]
Est_surv_12 <- KMdf[which(KMdf$time == 365), "Surv"]
Est_surv_24 <- KMdf[which(KMdf$time == 730), "Surv"]
Est_surv_30 <- KMdf[which(KMdf$time == 900), "Surv"]

SurvPro <- surv %>%
  t() %>%
  as.data.frame()
SurvPro$time <- time
rownames(SurvPro) <- time
SurvPro$se <- se.pred[1:length(se.pred)]
colnames(SurvPro) <- c("Surv", "time", "se")
SurvPro <- SurvPro %>% filter(time == c(180) |
  time == c(365) |
  time == c(730) |
  time == c(900) |
  time == c(Timepoint))
SurvPro$time <- factor(as.character(SurvPro$time), levels = c(180, 365, 730, 900, Timepoint))

Survplot <- ggplot(SurvPro, aes(y = Surv, x = time, fill = time)) +
  geom_col(
    position = position_dodge(0.8), width = .5, size = .4,
    color = "black"
  ) +
  geom_errorbar(aes(x = time, ymin = Surv - se, ymax = Surv + se),
    width = 0.2,
    size = .3, color = "black"
  ) +
  theme_bw(base_rect_size = .7) +
  labs(x = "Time Point (month)", y = "Survival Probability") +
  ggtitle("95% Credible Interval") +
  theme(
    plot.margin = ggplot2::margin(
      t = 0.0, 
      r = 0.1, 
      b = 0.0, 
      l = 0.0, 
      unit = "cm"
    ),
    axis.text.x = element_text(size = 3.4, colour = "grey15", face = "bold", vjust = 6),
    axis.text.y = element_text(size = 3.4, colour = "grey15", face = "bold", hjust = 3),
    axis.title.x = element_text(size = 4.2, colour = "grey15", face = "bold", vjust = 6),
    axis.title.y = element_text(size = 4.2, colour = "grey15", face = "bold", vjust = -1),
    plot.title = element_text(size = 4.2, hjust = 0.5, colour = "darkred", face = "bold", vjust = -3),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey95", size = .23, linetype = "solid"),
    panel.grid.minor = element_line(color = "white", size = .13, linetype = "solid"),
    panel.background = element_rect(fill = "white"),
    strip.text = element_text(size = 12),
    axis.ticks.y = element_line(size = .3),
    axis.ticks.length.y = unit(0.05, "cm"),
    axis.ticks.x = element_blank()
  ) +
  scale_fill_manual(values = c("#B0997F", "#93a8ac", "#ffc857", "#61a5c2", "#119da4", "#FF6666")) +
  scale_y_continuous(expand = c(0, 0.09)) +
  facet_stereo() +
  facet_zoom(ylim = c(sort(SurvPro$Surv)[1] - 0.02, .87), zoom.size = 2, shrink = T)

png(paste0(filename, "_Survplot.png"), width = 5.9, height = 4, res = 700, units = "cm")
Survplot
dev.off()


pdf(paste0(filename, "_Survplot.pdf"), width = 4.5, height = 3.55)
Survplot
dev.off()


print(Risk)
print(round(SurvPro["180", "Surv"], 4))
print(round(SurvPro["365", "Surv"], 4))
print(round(SurvPro["730", "Surv"], 4))
print(round(SurvPro["900", "Surv"], 4))


####################################'@end!!!!!!!!





















