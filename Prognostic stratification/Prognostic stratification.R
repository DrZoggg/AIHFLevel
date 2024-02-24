


if(T){
  options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) ## 设置镜像
  options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/") ## 设置镜像
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
  library(party)
  library(xgboost)
  library(uwot)
  library(umap)
  library(Rtsne)
  library(psych)
  library(rgl)
  library(vegan)
  library(trinROC)
  library(rgl)
  options(stringsAsFactors = FALSE) 
}






#############################################################################################
################################################'@Survival-decision-tree
#############################################################################################

load("val_dd_list.rda")
rm(fit, val_dd_list, task_test, task_train, task_meta)
load("HF_metadata.rda")
df <- dataset3 %>% dplyr::select(-c("Name", "ID"))
rm(dataset3)
load("val_dd_list_all.rda")
df$ID <- c(paste0("G", c(1:712)))
row.names(df) <- df$ID
df <- df %>% dplyr::select(-c("ID"))
df <- df %>% mutate(
  CKD_stage = case_when(
    CKD_stage %in% "I" ~ 1, CKD_stage %in% "II" ~ 2,
    CKD_stage %in% "IIIa" | CKD_stage %in% "IIIb" ~ 3,
    CKD_stage %in% "IV" ~ 4, CKD_stage %in% "V" ~ 5
  )
)
df <- df[, c("ACM", "ACM_Censor", rid)]
df <- apply(df, 2, as.numeric) %>% as.data.frame()

expr.surv <- df
lrn_xgboost <- lrn("surv.xgboost",
  max_depth = 3, eta = 0.01, gamma = 0, subsample = 0.75, colsample_bytree = 0.8,
  min_child_weight = 3, lambda = 1, alpha = 0, nrounds = 500
)
measures <- c("surv.cindex")
resampling_HD <- rsmp("holdout", ratio = 0.7)
set.seed(1771)
rr_HD <- resample(as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor"), lrn_xgboost, resampling_HD, store_models = T)
index_HD <- msrs(measures) %>% rr_HD$aggregate()
index_HD
fit2 <- rr_HD$learners[[1]]


rr_HD$resampling$train_set(1)
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task$select(rid)
task_train <- task$filter(rr_HD$resampling$train_set(1))
task_train
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task$select(rid)
task_test <- task$filter(rr_HD$resampling$test_set(1))
task_test
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task$select(rid)
task_meta <- task
task_meta

fit2 <- rr_HD$learners[[1]]
fit2$predict(task_test)$score(msrs(measures))
fit2$predict(task_train)$score(msrs(measures))
fit2$predict(task_meta)$score(msrs(measures))

val_dd_list2 <- list()
val_dd_list2[[1]] <- task_train$data()
val_dd_list2[[2]] <- task_test$data()
val_dd_list2[[3]] <- task_meta$data()
names(val_dd_list2) <- c("Train_1", "Test_2", "Meta_3")

rs2 <- lapply(val_dd_list2, function(x) {
  dat.mat <- x %>% as.matrix()
  dtrain <- list(
    data = dat.mat[, 3:ncol(dat.mat), drop = F],
    label = dat.mat[, "ACM"] * (-(-1)^(as.numeric(dat.mat[, "ACM_Censor"])))
  )
  Dtrain <- xgb.DMatrix(dtrain$data, label = dtrain$label)
  RS <- predict(fit2$model, newdata = Dtrain, type = "risk")
  cbind(x[, 1:2], RS = RS)
})

C_index <- sapply(rs2, function(x) {
  concordance.index(scale(x$RS),
    surv.time = x$ACM,
    surv.event = x$ACM_Censor,
    method = "noether"
  )$c.index
})
C_index

for (i in names(rs2)) {
  if (identical(val_dd_list_all[[i]]$ACM_Censor, rs2[[i]]$ACM_Censor)) {
    val_dd_list_all[[i]] <- cbind(val_dd_list_all[[i]], RS = rs2[[i]][, 3])
  } else {
    break
  }
  colnames(val_dd_list_all[[i]])[ncol(val_dd_list_all[[i]])] <- "RS"
}

pbc <- rs2$Meta_3
set.seed(5000)
tree <- party::ctree(Surv(ACM, ACM_Censor) ~ RS,
  data = pbc,
  controls = party::ctree_control(maxdepth = 2)
)
plot(tree, type = "extended")
pdf("tree_plot.pdf", width = 7, height = 5)
plot(tree, type = "extended")
dev.off()
#' @(1.548,0.435)
save(tree, file = "tree.rda")

lp.ctree_meta <- predict(tree, newdata = rs2$Meta_3, type = "node")
table(lp.ctree_meta)
res.cat_meta <- rs2$Meta_3 %>% cbind(., Risk = lp.ctree_meta)
lp.ctree_test <- predict(tree, newdata = rs2$Test_2, type = "node")
table(lp.ctree_meta)
lp.ctree_train <- predict(tree, newdata = rs2$Train_1, type = "node")
table(lp.ctree_train)


lp.ctree <- lapply(rs2, function(x) {
  lp.ctree_pred <- predict(tree, newdata = x, type = "node")
  lp.ctree <- cbind(x, group = lp.ctree_pred)
  lp.ctree <- lp.ctree %>% mutate(
    group = case_when(
      group %in% 3 ~ "Low risk",
      group %in% 4 ~ "Intermediate risk",
      group %in% 5 ~ "High risk"
    )
  )
  return(lp.ctree)
})

table(lp.ctree$Test_2$group)
# High risk Intermediate risk          Low risk
# 23                61               130
table(lp.ctree$Meta_3$group)
# High risk Intermediate risk          Low risk
# 108               177               427


#############'@===KM_plot
KMPLOT <- lapply(lp.ctree, function(x) {
  group <- factor(x$group, levels = c("Low risk", "Intermediate risk", "High risk"))
  my.surv <- Surv(x$ACM, x$ACM_Censor)
  data.survdiff <- survdiff(my.surv ~ group)
  p.val <- 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  p.val
  KMplot <- ggplot(x, aes(time = ACM, status = ACM_Censor, color = group)) +
    geom_km(size = 1.5) +
    geom_kmticks(size = 1.2) +
    theme_bw(base_rect_size = 2.3) +
    scale_color_manual(values = c("#2C8EAE", "#5AB4C2", "#D1D2D4")) +
    scale_y_continuous(limits = c(0.15, 1)) +
    labs(y = "ACM", x = "Time (months)") +
    ggtitle("AIHFLevel") +
    theme(
      legend.position = c(0.21, 0.13),
      legend.text = element_text(size = 12, colour = "grey30", face = "bold.italic"),
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.key = element_rect(fill = "transparent", colour = "transparent"),
      axis.text = element_text(size = 12, colour = "grey25"),
      axis.title = element_text(size = 14, colour = "darkred", face = "bold"),
      plot.title = element_blank(),
      panel.grid.major = element_line(color = "#cacfd2", size = .25, linetype = "dashed"),
      panel.background = element_rect(fill = scales::alpha("#F0F8FF", .7)),
      panel.grid.minor = element_line(color = "white", size = .7, linetype = "solid")
    ) +
    annotate("text",
      x = 0.07, y = 0.38,
      label = paste(pval = ifelse(p.val < 0.0001, "Log-rank\np < 0.0001",
        paste("p = ", round(p.val, 4), sep = "")
      )),
      hjust = 0, size = 4, color = "grey20", fontface = "bold.italic"
    ) +
    coord_cartesian(clip = "on")
  KMplot
  return(KMplot)
})


for (i in names(KMPLOT)) {
  ggsave(KMPLOT[[i]],
    filename = paste0("KMplot_ctree ", i, " .pdf"),
    width = 5.4, height = 4.6
  )
}

lp.ctree2 <- lapply(val_dd_list_all, function(x) {
  lp.ctree_pred <- predict(tree, newdata = x, type = "node")
  lp.ctree <- cbind(x, group = lp.ctree_pred)
  lp.ctree <- lp.ctree %>% mutate(
    group = case_when(
      group %in% 3 ~ "Low risk",
      group %in% 4 ~ "Intermediate risk",
      group %in% 5 ~ "High risk"
    )
  )
  return(lp.ctree)
})




#############'@===dimensionality-reduction-analysis

my <- lp.ctree2$Test_2 %>%
  dplyr::select(group, RS, rid)
str(my)
pc <- principal(my %>% dplyr::select(RS, rid) %>% as.matrix(),
  nfactors = 2,
  rotate = "cluster"
)
res <- pc$scores %>% as.data.frame()
res$group <- my[, "group"]

col <- c("#2C8EAE", scales::alpha("#5AB4C2", 0.5), "#D1D2D4")
ggplot(res, aes(RC1, RC2, color = group, fill = group)) +
  scale_fill_manual(values = col) +
  geom_point(size = 2.5, shape = 21, color = "grey90") +
  labs(x = "UMAP1", y = "UMAP2") +
  theme_bw(base_rect_size = 1.5) +
  theme(
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 14, colour = "darkred", face = "bold"),
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.grid.major = element_line(color = "grey95", size = .5, linetype = "solid"),
    panel.grid.minor = element_line(color = "white", size = .3, linetype = "solid"),
    panel.background = element_rect(fill = "white"),
    legend.position = "top"
  )

my <- lp.ctree2$Test_2 %>%
  dplyr::select(group, RS, rid)
str(my)
brca.dist <- vegdist(my %>% dplyr::select(RS, rid) %>% as.matrix(),
  method = "cao"
)
pc <- cmdscale(brca.dist, eig = FALSE)
res <- pc %>% as.data.frame()
res$group <- my$group
col <- c("#2C8EAE", scales::alpha("#5AB4C2", 0.6), scales::alpha("#D1D2D4", .6))
P1 <- ggplot(res, aes(V1, V2, color = group, fill = group)) +
  scale_fill_manual(values = col) +
  geom_point(size = 2.8, shape = 21, color = "grey88") +
  labs(x = "PCoA1", y = "PCoA2") +
  theme_bw(base_rect_size = 1.5) +
  theme(
    axis.text = element_text(size = 11, colour = "grey35"),
    axis.title = element_text(size = 14, colour = "darkred", face = "bold"),
    plot.title = element_text(size = 15, hjust = 0.5),
    panel.grid.major = element_line(color = "grey95", size = .5, linetype = "solid"),
    panel.grid.minor = element_line(color = "white", size = .3, linetype = "solid"),
    panel.background = element_rect(fill = "white"), 
    legend.position = c(0.19, 0.9),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 11, colour = "grey10")
  ) +
  scale_x_continuous(expand = c(0.03, 0.09)) +
  scale_y_continuous(expand = c(0.02, 0.09))
P1
dev.off()

ggsave(P1,
  filename = "pcoa_ctree.pdf",
  width = 5.2, height = 4.8
)


PCAplot <- lapply(lp.ctree2, function(x) {
  my <- x %>%
    dplyr::select(group, RS, rid)
  str(my)
  brca.dist <- vegdist(my %>% dplyr::select(RS, rid) %>% as.matrix(),
    method = "cao"
  )
  pc <- cmdscale(brca.dist, eig = FALSE)
  res <- pc %>% as.data.frame()
  res$group <- my$group
  col <- c("#2C8EAE", scales::alpha("#5AB4C2", 0.6), scales::alpha("#D1D2D4", .6))
  P1 <- ggplot(res, aes(V1, V2, color = group, fill = group)) +
    scale_fill_manual(values = col) +
    geom_point(
      size = 3.6, shape = 21,
      color = "grey80"
    ) +
    labs(x = "t-SNE1", y = "t-SNE2") +
    theme_bw(base_rect_size = 2.3) +
    theme(
      axis.text = element_text(size = 11, colour = "grey5"),
      axis.title = element_text(size = 14, colour = "darkred", face = "bold"),
      plot.title = element_text(size = 15, hjust = 0.5),
      panel.grid.major = element_line(color = "white", size = .9, linetype = "solid"),
      panel.grid.minor = element_line(color = "white", size = .5, linetype = "solid"),
      panel.background = element_rect(fill = scales::alpha("#F0F8FF", .7)),
      legend.position = c(0.19, 0.9),
      legend.background = element_rect(fill = "transparent"),
      legend.title = element_blank(),
      legend.key = element_rect(fill = "transparent", colour = "transparent"),
      legend.text = element_text(size = 12, colour = "grey20", face = "bold")
    ) +
    scale_x_continuous(expand = c(0.03, 0.09)) +
    scale_y_continuous(expand = c(0.02, 0.09))

  return(P1)
})

PCAplot[[3]]


for (i in names(PCAplot)) {
  ggsave(PCAplot[[i]],
    filename = paste0("pcoa_ctree ", i, " .pdf"),
    width = 5.4, height = 4.6
  )
}





#############'@===trinROC

my <- lp.ctree2$Meta_3 %>%
  dplyr::select(group, RS, rid)
str(my)
dat <- my %>% dplyr::select(group, RS)
dat$group <- factor(dat$group, levels = c(
  "Low risk",
  "Intermediate risk",
  "High risk"
))

green <- "#C7EAB2"
cyan <- "#5FC1C2"
blue <- "#1B90BE"
p_top <- ggplot(dat, aes(x = RS, color = group, fill = group)) +
  geom_density() +
  scale_color_manual(values = c(scales::alpha(green, 0.7), scales::alpha(cyan, 0.7), scales::alpha(blue, 0.7))) +
  scale_fill_manual(values = c(scales::alpha(green, 0.7), scales::alpha(cyan, 0.7), scales::alpha(blue, 0.7))) +
  theme_classic() +
  xlab(paste0("HF&RD-CPS")) +
  ylab(NULL) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 11, colour = "black", face = "bold.italic"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  geom_rug()
p_top

p.val <- kruskal.test(RS ~ group, data = dat)
p.lab <- paste0(
  "P",
  ifelse(p.val$p.value < 0.0001, " < 0.0001",
    paste0(" = ", round(p.val$p.value, 4))
  )
)

p_bot <- ggplot(dat, aes(group, RS, fill = group)) +
  geom_boxplot(aes(col = group)) +
  scale_fill_manual(values = c(green, cyan, blue)) +
  scale_color_manual(values = c(green, cyan, blue)) +
  xlab(NULL) +
  ylab("HF&RD-CPS") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10, color = "grey25", face = "bold"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = c(0.00, 0.05)) +
  coord_cartesian(clip = "off") +
  annotate(
    geom = "text",
    x = 1.5, hjust = 1, y = max(dat$RS), size = 4, angle = 270, fontface = "bold",
    label = p.lab
  ) +
  coord_flip()
p_bot

stat <- ggplot_build(p_bot)$data[[1]]
p_bot <- p_bot + geom_segment(
  data = stat, aes(x = xmin, xend = xmax, y = middle, yend = middle), size = .8,
  color = "white", inherit.aes = F
)
p_bot
library(aplot)
library(patchwork)
p <- p_top %>% insert_bottom(p_bot, height = 0.4)
p

pdf(file = "plotVUS_data.pdf", width = 6, height = 3.5)
p
invisible(dev.off())


my <- lp.ctree2$Test_2 %>%
  dplyr::select(group, RS, rid)
str(my)
dat <- my %>% dplyr::select(group, RS)
dat$group <- factor(dat$group, levels = c(
  "Low risk",
  "Intermediate risk",
  "High risk"
))
roc.eda(
  dat = dat, type = "trinormal",
  plotVUS = T, saveVUS = T, sep.dens = F, scatter = T
)
rgl.postscript(
  filename = "plotVUS Test_2.pdf",
  fmt = "pdf"
)

# data: Low risk, Intermediate risk and High risk
#
# ROC test statistic: 1063.675, ROC p.value: 0
# VUS test statistic:  14.759 ,  VUS p.value:  0
#
# trinormal VUS:  0.876
#
# Parameters:
#   a	b	c	d
# 3.1431	-6.5434	0.2597	1.316



my <- lp.ctree2$Meta_3 %>%
  dplyr::select(group, RS, rid)
str(my)
dat <- my %>% dplyr::select(group, RS)
dat$group <- factor(dat$group, levels = c(
  "Low risk",
  "Intermediate risk",
  "High risk"
))
roc.eda(
  dat = dat, type = "trinormal",
  plotVUS = T, saveVUS = T, sep.dens = F, scatter = T
)
rgl.postscript(
  filename = "plotVUS Meta_3.pdf",
  fmt = "pdf"
)

# data: Low risk, Intermediate risk and High risk
#
# ROC test statistic: 8682.598, ROC p.value: 0
# VUS test statistic:  29.331 ,  VUS p.value:  0
#
# trinormal VUS:  0.867
#
# Parameters:
#   a	b	c	d
# 3.1341	-6.558	0.1817	1.2395



plotVUS_PLOT <- lapply(lp.ctree2, function(x) {
  my <- x %>%
    dplyr::select(group, RS, rid)
  dat <- my %>% dplyr::select(group, RS)
  dat$group <- factor(dat$group, levels = c(
    "Low risk",
    "Intermediate risk",
    "High risk"
  ))
  green <- "#C7EAB2"
  cyan <- "#5FC1C2"
  blue <- "#1B90BE"
  p_top <- ggplot(dat, aes(x = RS, color = group, fill = group)) +
    geom_density() +
    scale_color_manual(values = c(scales::alpha(green, 0.7), scales::alpha(cyan, 0.7), scales::alpha(blue, 0.7))) +
    scale_fill_manual(values = c(scales::alpha(green, 0.7), scales::alpha(cyan, 0.7), scales::alpha(blue, 0.7))) +
    theme_classic() +
    xlab(paste0("AIHFLevel")) +
    ylab(NULL) +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      axis.text.x = element_text(size = 11, color = "black"),
      axis.title = element_text(size = 11, colour = "black", face = "bold.italic"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      rect = element_rect(fill = scales::alpha("#F7FCFF", 1)),
      panel.grid = element_line(color = "transparent")
    ) +
    geom_rug()

  p.val <- kruskal.test(RS ~ group, data = dat)
  p.lab <- paste0(
    "P",
    ifelse(p.val$p.value < 0.0001, " < 0.0001",
      paste0(" = ", round(p.val$p.value, 4))
    )
  )

  p_bot <- ggplot(dat, aes(group, RS, fill = group)) +
    geom_boxplot(aes(col = group)) +
    scale_fill_manual(values = c(green, cyan, blue)) +
    scale_color_manual(values = c(green, cyan, blue)) +
    xlab(NULL) +
    ylab("AIHFLevel") +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 10, color = "grey25", face = "bold"),
      panel.background = element_rect(fill = "transparent"),
      rect = element_rect(fill = scales::alpha("#F7FCFF", 1)),
      panel.grid = element_line(color = "transparent"),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    ) +
    scale_y_continuous(expand = c(0.00, 0.05)) +
    scale_x_discrete(position = "top") +
    coord_cartesian(clip = "off") +
    annotate(
      geom = "text",
      x = 1.5, hjust = 1, y = max(dat$RS), size = 4, angle = 270, fontface = "bold",
      label = p.lab
    ) +
    coord_flip()

  stat <- ggplot_build(p_bot)$data[[1]]
  p_bot <- p_bot + geom_segment(
    data = stat, aes(x = xmin, xend = xmax, y = middle, yend = middle), size = .8,
    color = "white", inherit.aes = F
  )
  library(aplot)
  library(patchwork)
  p <- p_top %>% insert_bottom(p_bot, height = 0.4)

  return(p)
})

plotVUS_PLOT[[3]]

for (i in names(plotVUS_PLOT)) {
  ggsave(plotVUS_PLOT[[i]],
    filename = paste0("plotVUS_data ", i, " .pdf"),
    width = 9.5, height = 3.5
  )
}







##########################' @pie!plot

p_fisher <- data.frame()
for (i in colnames(clin)[-c(1:3)]) {
  rt <- clin %>% dplyr::select(group, i)
  tab_classify <- as.data.frame.array(table(rt$group, rt[, 2]))
  tab_classify
  p.fisher <- fisher.test(tab_classify, workspace = 1e9)$p.value
  p_fisher <- rbind(p_fisher, data.frame(id = i, P = p.fisher))
}

clin_long <- clin %>%
  pivot_longer(ACM_Censor:HTN,
    names_to = "subgroup",
    values_to = "value"
  )

p <- ggplot(clin_long, aes(x = group, fill = value)) +
  geom_bar(width = 0.9, position = "fill") +
  scale_y_continuous(labels = percent) +
  labs(
    title = "mUC",
    y = "Patients Percentage"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  facet_wrap(~subgroup)
p


clin2 <- clin %>% dplyr::select(-ACM)
gname <- "group"
vname <- setdiff(colnames(clin2), gname)
vname
pie.high <- pie.low <- pie.mid <- list()
fisher.p <- c()
for (i in vname) {
  tmp <- as.data.frame.array(table(clin2[, gname], clin2[, i]))
  p <- fisher.test(tmp, workspace = 1e9)$p.value
  names(p) <- i
  fisher.p <- append(fisher.p, p)

  tmp <- table(clin2[, gname], clin2[, i])
  pie.dat <-
    tmp %>%
    as.data.frame() %>%
    group_by(Var1) %>%
    mutate(Pct = Freq / sum(Freq)) %>%
    as.data.frame()
  pie.high[[i]] <- pie.dat[which(pie.dat$Var1 == "High risk"), ]
  pie.mid[[i]] <- pie.dat[which(pie.dat$Var1 == "Intermediate risk"), ]
  pie.low[[i]] <- pie.dat[which(pie.dat$Var1 == "Low risk"), ]
}


gname <- "group"
vname <- setdiff(colnames(clin2), gname)
vname
pie.high <- pie.low <- pie.mid <- list()
fisher.p <- c()


i <- vname[12]

tmp <- as.data.frame.array(table(clin2[, gname], clin2[, i]))
p <- format(fisher.test(tmp, workspace = 1e9)$p.value, digits = 3)
p
names(p) <- i
fisher.p <- append(fisher.p, p)
fisher.p

tmp <- table(clin2[, gname], clin2[, i])
pie.dat <-
  tmp %>%
  as.data.frame() %>%
  group_by(Var1) %>%
  mutate(Pct = Freq / sum(Freq)) %>%
  as.data.frame()

pie.high[[i]] <- pie.dat[which(pie.dat$Var1 == "High risk"), ]
pie.mid[[i]] <- pie.dat[which(pie.dat$Var1 == "Intermediate risk"), ]
pie.low[[i]] <- pie.dat[which(pie.dat$Var1 == "Low risk"), ]
fisher.p

black <- "#1E1E1B"
blue <- "#3C4E98"
yellow <- "#E4DB36"
orange <- "#E19143"
green <- "#57A12B"
cherry <- "#8D3A86"
red <- "#E36258"
lblue <- "#5AB4C2"

black <- "#1E1E1B"
yellow <- "#a26769"
orange <- "#d5b9b2"
blue <- "#FF6666"
green <- "#119da4"
cherry <- "#61a5c2"
red <- "#ffc857"
lblue <- "#B0997F"

ACM_Censor.col <- c("grey80", black)
Age.col <- c(yellow, orange)
CKD_stage.col <- alpha(blue, c(0.1, 0.35, 0.6, 0.8, 1))
HF_subtype.col <- alpha(green, c(0.2, 0.6, 1))
NYHA.col <- alpha(cherry, c(0.2, 0.5, 0.75, 1))
Gender.col <- alpha(red, c(0.5, 1))
CAD.col <- alpha(lblue, c(1, 0.5))



if (T) {
  pdf("pie_ctree.pdf", width = 14, height = 4.5)
  showLayout <- F
  layout(matrix(c(
    1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8,
    9, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 16, 16, 16,
    9, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 16, 16, 16,
    17, 17, 17, 18, 18, 18, 19, 19, 19, 20, 20, 20, 21, 21, 21, 22, 22, 22, 23, 23, 23, 24, 24, 24,
    17, 17, 17, 18, 18, 18, 19, 19, 19, 20, 20, 20, 21, 21, 21, 22, 22, 22, 23, 23, 23, 24, 24, 24,
    25, 25, 25, 26, 26, 26, 27, 27, 27, 28, 28, 28, 29, 29, 29, 30, 30, 30, 31, 31, 31, 32, 32, 32,
    25, 25, 25, 26, 26, 26, 27, 27, 27, 28, 28, 28, 29, 29, 29, 30, 30, 30, 31, 31, 31, 32, 32, 32,
    33, 33, 33, 34, 34, 34, 35, 35, 35, 36, 36, 36, 37, 37, 37, 38, 38, 38, 39, 39, 39, 40, 40, 40,
    41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41,
    41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41,
    41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41,
    41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41,
    41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41
  ),
  byrow = T, nrow = 13
  ))

  if (showLayout) {
    layout.show(n = 41)
  }


  if (T) {
    par(bty = "n", mgp = c(0, 0, 0), mar = c(0, 0, 0, 0), lwd = 2)
    plot(1, 1,
      xlab = "", xaxt = "n",
      ylab = "", yaxt = "n"
    )
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "black")
    text((par("usr")[1] + par("usr")[2]) / 2,
      (par("usr")[3] + par("usr")[4]) / 2,
      "HF&RD",
      cex = 1.7, col = "white"
    )

    plot(1, 1,
      xlab = "", xaxt = "n",
      ylab = "", yaxt = "n"
    )
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "black")
    text((par("usr")[1] + par("usr")[2]) / 2,
      (par("usr")[3] + par("usr")[4]) / 2,
      "Status",
      cex = 1.7, col = "white"
    )

    plot(1, 1,
      xlab = "", xaxt = "n",
      ylab = "", yaxt = "n"
    )
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "black")
    text((par("usr")[1] + par("usr")[2]) / 2,
      (par("usr")[3] + par("usr")[4]) / 2,
      "Age",
      cex = 1.7, col = "white"
    )

    plot(1, 1,
      xlab = "", xaxt = "n",
      ylab = "", yaxt = "n"
    )
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "black")
    text((par("usr")[1] + par("usr")[2]) / 2,
      (par("usr")[3] + par("usr")[4]) / 2,
      "CKD stage",
      cex = 1.7, col = "white"
    )

    plot(1, 1,
      xlab = "", xaxt = "n",
      ylab = "", yaxt = "n"
    )
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "black")
    text((par("usr")[1] + par("usr")[2]) / 2,
      (par("usr")[3] + par("usr")[4]) / 2,
      "Subtype",
      cex = 1.7, col = "white"
    )

    plot(1, 1,
      xlab = "", xaxt = "n",
      ylab = "", yaxt = "n"
    )
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "black")
    text((par("usr")[1] + par("usr")[2]) / 2,
      (par("usr")[3] + par("usr")[4]) / 2,
      "NYHA",
      cex = 1.7, col = "white"
    )
    plot(1, 1,
      xlab = "", xaxt = "n",
      ylab = "", yaxt = "n"
    )
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "black")
    text((par("usr")[1] + par("usr")[2]) / 2,
      (par("usr")[3] + par("usr")[4]) / 2,
      "CAD",
      cex = 1.7, col = "white"
    )
    plot(1, 1,
      xlab = "", xaxt = "n",
      ylab = "", yaxt = "n"
    )
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "black")
    text((par("usr")[1] + par("usr")[2]) / 2,
      (par("usr")[3] + par("usr")[4]) / 2,
      "Gender",
      cex = 1.7, col = "white"
    )
  }

  table(clin$group)
  # High risk Intermediate risk          Low risk
  # 108               177               427

  if (T) {
    plot(1, 1,
      xlab = "", xaxt = "n",
      ylab = "", yaxt = "n"
    )
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "black")
    text((par("usr")[1] + par("usr")[2]) / 2,
      (par("usr")[3] + par("usr")[4]) / 2,
      "High risk\n(n = 108)",
      cex = 1.6, col = "white"
    )

    pie(pie.high$ACM_Censor$Pct,
      col = ACM_Censor.col,
      border = "white",
      radius = 1,
      labels = NA,
      init.angle = 90
    )
    symbols(0, 0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

    pie(pie.high$Age$Pct,
      col = Age.col,
      border = "white",
      radius = 1,
      labels = NA,
      init.angle = 90
    )
    symbols(0, 0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

    pie(pie.high$CKD_stage$Pct,
      col = CKD_stage.col,
      border = "white",
      radius = 1,
      labels = NA,
      init.angle = 90
    )
    symbols(0, 0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

    pie(pie.high$HF_subtype$Pct,
      col = HF_subtype.col,
      border = "white",
      radius = 1,
      labels = NA,
      init.angle = 90
    )
    symbols(0, 0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

    pie(pie.high$NYHA$Pct,
      col = NYHA.col,
      border = "white",
      radius = 1,
      labels = NA,
      init.angle = 90
    )
    symbols(0, 0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

    pie(pie.high$CAD$Pct,
      col = CAD.col,
      border = "white",
      radius = 1,
      labels = NA,
      init.angle = 90
    )
    symbols(0, 0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

    pie(pie.high$Gender$Pct,
      col = Gender.col,
      border = "white",
      radius = 1,
      labels = NA,
      init.angle = 90
    )
    symbols(0, 0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

    abline(v = par("usr")[2.5], col = "black")
  }


  if (T) {
    plot(1, 1,
      xlab = "", xaxt = "n",
      ylab = "", yaxt = "n"
    )
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "black")
    text((par("usr")[1] + par("usr")[2]) / 2,
      (par("usr")[3] + par("usr")[4]) / 2,
      "Intermediate\n(n = 177)",
      cex = 1.6, col = "white"
    )


    pie(pie.mid$ACM_Censor$Pct,
      col = ACM_Censor.col,
      border = "white",
      radius = 1,
      labels = NA,
      init.angle = 90
    )
    symbols(0, 0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

    pie(pie.mid$Age$Pct,
      col = Age.col,
      border = "white",
      radius = 1,
      labels = NA,
      init.angle = 90
    )
    symbols(0, 0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

    pie(pie.mid$CKD_stage$Pct,
      col = CKD_stage.col,
      border = "white",
      radius = 1,
      labels = NA,
      init.angle = 90
    )
    symbols(0, 0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

    pie(pie.mid$HF_subtype$Pct,
      col = HF_subtype.col,
      border = "white",
      radius = 1,
      labels = NA,
      init.angle = 90
    )
    symbols(0, 0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

    pie(pie.mid$NYHA$Pct,
      col = NYHA.col,
      border = "white",
      radius = 1,
      labels = NA,
      init.angle = 90
    )
    symbols(0, 0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

    pie(pie.mid$CAD$Pct,
      col = CAD.col,
      border = "white",
      radius = 1,
      labels = NA,
      init.angle = 90
    )
    symbols(0, 0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

    pie(pie.mid$Gender$Pct,
      col = Gender.col,
      border = "white",
      radius = 1,
      labels = NA,
      init.angle = 90
    )
    symbols(0, 0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)


    abline(v = par("usr")[2.5], col = "black")
  }


  table(clin$group)
  # High risk Intermediate risk          Low risk
  # 108               177               427


  if (T) {
    plot(1, 1,
      xlab = "", xaxt = "n",
      ylab = "", yaxt = "n"
    )
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "black")
    text((par("usr")[1] + par("usr")[2]) / 2,
      (par("usr")[3] + par("usr")[4]) / 2,
      "Low risk\n(n = 427)",
      cex = 1.6, col = "white"
    )

    # Low group
    pie(pie.low$ACM_Censor$Pct,
      col = ACM_Censor.col,
      border = "white",
      radius = 1,
      labels = NA,
      init.angle = 90
    )
    symbols(0, 0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

    pie(pie.low$Age$Pct,
      col = Age.col,
      border = "white",
      radius = 1,
      labels = NA,
      init.angle = 90
    )
    symbols(0, 0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

    pie(pie.low$CKD_stage$Pct,
      col = CKD_stage.col,
      border = "white",
      radius = 1,
      labels = NA,
      init.angle = 90
    )
    symbols(0, 0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

    pie(pie.low$HF_subtype$Pct,
      col = HF_subtype.col,
      border = "white",
      radius = 1,
      labels = NA,
      init.angle = 90
    )
    symbols(0, 0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

    pie(pie.low$NYHA$Pct,
      col = NYHA.col,
      border = "white",
      radius = 1,
      labels = NA,
      init.angle = 90
    )
    symbols(0, 0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

    pie(pie.low$CAD$Pct,
      col = CAD.col,
      border = "white",
      radius = 1,
      labels = NA,
      init.angle = 90
    )
    symbols(0, 0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

    pie(pie.low$Gender$Pct,
      col = Gender.col,
      border = "white",
      radius = 1,
      labels = NA,
      init.angle = 90
    )
    symbols(0, 0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)


    abline(v = par("usr")[2], col = "black")
  }



  plot(1, 1,
    xlab = "", xaxt = "n",
    ylab = "", yaxt = "n"
  )
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "black")

  plot(1, 1,
    col = "white",
    xlab = "", xaxt = "n",
    ylab = "", yaxt = "n"
  )
  text((par("usr")[1] + par("usr")[2]) / 2,
    (par("usr")[3] + par("usr")[4]) / 2,
    paste0("p = ", fisher.p["ACM_Censor"]),
    cex = 1.5, col = "black"
  )
  abline(h = par("usr")[3], col = "black")

  plot(1, 1,
    col = "white",
    xlab = "", xaxt = "n",
    ylab = "", yaxt = "n"
  )
  text((par("usr")[1] + par("usr")[2]) / 2,
    (par("usr")[3] + par("usr")[4]) / 2,
    paste0("p = ", fisher.p["Age"]),
    cex = 1.5, col = "black"
  )
  abline(h = par("usr")[3], col = "black")

  plot(1, 1,
    col = "white",
    xlab = "", xaxt = "n",
    ylab = "", yaxt = "n"
  )
  text((par("usr")[1] + par("usr")[2]) / 2,
    (par("usr")[3] + par("usr")[4]) / 2,
    paste0("p = ", fisher.p["CKD_stage"]),
    cex = 1.5, col = "black"
  )
  abline(h = par("usr")[3], col = "black")

  plot(1, 1,
    col = "white",
    xlab = "", xaxt = "n", #
    ylab = "", yaxt = "n"
  ) #
  text((par("usr")[1] + par("usr")[2]) / 2,
    (par("usr")[3] + par("usr")[4]) / 2,
    paste0("p = ", fisher.p["HF_subtype"]),
    cex = 1.5, col = "black"
  )
  abline(h = par("usr")[3], col = "black")

  plot(1, 1,
    col = "white",
    xlab = "", xaxt = "n",
    ylab = "", yaxt = "n"
  )
  text((par("usr")[1] + par("usr")[2]) / 2,
    (par("usr")[3] + par("usr")[4]) / 2,
    paste0("p = ", fisher.p["NYHA"]),
    cex = 1.5, col = "black"
  )

  plot(1, 1,
    col = "white",
    xlab = "", xaxt = "n",
    ylab = "", yaxt = "n"
  )
  text((par("usr")[1] + par("usr")[2]) / 2, #
    (par("usr")[3] + par("usr")[4]) / 2,
    paste0("p = ", format(fisher.p["CAD"], digits = 3)),
    cex = 1.5, col = "black"
  )

  plot(1, 1,
    col = "white",
    xlab = "", xaxt = "n",
    ylab = "", yaxt = "n"
  )
  text((par("usr")[1] + par("usr")[2]) / 2,
    (par("usr")[3] + par("usr")[4]) / 2,
    paste0("p = ", fisher.p["Gender"]),
    cex = 1.5, col = "black"
  )


  abline(h = par("usr")[3], col = "black")
  abline(v = par("usr")[2.5], col = "black")




  plot(0, 0,
    col = "white",
    xlab = "", xaxt = "n", #
    ylab = "", yaxt = "n"
  ) #
  legend("topleft",
    legend = c(
      "Alive", "Dead",
      ">65", "≤65",
      "I", "II", "III", "IV", "V",
      "HFmrEF", "HFpEF", "HFrEF",
      "I", "II", "III", "IV",
      "Concomitant", "NO",
      "Female", "Male"
    ),
    fill = c(
      ACM_Censor.col,
      Age.col,
      CKD_stage.col,
      HF_subtype.col,
      NYHA.col,
      CAD.col,
      Gender.col
    ),
    border = NA, #
    bty = "n", #
    cex = 1.2,
    x.intersp = 0.05,
    y.intersp = 2,
    ncol = 4,
    text.width = 0.075,
    horiz = F
  )
}
invisible(dev.off())
ACM_Censor.col <- c("grey80", black)
Age.col <- c(yellow, orange)
CKD_stage.col <- alpha(blue, c(0.05, 0.25, 0.45, 0.65, 0.85, 1))
HF_subtype.col <- alpha(green, c(0.5, 0.7, 1))
NYHA.col <- alpha(cherry, c(0.4, 0.6, 0.8, 1))


table(clin$group, clin$Gender)
