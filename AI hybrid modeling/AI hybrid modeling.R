



if(T){
  options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) ##
  options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/") ## 
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
  options(stringsAsFactors = FALSE)
}




####################### ####################### #######################
#######################    1 Baseline Data      #######################
####################### ####################### #######################
load("HF_metadata.rda")
df <- dataset %>% dplyr::select(-c("Name", "ID"))
rm(dataset)
if (F) {
  library(gtsummary)
  suppressPackageStartupMessages(library(tidyverse))
  tbl_summary1 <- df %>%
    tbl_summary(
      by = HF_subtype,
      statistic = list(
        all_continuous() ~ "{mean} ({sd})",
        all_categorical() ~ "{n} / {N} ({p}%)"
        # all_categorical() ~ "{n} ({p}%)"
      ),
      missing_text = "(Missing)"
    ) %>%
    add_p() %>%
    add_overall() %>%
    modify_spanning_header(c("stat_1", "stat_2", , "stat_3") ~ "**Subtype**") %>%
    modify_caption("***Patient Characteristics***")

  tbl_summary1
  tbl_summary1 %>%
    as_flex_table() %>%
    flextable::save_as_docx(tbl_summary1, path = "baseline_data.docx")
}


####################### ####################### ####################### 
################    2 latent survival features      ###################
####################### ####################### ####################### 

####################'@=====cox
rt <- df
##### ACM
covariates <- colnames(rt)[-which(colnames(rt) %in% c("ACM", "ACM_Censor", "MACE", "MACE_Censor"))]
if (T) {
  df <- data.frame()
  for (i in covariates) {
    fit <- survival::coxph(Surv(ACM, ACM_Censor) ~ get(i), data = rt)
    fit1 <- summary(fit)
    coef <- coef(fit)
    HR <- exp(coef(fit))
    HR_l <- exp(confint(fit))[, 1]
    HR_H <- exp(confint(fit))[, 2]
    p_val <- fit1$coefficients[, 5]
    df <- rbind(
      df,
      cbind(i, coef, HR, HR_l, HR_H, p_val)
    )
  }
  name <- df$i
  unicox_ACM <- df[, -1] %>%
    apply(2, as.numeric) %>%
    as.data.frame()
  rownames(unicox_ACM) <- name
}

covariates <- colnames(rt)[-which(colnames(rt) %in% c("ACM", "ACM_Censor", "MACE", "MACE_Censor"))]
if (T) {
  univ_formulas <- lapply(
    covariates,
    function(x) as.formula(paste("Surv(ACM,ACM_Censor)~", x))
  )
  names(univ_formulas) <- covariates
  univ_models <- lapply(univ_formulas, function(x) {
    coxph(x, data = rt)
  })
  univ_results <- lapply(
    univ_models,
    function(fit) {
      fit1 <- summary(fit)
      coef <- coef(fit)
      HR <- exp(coef(fit))
      HR_l <- exp(confint(fit))[, 1]
      HR_H <- exp(confint(fit))[, 2]
      p_val <- fit1$coefficients[, 5]
      res <- c(coef, HR, HR_l, HR_H, p_val)
      names(res) <- c("coef", "HR", "HR_l", "HR_H", "p_val")
      return(res)
    }
  )
  unicox_ACM <- t(as.data.frame(univ_results, check.names = FALSE)) %>% as.data.frame()
}

i <- rownames(unicox_ACM)[unicox_ACM$p_val < 0.05]
unicox_ACM_sig <- unicox_ACM[i, ]



precox_multi <- rt[, c(i, "ACM", "ACM_Censor")]
fit <- coxph(Surv(ACM, ACM_Censor) ~ ., data = precox_multi)
fit1 <- summary(fit)
multi_ACM <- data.frame(
  coef = coef(fit),
  HR = exp(coef(fit)),
  HR_l = exp(confint(fit))[, 1],
  HR_H = exp(confint(fit))[, 2],
  p_val = fit1$coefficients[, 5]
)



##### MACE
covariates <- colnames(rt)[-which(colnames(rt) %in% c("ACM", "ACM_Censor", "MACE", "MACE_Censor"))]
if (T) {
  df <- data.frame()
  for (i in covariates) {
    fit <- survival::coxph(Surv(MACE, MACE_Censor) ~ get(i), data = rt)
    fit1 <- summary(fit)
    coef <- coef(fit)
    HR <- exp(coef(fit))
    HR_l <- exp(confint(fit))[, 1]
    HR_H <- exp(confint(fit))[, 2]
    p_val <- fit1$coefficients[, 5]
    df <- rbind(
      df,
      cbind(i, coef, HR, HR_l, HR_H, p_val)
    )
  }
  name <- df$i
  unicox_MACE <- df[, -1] %>%
    apply(2, as.numeric) %>%
    as.data.frame()
  rownames(unicox_MACE) <- name
}



covariates <- colnames(rt)[-which(colnames(rt) %in% c("ACM", "ACM_Censor", "MACE", "MACE_Censor"))]
if (T) {
  univ_formulas <- lapply(
    covariates,
    function(x) as.formula(paste("Surv(MACE,MACE_Censor)~", x))
  )
  names(univ_formulas) <- covariates
  univ_models <- lapply(univ_formulas, function(x) {
    coxph(x, data = rt)
  })
  univ_results <- lapply(
    univ_models,
    function(fit) {
      fit1 <- summary(fit)
      coef <- coef(fit)
      HR <- exp(coef(fit))
      HR_l <- exp(confint(fit))[, 1]
      HR_H <- exp(confint(fit))[, 2]
      p_val <- fit1$coefficients[, 5]
      res <- c(coef, HR, HR_l, HR_H, p_val)
      names(res) <- c("coef", "HR", "HR_l", "HR_H", "p_val")
      return(res)
    }
  )
  unicox_MACE <- t(as.data.frame(univ_results, check.names = FALSE)) %>% as.data.frame()
}

i <- rownames(unicox_MACE)[unicox_MACE$p_val < 0.05]
unicox_MACE_sig <- unicox_MACE[i, ]


precox_multi <- rt[, c(i, "MACE", "MACE_Censor")]
fit <- coxph(Surv(MACE, MACE_Censor) ~ ., data = precox_multi)
fit1 <- summary(fit)
multi_MACE <- data.frame(
  coef = coef(fit),
  HR = exp(coef(fit)),
  HR_l = exp(confint(fit))[, 1],
  HR_H = exp(confint(fit))[, 2],
  p_val = fit1$coefficients[, 5]
)


write.csv(unicox_ACM, file = "unicox_ACM.csv")
write.csv(unicox_MACE, file = "unicox_MACE.csv")
write.csv(multi_ACM, file = "multi_ACM.csv")
write.csv(multi_MACE, file = "multi_MACE.csv")
save(unicox_ACM, unicox_ACM_sig,
  unicox_MACE, unicox_MACE_sig,
  multi_ACM, multi_MACE,
  unicoxgene_ACM, unicoxgene_MACE,
  file = "ACM MACE cox result.RData"
)
 




####################' @=====KM
desc_number <- stat.desc(rt)
list <- colnames(rt)[rowSums(is.na(t(desc_number))) > 1]
list <- list[-which(list %in% c("MACE_Censor", "ACM_Censor"))]
rt <- rt[, c("MACE_Censor", "MACE", "ACM_Censor", "ACM", list)]
list2 <- list[-which(list %in% c("Gender", "NYHA", "HTN", "CKD_stage", "HF_subtype"))]

for (i in list2) {
  rt[, i] <- ifelse(rt[, i] %in% 1, "Y", "N")
}
rt$HTN <- ifelse(rt$HTN %in% c("0"), "N", "Y")
rt$CKD_stage <- ifelse(rt$CKD_stage %in% c("I", "II"), "I|II", ifelse(rt$CKD_stage %in% c("IIIa", "IIIb"), "IIIa|IIIb", "IV|V"))
rt$MACE_Censor <- as.numeric(as.character(rt$MACE_Censor))
rt$ACM_Censor <- as.numeric(as.character(rt$ACM_Censor))

picDir <- "clinic_OS_KM_uncontinuous_picture2"
dir.create(picDir)


#######'@Uncontinuous
outTab <- data.frame()
setwd(picDir)
getwd()

my.surv <- Surv(rt$ACM, rt$ACM_Censor)
pl <- pl_sig <- list()
ptable_KM_ACM_uncon <- data.frame(matrix(NA, ncol(rt) - 4, 5))
colnames(ptable_KM_ACM_uncon) <- c("ID", "pvalue", "HR", "CIlow", "CIup")

n <- 0
for (i in colnames(rt)[5:ncol(rt)]) {
  group <- rt[, i]
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)

  group <- factor(group)
  data.survdiff <- survdiff(my.surv ~ group)
  p.val <- 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR <- (data.survdiff$obs[2] / data.survdiff$exp[2]) / (data.survdiff$obs[1] / data.survdiff$exp[1])
  up95 <- exp(log(HR) + qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
  low95 <- exp(log(HR) - qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
  n <- n + 1
  ptable_KM_ACM_uncon$ID[n] <- i
  ptable_KM_ACM_uncon$pvalue[n] <- p.val
  ptable_KM_ACM_uncon$HR[n] <- HR
  ptable_KM_ACM_uncon$CIlow[n] <- low95
  ptable_KM_ACM_uncon$CIup[n] <- up95
  HR <- paste("Hazard Ratio = ", round(HR, 2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95, 2), round(up95, 2), sep = " - "), sep = "")
  pl[[i]] <- ggplot(rt, aes(time = ACM, status = ACM_Censor, color = get(i))) +
    geom_km(size = 2.2) +
    geom_kmticks(size = 1.8) +
    theme_bw(base_rect_size = 1.5) +
    scale_color_manual(
      values =
        scales::alpha(rev(c("#ffc857", "#93a8ac", "#61a5c2", "#FF6666")), .7)
    ) +
    guides(color = guide_legend(override.aes = list(size = 6))) +
    scale_y_continuous(limits = c(0.56, 1)) +
    labs(y = "ACM", x = "Time (months)") +
    ggtitle(i) +
    theme(
      legend.position = c(0.72, 0.72),
      legend.text = element_text(size = 25, colour = "black"),
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      axis.text = element_text(size = 19, colour = "black"),
      axis.title = element_text(size = 25, colour = "black", face = "bold"),
      plot.title = element_text(size = 28, hjust = 0.5, colour = "darkred", face = "bold"),
      panel.grid.major = element_line(color = "#cacfd2", size = .4, linetype = "dashed"),
      panel.grid.minor = element_line(color = "white", size = .3, linetype = "solid"),
      panel.background = element_rect(fill = "#F0F8FF"),
      panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid"),
    ) +
    annotate("text",
      x = 0.25, y = 0.63,
      label = paste(pval = ifelse(p.val < 0.0001, "Log-rank\np < 0.0001",
        paste("Log-rank\np = ", round(p.val, 5), sep = "")
      )),
      hjust = 0, size = 8, fontface = "bold.italic"
    ) +
    coord_cartesian(clip = "on")

  if (p.val < 0.05) {
    pl_sig[[i]] <- pl[[i]]
    png(filename = paste0(i, "_ACM", ".png"), width = 2000, height = 2000, res = 500)
    print(pl[[i]])
    dev.off()
  }
}
pl2_ACM_uncon <- pl2
save(pl2_ACM_uncon, file = "pl2_ACM_uncon.rda")
ptable_KM_ACM_uncon_sig <- ptable_KM_ACM_uncon %>% filter(pvalue < 0.05)
ptable_KM_ACM_uncon_sig_gene <- ptable_KM_ACM_uncon %>%
  filter(pvalue < 0.05) %>%
  pull(ID)
write.csv(ptable_KM_ACM_uncon, file = "ptable_KM_ACM_uncon.csv")
write.csv(ptable_KM_ACM_uncon_sig, file = "ptable_KM_ACM_uncon_sig.csv")



#######' @Continuous
desc_number <- stat.desc(rt)
list <- colnames(rt)[rowSums(is.na(t(desc_number))) == 0]
list
list <- list[-which(list %in% c("Follow_up_time", "ACM", "MACE"))]
list
rt <- rt[, c("MACE_Censor", "MACE", "ACM_Censor", "ACM", list)]
rt$MACE_Censor <- as.numeric(as.character(rt$MACE_Censor))
rt$ACM_Censor <- as.numeric(as.character(rt$ACM_Censor))

picDir <- "clinic_OS_KM_continuous_picture2"
dir.create(picDir)


setwd(picDir)
getwd()
res.cut <- surv_cutpoint(rt,
  time = "ACM",
  event = "ACM_Censor",
  variables = names(rt)[5:ncol(rt)],
  minprop = 0.3
)
res.cat <- surv_categorize(res.cut)
my.surv <- Surv(rt$ACM, rt$ACM_Censor)
pl2 <- pl2_sig <- list()
ptable_KM_ACM_con <- data.frame(matrix(NA, ncol(rt) - 4, 5))
colnames(ptable_KM_ACM_con) <- c("ID", "pvalue", "HR", "CIlow", "CIup")

n <- 0
for (i in colnames(rt)[5:ncol(rt)]) {
  group <- res.cat[, i]
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val <- 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR <- (data.survdiff$obs[2] / data.survdiff$exp[2]) / (data.survdiff$obs[1] / data.survdiff$exp[1])
  up95 <- exp(log(HR) + qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
  low95 <- exp(log(HR) - qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
  n <- n + 1
  ptable_KM_ACM_con$ID[n] <- i
  ptable_KM_ACM_con$pvalue[n] <- p.val
  ptable_KM_ACM_con$HR[n] <- HR
  ptable_KM_ACM_con$CIlow[n] <- low95
  ptable_KM_ACM_con$CIup[n] <- up95
  HR <- paste("Hazard Ratio = ", round(HR, 2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95, 2), round(up95, 2), sep = " - "), sep = "")

  pl2[[i]] <- ggplot(res.cat, aes(time = ACM, status = ACM_Censor, color = get(i))) +
    geom_km(size = 2.8) +
    geom_kmticks(size = 1.8) +
    theme_bw(base_rect_size = 1.5) +
    scale_color_manual(
      values =
        scales::alpha(rev(c("#ffc857", "#93a8ac", "#61a5c2", "#FF6666")), .7)
    ) +
    scale_y_continuous(limits = c(0.55, 1)) +
    labs(y = "ACM", x = "Time (months)") +
    ggtitle(i) +
    theme(
      legend.position = c(0.77, 0.81),
      legend.text = element_text(size = 28, colour = "black"),
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      axis.text = element_text(size = 19, colour = "black"),
      axis.title = element_text(size = 21, colour = "black", face = "bold"),
      plot.title = element_text(size = 33, hjust = 0.5, colour = "darkred", face = "bold"),
      panel.grid.major = element_line(color = "#cacfd2", size = .4, linetype = "dashed"),
      panel.grid.minor = element_line(color = "white", size = .3, linetype = "solid"),
      panel.background = element_rect(fill = "#F0F8FF"),
      panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid"),
    ) +
    annotate("text",
      x = 0.3, y = 0.67,
      label = paste(pval = ifelse(p.val < 0.0001, "Log-rank\np < 0.0001",
        paste("Log-rank\np = ", round(p.val, 5), sep = "")
      )),
      hjust = 0, size = 10.7, fontface = "italic"
    ) +
    coord_cartesian(clip = "on")

  if (p.val < 0.05) {
    pl2_sig[[i]] <- pl2[[i]]
    png(filename = paste0(i, "_ACM", ".png"), width = 1800, height = 2000, res = 500)
    print(pl2[[i]])
    dev.off()
  }
}
pl2_ACM_con <- pl2
save(pl2_ACM_con, file = "pl2_ACM_con.rda")
ptable_KM_ACM_con_sig <- ptable_KM_ACM_con %>% filter(pvalue < 0.05)
ptable_KM_ACM_con_sig_gene <- ptable_KM_ACM_con %>%
  filter(pvalue < 0.05) %>%
  pull(ID)
write.csv(ptable_KM_ACM_con, file = "ptable_KM_ACM_con.csv")
write.csv(ptable_KM_ACM_con_sig, file = "ptable_KM_ACM_con_sig.csv")





####################' @=====shared-features
load("ACM MACE cox result.RData")
a <- ls()
rm(list = a[which(a != "unicox_ACM_sig" & a != "ZG_class_fordf")])
cox_gene <- row.names(unicox_ACM_sig)
sort(cox_gene)

ptable_KM_ACM_uncon_sig <- read.csv("clinic_OS_KM_uncontinuous_picture\\ptable_KM_ACM_uncon_sig.csv", row.names = 1)
ptable_KM_ACM_con_sig <- read.csv("clinic_OS_KM_continuous_picture\\ptable_KM_ACM_con_sig.csv", row.names = 1)
km_gene <- c(ptable_KM_ACM_uncon_sig %>% pull(1), ptable_KM_ACM_con_sig %>% pull(1))
sort(km_gene)
keygene <- intersect(cox_gene, km_gene)
sort(keygene)
save(keygene, file = "ACM_keygene.rda")
load("ACM_keygene.rda")


rt <- rt[, c("MACE_Censor", "MACE", "ACM_Censor", "ACM", keygene)]
rt$HF_subtype <- factor(rt$HF_subtype, levels = c("HFpEF", "HFmrEF", "HFrEF"))
rt$MACE_Censor <- as.numeric(as.character(rt$MACE_Censor))
rt$ACM_Censor <- as.numeric(as.character(rt$ACM_Censor))
desc_number <- stat.desc(rt)
list <- colnames(rt)[rowSums(is.na(t(desc_number))) != 0]
list <- list[-which(list %in% c("Follow_up_time", "ACM", "MACE"))]
if (T) {
  rt <- rt %>% mutate(
    NYHA = case_when(
      NYHA %in% "I" ~ 1,
      NYHA %in% "II" ~ 2,
      NYHA %in% "III" ~ 3,
      NYHA %in% "IV" ~ 4,
    )
  )
  rt <- rt %>% mutate(
    CKD_stage = case_when(
      CKD_stage %in% "I" ~ 1,
      CKD_stage %in% "II" ~ 2,
      CKD_stage %in% "IIIa" | CKD_stage %in% "IIIb" ~ 3,
      CKD_stage %in% "IV" ~ 4,
      CKD_stage %in% "V" ~ 5
    )
  )
  rt <- rt %>% mutate(
    HF_subtype = case_when(
      HF_subtype %in% "HFpEF" ~ 1,
      HF_subtype %in% "HFmrEF" | HF_subtype %in% "HFrEF" ~ 2
    )
  )
}
rt$Gender <- ifelse(rt$Gender %in% c("female"), 0, 1)
desc_number <- stat.desc(rt)
list <- colnames(rt)[rowSums(is.na(t(desc_number))) > 1]
list
rt[, list] <- lapply(rt[, list], as.numeric)
expr.surv <- rt[, -c(1, 2)]
if (T) {
  desc_number <- stat.desc(expr.surv)
  list2 <- colnames(expr.surv)[rowSums(is.na(t(desc_number))) == 0]
  list2 <- list2[-which(list2 %in% c("MACE_Censor", "MACE", "ACM", "ACM_Censor"))]
  list2 <- list2[-which(list2 %in% c("Gender", "NYHA", "HTN", "HF_subtype", "CKD_stage"))]
  list2 <- list2[which(list2 %in% keygene)]
  list2 <- list2[-which(list2 %in% list)]
  expr.surv[, list2] <- lapply(expr.surv[, list2], scale)
}
save(expr.surv, file = "ACM_expr.surv.rda")




####################' @=====Modeling
load(file = "ACM_expr.surv.rda");expr.surv
load(file = "ACM_keygene.rda");keygene
if (T) {
  library(mlr3)
  library(mlr3viz)
  library(mlr3learners)
  library(mlr3verse)
  library(mlr3tuning)
  library(data.table)
  library(mlr3proba)
  library("magrittr")
  library(mlr3extralearners)
  mlr_tasks
  mlr_learners
  mlr_filters
  mlr_measures
  list = list_mlr3learners(select = c("id", "required_packages"))
}
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor")
print(task)
result = data.frame()




##############################################################################
#######################################' @0-1.L1,L2(lambda.min)(lambda.1se)
##############################################################################
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")

for (alpha in c(0, 1)) {
  for (s in c("lambda.min", "lambda.1se")) {
    set.seed(5000)
    lrn_cv_glmnet <- lrn("surv.cv_glmnet", alpha = alpha, nfolds = 10, s = s, type.measure = "C")
    measures <- c("surv.cindex")
    resampling_HD <- rsmp("holdout", ratio = 0.7)
    print(resampling_HD)
    set.seed(5000)
    rr_HD <- resample(task, lrn_cv_glmnet, resampling_HD, store_models = T)
    index_HD <- msrs(measures) %>% rr_HD$aggregate()

    resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
    print(resampling_SS)
    set.seed(5000)
    rr_SS <- resample(task, lrn_cv_glmnet, resampling_SS, store_models = T)
    index_SS <- msrs(measures) %>% rr_SS$aggregate()

    resampling_CV <- rsmp("repeated_cv", folds = 10, repeats = 2)
    print(resampling_CV)
    set.seed(5000)
    rr_CV <- resample(task, lrn_cv_glmnet, resampling_CV, store_models = T)
    index_CV <- msrs(measures) %>% rr_CV$aggregate()

    resampling_BS <- rsmp("bootstrap", repeats = 10)
    print(resampling_BS)
    set.seed(5000)
    rr_BS <- resample(task, lrn_cv_glmnet, resampling_BS, store_models = T)
    index_BS <- msrs(measures) %>% rr_BS$aggregate()

    name <- ifelse(alpha == 0, "L2", "L1")
    res <- data.frame(
      holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
    ) %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("ID")
    res$mod <- paste0(name, " regularization ", "(", s, ")")
    result <- rbind(result, res)
  }
}

##' @0-1.L1(lambda.min)
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
set.seed(1000)
rid <- lrn("surv.cv_glmnet", alpha = 1, nfolds = 10, s = "lambda.min", type.measure = "C")$train(task)$selected_features()
rid
save(rid, file = "0-1.L1(lambda.min)_atf_rid.rda")

##' @0-1.L1(lambda.1se)
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
set.seed(1000)
rid <- lrn("surv.cv_glmnet", alpha = 1, nfolds = 10, s = "lambda.1se", type.measure = "C")$train(task)$selected_features()
rid
save(rid, file = "0-1.L1(lambda.1se)_atf_rid.rda")


##############################################################################
#######################################' @1-1.Enet(lambda.1se)(lambda.min)
##############################################################################
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")

##############' @(alpha)(---lambda.min)
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet", alpha = 0.2, nfolds = 10, s = "lambda.min", type.measure = "C")
lrn_cv_glmnet$param_set
search_space <- ps(
  alpha = p_dbl(lower = 0.1, upper = 0.9)
)
search_space

resolution <- 50
n <- as.data.frame(rbindlist(generate_design_grid(search_space, resolution)$transpose())) %>% nrow()
n
tuner <- tnr("grid_search", resolution = resolution)
terminator <- trm("evals", n_evals = n)

set.seed(5000)
at <- AutoTuner$new(
  learner = lrn_cv_glmnet, resampling = rsmp("cv", folds = 10),
  measure = msr("surv.cindex"), search_space = search_space, terminator = terminator, tuner = tuner
)
at$train(task)

###############' @(alpha)(---lambda.1se)
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet", alpha = 0.2, nfolds = 10, s = "lambda.1se", type.measure = "C")
lrn_cv_glmnet$param_set
search_space <- ps(
  alpha = p_dbl(lower = 0.1, upper = 0.9)
)
search_space

resolution <- 50
n <- as.data.frame(rbindlist(generate_design_grid(search_space, resolution)$transpose())) %>% nrow()
n
tuner <- tnr("grid_search", resolution = resolution)
terminator <- trm("evals", n_evals = n)

set.seed(5000)
at <- AutoTuner$new(
  learner = lrn_cv_glmnet, resampling = rsmp("cv", folds = 10),
  measure = msr("surv.cindex"), search_space = search_space, terminator = terminator, tuner = tuner
)
at$train(task)


##############' @(lambda.min)
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet", alpha = 0.1, nfolds = 10, s = "lambda.min", type.measure = "C")
measures <- c("surv.cindex")
resampling_HD <- rsmp("holdout", ratio = 0.7)
print(resampling_HD)
set.seed(5000)
rr_HD <- resample(task, lrn_cv_glmnet, resampling_HD, store_models = T)
index_HD <- msrs(measures) %>% rr_HD$aggregate()

resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
print(resampling_SS)
set.seed(5000)
rr_SS <- resample(task, lrn_cv_glmnet, resampling_SS, store_models = T)
index_SS <- msrs(measures) %>% rr_SS$aggregate()

resampling_CV <- rsmp("repeated_cv", folds = 10, repeats = 2)
print(resampling_CV)
set.seed(5000)
rr_CV <- resample(task, lrn_cv_glmnet, resampling_CV, store_models = T)
index_CV <- msrs(measures) %>% rr_CV$aggregate()

resampling_BS <- rsmp("bootstrap", repeats = 10)
print(resampling_BS)
set.seed(5000)
rr_BS <- resample(task, lrn_cv_glmnet, resampling_BS, store_models = T)
index_BS <- msrs(measures) %>% rr_BS$aggregate()

res <- data.frame(
  holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("ID")
res$mod <- paste0("ENet regularization (lambda.min)")

result <- rbind(result, res)


##############' @(lambda.1se)
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet", alpha = 0.3285714, nfolds = 10, s = "lambda.1se", type.measure = "C")
measures <- c("surv.cindex")
resampling_HD <- rsmp("holdout", ratio = 0.7)
print(resampling_HD)
set.seed(5000)
rr_HD <- resample(task, lrn_cv_glmnet, resampling_HD, store_models = T)
index_HD <- msrs(measures) %>% rr_HD$aggregate()

resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
print(resampling_SS)
set.seed(5000)
rr_SS <- resample(task, lrn_cv_glmnet, resampling_SS, store_models = T)
index_SS <- msrs(measures) %>% rr_SS$aggregate()

resampling_CV <- rsmp("repeated_cv", folds = 10, repeats = 2)
print(resampling_CV)
set.seed(5000)
rr_CV <- resample(task, lrn_cv_glmnet, resampling_CV, store_models = T)
index_CV <- msrs(measures) %>% rr_CV$aggregate()

resampling_BS <- rsmp("bootstrap", repeats = 10)
print(resampling_BS)
set.seed(5000)
rr_BS <- resample(task, lrn_cv_glmnet, resampling_BS, store_models = T)
index_BS <- msrs(measures) %>% rr_BS$aggregate()

res <- data.frame(
  holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("ID")
res$mod <- paste0("ENet regularization (lambda.1se)")

result <- rbind(result, res)

##' @1-1.L1(lambda.min)
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
set.seed(1000)
rid <- lrn("surv.cv_glmnet", alpha = 0.1, nfolds = 10, s = "lambda.min", type.measure = "C")$train(task)$selected_features()
rid
save(rid, file = "1-1.ENet(lambda.min)_atf_rid.rda")

##' @1-1.L1(lambda.1se)
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
set.seed(1000)
rid <- lrn("surv.cv_glmnet", alpha = 0.3285714, nfolds = 10, s = "lambda.1se", type.measure = "C")$train(task)$selected_features()
rid
save(rid, file = "1-1.ENet(lambda.1se)_atf_rid.rda")


#################################################
############' @2-1.rsf
#################################################

set.seed(5000)
lrn_rfsrc <- lrn("surv.rfsrc", ntree = 2000, importance = "TRUE")
lrn_rfsrc
lrn_rfsrc$param_set # 展示超参数、选择范围等等

search_space <- ps(
  mtry = p_int(lower = 4, upper = 50),
  nodesize = p_int(lower = 1, upper = 30)
)
search_space # 设置超参数调优空间

resolution <- 14
n <- as.data.frame(rbindlist(generate_design_grid(search_space, resolution)$transpose())) %>% nrow()
n
tuner <- tnr("grid_search", resolution = resolution) # 选择超参数搜索算法
terminator <- trm("evals", n_evals = n) # 选择停止指标

set.seed(5000)
at <- AutoTuner$new(
  learner = lrn_rfsrc, resampling = rsmp("cv", folds = 4),
  measure = msr("surv.cindex"), search_space = search_space, terminator = terminator, tuner = tuner
)
at$train(task)
save(at, file = "2-1.rsf_at.rda")
load("2-1.rsf_at.rda")

at$archive$benchmark_result$aggregate(msr("surv.cindex"))
at$learner$param_set$values # 获取最佳超参数配置
lrn_rfsrc <- at$learner
lrn_rfsrc

#####' @(ntree=2000)-Wrapper—Methods_for_feature_selection
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
lrn_rfsrc <- lrn("surv.rfsrc", ntree = 2000, importance = "TRUE")
lrn_rfsrc$param_set

atf <- AutoFSelector$new(
  learner = lrn_rfsrc,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("none"),
  fselector = fs("sequential")
)
atf$train(task)

save(atf, file = "E2-1.rsf_atf.rda")
load("2-1.rsf_atf.rda")

rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task$select(rid)
task

task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
lrn_rfsrc <- lrn("surv.rfsrc", ntree = 2000, importance = "TRUE")
if (T) {
  measures <- c("surv.cindex")

  resampling_HD <- rsmp("holdout", ratio = 0.7)
  print(resampling_HD)
  set.seed(5000)
  rr_HD <- resample(task, lrn_rfsrc, resampling_HD, store_models = T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate()
  index_HD

  resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
  print(resampling_SS)
  set.seed(5000)
  rr_SS <- resample(task, lrn_rfsrc, resampling_SS, store_models = T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate()
  index_SS

  resampling_CV <- rsmp("repeated_cv", folds = 10, repeats = 2)
  print(resampling_CV)
  set.seed(5000)
  rr_CV <- resample(task, lrn_rfsrc, resampling_CV, store_models = T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate()
  index_CV

  resampling_BS <- rsmp("bootstrap", repeats = 10)
  print(resampling_BS)
  set.seed(5000)
  rr_BS <- resample(task, lrn_rfsrc, resampling_BS, store_models = T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate()
  index_BS

  res <- data.frame(
    holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  res$mod <- paste0("Surv.RSF")
}
result <- rbind(result, res)



#################################################
############' @3-1.surv.svm
#################################################

#####' @Wrapper—Methods_for_feature_selection
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
lrn_surv_svm <- lrn("surv.svm",
  type = "vanbelle1", kernel = "add_kernel", opt.meth = "ipop",
  diff.meth = "makediff1", gamma.mu = 0.1
)

set.seed(5000)
atf <- AutoFSelector$new(
  learner = lrn_surv_svm,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("none"),
  fselector = fs("sequential")
)

atf$train(task)
save(atf, file = "3-1.svm_atf.rda")

rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task$select(rid)
task

task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
lrn_surv_svm <- lrn("surv.svm",
  type = "vanbelle1", kernel = "add_kernel", opt.meth = "ipop",
  diff.meth = "makediff1", gamma.mu = 0.1
)
if (T) {
  measures <- c("surv.cindex")

  resampling_HD <- rsmp("holdout", ratio = 0.7)
  print(resampling_HD)
  set.seed(5000)
  rr_HD <- resample(task, lrn_surv_svm, resampling_HD, store_models = T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate()
  index_HD

  resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
  print(resampling_SS)
  set.seed(5000)
  rr_SS <- resample(task, lrn_surv_svm, resampling_SS, store_models = T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate()
  index_SS

  resampling_CV <- rsmp("repeated_cv", folds = 10, repeats = 2)
  print(resampling_CV)
  set.seed(5000)
  rr_CV <- resample(task, lrn_surv_svm, resampling_CV, store_models = T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate()
  index_CV

  resampling_BS <- rsmp("bootstrap", repeats = 10)
  print(resampling_BS)
  set.seed(5000)
  rr_BS <- resample(task, lrn_surv_svm, resampling_BS, store_models = T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate()
  index_BS

  res <- data.frame(
    holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  res$mod <- paste0("Surv.SVM")
}
result <- rbind(result, res)


#################################################
############' @4-1.surv.rpart
#################################################
set.seed(5000)
lrn_rpart_pre <- lrn("surv.rpart")
search_space <- ps(
  cp = p_dbl(lower = 0.001, upper = 0.1),
  minsplit = p_int(lower = 1, upper = 10),
  maxdepth = p_int(lower = 3, upper = 10),
  minbucket = p_int(lower = 3, upper = 10)
)
resolution <- 5
n <- as.data.frame(rbindlist(generate_design_grid(search_space, resolution)$transpose())) %>% nrow()
n
tuner <- tnr("grid_search", resolution = resolution)
terminator <- trm("evals", n_evals = n)

set.seed(5000)
at <- AutoTuner$new(
  learner = lrn_rpart_pre, resampling = rsmp("cv", folds = 3),
  measure = msr("surv.cindex"), search_space = search_space, terminator = terminator, tuner = tuner
)
at$train(task)
save(at, file = "4-1.rpart_at.rda")
load("4-1.rpart_at.rda")
lrn_rpart <- at$learner
lrn_rpart


task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
lrn_rpart <- lrn("surv.rpart")
lrn_rpart$param_set
set.seed(5000)
atf <- AutoFSelector$new(
  learner = lrn_rpart,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("none"),
  fselector = fs("sequential")
)
atf$train(task)

save(atf, file = "4-1.rpart_atf.rda")
load("4-1.rpart_atf.rda")

rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task$select(rid)
task


task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
lrn_rpart <- lrn("surv.rpart")
if (T) {
  measures <- c("surv.cindex")

  resampling_HD <- rsmp("holdout", ratio = 0.7)
  print(resampling_HD)
  set.seed(5000)
  rr_HD <- resample(task, lrn_rpart, resampling_HD, store_models = T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate()
  index_HD

  resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
  print(resampling_SS)
  set.seed(5000)
  rr_SS <- resample(task, lrn_rpart, resampling_SS, store_models = T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate()
  index_SS

  resampling_CV <- rsmp("repeated_cv", folds = 10, repeats = 2)
  print(resampling_CV)
  set.seed(5000)
  rr_CV <- resample(task, lrn_rpart, resampling_CV, store_models = T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate()
  index_CV

  resampling_BS <- rsmp("bootstrap", repeats = 10)
  print(resampling_BS)
  set.seed(5000)
  rr_BS <- resample(task, lrn_rpart, resampling_BS, store_models = T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate()
  index_BS

  res <- data.frame(
    holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  res$mod <- paste0("Surv.rpart")
}
result <- rbind(result, res)

##############################################################################
#######################################'@5-1.surv.xgboost
##############################################################################

set.seed(5000)
lrn_xgboost_pre <- lrn("surv.xgboost")
search_space <- ps(
  eta = p_int(lower = 0, upper = 1),
  gamma = p_int(lower = 0, upper = 5),
  max_depth = p_int(lower = 1, upper = 5),
  nrounds = p_int(lower = 500, upper = 800),
  subsample = p_dbl(lower = 0.5, upper = 1)
)
resolution <- 4
n <- as.data.frame(rbindlist(generate_design_grid(search_space, resolution)$transpose())) %>% nrow()
n
tuner <- tnr("grid_search", resolution = resolution) # 选择超参数搜索算法
terminator <- trm("evals", n_evals = n) # 选择停止指标

set.seed(5000)
at <- AutoTuner$new(
  learner = lrn_xgboost_pre, resampling = rsmp("cv", folds = 3),
  measure = msr("surv.cindex"), search_space = search_space, terminator = terminator, tuner = tuner
)
at$train(task)
save(at, file = "5-1.xgboost_at.rda")
load("5-1.xgboost_at.rda")
lrn_xgboost <- at$learner
lrn_xgboost


task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
lrn_xgboost <- lrn("surv.xgboost",
  max_depth = 3, eta = 0.01, gamma = 0, subsample = 0.8, colsample_bytree = 0.8,
  min_child_weight = 2, lambda = 1, alpha = 0, nrounds = 200
)
lrn_xgboost$param_set

set.seed(5000)
atf <- AutoFSelector$new(
  learner = lrn_xgboost,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("none"),
  fselector = fs("sequential")
)
atf$train(task)

rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task$select(rid)
task


task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
lrn_xgboost <- lrn("surv.xgboost",
  max_depth = 3, eta = 0.01, gamma = 0, subsample = 0.8, colsample_bytree = 0.8,
  min_child_weight = 2, lambda = 1, alpha = 0, nrounds = 200
)
if (T) {
  measures <- c("surv.cindex")

  resampling_HD <- rsmp("holdout", ratio = 0.7)
  print(resampling_HD)
  set.seed(5000)
  rr_HD <- resample(task, lrn_xgboost, resampling_HD, store_models = T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate()
  index_HD

  resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
  print(resampling_SS)
  set.seed(5000)
  rr_SS <- resample(task, lrn_xgboost, resampling_SS, store_models = T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate()
  index_SS

  resampling_CV <- rsmp("repeated_cv", folds = 10, repeats = 2)
  print(resampling_CV)
  set.seed(5000)
  rr_CV <- resample(task, lrn_xgboost, resampling_CV, store_models = T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate()
  index_CV

  resampling_BS <- rsmp("bootstrap", repeats = 10)
  print(resampling_BS)
  set.seed(5000)
  rr_BS <- resample(task, lrn_xgboost, resampling_BS, store_models = T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate()
  index_BS

  res <- data.frame(
    holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  res$mod <- paste0("Surv.Xgboost")
}
result <- rbind(result, res)


##############################################################################
#######################################' @6-1.surv.cv_coxboost
##############################################################################
library(CoxBoost)

task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
set.seed(5000)
pen <- suppressWarnings(suppressMessages(optimCoxBoostPenalty(expr.surv[, "ACM"], expr.surv[, "ACM_Censor"], as.matrix(expr.surv[, -c(1, 2)]),
  trace = F, start.penalty = 1000, parallel = T
)))
set.seed(5000)
cv.res <- cv.CoxBoost(expr.surv[, "ACM"], expr.surv[, "ACM_Censor"], as.matrix(expr.surv[, -c(1, 2)]),
  maxstepno = 1000, K = 10, type = "verweij", penalty = pen$penalty
)
lrn_coxboost <- lrn("surv.coxboost", stepno = cv.res$optimal.step)
lrn_coxboost$param_set

if (T) {
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout", ratio = 0.7)
  print(resampling_HD)
  set.seed(5000)
  rr_HD <- resample(task, lrn_coxboost, resampling_HD, store_models = T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate()
  index_HD

  resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
  print(resampling_SS)
  set.seed(5000)
  rr_SS <- resample(task, lrn_coxboost, resampling_SS, store_models = T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate()
  index_SS

  resampling_CV <- rsmp("repeated_cv", folds = 10, repeats = 2)
  print(resampling_CV)
  set.seed(5000)
  rr_CV <- resample(task, lrn_coxboost, resampling_CV, store_models = T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate()
  index_CV

  resampling_BS <- rsmp("bootstrap", repeats = 10)
  print(resampling_BS)
  set.seed(5000)
  rr_BS <- resample(task, lrn_coxboost, resampling_BS, store_models = T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate()
  index_BS

  res <- data.frame(
    holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  res$mod <- paste0("Surv.Coxboost")
}
result <- rbind(result, res)


set.seed(1000)
fit <- lrn_coxboost$train(task)$model
rid <- names(coef(fit)[coef(fit) != 0])
rid


##############################################################################
#######################################' @7-1.surv.glmboost
##############################################################################

task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
lrn_glmboost <- lrn("surv.glmboost", mstop = 500)

if (T) {
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout", ratio = 0.7)
  print(resampling_HD)
  set.seed(5000)
  rr_HD <- resample(task, lrn_glmboost, resampling_HD, store_models = T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate()
  index_HD

  resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
  print(resampling_SS)
  set.seed(5000)
  rr_SS <- resample(task, lrn_glmboost, resampling_SS, store_models = T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate()
  index_SS

  resampling_CV <- rsmp("repeated_cv", folds = 10, repeats = 2)
  print(resampling_CV)
  set.seed(5000)
  rr_CV <- resample(task, lrn_glmboost, resampling_CV, store_models = T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate()
  index_CV

  resampling_BS <- rsmp("bootstrap", repeats = 10)
  print(resampling_BS)
  set.seed(5000)
  rr_BS <- resample(task, lrn_glmboost, resampling_BS, store_models = T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate()
  index_BS

  res <- data.frame(
    holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  res$mod <- paste0("Surv.glmboost")
}
result <- rbind(result, res)


set.seed(5000)
fit <- lrn_glmboost$train(task)$model
rid <- names(coef(fit)[coef(fit) != 0])
rid


##############################################################################
#######################################' @8-1.surv.gbm
##############################################################################
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
lrn_gbm <- lrn("surv.gbm",
  distribution = "coxph", train.fraction = 0.7,
  n.trees = 1000
)
lrn_gbm$param_set

task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
lrn_gbm <- lrn("surv.gbm",
  distribution = "coxph", train.fraction = 0.7,
  n.trees = 1000
)
set.seed(5000)
atf <- AutoFSelector$new(
  learner = lrn_gbm,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("none"),
  fselector = fs("sequential")
)
atf$train(task)


rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task$select(rid)
task


task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
lrn_gbm <- lrn("surv.gbm", distribution = "coxph", train.fraction = 0.7, n.trees = 1000)
if (T) {
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout", ratio = 0.7)
  print(resampling_HD)
  set.seed(5000)
  rr_HD <- resample(task, lrn_gbm, resampling_HD, store_models = T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate()
  index_HD

  resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
  print(resampling_SS)
  set.seed(5000)
  rr_SS <- resample(task, lrn_gbm, resampling_SS, store_models = T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate()
  index_SS

  resampling_CV <- rsmp("repeated_cv", folds = 10, repeats = 2)
  print(resampling_CV)
  set.seed(5000)
  rr_CV <- resample(task, lrn_gbm, resampling_CV, store_models = T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate()
  index_CV

  resampling_BS <- rsmp("bootstrap", repeats = 10)
  print(resampling_BS)
  set.seed(5000)
  rr_BS <- resample(task, lrn_gbm, resampling_BS, store_models = T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate()
  index_BS

  res <- data.frame(
    holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  res$mod <- paste0("Surv.gbm")
}
result <- rbind(result, res)


##############################################################################
#######################################' @9-1.surv.parametric
##############################################################################
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
lrn_parametric <- lrn("surv.parametric")
lrn_parametric$param_set

task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
lrn_parametric <- lrn("surv.parametric")

set.seed(5000)
atf <- AutoFSelector$new(
  learner = lrn_parametric,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("none"),
  fselector = fs("sequential")
)
atf$train(task)


rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task$select(rid)
task


task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
lrn_parametric <- lrn("surv.parametric")
if (T) {
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout", ratio = 0.7)
  print(resampling_HD)
  set.seed(5000)
  rr_HD <- resample(task, lrn_parametric, resampling_HD, store_models = T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate()
  index_HD

  resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
  print(resampling_SS)
  set.seed(5000)
  rr_SS <- resample(task, lrn_parametric, resampling_SS, store_models = T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate()
  index_SS

  resampling_CV <- rsmp("repeated_cv", folds = 10, repeats = 2)
  print(resampling_CV)
  set.seed(5000)
  rr_CV <- resample(task, lrn_parametric, resampling_CV, store_models = T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate()
  index_CV

  resampling_BS <- rsmp("bootstrap", repeats = 10)
  print(resampling_BS)
  set.seed(5000)
  rr_BS <- resample(task, lrn_parametric, resampling_BS, store_models = T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate()
  index_BS

  res <- data.frame(
    holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  res$mod <- paste0("Surv.Fully-Parametric")
}
result <- rbind(result, res)




##############################################################################
#######################################' @10-1.surv.aorsf
##############################################################################
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
lrn_aorsf <- lrn("surv.aorsf", n_tree = 1000)
lrn_aorsf$param_set

task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
lrn_aorsf <- lrn("surv.aorsf", n_tree = 1000)

set.seed(5000)
atf <- AutoFSelector$new(
  learner = lrn_aorsf,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("none"),
  fselector = fs("sequential", min_features = 2)
)
atf$train(task)

rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task$select(rid)
task


task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
lrn_aorsf <- lrn("surv.aorsf", n_tree = 1000)
if (T) {
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout", ratio = 0.7)
  print(resampling_HD)
  set.seed(5000)
  rr_HD <- resample(task, lrn_aorsf, resampling_HD, store_models = T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate()
  index_HD

  resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
  print(resampling_SS)
  set.seed(5000)
  rr_SS <- resample(task, lrn_aorsf, resampling_SS, store_models = T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate()
  index_SS

  resampling_CV <- rsmp("repeated_cv", folds = 10, repeats = 2)
  print(resampling_CV)
  set.seed(5000)
  rr_CV <- resample(task, lrn_aorsf, resampling_CV, store_models = T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate()
  index_CV

  resampling_BS <- rsmp("bootstrap", repeats = 10)
  print(resampling_BS)
  set.seed(5000)
  rr_BS <- resample(task, lrn_aorsf, resampling_BS, store_models = T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate()
  index_BS

  res <- data.frame(
    holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  res$mod <- paste0("Surv.aorsf")
}
result <- rbind(result, res)


##############################################################################
#######################################'@0-2.L1(lambda.1se)&Surv.RSF
##############################################################################
load(file = '0-1.L1(lambda.1se)_atf_rid.rda')
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
lrn_rfsrc <- lrn("surv.rfsrc", ntree = 2000,importance = 'TRUE')
load('0-2.L1(lambda.1se)&Surv.RSF_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task
if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_rfsrc,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_rfsrc,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_rfsrc,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_rfsrc,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('L1 regularization & Surv.RSF')
}
result <- rbind(result,res)



##############################################################################
#######################################' @0-3.L1(lambda.1se)&Surv.SVM
##############################################################################
load(file = "0-1.L1(lambda.1se)_atf_rid.rda")
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
task$select(rid)
task
lrn_surv_svm <- lrn("surv.svm",
  type = "vanbelle1", kernel = "add_kernel", opt.meth = "ipop",
  diff.meth = "makediff1", gamma.mu = 0.1
)
load("0-3.L1(lambda.1se)&Surv.SVM_atf.rda")
rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task$select(rid)
task

if (T) {
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout", ratio = 0.7)
  print(resampling_HD)
  set.seed(5000)
  rr_HD <- resample(task, lrn_surv_svm, resampling_HD, store_models = T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate()
  index_HD

  resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
  print(resampling_SS)
  set.seed(5000)
  rr_SS <- resample(task, lrn_surv_svm, resampling_SS, store_models = T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate()
  index_SS

  resampling_CV <- rsmp("repeated_cv", folds = 5, repeats = 1)
  print(resampling_CV)
  set.seed(5000)
  rr_CV <- resample(task, lrn_surv_svm, resampling_CV, store_models = T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate()
  index_CV

  resampling_BS <- rsmp("bootstrap", repeats = 10)
  print(resampling_BS)
  set.seed(5000)
  rr_BS <- resample(task, lrn_surv_svm, resampling_BS, store_models = T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate()
  index_BS

  res <- data.frame(
    holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  res$mod <- paste0("L1 regularization & Surv.SVM")
}
result <- rbind(result, res)



##############################################################################
#######################################' @0-4.L1(lambda.1se)&surv.rpart
##############################################################################
load(file = "0-1.L1(lambda.1se)_atf_rid.rda")
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
task$select(rid)
task
lrn_rpart <- lrn("surv.rpart")
load("0-4.L1(lambda.1se)&Surv.rpart_atf.rda")
rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task$select(rid)
task

if (T) {
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout", ratio = 0.7)
  print(resampling_HD)
  set.seed(5000)
  rr_HD <- resample(task, lrn_rpart, resampling_HD, store_models = T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate()
  index_HD

  resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
  print(resampling_SS)
  set.seed(5000)
  rr_SS <- resample(task, lrn_rpart, resampling_SS, store_models = T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate()
  index_SS

  resampling_CV <- rsmp("repeated_cv", folds = 5, repeats = 1)
  print(resampling_CV)
  set.seed(5000)
  rr_CV <- resample(task, lrn_rpart, resampling_CV, store_models = T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate()
  index_CV

  resampling_BS <- rsmp("bootstrap", repeats = 10, ratio = 0.8)
  print(resampling_BS)
  set.seed(5000)
  rr_BS <- resample(task, lrn_rpart, resampling_BS, store_models = T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate()
  index_BS

  res <- data.frame(
    holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  res$mod <- paste0("L1 regularization & Surv.rpart")
}
result <- rbind(result, res)




##############################################################################
#######################################' @0-5.L1(lambda.1se)&surv.xgboost
##############################################################################
load(file = "0-1.L1(lambda.1se)_atf_rid.rda")
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
task$select(rid)
task
lrn_xgboost <- lrn("surv.xgboost",
  max_depth = 3, eta = 0.01, gamma = 0, subsample = 0.8, colsample_bytree = 0.8,
  min_child_weight = 2, lambda = 1, alpha = 0, nrounds = 200
)

load("0-5.L1(lambda.1se)&Surv.xgboost_atf.rda")
rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task$select(rid)
task

if (T) {
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout", ratio = 0.7)
  print(resampling_HD)
  set.seed(5000)
  rr_HD <- resample(task, lrn_xgboost, resampling_HD, store_models = T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate()
  index_HD

  resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
  print(resampling_SS)
  set.seed(5000)
  rr_SS <- resample(task, lrn_xgboost, resampling_SS, store_models = T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate()
  index_SS

  resampling_CV <- rsmp("repeated_cv", folds = 5, repeats = 1)
  print(resampling_CV)
  set.seed(5000)
  rr_CV <- resample(task, lrn_xgboost, resampling_CV, store_models = T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate()
  index_CV

  resampling_BS <- rsmp("bootstrap", repeats = 10, ratio = 0.8)
  print(resampling_BS)
  set.seed(5000)
  rr_BS <- resample(task, lrn_xgboost, resampling_BS, store_models = T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate()
  index_BS

  res <- data.frame(
    holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  res$mod <- paste0("L1 regularization & Surv.xgboost")
}
result <- rbind(result, res)



##############################################################################
#######################################' @0-6.L1(lambda.1se)&surv.coxboost
##############################################################################
load(file = "0-1.L1(lambda.1se)_atf_rid.rda")
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
task$select(rid)
task
set.seed(5000)
library(CoxBoost)
pen <- suppressWarnings(suppressMessages(optimCoxBoostPenalty(expr.surv[, "ACM"], expr.surv[, "ACM_Censor"], as.matrix(expr.surv[, -c(1, 2)]),
  trace = F, start.penalty = 1000, parallel = T
)))
set.seed(5000)
cv.res <- cv.CoxBoost(expr.surv[, "ACM"], expr.surv[, "ACM_Censor"], as.matrix(expr.surv[, -c(1, 2)]),
  maxstepno = 1000, K = 10, type = "verweij", penalty = pen$penalty
)
lrn_coxboost <- lrn("surv.coxboost", stepno = cv.res$optimal.step)
load("0-6.L1(lambda.1se)&Surv.coxboost_atf.rda")
rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task$select(rid)
task

if (T) {
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout", ratio = 0.7)
  print(resampling_HD)
  set.seed(5000)
  rr_HD <- resample(task, lrn_coxboost, resampling_HD, store_models = T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate()
  index_HD

  resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
  print(resampling_SS)
  set.seed(5000)
  rr_SS <- resample(task, lrn_coxboost, resampling_SS, store_models = T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate()
  index_SS

  resampling_CV <- rsmp("repeated_cv", folds = 5, repeats = 1)
  print(resampling_CV)
  set.seed(5000)
  rr_CV <- resample(task, lrn_coxboost, resampling_CV, store_models = T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate()
  index_CV

  resampling_BS <- rsmp("bootstrap", repeats = 10, ratio = 0.8)
  print(resampling_BS)
  set.seed(5000)
  rr_BS <- resample(task, lrn_coxboost, resampling_BS, store_models = T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate()
  index_BS

  res <- data.frame(
    holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  res$mod <- paste0("L1 regularization & Surv.coxboost")
}
result <- rbind(result, res)



##############################################################################
#######################################' @0-7.L1(lambda.1se)&surv.glmboost
##############################################################################
load(file = "0-1.L1(lambda.1se)_atf_rid.rda")
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
task$select(rid)
task
set.seed(5000)
lrn_glmboost <- lrn("surv.glmboost", mstop = 500)
load("0-7.L1(lambda.1se)&Surv.glmboost_atf.rda")
rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task$select(rid)
task

if (T) {
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout", ratio = 0.7)
  print(resampling_HD)
  set.seed(5000)
  rr_HD <- resample(task, lrn_glmboost, resampling_HD, store_models = T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate()
  index_HD

  resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
  print(resampling_SS)
  set.seed(5000)
  rr_SS <- resample(task, lrn_glmboost, resampling_SS, store_models = T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate()
  index_SS

  resampling_CV <- rsmp("repeated_cv", folds = 5, repeats = 1)
  print(resampling_CV)
  set.seed(5000)
  rr_CV <- resample(task, lrn_glmboost, resampling_CV, store_models = T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate()
  index_CV

  resampling_BS <- rsmp("bootstrap", repeats = 10, ratio = 0.8)
  print(resampling_BS)
  set.seed(5000)
  rr_BS <- resample(task, lrn_glmboost, resampling_BS, store_models = T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate()
  index_BS

  res <- data.frame(
    holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  res$mod <- paste0("L1 regularization & Surv.glmboost")
}
result <- rbind(result, res)





##############################################################################
#######################################' @0-8.L1(lambda.1se)&surv.gbm
##############################################################################
load(file = "0-1.L1(lambda.1se)_atf_rid.rda")
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
task$select(rid)
task
set.seed(5000)
lrn_gbm <- lrn("surv.gbm",
  distribution = "coxph", train.fraction = 0.7,
  n.trees = 1000
)
load("0-8.L1(lambda.1se)&Surv.gbm_atf.rda")
rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task$select(rid)
task

if (T) {
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout", ratio = 0.7)
  print(resampling_HD)
  set.seed(5000)
  rr_HD <- resample(task, lrn_gbm, resampling_HD, store_models = T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate()
  index_HD

  resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
  print(resampling_SS)
  set.seed(5000)
  rr_SS <- resample(task, lrn_gbm, resampling_SS, store_models = T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate()
  index_SS

  resampling_CV <- rsmp("repeated_cv", folds = 5, repeats = 1)
  print(resampling_CV)
  set.seed(5000)
  rr_CV <- resample(task, lrn_gbm, resampling_CV, store_models = T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate()
  index_CV

  resampling_BS <- rsmp("bootstrap", repeats = 10, ratio = 0.8)
  print(resampling_BS)
  set.seed(5000)
  rr_BS <- resample(task, lrn_gbm, resampling_BS, store_models = T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate()
  index_BS

  res <- data.frame(
    holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  res$mod <- paste0("L1 regularization & Surv.gbm")
}
result <- rbind(result, res)



##############################################################################
#######################################' @0-9.L1(lambda.1se)&surv.parametric
##############################################################################
load(file = "0-1.L1(lambda.1se)_atf_rid.rda")
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
task$select(rid)
task
set.seed(5000)
lrn_parametric <- lrn("surv.parametric")

load("0-9.L1(lambda.1se)&Surv.parametric_atf.rda")
rid
task$select(rid)
task

if (T) {
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout", ratio = 0.7)
  print(resampling_HD)
  set.seed(5000)
  rr_HD <- resample(task, lrn_parametric, resampling_HD, store_models = T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate()
  index_HD

  resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
  print(resampling_SS)
  set.seed(5000)
  rr_SS <- resample(task, lrn_parametric, resampling_SS, store_models = T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate()
  index_SS

  resampling_CV <- rsmp("repeated_cv", folds = 5, repeats = 1)
  print(resampling_CV)
  set.seed(5000)
  rr_CV <- resample(task, lrn_parametric, resampling_CV, store_models = T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate()
  index_CV

  resampling_BS <- rsmp("bootstrap", repeats = 10, ratio = 0.8)
  print(resampling_BS)
  set.seed(5000)
  rr_BS <- resample(task, lrn_parametric, resampling_BS, store_models = T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate()
  index_BS

  res <- data.frame(
    holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  res$mod <- paste0("L1 regularization & Surv.Fully-Parametric")
}
result <- rbind(result, res)


##############################################################################
#######################################' @0-10.L1(lambda.1se)&surv.aorsf
##############################################################################
load(file = "0-1.L1(lambda.1se)_atf_rid.rda")
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
task$select(rid)
task
set.seed(5000)
lrn_aorsf <- lrn("surv.aorsf", n_tree = 1000)

load("0-10.L1(lambda.1se)&Surv.aorsf_atf.rda")
rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task$select(rid)
task


if (T) {
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout", ratio = 0.7)
  print(resampling_HD)
  set.seed(5000)
  rr_HD <- resample(task, lrn_aorsf, resampling_HD, store_models = T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate()
  index_HD

  resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
  print(resampling_SS)
  set.seed(5000)
  rr_SS <- resample(task, lrn_aorsf, resampling_SS, store_models = T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate()
  index_SS

  resampling_CV <- rsmp("repeated_cv", folds = 5, repeats = 1)
  print(resampling_CV)
  set.seed(5000)
  rr_CV <- resample(task, lrn_aorsf, resampling_CV, store_models = T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate()
  index_CV

  resampling_BS <- rsmp("bootstrap", repeats = 10, ratio = 0.8)
  print(resampling_BS)
  set.seed(5000)
  rr_BS <- resample(task, lrn_aorsf, resampling_BS, store_models = T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate()
  index_BS

  res <- data.frame(
    holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  res$mod <- paste0("L1 regularization & Surv.aorsf")
}
result <- rbind(result, res)



#############################################################################
#######################################' @1-2.ENet(lambda.min)&Surv.RSF
##############################################################################
load(file = "1-1.ENet(lambda.1se)_atf_rid.rda")
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
task$select(rid)
task
lrn_rfsrc <- lrn("surv.rfsrc", ntree = 2000, importance = "TRUE")

load("1-2.ENet(lambda.1se)&Surv.RSF_atf.rda")
rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task$select(rid)
task


if (T) {
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout", ratio = 0.7)
  print(resampling_HD)
  set.seed(1000)
  rr_HD <- resample(task, lrn_rfsrc, resampling_HD, store_models = T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate()
  index_HD

  resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
  print(resampling_SS)
  set.seed(1000)
  rr_SS <- resample(task, lrn_rfsrc, resampling_SS, store_models = T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate()
  index_SS

  resampling_CV <- rsmp("repeated_cv", folds = 5, repeats = 1)
  print(resampling_CV)
  set.seed(5000)
  rr_CV <- resample(task, lrn_rfsrc, resampling_CV, store_models = T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate()
  index_CV

  resampling_BS <- rsmp("bootstrap", repeats = 10)
  print(resampling_BS)
  set.seed(5000)
  rr_BS <- resample(task, lrn_rfsrc, resampling_BS, store_models = T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate()
  index_BS

  res <- data.frame(
    holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  res$mod <- paste0("ENet regularization & Surv.RSF")
}
result <- rbind(result, res)



##############################################################################
#######################################'@1-3.ENet(lambda.1se)&Surv.SVM
##############################################################################
load(file = '1-1.ENet(lambda.1se)_atf_rid.rda')
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
lrn_surv_svm <- lrn("surv.svm",type="vanbelle1",kernel="add_kernel",opt.meth="ipop",
                    diff.meth="makediff1",gamma.mu = 0.1
)

load('1-3.ENet(lambda.1se)&Surv.SVM_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(5000)
  rr_HD <- resample(task,lrn_surv_svm,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(5000)
  rr_SS <- resample(task,lrn_surv_svm,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_surv_svm,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_surv_svm,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('ENet regularization & Surv.SVM')
}
result <- rbind(result,res)



##############################################################################
#######################################'@1-4.ENet(lambda.1se)&surv.rpart
##############################################################################
load(file = '1-1.ENet(lambda.1se)_atf_rid.rda')
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
lrn_rpart <- lrn("surv.rpart")

load('1-4.ENet(lambda.1se)&Surv.rpart_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(5000)
  rr_HD <- resample(task,lrn_rpart,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(5000)
  rr_SS <- resample(task,lrn_rpart,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_rpart,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10,ratio = 0.8);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_rpart,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('ENet regularization & Surv.rpart')
}
result <- rbind(result,res)





##############################################################################
#######################################'@1-5.ENet(lambda.1se)&surv.xgboost
##############################################################################
load(file = '1-1.ENet(lambda.1se)_atf_rid.rda')
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
lrn_xgboost <- lrn("surv.xgboost",
                   max_depth = 3,eta = 0.01,gamma = 0,subsample = 0.8,colsample_bytree = 0.8,
                   min_child_weight = 2,lambda=1,alpha=0,nrounds = 200
)
load('1-5.ENet(lambda.1se)&Surv.xgboost_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(5000)
  rr_HD <- resample(task,lrn_xgboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(5000)
  rr_SS <- resample(task,lrn_xgboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_xgboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10,ratio = 0.8);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_xgboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('ENet regularization & Surv.xgboost')
}
result <- rbind(result,res)



##############################################################################
#######################################'@1-6.ENet(lambda.1se)&surv.coxboost
##############################################################################
load(file = '1-1.ENet(lambda.1se)_atf_rid.rda')
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
library(CoxBoost)
pen <- suppressWarnings(suppressMessages(optimCoxBoostPenalty(expr.surv[,'ACM'],expr.surv[,'ACM_Censor'],as.matrix(expr.surv[,-c(1,2)]),
                                                              trace=F,start.penalty=1000,parallel = T)));set.seed(5000)
cv.res <- cv.CoxBoost(expr.surv[,'ACM'],expr.surv[,'ACM_Censor'],as.matrix(expr.surv[,-c(1,2)]),
                      maxstepno=1000,K=10,type="verweij",penalty=pen$penalty)
lrn_coxboost <- lrn("surv.coxboost",stepno=cv.res$optimal.step
)

load('1-6.ENet(lambda.1se)&Surv.coxboost_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(5000)
  rr_HD <- resample(task,lrn_coxboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(5000)
  rr_SS <- resample(task,lrn_coxboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_coxboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10,ratio = 0.8);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_coxboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('ENet regularization & Surv.coxboost')
}
result <- rbind(result,res)



##############################################################################
#######################################'@1-7.ENet(lambda.1se)&surv.glmboost
##############################################################################
load(file = '1-1.ENet(lambda.1se)_atf_rid.rda')
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_glmboost <- lrn("surv.glmboost",mstop = 500
)

load('1-7.ENet(lambda.1se)&Surv.glmboost_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(5000)
  rr_HD <- resample(task,lrn_glmboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(5000)
  rr_SS <- resample(task,lrn_glmboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_glmboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10,ratio = 0.8);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_glmboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('ENet regularization & Surv.glmboost')
}
result <- rbind(result,res)



##############################################################################
#######################################'@1-8.ENet(lambda.1se)&surv.gbm
##############################################################################
load(file = 'D:\\KEYAN2\\临床模型\\心衰\\HF\\4 NB model\\1-1.ENet(lambda.1se)_atf_rid.rda')
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_gbm <- lrn("surv.gbm",distribution = 'coxph',train.fraction = 0.7,
               n.trees = 1000
)

load('1-8.ENet(lambda.1se)&Surv.gbm_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(5000)
  rr_HD <- resample(task,lrn_gbm,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(5000)
  rr_SS <- resample(task,lrn_gbm,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_gbm,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10,ratio = 0.8);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_gbm,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('ENet regularization & Surv.gbm')
}
result <- rbind(result,res)



##############################################################################
#######################################'@1-9.ENet(lambda.1se)&surv.parametric
##############################################################################
load(file = '1-1.ENet(lambda.1se)_atf_rid.rda')
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_parametric <- lrn("surv.parametric"
)

load('1-9.ENet(lambda.1se)&Surv.parametric_atf.rda');rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(5000)
  rr_HD <- resample(task,lrn_parametric,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(5000)
  rr_SS <- resample(task,lrn_parametric,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_parametric,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10,ratio = 0.8);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_parametric,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('ENet regularization & Surv.Fully-Parametric')
}
result <- rbind(result,res)




##############################################################################
#######################################'@1-10.ENet(lambda.1se)&surv.aorsf
##############################################################################
load(file = '1-1.ENet(lambda.1se)_atf_rid.rda')
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_aorsf <- lrn("surv.aorsf",n_tree=1000
)

load('1-10.ENet(lambda.1se)&Surv.aorsf_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(5000)
  rr_HD <- resample(task,lrn_aorsf,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(5000)
  rr_SS <- resample(task,lrn_aorsf,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_aorsf,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10,ratio = 0.8);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_aorsf,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('ENet regularization & Surv.aorsf')
}
result <- rbind(result,res)



##############################################################################
#######################################'@2-01.Surv.RSF&L2(lambda.min)
##############################################################################
load(file = '2-1.rsf_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 0,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('2-01.Surv.RSF&L2(lambda.min)_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.RSF & L2 regularization')
}
result <- rbind(result,res)


##############################################################################
#######################################'@2-02.Surv.RSF&L1(lambda.min)
##############################################################################
load(file = '2-1.rsf_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid ###16个特征！！！！
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 1,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('2-02.Surv.RSF&L1(lambda.min)_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.RSF & L1 regularization')
}
result <- rbind(result,res)




##############################################################################
#######################################'@2-1.Surv.RSF&Enet(lambda.min)
##############################################################################
load(file = '2-1.rsf_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 0.1,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('2-1.Surv.RSF&Enet(lambda.min)_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.RSF & Enet regularization')
}
result <- rbind(result,res)




##############################################################################
#######################################'@2-3.Surv.RSF&Surv.SVM
##############################################################################
load(file = '2-1.rsf_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_surv_svm <- lrn("surv.svm",type="vanbelle1",kernel="add_kernel",opt.meth="ipop",
                    diff.meth="makediff1",gamma.mu = 0.1
)

load('2-3.Surv.RSF&Surv.SVM_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_surv_svm,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_surv_svm,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_surv_svm,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_surv_svm,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.RSF & Surv.SVM')
}
result <- rbind(result,res)


##############################################################################
#######################################'@2-4.Surv.RSF&Surv.rpart
##############################################################################
load(file = '2-1.rsf_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_rpart <- lrn("surv.rpart")

load('2-4.Surv.RSF&Surv.rpart_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_rpart,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_rpart,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_rpart,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_rpart,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.RSF & Surv.rpart')
}
result <- rbind(result,res)




##############################################################################
#######################################'@2-5.Surv.RSF&Surv.xgboost
##############################################################################
load(file = '2-1.rsf_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_xgboost <- lrn("surv.xgboost",
                   max_depth = 3,eta = 0.01,gamma = 0,subsample = 0.8,colsample_bytree = 0.8,
                   min_child_weight = 2,lambda=1,alpha=0,nrounds = 200
)

load('2-5.Surv.RSF&Surv.xgboost_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_xgboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_xgboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_xgboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_xgboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.RSF & Surv.xgboost')
}
result <- rbind(result,res)




##############################################################################
#######################################'@2-6.Surv.RSF&Surv.coxboost
##############################################################################
library(CoxBoost)
load(file = '2-1.rsf_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task

load('2-6.Surv.RSF&Surv.coxboost_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task#新特征的task!!

set.seed(5000)
pen <- suppressWarnings(suppressMessages(optimCoxBoostPenalty(expr.surv[,'ACM'],expr.surv[,'ACM_Censor'],as.matrix(expr.surv[,-c(1,2)]),
                                                              trace=F,start.penalty=1000,parallel = T)));set.seed(5000)
cv.res <- cv.CoxBoost(expr.surv[,'ACM'],expr.surv[,'ACM_Censor'],as.matrix(expr.surv[,-c(1,2)]),
                      maxstepno=1000,K=10,type="verweij",penalty=pen$penalty)
lrn_coxboost <- lrn("surv.coxboost",stepno=cv.res$optimal.step
)

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_coxboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_coxboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_coxboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_coxboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.RSF & Surv.coxboost')
}
result <- rbind(result,res)


##############################################################################
#######################################'@2-7.Surv.RSF&Surv.glmboost
##############################################################################
load(file = '2-1.rsf_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid ###16个特征！！！！
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_glmboost <- lrn("surv.glmboost",mstop = 500
)

load('2-7.Surv.RSF&Surv.glmboost_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_glmboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_glmboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_glmboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_glmboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.RSF & Surv.glmboost')
}
result <- rbind(result,res)



##############################################################################
#######################################'@2-8.Surv.RSF&Surv.gbm
##############################################################################
load(file = '2-1.rsf_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_gbm <- lrn("surv.gbm",distribution = 'coxph',train.fraction = 0.7,
               n.trees = 1000
)
load('2-8.Surv.RSF&Surv.gbm_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_gbm,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_gbm,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_gbm,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_gbm,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.RSF & Surv.gbm')
}
result <- rbind(result,res)



##############################################################################
#######################################'@2-9.Surv.RSF&Surv.parametric-----15gene
##############################################################################
load(file = '2-1.rsf_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid ###16个特征！！！！
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_parametric <- lrn("surv.parametric"
)

load('2-9.Surv.RSF&Surv.parametric_atf.rda');rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_parametric,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_parametric,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_parametric,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_parametric,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.RSF & Surv.Fully-Parametric')
}
result <- rbind(result,res)


##############################################################################
#######################################'@2-10.Surv.RSF&Surv.aorsf
##############################################################################
load(file = '2-1.rsf_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task

set.seed(5000)
lrn_aorsf <- lrn("surv.aorsf",n_tree=1000
)

load('2-10.Surv.RSF&Surv.aorsf_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_aorsf,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_aorsf,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_aorsf,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_aorsf,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.RSF & Surv.aorsf')
}
result <- rbind(result,res)




##############################################################################
#######################################'@3-01.Surv.SVM&L2(lambda.min)
##############################################################################
load(file = '3-1.svm_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 0,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('3-01.Surv.SVM&L2(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.SVM & L2 regularization')
}
result <- rbind(result,res)


##############################################################################
#######################################'@3-02.Surv.SVM&L1(lambda.min)
##############################################################################
load(file = '3-1.svm_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 1,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('3-02.Surv.SVM&L1(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.SVM & L1 regularization')
}
result <- rbind(result,res)



##############################################################################
#######################################'@3-1.Surv.SVM&Enet(lambda.min)
##############################################################################
load(file = '3-1.svm_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 0.1,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('3-1.Surv.SVM&Enet(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.SVM & Enet regularization')
}
result <- rbind(result,res)


##############################################################################
#######################################'@3-2.Surv.SVM&Surv.RSF
##############################################################################
load(file = '3-1.svm_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_rfsrc <- lrn("surv.rfsrc", ntree = 2000,importance = 'TRUE')

load('3-2.Surv.SVM&Surv.RSF_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_rfsrc,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_rfsrc,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_rfsrc,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_rfsrc,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.SVM & Surv.RSF')
}
result <- rbind(result,res)


##############################################################################
#######################################'@3-4.Surv.SVM&Surv.rpart
##############################################################################
load(file = '3-1.svm_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_rpart <- lrn("surv.rpart")

load('3-4.Surv.SVM&Surv.rpart_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_rpart,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_rpart,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_rpart,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_rpart,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.SVM & Surv.rpart')
}
result <- rbind(result,res)


##############################################################################
#######################################'@3-5.Surv.SVM&Surv.xgboost
##############################################################################
load(file = '3-1.svm_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_xgboost <- lrn("surv.xgboost",
                   max_depth = 3,eta = 0.01,gamma = 0,subsample = 0.8,colsample_bytree = 0.8,
                   min_child_weight = 2,lambda=1,alpha=0,nrounds = 200
)

load('3-5.Surv.SVM&Surv.xgboost_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task



if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_xgboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_xgboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_xgboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_xgboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.SVM & Surv.xgboost')
}
result <- rbind(result,res)



##############################################################################
#######################################'@3-6.Surv.SVM&Surv.coxboost
##############################################################################
load(file = '3-1.svm_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)

load('3-6.Surv.SVM&Surv.coxboost_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

set.seed(5000)
pen <- suppressWarnings(suppressMessages(optimCoxBoostPenalty(expr.surv[,'ACM'],expr.surv[,'ACM_Censor'],as.matrix(expr.surv[,-c(1,2)]),
                                                              trace=F,start.penalty=1000,parallel = T)));set.seed(5000)
cv.res <- cv.CoxBoost(expr.surv[,'ACM'],expr.surv[,'ACM_Censor'],as.matrix(expr.surv[,-c(1,2)]),
                      maxstepno=1000,K=10,type="verweij",penalty=pen$penalty)
lrn_coxboost <- lrn("surv.coxboost",stepno=cv.res$optimal.step
)

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_coxboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_coxboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_coxboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_coxboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.SVM & Surv.coxboost')
}
result <- rbind(result,res)



##############################################################################
#######################################'@3-7.Surv.SVM&Surv.glmboost
##############################################################################
load(file = '3-1.svm_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_glmboost <- lrn("surv.glmboost",mstop = 500
)

load('3-7.Surv.SVM&Surv.glmboost_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_glmboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_glmboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_glmboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_glmboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.SVM & Surv.glmboost')
}
result <- rbind(result,res)



##############################################################################
#######################################'@3-8.Surv.SVM&Surv.gbm
##############################################################################
load(file = '3-1.svm_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_gbm <- lrn("surv.gbm",distribution = 'coxph',train.fraction = 0.7,
               n.trees = 1000
)

load('3-8.Surv.SVM&Surv.gbm_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task



if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_gbm,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_gbm,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_gbm,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_gbm,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.SVM & Surv.gbm')
}
result <- rbind(result,res)


##############################################################################
#######################################'@3-9.Surv.SVM&Surv.parametric
##############################################################################
load(file = '3-1.svm_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_parametric <- lrn("surv.parametric"
)

load('3-9.Surv.SVM&Surv.parametric_atf.rda');rid
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task



if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_parametric,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_parametric,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_parametric,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_parametric,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.SVM & Surv.Fully-Parametric')
}
result <- rbind(result,res)



##############################################################################
#######################################'@3-10.Surv.SVM&Surv.aorsf
##############################################################################
load(file = '3-1.svm_atf.rda')
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_aorsf <- lrn("surv.aorsf",n_tree=1000
)

load('3-10.Surv.SVM&Surv.aorsf_atf.rda');rid
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_aorsf,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_aorsf,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_aorsf,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_aorsf,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.SVM & Surv.aorsf')
}
result <- rbind(result,res)


##############################################################################
#######################################'@4-01.Surv.rpart&L2(lambda.min)
##############################################################################
load(file = '4-1.rpart_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 0,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('4-01.Surv.rpart&L2(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.rpart & L2 regularization')
}
result <- rbind(result,res)



##############################################################################
#######################################'@4-02.Surv.rpart&L1(lambda.min)
##############################################################################
load(file = '4-1.rpart_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 1,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('4-02.Surv.rpart&L1(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task



if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.rpart & L1 regularization')
}
result <- rbind(result,res)


##############################################################################
#######################################'@4-1.Surv.rpart&Enet(lambda.min)
##############################################################################
load(file = '4-1.rpart_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 0.1,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('4-1.Surv.rpart&Enet(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task



if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.rpart & Enet regularization')
}
result <- rbind(result,res)



##############################################################################
#######################################'@4-2.Surv.rpart&Surv.RSF
##############################################################################
load(file = '4-1.rpart_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_rfsrc <- lrn("surv.rfsrc", ntree = 2000,importance = 'TRUE')


load('4-2.Surv.rpart&Surv.RSF_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_rfsrc,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_rfsrc,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_rfsrc,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_rfsrc,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.rpart & Surv.RSF')
}
result <- rbind(result,res)


##############################################################################
#######################################'@4-3.Surv.rpart&Surv.SVM
##############################################################################
load(file = '4-1.rpart_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_surv_svm <- lrn("surv.svm",type="vanbelle1",kernel="add_kernel",opt.meth="ipop",
                    diff.meth="makediff1",gamma.mu = 0.1
)

load('4-3.Surv.rpart&Surv.SVM_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_surv_svm,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_surv_svm,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_surv_svm,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_surv_svm,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.rpart & Surv.SVM')
}
result <- rbind(result,res)


##############################################################################
#######################################'@4-5.Surv.rpart&Surv.xgboost
##############################################################################
load(file = '4-1.rpart_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_xgboost <- lrn("surv.xgboost",
                   max_depth = 3,eta = 0.01,gamma = 0,subsample = 0.8,colsample_bytree = 0.8,
                   min_child_weight = 2,lambda=1,alpha=0,nrounds = 200
)

load('4-5.Surv.rpart&Surv.xgboost_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_xgboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_xgboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_xgboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_xgboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.rpart & Surv.xgboost')
}
result <- rbind(result,res)

##############################################################################
#######################################'@4-6.Surv.rpart&Surv.coxboost
##############################################################################
load(file = '4-1.rpart_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task

load('4-6.Surv.rpart&Surv.coxboost_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task
library(CoxBoost)
set.seed(5000)
pen <- suppressWarnings(suppressMessages(optimCoxBoostPenalty(expr.surv[,'ACM'],expr.surv[,'ACM_Censor'],as.matrix(expr.surv[,-c(1,2)]),
                                                              trace=F,start.penalty=1000,parallel = T)));set.seed(5000)
cv.res <- cv.CoxBoost(expr.surv[,'ACM'],expr.surv[,'ACM_Censor'],as.matrix(expr.surv[,-c(1,2)]),
                      maxstepno=1000,K=10,type="verweij",penalty=pen$penalty)
lrn_coxboost <- lrn("surv.coxboost",stepno=cv.res$optimal.step
)

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_coxboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_coxboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_coxboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_coxboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.rpart & Surv.coxboost')
}
result <- rbind(result,res)


##############################################################################
#######################################'@4-7.Surv.rpart&Surv.glmboost
##############################################################################
load(file = '4-1.rpart_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task

set.seed(5000)
lrn_glmboost <- lrn("surv.glmboost",mstop = 500
)

load('4-7.Surv.rpart&Surv.glmboost_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_glmboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_glmboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_glmboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_glmboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.rpart & Surv.glmboost')
}
result <- rbind(result,res)



##############################################################################
#######################################'@4-8.Surv.rpart&Surv.gbm
##############################################################################
load(file = '4-1.rpart_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task

set.seed(5000)
lrn_gbm <- lrn("surv.gbm",distribution = 'coxph',train.fraction = 0.7,
               n.trees = 1000
)

load('4-8.Surv.rpart&Surv.gbm_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task#新特征的task!!


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_gbm,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_gbm,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_gbm,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_gbm,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.rpart & Surv.gbm')
}
result <- rbind(result,res)



##############################################################################
#######################################'@4-9.Surv.rpart&Surv.parametric
##############################################################################
load(file = '4-1.rpart_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task

set.seed(5000)
lrn_parametric <- lrn("surv.parametric"
)

load('4-9.Surv.rpart&Surv.parametric_atf.rda');rid
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_parametric,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_parametric,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_parametric,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_parametric,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.rpart & Surv.Fully-Parametric')
}
result <- rbind(result,res)


##############################################################################
#######################################'@4-10.Surv.rpart&Surv.aorsf
##############################################################################
load(file = '4-1.rpart_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task

set.seed(5000)
lrn_aorsf <- lrn("surv.aorsf",n_tree=1000
)

load('4-10.Surv.rpart&Surv.aorsf_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_aorsf,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_aorsf,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_aorsf,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_aorsf,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.rpart & Surv.aorsf')
}
result <- rbind(result,res)


##############################################################################
#######################################'@5-01.Surv.xgboost&L2(lambda.min)
##############################################################################
load(file = '5-1.xgboost_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 0,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('5-01.Surv.xgboost&L2(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.xgboost & L2 regularization')
}
result <- rbind(result,res)



##############################################################################
#######################################'@5-02.Surv.xgboost&L1(lambda.min)
##############################################################################
load(file = '5-1.xgboost_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 1,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('5-02.Surv.xgboost&L1(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.xgboost & L1 regularization')
}
result <- rbind(result,res)


##############################################################################
#######################################'@5-1.Surv.xgboost&Enet(lambda.min)
##############################################################################
load(file = '5-1.xgboost_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 0.1,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('5-1.Surv.xgboost&Enet(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.xgboost & Enet regularization')
}
result <- rbind(result,res)



##############################################################################
#######################################'@5-2.Surv.xgboost&Surv.RSF
##############################################################################
load(file = '5-1.xgboost_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_rfsrc <- lrn("surv.rfsrc", ntree = 2000,importance = 'TRUE')


load('5-2.Surv.xgboost&Surv.RSF_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_rfsrc,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_rfsrc,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_rfsrc,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_rfsrc,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.xgboost & Surv.RSF')
}
result <- rbind(result,res)



##############################################################################
#######################################'@5-3.Surv.xgboost&Surv.SVM
##############################################################################
load(file = '5-1.xgboost_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_surv_svm <- lrn("surv.svm",type="vanbelle1",kernel="add_kernel",opt.meth="ipop",
                    diff.meth="makediff1",gamma.mu = 0.1
)

load('5-3.Surv.xgboost&Surv.SVM_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_surv_svm,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_surv_svm,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_surv_svm,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_surv_svm,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.xgboost & Surv.SVM')
}
result <- rbind(result,res)


##############################################################################
#######################################'@5-4.Surv.xgboost&Surv.rpart
##############################################################################
load(file = '5-1.xgboost_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_rpart <- lrn("surv.rpart")


load('5-4.Surv.xgboost&Surv.rpart_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_rpart,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_rpart,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_rpart,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_rpart,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.xgboost & Surv.rpart')
}
result <- rbind(result,res)



##############################################################################
#######################################'@5-6.Surv.xgboost&Surv.coxboost
##############################################################################
load(file = '5-1.xgboost_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
library(CoxBoost)

load('5-6.Surv.xgboost&Surv.coxboost_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

set.seed(5000)
pen <- suppressWarnings(suppressMessages(optimCoxBoostPenalty(expr.surv[,'ACM'],expr.surv[,'ACM_Censor'],as.matrix(expr.surv[,-c(1,2)]),
                                                              trace=F,start.penalty=1000,parallel = T)));set.seed(5000)
cv.res <- cv.CoxBoost(expr.surv[,'ACM'],expr.surv[,'ACM_Censor'],as.matrix(expr.surv[,-c(1,2)]),
                      maxstepno=1000,K=10,type="verweij",penalty=pen$penalty)
lrn_coxboost <- lrn("surv.coxboost",stepno=cv.res$optimal.step
)

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_coxboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_coxboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_coxboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_coxboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.xgboost & Surv.coxboost')
}
result <- rbind(result,res)



##############################################################################
#######################################'@5-7.Surv.xgboost&Surv.glmboost
##############################################################################
load(file = '5-1.xgboost_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_glmboost <- lrn("surv.glmboost",mstop = 500
)

load('5-7.Surv.xgboost&Surv.glmboost_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_glmboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_glmboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_glmboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_glmboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.xgboost & Surv.glmboost')
}
result <- rbind(result,res)



##############################################################################
#######################################'@5-8.Surv.xgboost&Surv.gbm
##############################################################################
load(file = '5-1.xgboost_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_gbm <- lrn("surv.gbm",distribution = 'coxph',train.fraction = 0.7,
               n.trees = 1000
)

load('5-8.Surv.xgboost&Surv.gbm_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_gbm,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_gbm,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_gbm,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_gbm,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.xgboost & Surv.gbm')
}
result <- rbind(result,res)



##############################################################################
#######################################'@5-9.Surv.xgboost&Surv.parametric
##############################################################################
load(file = '5-1.xgboost_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_parametric <- lrn("surv.parametric"
)

load('\5-9.Surv.xgboost&Surv.parametric_atf.rda');rid
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_parametric,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_parametric,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_parametric,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_parametric,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.xgboost & Surv.Fully-Parametric')
}
result <- rbind(result,res)



##############################################################################
#######################################'@5-10.Surv.xgboost&Surv.aorsf
##############################################################################
load(file = '5-1.xgboost_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_aorsf <- lrn("surv.aorsf",n_tree=1000
)

load('5-10.Surv.xgboost&Surv.aorsf_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_aorsf,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_aorsf,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_aorsf,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_aorsf,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.xgboost & Surv.aorsf')
}
result <- rbind(result,res)





##############################################################################
#######################################'@6-01.Surv.coxboost&L2(lambda.min)
##############################################################################
load(file = '6-1.coxboost_atf_rid.rda');rid
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 0,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('6-01.Surv.coxboost&L2(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.coxboost & L2 regularization')
}
result <- rbind(result,res)



##############################################################################
#######################################'@6-02.Surv.coxboost&L1(lambda.min)
##############################################################################
load(file = '6-1.coxboost_atf_rid.rda');rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 1,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('6-02.Surv.coxboost&L1(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.coxboost & L1 regularization')
}
result <- rbind(result,res)


##############################################################################
#######################################'@6-1.Surv.coxboost&Enet(lambda.min)
##############################################################################
load(file = '6-1.coxboost_atf_rid.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 0.1,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('6-1.Surv.coxboost&Enet(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.coxboost & Enet regularization')
}
result <- rbind(result,res)



##############################################################################
#######################################'@6-2.Surv.coxboost&Surv.RSF
##############################################################################
load(file = '6-1.coxboost_atf_rid.rda');rid
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_rfsrc <- lrn("surv.rfsrc", ntree = 2000,importance = 'TRUE')

load('6-2.Surv.coxboost&Surv.RSF_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_rfsrc,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_rfsrc,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_rfsrc,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_rfsrc,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.coxboost & Surv.RSF')
}
result <- rbind(result,res)




##############################################################################
#######################################'@6-3.Surv.coxboost&Surv.SVM-----19gene
##############################################################################
load(file = '6-1.coxboost_atf_rid.rda');rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_surv_svm <- lrn("surv.svm",type="vanbelle1",kernel="add_kernel",opt.meth="ipop",
                    diff.meth="makediff1",gamma.mu = 0.1
)

load('6-3.Surv.coxboost&Surv.SVM_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(5000)
  rr_HD <- resample(task,lrn_surv_svm,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(5000)
  rr_SS <- resample(task,lrn_surv_svm,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_surv_svm,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_surv_svm,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.coxboost & Surv.SVM')
}
result <- rbind(result,res)


##############################################################################
#######################################'@6-4.Surv.coxboost&Surv.rpart
##############################################################################
load(file = '6-1.coxboost_atf_rid.rda');rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_rpart <- lrn("surv.rpart")

load('6-4.Surv.coxboost&Surv.rpart_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_rpart,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_rpart,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_rpart,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_rpart,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.coxboost & Surv.rpart')
}
result <- rbind(result,res)


##############################################################################
#######################################'@6-4.Surv.coxboost&Surv.rpart
##############################################################################
load(file = '6-1.coxboost_atf_rid.rda');rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_rpart <- lrn("surv.rpart")


load('6-4.Surv.coxboost&Surv.rpart_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_rpart,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_rpart,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_rpart,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_rpart,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.coxboost & Surv.rpart')
}
result <- rbind(result,res)



##############################################################################
#######################################'@6-5.Surv.coxboost&Surv.xgboost
##############################################################################
load(file = '6-1.coxboost_atf_rid.rda');rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_xgboost <- lrn("surv.xgboost",
                   max_depth = 3,eta = 0.01,gamma = 0,subsample = 0.8,colsample_bytree = 0.8,
                   min_child_weight = 2,lambda=1,alpha=0,nrounds = 200
)

load('6-5.Surv.coxboost&Surv.xgboost_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_xgboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_xgboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_xgboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_xgboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.coxboost & Surv.xgboost')
}
result <- rbind(result,res)




##############################################################################
#######################################'@6-7.Surv.coxboost&Surv.glmboost
##############################################################################
load(file = '6-1.coxboost_atf_rid.rda');rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_glmboost <- lrn("surv.glmboost",mstop = 500
)

load('6-7.Surv.coxboost&Surv.glmboost_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_glmboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_glmboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_glmboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_glmboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.coxboost & Surv.glmboost')
}
result <- rbind(result,res)



##############################################################################
#######################################'@6-8.Surv.coxboost&Surv.gbm
##############################################################################
load(file = '6-1.coxboost_atf_rid.rda');rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_gbm <- lrn("surv.gbm",distribution = 'coxph',train.fraction = 0.7,
               n.trees = 1000
)

load('6-8.Surv.coxboost&Surv.gbm_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_gbm,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_gbm,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_gbm,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_gbm,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.coxboost & Surv.gbm')
}
result <- rbind(result,res)




##############################################################################
#######################################'@6-9.Surv.coxboost&Surv.parametric
##############################################################################
load(file = '6-1.coxboost_atf_rid.rda');rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid ###16个特征！！！！
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_parametric <- lrn("surv.parametric"
)

load('6-9.Surv.coxboost&Surv.parametric_atf.rda');rid
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_parametric,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_parametric,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_parametric,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_parametric,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.coxboost & Surv.parametric')
}
result <- rbind(result,res)




##############################################################################
#######################################'@6-10.Surv.coxboost&Surv.aorsf
##############################################################################
load(file = '6-1.coxboost_atf_rid.rda');rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_aorsf <- lrn("surv.aorsf",n_tree=1000
)

load('6-10.Surv.coxboost&Surv.aorsf_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_aorsf,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_aorsf,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_aorsf,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_aorsf,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.coxboost & Surv.aorsf')
}
result <- rbind(result,res)



##############################################################################
#######################################'@7-01.Surv.glmboost&L2(lambda.min)
##############################################################################
load(file = '7-1.glmboost_atf_rid.rda');rid<-rid[-1];rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 0,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('7-01.Surv.glmboost&L2(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.glmboost & L2 regularization')
}
result <- rbind(result,res)


##############################################################################
#######################################'@7-02.Surv.glmboost&L1(lambda.min)
##############################################################################
load(file = '7-1.glmboost_atf_rid.rda');rid<-rid[-1];rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 1,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('7-02.Surv.glmboost&L1(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.glmboost & L1 regularization')
}
result <- rbind(result,res)



##############################################################################
#######################################'@7-1.Surv.glmboost&Enet(lambda.min)-
##############################################################################
load(file = '7-1.glmboost_atf_rid.rda');rid<-rid[-1];rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 0.1,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('7-1.Surv.glmboost&Enet(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.glmboost & Enet regularization')
}
result <- rbind(result,res)




##############################################################################
#######################################'@7-2.Surv.glmboost&Surv.RSF
##############################################################################
load(file = '7-1.glmboost_atf_rid.rda');rid<-rid[-1];rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_rfsrc <- lrn("surv.rfsrc", ntree = 2000,importance = 'TRUE')


load('7-2.Surv.glmboost&Surv.RSF_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_rfsrc,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_rfsrc,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_rfsrc,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_rfsrc,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.glmboost & Surv.RSF')
}
result <- rbind(result,res)


##############################################################################
#######################################'@7-3.Surv.glmboost&Surv.SVM
##############################################################################
load(file = '7-1.glmboost_atf_rid.rda');rid<-rid[-1];rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_surv_svm <- lrn("surv.svm",type="vanbelle1",kernel="add_kernel",opt.meth="ipop",
                    diff.meth="makediff1",gamma.mu = 0.1
)

load('7-3.Surv.glmboost&Surv.SVM_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_surv_svm,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_surv_svm,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_surv_svm,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_surv_svm,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.glmboost & Surv.SVM')
}
result <- rbind(result,res)


##############################################################################
#######################################'@7-4.Surv.glmboost&Surv.rpart
##############################################################################
load(file = '7-1.glmboost_atf_rid.rda');rid<-rid[-1];rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_rpart <- lrn("surv.rpart")


load('7-4.Surv.glmboost&Surv.rpart_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_rpart,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_rpart,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_rpart,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_rpart,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.glmboost & Surv.rpart')
}
result <- rbind(result,res)


##############################################################################
#######################################'@7-5.Surv.glmboost&Surv.xgboost
##############################################################################
load(file = '7-1.glmboost_atf_rid.rda');rid<-rid[-1];rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_xgboost <- lrn("surv.xgboost",
                   max_depth = 3,eta = 0.01,gamma = 0,subsample = 0.8,colsample_bytree = 0.8,
                   min_child_weight = 2,lambda=1,alpha=0,nrounds = 200
)

load('7-5.Surv.glmboost&Surv.xgboost_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_xgboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_xgboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_xgboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_xgboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.glmboost & Surv.xgboost')
}
result <- rbind(result,res)



##############################################################################
#######################################'@7-6.Surv.glmboost&Surv.coxboost
##############################################################################
library(CoxBoost)
load(file = '7-1.glmboost_atf_rid.rda');rid<-rid[-1];rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task


load('7-6.Surv.glmboost&Surv.coxboost_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

set.seed(5000)
pen <- suppressWarnings(suppressMessages(optimCoxBoostPenalty(expr.surv[,'ACM'],expr.surv[,'ACM_Censor'],as.matrix(expr.surv[,-c(1,2)]),
                                                              trace=F,start.penalty=1000,parallel = T)));set.seed(5000)
cv.res <- cv.CoxBoost(expr.surv[,'ACM'],expr.surv[,'ACM_Censor'],as.matrix(expr.surv[,-c(1,2)]),
                      maxstepno=1000,K=10,type="verweij",penalty=pen$penalty)
lrn_coxboost <- lrn("surv.coxboost",stepno=cv.res$optimal.step
)

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_coxboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_coxboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_coxboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_coxboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.glmboost & Surv.coxboost')
}
result <- rbind(result,res)



##############################################################################
#######################################'@7-8.Surv.glmboost&Surv.gbm
##############################################################################
library(CoxBoost)
load(file = '7-1.glmboost_atf_rid.rda');rid<-rid[-1];rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_gbm <- lrn("surv.gbm",distribution = 'coxph',train.fraction = 0.7,
               n.trees = 1000
)


load('7-8.Surv.glmboost&Surv.gbm_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_gbm,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_gbm,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_gbm,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_gbm,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.glmboost & Surv.gbm')
}
result <- rbind(result,res)


##############################################################################
#######################################'@7-9.Surv.glmboost&Surv.parametric
##############################################################################
library(CoxBoost)
load(file = '7-1.glmboost_atf_rid.rda');rid<-rid[-1];rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_parametric <- lrn("surv.parametric"
)


load('7-9.Surv.glmboost&Surv.parametric_atf.rda');rid
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_parametric,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_parametric,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_parametric,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_parametric,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.glmboost & Surv.Fully-Parametric')
}
result <- rbind(result,res)



##############################################################################
#######################################'@7-10.Surv.glmboost&Surv.aorsf
##############################################################################
library(CoxBoost)
load(file = '7-1.glmboost_atf_rid.rda');rid<-rid[-1];rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_aorsf <- lrn("surv.aorsf",n_tree=1000
)


load('7-10.Surv.glmboost&Surv.aorsf_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_aorsf,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_aorsf,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_aorsf,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_aorsf,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.glmboost & Surv.aorsf')
}
result <- rbind(result,res)



##############################################################################
#######################################'@8-01.Surv.gbm&L2(lambda.min)
##############################################################################
load(file = '8-1.gbm_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 0,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('8-01.Surv.gbm&L2(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.gbm & L2 regularization')
}
result <- rbind(result,res)



##############################################################################
#######################################'@8-02.Surv.gbm&L1(lambda.min)
##############################################################################
load(file = '8-1.gbm_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 1,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('8-02.Surv.gbm&L1(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.gbm & L1 regularization')
}
result <- rbind(result,res)


##############################################################################
#######################################'@8-1.Surv.gbm&Enet(lambda.min)
##############################################################################
load(file = '8-1.gbm_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 0.1,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('8-1.Surv.gbm&Enet(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.gbm & Enet regularization')
}
result <- rbind(result,res)


##############################################################################
#######################################'@8-2.Surv.gbm&Surv.RSF
##############################################################################
load(file = '8-1.gbm_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_rfsrc <- lrn("surv.rfsrc", ntree = 2000,importance = 'TRUE')


load('8-2.Surv.gbm&Surv.RSF_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_rfsrc,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_rfsrc,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_rfsrc,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_rfsrc,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.gbm & Surv.RSF')
}
result <- rbind(result,res)



##############################################################################
#######################################'@8-3.Surv.gbm&Surv.SVM
##############################################################################
load(file = '8-1.gbm_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_surv_svm <- lrn("surv.svm",type="vanbelle1",kernel="add_kernel",opt.meth="ipop",
                    diff.meth="makediff1",gamma.mu = 0.1
)

load('8-3.Surv.gbm&Surv.SVM_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_surv_svm,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_surv_svm,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_surv_svm,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_surv_svm,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.gbm & Surv.SVM')
}
result <- rbind(result,res)




##############################################################################
#######################################'@8-4.Surv.gbm&Surv.rpart
##############################################################################
load(file = '8-1.gbm_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_rpart <- lrn("surv.rpart")


load('8-4.Surv.gbm&Surv.rpart_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_rpart,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_rpart,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_rpart,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_rpart,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.gbm & Surv.rpart')
}
result <- rbind(result,res)


##############################################################################
#######################################'@8-5.Surv.gbm&Surv.xgboost
##############################################################################
load(file = '8-1.gbm_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_xgboost <- lrn("surv.xgboost",
                   max_depth = 3,eta = 0.01,gamma = 0,subsample = 0.8,colsample_bytree = 0.8,
                   min_child_weight = 2,lambda=1,alpha=0,nrounds = 200
)

load('8-5.Surv.gbm&Surv.xgboost_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_xgboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_xgboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_xgboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_xgboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.gbm & Surv.xgboost')
}
result <- rbind(result,res)


##############################################################################
#######################################'@8-6.Surv.gbm&Surv.coxboost
##############################################################################
library(CoxBoost)
load(file = '8-1.gbm_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task


load('8-6.Surv.gbm&Surv.coxboost_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

set.seed(5000)
pen <- suppressWarnings(suppressMessages(optimCoxBoostPenalty(expr.surv[,'ACM'],expr.surv[,'ACM_Censor'],as.matrix(expr.surv[,-c(1,2)]),
                                                              trace=F,start.penalty=1000,parallel = T)));set.seed(5000)
cv.res <- cv.CoxBoost(expr.surv[,'ACM'],expr.surv[,'ACM_Censor'],as.matrix(expr.surv[,-c(1,2)]),
                      maxstepno=1000,K=10,type="verweij",penalty=pen$penalty)
lrn_coxboost <- lrn("surv.coxboost",stepno=cv.res$optimal.step
)

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_coxboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_coxboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_coxboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_coxboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.gbm & Surv.coxboost')
}
result <- rbind(result,res)



##############################################################################
#######################################'@8-7.Surv.gbm&Surv.glmboost
##############################################################################
load(file = '8-1.gbm_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_glmboost <- lrn("surv.glmboost",mstop = 500
)

load('8-7.Surv.gbm&Surv.glmboost_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_glmboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_glmboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_glmboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_glmboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.gbm & Surv.glmboost')
}
result <- rbind(result,res)


##############################################################################
#######################################'@8-9.Surv.gbm&Surv.parametric
##############################################################################
load(file = '8-1.gbm_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_parametric <- lrn("surv.parametric"
)

load('8-9.Surv.gbm&Surv.parametric_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_parametric,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_parametric,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_parametric,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_parametric,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.gbm & Surv.Fully-Parametric')
}
result <- rbind(result,res)




##############################################################################
#######################################'@8-10.Surv.gbm&Surv.aorsf
##############################################################################
load(file = '8-1.gbm_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_aorsf <- lrn("surv.aorsf",n_tree=1000
)

load('8-10.Surv.gbm&Surv.aorsf_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_aorsf,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_aorsf,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_aorsf,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_aorsf,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.gbm & Surv.aorsf')
}
result <- rbind(result,res)




##############################################################################
#######################################'@9-01.Surv.parametric&L2(lambda.min)
##############################################################################
load(file = '9-1.parametric_atf_rid.rda');rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 0,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('9-01.Surv.parametric&L2(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.Fully-Parametric & L2 regularization')
}
result <- rbind(result,res)




##############################################################################
#######################################'@9-02.Surv.parametric&L1(lambda.min)
##############################################################################
load(file = '9-1.parametric_atf_rid.rda');rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 1,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('9-02.Surv.parametric&L1(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.Fully-Parametric & L1 regularization')
}
result <- rbind(result,res)


##############################################################################
#######################################'@9-1.Surv.parametric&Enet(lambda.min)
##############################################################################
load(file = '9-1.parametric_atf_rid.rda');rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 0.1,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('9-02.Surv.parametric&L1(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.Fully-Parametric & Enet regularization')
}
result <- rbind(result,res)




##############################################################################
#######################################'@9-2.Surv.parametric&Surv.RSF
##############################################################################
load(file = '9-1.parametric_atf_rid.rda');rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_rfsrc <- lrn("surv.rfsrc", ntree = 2000,importance = 'TRUE')


load('9-2.Surv.parametric&Surv.RSF_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_rfsrc,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_rfsrc,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_rfsrc,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_rfsrc,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.Fully-Parametric & Surv.RSF')
}
result <- rbind(result,res)


##############################################################################
#######################################'@9-3.Surv.parametric&Surv.SVM
##############################################################################
load(file = '9-1.parametric_atf_rid.rda');rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_surv_svm <- lrn("surv.svm",type="vanbelle1",kernel="add_kernel",opt.meth="ipop",
                    diff.meth="makediff1",gamma.mu = 0.1
)

load('9-3.Surv.parametric&Surv.SVM_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_surv_svm,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_surv_svm,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_surv_svm,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_surv_svm,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.Fully-Parametric & Surv.SVM')
}
result <- rbind(result,res)


##############################################################################
#######################################'@9-4.Surv.parametric&Surv.rpart
##############################################################################
load(file = '9-1.parametric_atf_rid.rda');rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_rpart <- lrn("surv.rpart")


load('9-4.Surv.parametric&Surv.rpart_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_rpart,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_rpart,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_rpart,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_rpart,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.Fully-Parametric & Surv.rpart')
}
result <- rbind(result,res)


##############################################################################
#######################################'@9-5.Surv.parametric&Surv.xgboost
##############################################################################
load(file = '9-1.parametric_atf_rid.rda');rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_xgboost <- lrn("surv.xgboost",
                   max_depth = 3,eta = 0.01,gamma = 0,subsample = 0.8,colsample_bytree = 0.8,
                   min_child_weight = 2,lambda=1,alpha=0,nrounds = 200
)

load('9-5.Surv.parametric&Surv.xgboost_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_xgboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_xgboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_xgboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_xgboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.Fully-Parametric & Surv.xgboost')
}
result <- rbind(result,res)



##############################################################################
#######################################'@9-6.Surv.parametric&Surv.coxboost
##############################################################################
library(CoxBoost)
load(file = '9-1.parametric_atf_rid.rda');rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task

load('9-6.Surv.parametric&Surv.coxboost_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

set.seed(5000)
pen <- suppressWarnings(suppressMessages(optimCoxBoostPenalty(expr.surv[,'ACM'],expr.surv[,'ACM_Censor'],as.matrix(expr.surv[,-c(1,2)]),
                                                              trace=F,start.penalty=1000,parallel = T)));set.seed(5000)
cv.res <- cv.CoxBoost(expr.surv[,'ACM'],expr.surv[,'ACM_Censor'],as.matrix(expr.surv[,-c(1,2)]),
                      maxstepno=1000,K=10,type="verweij",penalty=pen$penalty)
lrn_coxboost <- lrn("surv.coxboost",stepno=cv.res$optimal.step
)

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_coxboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_coxboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_coxboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_coxboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.Fully-Parametric & Surv.coxboost')
}
result <- rbind(result,res)



##############################################################################
#######################################'@9-7.Surv.parametric&Surv.glmboost
##############################################################################
load(file = '9-1.parametric_atf_rid.rda');rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_glmboost <- lrn("surv.glmboost",mstop = 500
)

load('9-7.Surv.parametric&Surv.glmboost_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_glmboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_glmboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_glmboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_glmboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.Fully-Parametric & Surv.glmboost')
}
result <- rbind(result,res)




##############################################################################
#######################################'@9-8.Surv.parametric&Surv.gbm
##############################################################################
load(file = '9-1.parametric_atf_rid.rda');rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_gbm <- lrn("surv.gbm",distribution = 'coxph',train.fraction = 0.7,
               n.trees = 1000
)

load('9-8.Surv.parametric&Surv.gbm_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_gbm,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_gbm,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_gbm,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_gbm,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.Fully-Parametric & Surv.gbm')
}
result <- rbind(result,res)




##############################################################################
#######################################'@9-10.Surv.parametric&Surv.aorsf
##############################################################################
load(file = '9-1.parametric_atf_rid.rda');rid
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_aorsf <- lrn("surv.aorsf",n_tree=1000
)

load('9-10.Surv.parametric&Surv.aorsf_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_aorsf,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_aorsf,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_aorsf,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_aorsf,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.Fully-Parametric & Surv.aorsf')
}
result <- rbind(result,res)



##############################################################################
#######################################'@10-01.Surv.aorsf&L2(lambda.min)
##############################################################################
load(file = '10-1.aorsf_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 0,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('10-01.Surv.aorsf&L2(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.aorsf & L2 regularization')
}
result <- rbind(result,res)



##############################################################################
#######################################'@10-02.Surv.aorsf&L1(lambda.min)
##############################################################################
load(file = '10-1.aorsf_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 1,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('10-02.Surv.aorsf&L1(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task

if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.aorsf & L1 regularization')
}
result <- rbind(result,res)




##############################################################################
#######################################'@10-1.Surv.aorsf&Enet(lambda.min)
##############################################################################
load(file = '10-1.aorsf_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_cv_glmnet <- lrn("surv.cv_glmnet",alpha = 0.1,nfolds = 10, s = 'lambda.min',type.measure = "C"
)

load('10-1.Surv.aorsf&Enet(lambda.min)_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task#新特征的task!!


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_cv_glmnet,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_cv_glmnet,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_cv_glmnet,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_cv_glmnet,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.aorsf & Enet regularization')
}
result <- rbind(result,res)




##############################################################################
#######################################'@10-2.Surv.aorsf&Surv.RSF
##############################################################################
load(file = '10-1.aorsf_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_rfsrc <- lrn("surv.rfsrc", ntree = 2000,importance = 'TRUE')


load('10-2.Surv.aorsf&Surv.RSF_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_rfsrc,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_rfsrc,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_rfsrc,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_rfsrc,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.aorsf & Surv.RSF')
}
result <- rbind(result,res)




##############################################################################
#######################################'@10-3.Surv.aorsf&Surv.SVM
##############################################################################
load(file = '10-1.aorsf_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_surv_svm <- lrn("surv.svm",type="vanbelle1",kernel="add_kernel",opt.meth="ipop",
                    diff.meth="makediff1",gamma.mu = 0.1
)

load('10-3.Surv.aorsf&Surv.SVM_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_surv_svm,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_surv_svm,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_surv_svm,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_surv_svm,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.aorsf & Surv.SVM')
}
result <- rbind(result,res)




##############################################################################
#######################################'@10-4.Surv.aorsf&Surv.rpart
##############################################################################
load(file = '10-1.aorsf_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_rpart <- lrn("surv.rpart")


load('10-4.Surv.aorsf&Surv.rpart_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_rpart,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_rpart,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_rpart,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_rpart,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.aorsf & Surv.rpart')
}
result <- rbind(result,res)



##############################################################################
#######################################'@10-5.Surv.aorsf&Surv.xgboost
##############################################################################
load(file = '10-1.aorsf_atf.rda')
atf$base_learner()
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid 
task=as_task_surv(expr.surv,time="ACM",event="ACM_Censor");task
task$select(rid);task
set.seed(5000)
lrn_xgboost <- lrn("surv.xgboost",
                   max_depth = 3,eta = 0.01,gamma = 0,subsample = 0.8,colsample_bytree = 0.8,
                   min_child_weight = 2,lambda=1,alpha=0,nrounds = 200
)

load('10-5.Surv.aorsf&Surv.xgboost_atf.rda')
atf$base_learner();atf$fselect_result
rid<-atf$fselect_result%>%t()%>%.[ncol(atf$fselect_result)-1]%>%.[[1]];rid
task$select(rid);task


if(T){
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout",ratio = 0.7);print(resampling_HD);set.seed(1000)
  rr_HD <- resample(task,lrn_xgboost,resampling_HD,store_models=T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate();index_HD
  
  resampling_SS <- rsmp("subsampling",repeats = 10,ratio = 0.7);print(resampling_SS);set.seed(1000)
  rr_SS <- resample(task,lrn_xgboost,resampling_SS,store_models=T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate();index_SS
  
  resampling_CV <- rsmp("repeated_cv",folds = 5,repeats= 1);print(resampling_CV);set.seed(5000)
  rr_CV <- resample(task,lrn_xgboost,resampling_CV,store_models=T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate();index_CV
  
  resampling_BS <- rsmp("bootstrap",repeats= 10);print(resampling_BS);set.seed(5000)
  rr_BS <- resample(task,lrn_xgboost,resampling_BS,store_models=T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate();index_BS
  
  res = data.frame(
    holdout = index_HD,subsampling = index_SS,repeated_cv = index_CV,bootstrap = index_BS
  )%>%t()%>%as.data.frame()%>% rownames_to_column('ID');res$mod <- paste0('Surv.aorsf & Surv.xgboost')
}
result <- rbind(result,res)




##############################################################################
#######################################' @10-6.Surv.aorsf&Surv.coxboost
##############################################################################
library(CoxBoost)
load(file = "10-1.aorsf_atf.rda")
atf$base_learner()
rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
task$select(rid)
task

load("10-6.Surv.aorsf&Surv.coxboost_atf.rda")
atf$base_learner()
atf$fselect_result
rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task$select(rid)
task

set.seed(5000)
pen <- suppressWarnings(suppressMessages(optimCoxBoostPenalty(expr.surv[, "ACM"], expr.surv[, "ACM_Censor"], as.matrix(expr.surv[, -c(1, 2)]),
  trace = F, start.penalty = 1000, parallel = T
)))
set.seed(5000)
cv.res <- cv.CoxBoost(expr.surv[, "ACM"], expr.surv[, "ACM_Censor"], as.matrix(expr.surv[, -c(1, 2)]),
  maxstepno = 1000, K = 10, type = "verweij", penalty = pen$penalty
)
lrn_coxboost <- lrn("surv.coxboost", stepno = cv.res$optimal.step)

if (T) {
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout", ratio = 0.7)
  print(resampling_HD)
  set.seed(1000)
  rr_HD <- resample(task, lrn_coxboost, resampling_HD, store_models = T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate()
  index_HD

  resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
  print(resampling_SS)
  set.seed(1000)
  rr_SS <- resample(task, lrn_coxboost, resampling_SS, store_models = T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate()
  index_SS

  resampling_CV <- rsmp("repeated_cv", folds = 5, repeats = 1)
  print(resampling_CV)
  set.seed(5000)
  rr_CV <- resample(task, lrn_coxboost, resampling_CV, store_models = T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate()
  index_CV

  resampling_BS <- rsmp("bootstrap", repeats = 10)
  print(resampling_BS)
  set.seed(5000)
  rr_BS <- resample(task, lrn_coxboost, resampling_BS, store_models = T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate()
  index_BS

  res <- data.frame(
    holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  res$mod <- paste0("Surv.aorsf & Surv.coxboost")
}
result <- rbind(result, res)



##############################################################################
#######################################' @10-7.Surv.aorsf&Surv.glmboost
##############################################################################
library(CoxBoost)
load(file = "10-1.aorsf_atf.rda")
atf$base_learner()
rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
task$select(rid)
task

load("10-7.Surv.aorsf&Surv.glmboost_atf.rda")
atf$base_learner()
atf$fselect_result
rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task$select(rid)
task
set.seed(5000)
lrn_glmboost <- lrn("surv.glmboost", mstop = 500)

if (T) {
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout", ratio = 0.7)
  print(resampling_HD)
  set.seed(1000)
  rr_HD <- resample(task, lrn_glmboost, resampling_HD, store_models = T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate()
  index_HD

  resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
  print(resampling_SS)
  set.seed(1000)
  rr_SS <- resample(task, lrn_glmboost, resampling_SS, store_models = T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate()
  index_SS

  resampling_CV <- rsmp("repeated_cv", folds = 5, repeats = 1)
  print(resampling_CV)
  set.seed(5000)
  rr_CV <- resample(task, lrn_glmboost, resampling_CV, store_models = T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate()
  index_CV

  resampling_BS <- rsmp("bootstrap", repeats = 10)
  print(resampling_BS)
  set.seed(5000)
  rr_BS <- resample(task, lrn_glmboost, resampling_BS, store_models = T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate()
  index_BS

  res <- data.frame(
    holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  res$mod <- paste0("Surv.aorsf & Surv.glmboost")
}
result <- rbind(result, res)


##############################################################################
#######################################' @10-8.Surv.aorsf&Surv.gbm
##############################################################################
library(CoxBoost)
load(file = "10-1.aorsf_atf.rda")
atf$base_learner()
rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
task$select(rid)
task

load("10-8.Surv.aorsf&Surv.gbm_atf.rda")
atf$base_learner()
atf$fselect_result
rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task$select(rid)
task

set.seed(5000)
lrn_gbm <- lrn("surv.gbm",
  distribution = "coxph", train.fraction = 0.7,
  n.trees = 1000
)

if (T) {
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout", ratio = 0.7)
  print(resampling_HD)
  set.seed(1000)
  rr_HD <- resample(task, lrn_gbm, resampling_HD, store_models = T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate()
  index_HD

  resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
  print(resampling_SS)
  set.seed(1000)
  rr_SS <- resample(task, lrn_gbm, resampling_SS, store_models = T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate()
  index_SS

  resampling_CV <- rsmp("repeated_cv", folds = 5, repeats = 1)
  print(resampling_CV)
  set.seed(5000)
  rr_CV <- resample(task, lrn_gbm, resampling_CV, store_models = T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate()
  index_CV

  resampling_BS <- rsmp("bootstrap", repeats = 10)
  print(resampling_BS)
  set.seed(5000)
  rr_BS <- resample(task, lrn_gbm, resampling_BS, store_models = T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate()
  index_BS

  res <- data.frame(
    holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  res$mod <- paste0("Surv.aorsf & Surv.gbm")
}
result <- rbind(result, res)




##############################################################################
#######################################' @10-9.Surv.aorsf&Surv.parametric
##############################################################################
load(file = "10-1.aorsf_atf.rda")
atf$base_learner()
rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task <- as_task_surv(expr.surv, time = "ACM", event = "ACM_Censor")
task
task$select(rid)
task

load("10-9.Surv.aorsf&Surv.parametric_atf.rda")
rid
atf$base_learner()
atf$fselect_result
rid <- atf$fselect_result %>%
  t() %>%
  .[ncol(atf$fselect_result) - 1] %>%
  .[[1]]
rid
task$select(rid)
task

set.seed(5000)
lrn_parametric <- lrn("surv.parametric")

if (T) {
  measures <- c("surv.cindex")
  resampling_HD <- rsmp("holdout", ratio = 0.7)
  print(resampling_HD)
  set.seed(1000)
  rr_HD <- resample(task, lrn_parametric, resampling_HD, store_models = T)
  index_HD <- msrs(measures) %>% rr_HD$aggregate()
  index_HD

  resampling_SS <- rsmp("subsampling", repeats = 10, ratio = 0.7)
  print(resampling_SS)
  set.seed(1000)
  rr_SS <- resample(task, lrn_parametric, resampling_SS, store_models = T)
  index_SS <- msrs(measures) %>% rr_SS$aggregate()
  index_SS

  resampling_CV <- rsmp("repeated_cv", folds = 5, repeats = 1)
  print(resampling_CV)
  set.seed(5000)
  rr_CV <- resample(task, lrn_parametric, resampling_CV, store_models = T)
  index_CV <- msrs(measures) %>% rr_CV$aggregate()
  index_CV

  resampling_BS <- rsmp("bootstrap", repeats = 10)
  print(resampling_BS)
  set.seed(5000)
  rr_BS <- resample(task, lrn_parametric, resampling_BS, store_models = T)
  index_BS <- msrs(measures) %>% rr_BS$aggregate()
  index_BS

  res <- data.frame(
    holdout = index_HD, subsampling = index_SS, repeated_cv = index_CV, bootstrap = index_BS
  ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  res$mod <- paste0("Surv.aorsf & Surv.Fully-Parametric")
}
result <- rbind(result, res)







