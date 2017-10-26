#### Analysis of field data
# Sam Flake, 12 Sep 2017
# sflake@gmail.com or swflake@ncsu.edu


# load libraries
library(plyr)
library(boot)
library(lme4)
library(MuMIn)
library(car)

set.seed(665224291)
options(na.action = na.fail)
pardefault <- par(no.readonly = T)

#Source some useful functions
source("pj_mortality_functions_090617.R")

#read prepared data
all_tree_vars <- read.csv("./Clean data/all_tree_Vars_scaled_and_transformed.csv")

#--------------------------------------------------------------------------------------------
# tree-level models

#global model with tree-level vars and abiotic vars
pimo_global <- lmer(formula = Delta_pdc ~ Neighbor_larger_log + ENN_dist_log + Diam_log + Height + BA_4m_log +  
                          cwd_normal_cum + fdsi_anom +
                          cwd_peak_anom + tmean + ppt +
                          vpd_normal_annual_max + vpd_peak_anom + avg_summer_tmax +
                          peak_tmax + min_ann_ppt +
                          Pndjfm + Avg_depth + AWC +
                          (1|Cluster/site), 
                      data = all_tree_vars[all_tree_vars$Spp == "PIMO", ],
                    na.action = na.fail)

#only use combos that aren't pairwise correlated
cormat <- abs(cor(all_tree_vars[all_tree_vars$Spp == "PIMO", c('Neighbor_larger_log', 'ENN_dist_log', 'Diam_log', 'Height', 'BA_4m_log', 
                                                               'cwd_normal_cum', 'fdsi_anom',
                                                                   'cwd_peak_anom' , 'tmean', 'ppt',
                                                                   'vpd_normal_annual_max' , 'vpd_peak_anom' , 'avg_summer_tmax', 'peak_tmax',
                                                                   'min_ann_ppt',  'Pndjfm' , 'Avg_depth' , 'AWC')])) <=.5

#use all tree-level vars
cormat[1:5,1:5] <- TRUE

#use only one side of matrix (avoid duplicated models)
cormat[!lower.tri(cormat)] <- NA

#find all combinations with a few given conditions: 
# abiotic variables can change but stand structure variables are fixed;
# abiotic variables with correlations greater than .5 are not allowed;
# 9 maximum variables are allowed;
# progress is displayed in console

# this only takes a minute or two
pimo_dredge <- dredge(pimo_global, 
                      fixed = c('Neighbor_larger_log', 'ENN_dist_log', 'Diam_log', 'Height', 'BA_4m_log'),
                      subset = cormat, 
                      m.lim = c(5, 8), 
                      trace = 2)

#get call for top model and evaluate it
best_pimo_call <- as.character(unlist(attr((pimo_dredge)[1], "model.calls")))
pimo_model_dredge <- eval(parse(text = best_pimo_call))
summary(pimo_model_dredge)

pimo_dredge[1:10] #top 10 models

# #----------------------------------------------
# 
#global model with tree-level vars and abiotic vars
juos_global <- lmer(formula = Delta_pdc ~ Neighbor_larger + ENN_dist + Height + Diam + BA_4m +  
                      cwd_normal_cum + fdsi_anom +
                      cwd_peak_anom + tmean + ppt +
                      vpd_normal_annual_max + vpd_peak_anom + avg_summer_tmax +
                      peak_tmax + min_ann_ppt +
                      Pndjfm + Avg_depth + AWC +
                      (1|Cluster/site), 
                    data = all_tree_vars[all_tree_vars$Spp == "JUOS", ],
                    na.action = na.fail)

#only use combos that aren't pairwise correlated
cormat <- abs(cor(all_tree_vars[all_tree_vars$Spp == "JUOS", c('Neighbor_larger', 'ENN_dist', 'Height', 'Diam', 'BA_4m', 
                                                               'cwd_normal_cum', 'fdsi_anom',
                                                               'cwd_peak_anom' , 'tmean', 'ppt',
                                                               'vpd_normal_annual_max' , 'vpd_peak_anom' , 'avg_summer_tmax', 'peak_tmax',
                                                               'min_ann_ppt',  'Pndjfm' , 'Avg_depth' , 'AWC')])) <=.5

#use all tree-level vars
cormat[1:5,1:5] <- TRUE

#use only one side of matrix (avoid duplicated models)
cormat[!lower.tri(cormat)] <- NA

#find all combinations
juos_dredge <- dredge(juos_global, 
                      fixed = c('Neighbor_larger', 'ENN_dist', 'Height', 'Diam', 'BA_4m'),
                      subset = cormat, 
                      m.lim = c(5, 8), 
                      trace = 2)


#get call for top model and evaluate it
best_juos_call <- as.character(unlist(attr((juos_dredge)[1], "model.calls")))
juos_model_dredge <- eval(parse(text = best_juos_call))


# 
# ###################################
# # PIMO mortality
# ###################################
pmort_global <- glmer(formula = Died ~ Neighbor_larger_log + ENN_dist_log + Diam_log + Height + BA_4m_log +  
                      cwd_normal_cum + fdsi_anom +
                      cwd_peak_anom + tmean + ppt +
                      vpd_normal_annual_max + vpd_peak_anom + avg_summer_tmax +
                      peak_tmax + min_ann_ppt +
                      Pndjfm + Avg_depth + AWC +
                      (1|site), family = "binomial",
                      data = all_tree_vars[all_tree_vars$Spp == "PIMO", ], na.action = na.fail)

#only use combos that aren't pairwise correlated
cormat <- abs(cor(all_tree_vars[all_tree_vars$Spp == "PIMO", c('Neighbor_larger_log', 'ENN_dist_log', 'Diam_log', 'Height', 'BA_4m_log', 
                                                               'cwd_normal_cum', 'fdsi_anom',
                                                               'cwd_peak_anom' , 'tmean', 'ppt',
                                                               'vpd_normal_annual_max' , 'vpd_peak_anom' , 'avg_summer_tmax', 'peak_tmax',
                                                               'min_ann_ppt',  'Pndjfm' , 'Avg_depth' , 'AWC')])) <=.5

#use all tree-level vars
cormat[1:5,1:5] <- TRUE

#use only one side of matrix (avoid duplicated models)
cormat[!lower.tri(cormat)] <- NA

# find all combinations with a few given conditions: 
# abiotic variables can change but stand structure variables are fixed;
# abiotic variables with correlations greater than .5 are not allowed;
# 8 maximum variables are allowed;
# progress is displayed in console

# this only takes a half an hour or so
start.time <- Sys.time()
pmort_dredge <- dredge(pmort_global, 
                      fixed = c('Neighbor_larger_log', 'ENN_dist_log', 'Diam_log', 'Height', 'BA_4m_log'),
                      subset = cormat, 
                      m.lim = c(5, 8), 
                      trace = 2)
end.time <- Sys.time()

time.taken_pmort_dredge <- end.time - start.time
time.taken_pmort_dredge #30 minutes
(pmort_dredge)[c(1:10)]

#get call for top model and evaluate it
best_pmort_call <- as.character(unlist(attr((pmort_dredge)[1], "model.calls")))
pmort_model_dredge <- eval(parse(text = best_pmort_call))


#############################################################################################
#Tree-level models
#############################################################################################

#PIMO dieback

#add interaction term to the PIMO dieback model
pimo_int_ba <- update(pimo_model_dredge, . ~ . + BA_4m*cwd_normal_cum)

summary(pimo_int_ba)

AICc(pimo_model_dredge) - AICc(pimo_int_ba) #improves AICc by 3.8

#nonparametric bootstrap cis
bs <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- lmer(formula, data=d)
  return(fixef(fit)) 
} 

bs(formula = pimo_int_ba@call,
   data = all_tree_vars[all_tree_vars$Spp == "PIMO", ])

# bootstrapping with 5000 replications
#this takes a while (~6 hrs on an i7 3.4 Ghz processor. Would be faster with multithreading)
start.time <- Sys.time()
bootstrap_pimo_int <- boot(data=all_tree_vars[all_tree_vars$Spp == "PIMO", ], 
                                statistic=bs, 
                                R=5000, 
                                formula= pimo_int_ba@call,
                                parallel = "no")
end.time <- Sys.time()
time.taken_pimo <- end.time - start.time
time.taken_pimo

saveRDS(bootstrap_pimo_int, file = paste0("./model output/bootstrap_pimo_int", Sys.Date(), ".rds"))

# view results
bootstrap_pimo <- bootstrap_pimo_int
summary(bootstrap_pimo)
plot(bootstrap_pimo, index=1) # intercept
plot(bootstrap_pimo, index=2) #
plot(bootstrap_pimo, index=3) #

# get 95% confidence intervals 
ci_pimo <- cbind(attr(bootstrap_pimo_int$t0, "names"), summary(bootstrap_pimo_int))
ci_pimo$lo <- NA
ci_pimo$hi <- NA

for( i in 1:nrow(ci_pimo)){
  ci_boot <- boot.ci(bootstrap_pimo_int, type = "bca", index = i)
  #extract low and hi, add to ci_pimo df
  ci_pimo$lo[i] <- ci_boot$bca[1,4]
  ci_pimo$hi[i] <- ci_boot$bca[1,5]
  ci_pimo$mean[i] <- mean(bootstrap_pimo_int$t[, i])
}

saveRDS(ci_pimo, file = "./model output/ci_pimo.rds")
write.csv(ci_pimo, "./model output/ci_pimo.csv")

# plot(allEffects(pimo_model))
#    plot(Effect(focal.predictors = c("ENN_dist", "Neighbor_larger"), mod = pimo_model, 
#                ylim = c(-40, 0), partial.residuals = TRUE))
#    plot(Effect(focal.predictors = c("BA_4m", "cwd_normal_cum"), mod = pimo_int_ba, 
#                ylim = c(-40, 0), partial.residuals = TRUE))
# coefplot2::coefplot2(pimo_model)
# AICc(pimo_model)
# vif.mer(pimo_model)
# r.squaredGLMM(pimo_int_ba)
# 
# plot(resid(pimo_int_ba) ~ predict(pimo_int_ba))
# abline(h=0)
# hist(resid(pimo_int_ba))
# mse <- sum(resid(pimo_int_ba)^2)*(1 / length(resid(pimo_int_ba)))
# 
# 
# qqmath(ranef(pimo_model, postVar = TRUE))

#--------------------------------------------------------------

# # AICc(juos_int_ba)
juos_int_enn <- update(juos_model_dredge, . ~ . + ENN_dist*Neighbor_larger)

bs(formula = juos_int_enn@call,
   data = all_tree_vars[all_tree_vars$Spp == "JUOS", ])

# bootstrapping with 5000 replications
#this takes a while (~3.5 hrs on an i7 3.4 Ghz processor. Would be faster with multithreading)
start.time <- Sys.time()
bootstrap_juos <- boot(data=all_tree_vars[all_tree_vars$Spp == "JUOS", ], 
                       statistic=bs, 
                       R=5000, 
                       formula=juos_int_enn@call)


end.time <- Sys.time()
time.taken_juos <- end.time - start.time
time.taken_juos

saveRDS(bootstrap_juos, file = "./model output/bootstrap_juos.rds")

# view results
# bootstrap_juos
# summary(bootstrap_juos)
# plot(bootstrap_juos, index=1) # intercept 
# plot(bootstrap_juos, index=2) # 
# plot(bootstrap_juos, index=3) # 
# plot(bootstrap_juos, index=4) # 
# plot(bootstrap_juos, index=5) # 
# plot(bootstrap_juos, index=6) # 
# plot(bootstrap_juos, index=7) # 
# plot(bootstrap_juos, index=8) # 
# plot(bootstrap_juos, index=9) # 
# plot(bootstrap_juos, index=10) # 

# get 95% confidence intervals 
ci_juos <- cbind(attr(bootstrap_juos$t0, "names"), summary(bootstrap_juos))
ci_juos$lo <- NA
ci_juos$hi <- NA

for( i in 1:nrow(ci_juos)){
  ci_boot <- boot.ci(bootstrap_juos, type = "bca", index = i)
  #extract low and hi, add to ci_pimo df
  ci_juos$lo[i] <- ci_boot$bca[1,4]
  ci_juos$hi[i] <- ci_boot$bca[1,5]
  ci_juos$mean[i] <- mean(bootstrap_juos$t[, i])
}

saveRDS(ci_juos, file = "./model output/ci_juos.rds")

write.csv(ci_juos, "./model output/ci_juos.csv")

# plot(allEffects(juos_model))
#   plot(Effect(focal.predictors = c("ENN_dist", "Neighbor_larger"), mod = juos_model, partial.residuals = TRUE), ylim = c(-20, 0))
#   coefplot2::coefplot2(juos_model)
#   AICc(juos_model)
#   vif.mer(juos_model)
#   r.squaredGLMM(juos_int_enn)
#   
#   plot(resid(juos_int_enn) ~ predict(juos_int_enn))
#   abline(h = 0)
#   hist(resid(juos_int_enn))
#   mse <- sum(resid(juos_int_enn)^2)*(1 / length(resid(juos_int_enn)))

#####################################
# Below here is PIMO mortality stuff
## PIMO mortality

#full model
p_mort <- pmort_model_dredge
  ## scale the gradient by the hessian to get relative gradient, 
  # see Bolker comments at https://github.com/lme4/lme4/issues/120
  relgrad <- with(p_mort@optinfo$derivs,solve(Hessian,gradient))
  max(abs(relgrad))
  #looks good!

# plot(inv.logit(resid(p_mort)) ~ inv.logit(predict(p_mort)))


# bootstrapping with 5000 replications
#this takes a while (1.1 days on an i7 3.4 Ghz processor. Would be faster with multithreading)

bs <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- glmer(formula, data=d, nAGQ = 1, family = "binomial")
  return(fixef(fit)) 
} 
start.time <- Sys.time()
bootstrap_pmort <- boot(data=all_tree_vars[all_tree_vars$Spp == "PIMO", ], 
                       statistic=bs, 
                       R=5000, 
                       formula=Died ~ Neighbor_larger + ENN_dist + Height + Diam + BA_4m +
                         Avg_depth + AWC + cwd_normal_cum + fdsi_anom + (1|Cluster/site))
end.time <- Sys.time()
time.taken_pmort <- end.time - start.time
time.taken_pmort

saveRDS(bootstrap_pmort, file = "./model output/bootstrap_pmort.rds")

# view results
# # view results
# bootstrap_pmort
# summary(bootstrap_pmort)
# plot(bootstrap_pmort, index=1) # intercept x
# plot(bootstrap_pmort, index=2) #
# plot(bootstrap_pmort, index=3) #
# plot(bootstrap_pmort, index=4) #
# plot(bootstrap_pmort, index=5) #
# plot(bootstrap_pmort, index=6) #
# plot(bootstrap_pmort, index=7) #
# plot(bootstrap_pmort, index=8) #
# plot(bootstrap_pmort, index=9) #
# plot(bootstrap_pmort, index=10) #

# get 95% confidence intervals 
ci_pmort <- cbind(attr(bootstrap_pmort$t0, "names"), summary(bootstrap_pmort))
ci_pmort$lo <- NA
ci_pmort$hi <- NA

for( i in 1:nrow(ci_pmort)){
  ci_boot <- boot.ci(bootstrap_pmort, type = "bca", index = i)
  #extract low and hi, add to ci_pimo df
  ci_pmort$lo[i] <- ci_boot$bca[1,4]
  ci_pmort$hi[i] <- ci_boot$bca[1,5]
  ci_pmort$mean[i] <- mean(bootstrap_pmort$t[, i])
}

saveRDS(ci_pmort, file = "./model output/ci_pmort.rds")

write.csv(ci_pmort, "./model output/ci_pmort.csv")


##############################################
# Does PDC predict mortality?
# hist(all_tree_vars$PDC05)

# #only climate vars
# p_mort_clim <- glmer(Died ~
#                        Avg_depth + AWC + cwd_normal_cum + fdsi_anom + (1|Cluster/site), nAGQ = 1, family = "binomial", 
#                      data = all_tree_vars[all_tree_vars$Spp == "PIMO", ])
# summary(p_mort_clim)
# plot(allEffects(p_mort_clim))
# coefplot2::coefplot2(p_mort_clim)
# AICc(p_mort_clim)
# vif.mer(p_mort_clim)
# r.squaredGLMM(p_mort_clim)
# 
# #only stand structure vars
# p_mort_stand <- glmer(Died ~ ENN_dist*Neighbor_larger + BA_4m +
#                   (1|Cluster/site), nAGQ = 1, family = "binomial", 
#                 data = all_tree_vars[all_tree_vars$Spp == "PIMO", ])
# summary(p_mort_stand)
# plot(allEffects(p_mort_stand))
# coefplot2::coefplot2(p_mort_stand)
# AICc(p_mort_stand)
# vif.mer(p_mort_stand)
# r.squaredGLMM(p_mort_stand)
# 
# 
# 
# #only tree-level vars
# p_mort_tree <- glmer(Died ~ Height + BA_cm + 
#                         (1|Cluster/site), nAGQ = 1, family = "binomial", 
#                       data = all_tree_vars[all_tree_vars$Spp == "PIMO", ])
# summary(p_mort_tree)
# plot(allEffects(p_mort_tree))
# coefplot2::coefplot2(p_mort_tree)
# AICc(p_mort_tree)
# vif.mer(p_mort_tree)
# r.squaredGLMM(p_mort_tree)
# 
# 
# 
# #intercept-only model
# p_mort_int <- glm(Died ~ 1, family = "binomial",
#                   data = all_tree_vars[all_tree_vars$Spp == "PIMO", ])
# 
#random effects only model
# p_mort_int_mix <- glmer(Died ~ 1 + (1|Cluster/site), nAGQ = 1, family = "binomial",
#                         data = all_tree_vars[all_tree_vars$Spp == "PIMO", ])

# 
# pdc_mort <- glmer(Died ~ PDC05 + ENN_dist + Neighbor_larger + Height + BA_cm + Avg_depth + BA_4m +
#                     Avg_depth + AWC + cwd_normal_cum + fdsi_anom + (1|Cluster/site), nAGQ = 1, family = "binomial", 
#                   data = all_tree_vars[all_tree_vars$Spp == "PIMO", ])
# 
# 
# pdc_mort_poly <- glmer(Died ~ poly(PDC05, 2) + ENN_dist + Neighbor_larger + Height + BA_cm + Avg_depth + BA_4m +
#                     Avg_depth + AWC + cwd_normal_cum + fdsi_anom + (1|Cluster/site), nAGQ = 1, family = "binomial", 
#                     data = all_tree_vars[all_tree_vars$Spp == "PIMO", ])
# 
# pdc_mort_exp <- glmer(Died ~ exp(PDC05) + ENN_dist + Neighbor_larger + Height + BA_cm + Avg_depth + BA_4m +
#                     Avg_depth + AWC + cwd_normal_cum + fdsi_anom + (1|Cluster/site), nAGQ = 1, family = "binomial", 
#                     data = all_tree_vars[all_tree_vars$Spp == "PIMO", ])

# ## Model comparisons
# AICc(p_mort)
# AICc(p_mort_clim)
# AICc(p_mort_stand)
# AICc(p_mort_tree)
# AICc(p_mort_int)
# AICc(p_mort_int_mix)
# 
# AICc(pdc_mort)
# AICc(pdc_mort_poly)
# AICc(pdc_mort_exp)
# 
# r.squaredGLMM(p_mort)
# r.squaredGLMM(p_mort_clim)
# r.squaredGLMM(p_mort_stand)
# r.squaredGLMM(p_mort_tree)
# pR2(p_mort_int)
# r.squaredGLMM(p_mort_int_mix)
# 
# 
# r.squaredGLMM(pdc_mort)
# r.squaredGLMM(pdc_mort_poly)
# r.squaredGLMM(pdc_mort_exp)
# 
# # ROC plots for predicting mortality
# library("pROC")

# # hists and confusion matrices; not finished
# hist(fitted(p_mort)[all_tree_vars$Spp == "PIMO" & all_tree_vars$Died == "Y"])
# hist(fitted(p_mort)[all_tree_vars$Spp == "PIMO" & all_tree_vars$Died == "N"])
# hist(fitted(p_mort)[which(all_tree_vars$Spp == "PIMO" & all_tree_vars$Died == "Y")])
# hist(fitted(p_mort)[which(all_tree_vars$Spp == "PIMO" & all_tree_vars$Died == "N")])

######################################################
#Plot the ROCs, adding the AUC for each
######################################################
tiff(filename="./plots/roc_plots.tiff", 
    type="cairo",
    units="in", 
    width = 8, 
    height=8, 
    pointsize=15, 
    res=160)

par(mar = c(5.1, 4.1, 2.1, 2.1))

plot.roc(all_tree_vars[all_tree_vars$Spp == "PIMO", "Died"],fitted(p_mort_int),print.auc = TRUE, 
         col = "green", lty = 2)
plot.roc(all_tree_vars[all_tree_vars$Spp == "PIMO", "Died"],fitted(p_mort_int_mix),print.auc = TRUE,
         add = TRUE, col = "violet", lty = 4, print.auc.y = 0.45)
plot.roc(all_tree_vars[all_tree_vars$Spp == "PIMO", "Died"],fitted(p_mort),print.auc = TRUE,
         add = TRUE, col = "black", print.auc.y = 0.4)
plot.roc(all_tree_vars[all_tree_vars$Spp == "PIMO", "Died"],fitted(pdc_mort_poly),add = TRUE,
         print.auc = TRUE, lty = 5, col = "blue", print.auc.y = 0.35)
legend("bottomright", legend = c("Fixed Intercept", "Random Intercepts", "Full Model", "Full Model with PDC"),
       col = c("green", "violet", "black", "blue"), 
       lty = c(2, 4, 1, 5),
       lwd = 2,
       cex = 1)

dev.off()


#########################
## Classification accuracy
########################

#calculate classification accuracy for a given threshold
class_acc <- function(mod, thresh = .5){
  con <- data.frame(fitted(mod), as.factor(all_tree_vars$Died[all_tree_vars$Spp == "PIMO"]))
  nrow(all_tree_vars[all_tree_vars$Spp == "PIMO", ])
  true_pos <- nrow((con[con[, 1] > thresh & con[, 2] == "Y", ]))
  false_pos <- nrow((con[con[, 1] > thresh & con[, 2] == "N", ]))
  true_neg <- nrow((con[con[, 1] < thresh & con[, 2] == "N", ]))
  false_neg <- nrow((con[con[, 1] < thresh & con[, 2] == "Y", ]))
  accuracy <- (true_pos + true_neg) / length(fitted(mod))
  return(accuracy)
}


# 
# class_acc(p_mort, .5)
# class_acc(p_mort_clim, .5)
# class_acc(p_mort_stand, .5)
# class_acc(p_mort_tree, .5)
# class_acc(p_mort_int, .5)
# class_acc(p_mort_int_mix, .5)
# 
# class_acc(pdc_mort, .5)
# class_acc(pdc_mort_poly, .5)
# class_acc(pdc_mort_exp, .5)


ca_pmort <- NA
for(i in 1:100){
  ca_pmort[i] <- class_acc(p_mort, thresh = i / 100)
}
ca_pmort_clim <- NA
for(i in 1:100){
  ca_pmort_clim[i] <- class_acc(p_mort_clim, thresh = i / 100)
}
ca_pmort_pdc_poly <- NA
for(i in 1:100){
  ca_pmort_pdc_poly[i] <- class_acc(pdc_mort_poly, thresh = i / 100)
}
ca_pmort_int <- NA
for(i in 1:100){
  ca_pmort_int[i] <- class_acc(p_mort_int, thresh = i / 100)
}
ca_pmort_int_mix <- NA
for(i in 1:100){
  ca_pmort_int_mix[i] <- class_acc(p_mort_int_mix, thresh = i / 100)
}

## Plot of classification accuracy ~ threshold
tiff(filename="classification_accuracy_111916.tif", 
    type = "cairo",
    antialias = "gray",
    compression = "lzw",
    units="in", 
    width = 8, 
    height=7, 
    pointsize=15, 
    res=600)

par(mar = c(5.1, 4.1, 2.1, 2.1))

plot(NA,
     xlim = c(0, 100),
     ylim = c(0, 1),
     xlab = "Threshold",
     ylab = "Classification Accuracy",
     xaxt = 'n',
     yaxt = 'n',
     cex = 1.7)
lines(ca_pmort_int, lty = 2, lwd = 2, col = "green")
lines(ca_pmort_int_mix, lty = 4, lwd = 2, col = "violet")
lines(ca_pmort, lwd = 2, col = "black")
lines(ca_pmort_pdc_poly, lwd = 2, lty = 5, col = "blue")
axis(side = 1, at = c(0, 20, 40, 60, 80, 100), tick = TRUE,
     labels = c("0.00", "0.20", "0.40", "0.60", "0.80", "1.00"))
axis(side = 2, at = c(0, .2, .4, .6, .8, 1), tick = TRUE,
     labels = c("0.00", "0.20", "0.40", "0.60", "0.80", "1.00"))
legend("bottomright", legend = c("Fixed Intercept", "Random Intercepts", "Full Model", "Full Model with PDC"),
       col = c("green", "violet", "black", "blue"), 
       lty = c(2, 4, 1, 5),
       lwd = 2,
       cex = 1.2)

dev.off()


###############################################################
#### Generate figure 5
###############################################################

####partial effect plot of 2005 PDC on mortality risk 2005-2015

minx <- -mean(all_tree_vars_orig$PDC05)/sd(all_tree_vars_orig$PDC05)
maxx <- (100-mean(all_tree_vars_orig$PDC05))/sd(all_tree_vars_orig$PDC05)

#PDC05 effect from model
mort_eff <- Effect(pdc_mort_poly, focal.predictors = c("PDC05"), xlevels = list(PDC05 = seq(minx,maxx,length.out = 100)))

#prediction and 1-se envelope
x <- unlist(mort_eff$x)*sd(all_tree_vars_orig$PDC05) + mean(all_tree_vars_orig$PDC05) #unscale x variable
up <- exp(unlist(mort_eff$fit) + (mort_eff$se))
low <- exp(unlist(mort_eff$fit) - (mort_eff$se))
y <- exp(unlist(mort_eff$fit))

tiff(filename="pdc_mort_effects_102916.tif", 
    type = "cairo",
    antialias = "gray",
    compression = "lzw",
    units="in", 
    width = 3, 
    height=3, 
    pointsize=15, 
    res=600)

par(mar = c(2.5,2.5,0.2,0.2),
    oma = c(1,1,0.2,0.2),
    family = "serif")

plot(y ~ x, 
     ylim = c(0, 0.2),
     xlim 
     type = "l",
     bty = "n",
     xlab = "",
     ylab = "",
     xaxt = "n",
     yaxt = "n",
     lwd = 2)
lines(up ~ x)
lines(low ~ x)
polygon(c(x, rev(x)), c(up, rev(low)),
        col = addTrans("grey30",30), border = NA) #fills in the area between high and low confidence intervals
axis(1, cex.axis = 0.7, at= c(0, 50, 100))
axis(2, cex.axis = 0.7, at = c(0, 0.1, 0.2))
mtext(text = "Percent dead canopy in 2005", side = 1, cex = 0.8, line = 2.3)
mtext(text = "p(mortality), 2005-2015", side = 2, cex = 0.8, line = 2.3)


dev.off()


#------------------------------------------------------------
## Plots for bootstraps of variable values
#-------------------------------------------------
#Generate Figure 3
var_names <- as.character(ci_pimo[, 1])
order <- c(1, 4, 5, 3, 6, 2, 7, 8)
var_names[order] 
var_names[1] <- "Intercept"
var_names[2] <- "Neighbor Ratio"
var_names[3] <- "ENN Dist"
var_names[4] <- "Height"
var_names[5] <- "Diameter"
var_names[6] <- "4m BA"
var_names[7] <- "4m BA : CWD"
var_names[8] <- "ENN Dist : Neigh. Rat."

opar <- par(no.readonly = TRUE)

par(opar)

tiff(filename="coefficient_density_plots_082917.tif", 
    type = "cairo",
    antialias = "gray",
    compression = "lzw",
    units="in", 
    width = 6, 
    height=4, 
    pointsize=15, 
    res=600)

par(mar = c(1,1,1,0.6),
    oma = c(1,5,0.2,0),
    family = "serif")

layout(matrix(c(0,1,2,0,3,4,0,5,6,0,7,8,0,9,10,0,11,12,0,13,0,0,14,0,15,0,0), 9, 3, byrow = TRUE),
       widths = rep(c(.3,2,1), times = 9),
       heights = c(1.3,1,1,1,1,1,1,1))

for (i in 1:length(order)){
j <- ifelse(i==7, 11, 
            ifelse(i ==8, 11, order[i]))
d_p <- density(bootstrap_pimo_int$t[, j])
d_j <- density(bootstrap_juos$t[, j])
if(j<8){d_pm <- density(bootstrap_pmort$t[, j])}

if(i == 1){ #different axes for intercept
plot(NA,
     xlim = c(min(c(min(d_p$x), min(d_j$x))), max(c(max(d_p$x), max(d_j$x)))),
     ylim = c(0,max(c(max(d_p$y), max(d_j$y)))),
     xaxt = 'n', 
     yaxt = 'n',
     bty = 'n',
     xlab = "",
     ylab = "")
  lines(d_p, col = "blue")
  polygon(d_p, col = addTrans("grey50", 90))
  lines(d_j, col = "red", new = FALSE)
  polygon(d_j, col = addTrans("grey90", 60))
  axis(1)

  mtext(var_names[j],side=3,line=-1.8, 
        at=par("usr")[1] - 3,
        cex=.75,
        adj = 1)
  
  mtext("Canopy dieback", side = 3, cex = 0.9)
  
  plot(NA,
       xlim = c(min(d_pm$x), max(d_pm$x)),
       ylim = c(min(d_pm$y), max(d_pm$y)),
       xaxt = 'n', 
       yaxt = 'n',
       bty = 'n',
       xlab = "",
       ylab = "")
  lines(d_pm, col = "blue")
  polygon(d_pm, col = addTrans("grey50", 90))
  axis(1)
#  axis(2)
  
  mtext("Pinyon mortality", side = 3, cex = 0.9)
  
} 

if(i > 1 & i <7){ #all other variables
  
  par(mar = c(0,0,1,0))
  
  plot(NA,
       xlim = (c(-6,4)),
       ylim = c(0,max(c(max(d_p$y), max(d_j$y)))),
       xaxt = 'n', 
       yaxt = 'n',
       bty = 'n',
       xlab = "",
       ylab = "")
  lines(d_p, col = "black")
  polygon(d_p, col = addTrans("grey50", 90))
  lines(d_j, col = "black", new = FALSE)
  polygon(d_j, col = addTrans("grey90", 60))
  axis(1, labels = FALSE)
 # axis(2)
  abline(v = 0, lwd = 2, lty = 2)
  
  mtext(var_names[j],side=3,line=-1.5, 
        at=par("usr")[1] - 1,
        cex=.75,
        adj = 1)
  
  
  
  par(mar = c(0,0.6,1,0))
  
plot(NA,
       xlim = c(-1,2),
       ylim = c(min(d_pm$y), max(d_pm$y)),
       xaxt = 'n', 
       yaxt = 'n',
       bty = 'n',
     xlab = "",
     ylab = "")
  lines(d_pm, col = "blue")
  polygon(d_pm, col = addTrans("grey50", 90))
  axis(1, at = c(-1, 0, 1, 2), labels = FALSE)
 # axis(2)
  abline(v = 0, lwd = 2, lty = 2)
  
  if (i == 6){axis(1, at = c(-1, 0, 1, 2))}
}

if(i == 7){
  par(mar = c(0,0,1,0))
  
  plot(NA,
       xlim = (c(-6,4)),
       ylim = c(0,(max(d_p$y))),
       xaxt = 'n', 
       yaxt = 'n',
       bty = 'n',
       xlab = "",
       ylab = "")
  lines(d_p, col = "black")
  polygon(d_p, col = addTrans("grey50", 90))
  axis(1, labels = FALSE)
  abline(v = 0, lwd = 2, lty = 2)
  
  mtext(var_names[7],side=3,line=-1.5, 
        at=par("usr")[1] - 1,
        cex=.75,
        adj = 1)
  
}

if(i == 8){
  par(mar = c(0,0,1,0))
  
  plot(NA,
       xlim = (c(-6,4)),
       ylim = c(0,(max(d_j$y))),
       xaxt = 'n', 
       yaxt = 'n',
       bty = 'n',
       xlab = "",
       ylab = "")
  lines(d_j, col = "black")
  polygon(d_j, col = addTrans("grey90", 60))
  axis(1, labels = FALSE)
  abline(v = 0, lwd = 2, lty = 2)
  
  mtext("ENN Dist :\nNeigh. Rat.",side=3,line=-2.5, 
        at=par("usr")[1] - 1,
        cex=.75,
        adj = 1)
  
  axis(1)
}

}




mtext(text = "Bootstrap coefficient estimates (scaled)", side = 1,
      line = 2.3, at = 0, cex = 0.9)

dev.off()


#####################################################################
### Generate figure 4
### multipanel figure
#####################################################################
## this is the ugliest, WETtest piece of code you've ever seen but it should work.
## run this whole chunk until dev.off() way down there

tiff(filename="test_figure.tif",
     type="cairo",
     units="in",
     width = 6,
     height = 6,
     pointsize=15,
     res=600,
     compression = "lzw")



layout(matrix(c(1,2,3,4,5,6,7,8,9), nrow=3, ncol=3, byrow = TRUE))
par(mar = c(3,.5,1.2,0), oma = c(2,4,1,1), family = "serif", bty = 'n', 
    cex = 0.7, cex.lab = 0.7, cex.axis = 0.7)


# function to calculate standard errors for a given vcov matrix and row of x values, to get
# point estimates of prediction error for a given prediction. See http://www.ats.ucla.edu/stat/r/faq/deltamethod.htm
# This gets used for the confidence interval of the effects plot

calc.se <- function(x){
  sqrt(as.matrix(x) %*% vcov(model) %*% t(x))
} 

##############################################
## Creating new predictions
##############################################

# Calculate residuals from predicted values and observed values
# resid <- residuals(model)

# for unscaling variables -- these things come from a full model in the model set, you'll have to change
# the number or the path depending upon what kind of model object and the order of your models, etc.
# You can also just get them from the raw data.

alldata_scaled <- read.csv("pj_resample_plot_level_data_scaled.csv")

alldata_unscaled <- read.csv( "pj_resample_plot_level_data_unscaled.csv")

avg_BA_4m_mean <-  mean(alldata_unscaled$avg_BA_4m)
avg_BA_4m_sd <- sd(alldata_unscaled$avg_BA_4m)

vpd_peak_anom_mean <- mean(alldata_unscaled$vpd_peak_anom)
vpd_peak_anom_sd <- sd(alldata_unscaled$vpd_peak_anom)

BA_4m_mean <- mean(log(all_tree_vars_orig[all_tree_vars_orig$Spp == "PIMO", "BA_4m"]+0.01))
BA_4m_sd <- sd(log(all_tree_vars_orig[all_tree_vars_orig$Spp == "PIMO", "BA_4m"]+0.01))

ENN_dist_mean <- mean(log(all_tree_vars_orig[all_tree_vars_orig$Spp == "JUOS", "ENN_dist"]+0.01))
ENN_dist_sd <- sd(log(all_tree_vars_orig[all_tree_vars_orig$Spp == "JUOS", "ENN_dist"]+0.01))


################################
#generates figure 4a
################################
pardefault <- par(no.readonly = T)
source('addTrans.R', echo=FALSE)

# setwd("")

### set some parameters to be use throughout
ymax <- 50
ymin <- -500

# model used for effects plot
model <- lmer(Weighted_pdc ~ Avg_depth + avg_BA_4m + avg_summer_tmax + cwd_normal_cum + 
                vpd_peak_anom + avg_BA + avg_ENN + Pct_pimo + avg_BA_4m *vpd_peak_anom + (1|Cluster), data = alldata)


ranefs <- vector(mode = "numeric", length = 98) #initialize a vector for the ranefs

#pull out the random intercept for each plot
for (i in (1:98)){
  ranefs[i] <- ranef(model)$Cluster[rownames(ranef(model)$Cluster) == model@frame$Cluster[i], ]
}

calculate_preds_avg_BA_4m <- function(quant){
  C <- data.frame(
    "(Intercept)" = 1,  
    Avg_depth = 0,
    avg_BA_4m = seq(min(alldata$avg_BA_4m), 
                    max(alldata$avg_BA_4m), 
                    length.out=nrow(alldata)),
    avg_summer_tmax = 0,
    cwd_normal_cum=0,
    vpd_peak_anom = unname(quantile(alldata$vpd_peak_anom, quant)),
    avg_BA = 0,
    avg_ENN=0,
    Pct_pimo = 0,
    "avg_BA_4m:vpd_peak_anom" = seq(unname(quantile(alldata$vpd_peak_anom, quant))*min(alldata$avg_BA_4m), 
                                    unname(quantile(alldata$vpd_peak_anom, quant))*max(alldata$avg_BA_4m), 
                                    length.out=nrow(alldata))
  )
  
  se <- vector(mode='numeric', length=nrow(C))
  
  for (i in 1:nrow(C)){
    se[i] <- as.numeric(calc.se(C[i, ]))
  }
  
  preds <- data.frame(
    predictions = as.matrix(C) %*% as.matrix(fixef(model)),
    lo = as.matrix(C) %*% as.matrix(fixef(model)) - 1.96*se,
    hi = as.matrix(C) %*% as.matrix(fixef(model)) + 1.96*se,
    avg_BA_4m_orig = model@frame$avg_BA_4m,
    avg_BA_4m_sim = C$avg_BA_4m #,
    # part_resids = part_resids
  )
  
  return(preds)
  
}


for(i in 1:3){
  values = c(0.1,0.5,0.9)
  # Use the function to get all the stuff we need (predictions, partial residuals, and SEs)
  pred <- calculate_preds_avg_BA_4m(values[i])
  
  #set up the empty plot window
  plot(NA,
       ylim = c(-400, 100),
       xlim = c(0, max(I(((pred$avg_BA_4m_sim * avg_BA_4m_sd) + avg_BA_4m_mean)/10000))),
       xlab = "",
       ylab = "",
       # log = "x",
       cex.lab = 1.5,
       xaxt = "n",
       yaxt = "n"
  )
  
  #draw prediction line and confidence envelope of predictions
  lines(pred$predictions ~ I(((pred$avg_BA_4m_sim * avg_BA_4m_sd) + avg_BA_4m_mean)/10000), lwd = 2)
  lines((pred$lo) ~ I(((pred$avg_BA_4m_sim * avg_BA_4m_sd) + avg_BA_4m_mean)/10000), lwd = 1.3) 
  lines((pred$hi) ~ I(((pred$avg_BA_4m_sim * avg_BA_4m_sd) + avg_BA_4m_mean)/10000), lwd = 1.3)
  
  #fill in the confidence envelope
  polygon(c(I(((pred$avg_BA_4m_sim * avg_BA_4m_sd) + avg_BA_4m_mean)/10000), rev(I(((pred$avg_BA_4m_sim * avg_BA_4m_sd) + avg_BA_4m_mean)/10000))), c((pred$hi), rev((pred$lo))),
          col = addTrans("grey",30), border = NA)
  
  #label the levels of the interacting variable
  text(x = 0.2, y = -0, labels = paste0(values[i]*100, "% Peak VPD \n Anomaly = ", 
                                        round(((unname(quantile(alldata$vpd_peak_anom, values[i])) * vpd_peak_anom_sd + vpd_peak_anom_mean)), 3)),
       cex = 0.7)
  
  axis(1, at = c(0, .1, .2, .3)) #draw the x-axis
  
  #draw the y-axis just for the leftmost panel
  
  if(i == 1){
    mtext(side = 3, at = 0, text= "(a)")
    axis(2)
    mtext(side = 2, text = "Predicted Change in \nPlot Weighted Canopy", cex = 0.7, outer = FALSE, line = 1.7)}
  
  #label x axis
  if(i == 2){mtext(side = 1, text = expression(paste("Avg 4m BA (", m^2, ")")), cex= 0.7, outer = FALSE, line = 2)}
  
}


#label y axis





#############################################################
### Pimo part -- figure 4b
#############################################################

model <- pimo_int_ba

ymax <- 20
ymin <- -40

num_tree <- (nrow(all_tree_vars[all_tree_vars$Spp == "PIMO", ]))
ranefs <- vector(mode = "numeric", length = num_tree) 

###BA_4meter
calculate_preds_BA_4m <- function(quant){
  C <- data.frame(
    '(Intercept)' = 1,
    Neighbor_larger = 0,
    ENN_dist = 0,
    Height = 0,
    Diam = 0,
    BA_4m = seq(min(all_tree_vars[all_tree_vars$Spp == "PIMO", ]$BA_4m), 
                max(all_tree_vars[all_tree_vars$Spp == "PIMO", ]$BA_4m), 
                length.out=nrow(all_tree_vars[all_tree_vars$Spp == "PIMO", ])),
    Avg_depth = 0,
    cwd_normal_cum = unname(quantile(all_tree_vars$cwd_normal_cum, quant)),
    vpd_normal_annual_max = 0,
    vpd_peak_anom = 0,
    'BA_4m:cwd_normal_cum' = seq(unname(quantile(all_tree_vars$cwd_normal_cum, quant))*
                                   min(all_tree_vars[all_tree_vars$Spp == "PIMO", ]$BA_4m), 
                                 unname(quantile(all_tree_vars$cwd_normal_cum, quant))*
                                   max(all_tree_vars[all_tree_vars$Spp == "PIMO", ]$BA_4m), 
                                 length.out=nrow(all_tree_vars[all_tree_vars$Spp == "PIMO", ]))
  )
  
  se <- vector(mode='numeric', length=nrow(C))
  
  #get prediction ses from above function
  for (i in 1:nrow(C)){
    se[i] <- as.numeric(calc.se(C[i, ]))
  }
  
  preds <- data.frame(
    predictions = as.matrix(C) %*% as.matrix(fixef(model)),
    lo = as.matrix(C) %*% as.matrix(fixef(model)) - 1.96*se,
    hi = as.matrix(C) %*% as.matrix(fixef(model)) + 1.96*se,
    BA_4m = model@frame$`BA_4m`,
    BA_4m_sim = C$BA_4m
  )
  
  return(preds)
  
}

# setwd("E:\\IALE analysis\\plots") #where the plots go
# 

for(i in 1:3){
  values = c(0.1,0.5,0.9)
  # Use the function to get all the stuff we need (predictions, partial residuals, and SEs)
  pred <- calculate_preds_BA_4m(values[i])
  
  
  plot(NA,
       ylim = c(-50, 0),
       xlim = c(min(I(exp((pred$BA_4m_sim * BA_4m_sd) + BA_4m_mean))/10000), max(I(exp((pred$BA_4m_sim * BA_4m_sd) + BA_4m_mean))/10000)),
       xlab = "",
       ylab = "",
       # log = "x",
       cex.lab = 1.5,
       xaxt = "n",
       yaxt = "n"
  )
  
  lines(pred$predictions ~ I(exp((pred$BA_4m_sim * BA_4m_sd) + BA_4m_mean)/10000), lwd = 2)
  lines((pred$lo) ~ I(exp((pred$BA_4m_sim * BA_4m_sd) + BA_4m_mean)/10000), lwd = 1.3) 
  lines((pred$hi) ~ I(exp((pred$BA_4m_sim * BA_4m_sd) + BA_4m_mean)/10000), lwd = 1.3)
  polygon(c(I(exp((pred$BA_4m_sim * BA_4m_sd) + BA_4m_mean)/10000), rev(I(exp((pred$BA_4m_sim * BA_4m_sd) + BA_4m_mean))/10000)), c((pred$hi), rev((pred$lo))),
          col = addTrans("grey",30), border = NA) #fills in the area between high and low confidence intervals
  
  
  text(x = 0.7, y = -3, labels = paste0(values[i]*100, "% CWD = \n", 
                                        round(((unname(quantile(all_tree_vars$cwd_normal_cum, values[i])) * cwd_normal_cum_sd + cwd_normal_cum_mean)), 0), " mm"),
       cex = 0.7)
  
  axis(1)
  
  if(i == 1){
    mtext(side = 3, at = 0, text= "(b)")
    axis(2)
    mtext(side = 2, text = "Predicted Change in Crown", cex = 0.7, outer = FALSE, line = 1.7)}
  
  #label x axis
  if(i == 2){mtext(side = 1, text = expression(paste("BA, 4m buffer (", m^2, ")")), cex= 0.7, outer = FALSE, line = 2)}
  
}


#######################################################
### juos part -- figure 4c
#######################################################


model <- juos_int_enn

num_tree <- (nrow(all_tree_vars[all_tree_vars$Spp == "PIMO", ]))
ranefs <- vector(mode = "numeric", length = num_tree) 

###BA_4meter
calculate_preds_Neighbor_larger <- function(quant){
  C <- data.frame(
    '(Intercept)' = 1,
    Neighbor_larger = seq(min(all_tree_vars[all_tree_vars$Spp == "JUOS", ]$Neighbor_larger), 
                          max(all_tree_vars[all_tree_vars$Spp == "JUOS", ]$Neighbor_larger), 
                          length.out=nrow(all_tree_vars[all_tree_vars$Spp == "JUOS", ])),
    ENN_dist = unname(quantile(all_tree_vars$ENN_dist, quant)),
    Height = 0,
    Diam = 0,
    BA_4m = 0,
    fdsi_anom = 0,
    avg_summer_tmax = 0,
    Avg_depth = 0,
    AWC = 0,
    'Neighbor_larger:ENN_dist' = seq(unname(quantile(all_tree_vars$ENN_dist, quant))*
                                       min(all_tree_vars[all_tree_vars$Spp == "JUOS", ]$Neighbor_larger), 
                                     unname(quantile(all_tree_vars$ENN_dist, quant))*
                                       max(all_tree_vars[all_tree_vars$Spp == "JUOS", ]$Neighbor_larger), 
                                     length.out=nrow(all_tree_vars[all_tree_vars$Spp == "JUOS", ]))
  )
  
  se <- vector(mode='numeric', length=nrow(C))
  
  #get prediction ses from above function
  for (i in 1:nrow(C)){
    se[i] <- as.numeric(calc.se(C[i, ]))
  }
  
  preds <- data.frame(
    predictions = as.matrix(C) %*% as.matrix(fixef(model)),
    lo = as.matrix(C) %*% as.matrix(fixef(model)) - 1.96*se,
    hi = as.matrix(C) %*% as.matrix(fixef(model)) + 1.96*se,
    Neighbor_larger = model@frame$`Neighbor_larger`,
    Neighbor_larger_sim = C$Neighbor_larger
  )
  
  return(preds)
  
}


for(i in 1:3){
  values = c(0.1,0.5,0.9)
  # Use the function to get all the stuff we need (predictions, partial residuals, and SEs)
  pred <- calculate_preds_Neighbor_larger(values[i])
  
  
  plot(NA,
       ylim = c(-25, 0),
       xlim = c(min(I(exp((pred$Neighbor_larger * Neighbor_larger_sd) + Neighbor_larger_mean))), max(I(exp((pred$Neighbor_larger_sim * Neighbor_larger_sd) + Neighbor_larger_mean)))),
       xlab = "",
       ylab = "",
       log = "x",
       cex.lab = 1.5,
       xaxt = "n",
       yaxt = "n"
  )
  
  lines(pred$predictions ~ I(exp((pred$Neighbor_larger_sim * Neighbor_larger_sd) + Neighbor_larger_mean)), lwd = 2)
  lines((pred$lo) ~ I(exp((pred$Neighbor_larger_sim * Neighbor_larger_sd) + Neighbor_larger_mean)), lwd = 1.3) 
  lines((pred$hi) ~ I(exp((pred$Neighbor_larger_sim * Neighbor_larger_sd) + Neighbor_larger_mean)), lwd = 1.3)
  polygon(c(I(exp((pred$Neighbor_larger_sim * Neighbor_larger_sd) + Neighbor_larger_mean)), rev(I(exp((pred$Neighbor_larger_sim * Neighbor_larger_sd) + Neighbor_larger_mean)))), c((pred$hi), rev((pred$lo))),
          col = addTrans("grey",30), border = NA) #fills in the area between high and low confidence intervals
  
  ticks <- c(-2, 0, 2)
  labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
  axis(1, at=c(0.01, 1, 100), labels=labels)
  
  text(x = 0.9, y = -3, labels = paste0(values[i]*100, "% ENN Distance = \n", 
                                        round(exp((unname(quantile(all_tree_vars$ENN_dist, values[i])) * ENN_dist_sd + ENN_dist_mean)), 1), " m"),
       cex = 0.7)
  
  
  if(i == 1){
    mtext(side = 3, at = 0.003, text= "(c)")
    axis(2)
    mtext(side = 2, text = "Predicted Change in Crown", cex = 0.7, outer = FALSE, line = 1.7)}
  
  #label x axis
  if(i == 2){mtext(side = 1, text = "Neighbor Size Ratio (log)", cex= 0.7, outer = FALSE, line = 2)}
  
}

dev.off()


