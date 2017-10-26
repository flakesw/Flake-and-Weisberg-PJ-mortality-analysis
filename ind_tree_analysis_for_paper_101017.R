#### Analysis of field data
# Sam Flake, 12 Sep 2017
# sflake@gmail.com or swflake@ncsu.edu


# load libraries
library(plyr)
library(boot)
library(lme4)
library(MuMIn)
library(car)
library(effects)
library(pscl)
library("pROC")

## set up some options
set.seed(665224291)
options(na.action = na.fail)
pardefault <- par(no.readonly = T)

#Source some useful functions
source("pj_mortality_functions_090617.R")

#read prepared data
all_tree_vars <- read.csv("./Clean data/all_tree_Vars_scaled_and_transformed.csv")
all_tree_vars_orig <- read.csv("./Clean data/all_tree_Vars_unscaled_untrans.csv")

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

saveRDS(pimo_dredge, "./model output/pimo_dredge_outputs.rds")
#get call for top model and evaluate it
best_pimo_call <- as.character(unlist(attr((pimo_dredge)[1], "model.calls")))
pimo_model_dredge <- eval(parse(text = best_pimo_call))
summary(pimo_model_dredge)

pimo_dredge[1:10] #top 10 models

saveRDS(pimo_model_dredge, "./model output/pimo_model_dredge.rds")

# #----------------------------------------------
# 

#global model with tree-level vars and abiotic vars
juos_global <- lmer(formula = Delta_pdc ~ Neighbor_larger_log + ENN_dist_log + Diam_log + Height + BA_4m_log +  
                      cwd_normal_cum + fdsi_anom +
                      cwd_peak_anom + tmean + ppt +
                      vpd_normal_annual_max + vpd_peak_anom + avg_summer_tmax +
                      peak_tmax + min_ann_ppt +
                      Pndjfm + Avg_depth + AWC +
                      (1|Cluster/site), 
                    data = all_tree_vars[all_tree_vars$Spp == "JUOS", ],
                    na.action = na.fail)

#only use combos that aren't pairwise correlated
cormat <- abs(cor(all_tree_vars[all_tree_vars$Spp == "JUOS", c('Neighbor_larger_log', 'ENN_dist_log', 'Height', 'Diam_log', 'BA_4m_log', 
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
                      fixed = c('Neighbor_larger_log', 'ENN_dist_log', 'Height', 'Diam_log', 'BA_4m_log'),
                      subset = cormat, 
                      m.lim = c(5, 8), 
                      trace = 2)

saveRDS(juos_dredge, "./model output/juos_dredge_outputs.rds")


#get call for top model and evaluate it
best_juos_call <- as.character(unlist(attr((juos_dredge)[1], "model.calls")))
juos_model_dredge <- eval(parse(text = best_juos_call))

saveRDS(juos_model_dredge, "./model output/juos_model_dredge.rds")

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

saveRDS(pmort_dredge, "./model output/pmort_dredge_outputs.rds")
end.time <- Sys.time()

time.taken_pmort_dredge <- end.time - start.time
time.taken_pmort_dredge #30 minutes
(pmort_dredge)[c(1:10)]

#get call for top model and evaluate it
best_pmort_call <- as.character(unlist(attr((pmort_dredge)[1], "model.calls")))
pmort_model_dredge <- eval(parse(text = best_pmort_call))

saveRDS(pmort_model_dredge, "./model output/pmort_model_dredge.rds")

#############################################################################################
#Tree-level models
#############################################################################################

#PIMO dieback
pimo_model_dredge <- readRDS("./model output/pimo_model.rds")
summary(pimo_model)
#add interaction term to the PIMO dieback model
pimo_model <- update(pimo_model_dredge, . ~ . - BA_4m_log:cwd_normal_cum)
summary(pimo_model)


pimo_int_ba_cwd <- update(pimo_model, . ~ . + cwd_normal_cum:BA_4m_log)
pimo_int_ba_vpd <- update(pimo_model, . ~ . + BA_4m_log:vpd_normal_annual_max)
pimo_int_ba_depth <- update(pimo_model, . ~ . + BA_4m_log:Avg_depth)
pimo_int_enn_nl <- update(pimo_model, . ~ . + ENN_dist_log:Neighbor_larger_log)
AICc(pimo_model)
#only ba:cwd interaction decreases AICc; the others increase it
AICc(pimo_int_ba_cwd)
AICc(pimo_int_ba_vpd)
AICc(pimo_int_ba_depth)
AICc(pimo_int_enn_nl)

plot(allEffects(pimo_int_ba_cwd))

summary(pimo_int_ba)
saveRDS(pimo_int_ba, "./model output/pimo_model.rds")

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

saveRDS(bootstrap_pimo_int, file = paste0("./model output/bootstrap_pimo_int.rds"))

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

# add interaction to the selected model
juos_int_enn <- update(juos_model_dredge, . ~ . + ENN_dist_log:Neighbor_larger_log)

summary(juos_int_enn)
saveRDS(juos_int_enn, "./model output/juos_model.rds")


#nonparametric bootstrap cis
bs <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- lmer(formula, data=d)
  return(fixef(fit)) 
} 


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

summary(p_mort)

saveRDS(p_mort, "./model output/pimo_mort_model.rds")

  ## scale the gradient by the hessian to get relative gradient, 
  # see Bolker comments at https://github.com/lme4/lme4/issues/120
  relgrad <- with(p_mort@optinfo$derivs,solve(Hessian,gradient))
  max(abs(relgrad))
  #looks good!

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
                       formula = p_mort@call)
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


