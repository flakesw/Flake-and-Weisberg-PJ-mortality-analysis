#### Analysis of field data
library(plyr)
library(lme4)
library(MuMIn)
library(boot)
library(car)
source("pj_mortality_functions_090617.R")
options(na.action = na.fail)


plot_data <- read.csv("./Clean data/plot_level_vars.csv")
plot_data_scaled <- read.csv("./Clean data/plot_level_vars_scaled.csv")

#---------------------------------------------------------------------------------
####################### Plot analysis
###### Select climate variables
######################################

#global model with all abiotic variables
plot_climate_global <- lmer(Avg_delta_pdc ~ avg_height + 
                              avg_Diam + avg_ENN + BA.Plot + Pct_pimo + 
                              cwd_normal_cum + fdsi_anom + 
                              cwd_peak_anom + tmean + ppt + 
                              vpd_normal_annual_max + vpd_peak_anom + avg_summer_tmax + 
                              peak_tmax + min_ann_ppt +  
                              Pndjfm + Avg_depth + AWC + (1|Cluster), data = plot_data_scaled, na.action = na.fail)


#correlation matrix to select only variables pairwise correlated with r <=0.5
cormat <- abs(cor(plot_data_scaled[, c('avg_height', 
                                       'avg_Diam', 'avg_ENN', 'BA.Plot', 'Pct_pimo',
                                       'cwd_normal_cum', 'fdsi_anom',
                                       'cwd_peak_anom' , 'tmean', 'ppt',
                                       'vpd_normal_annual_max' , 'vpd_peak_anom' , 'avg_summer_tmax', 'peak_tmax', 
                                       'min_ann_ppt',  'Pndjfm' , 'Avg_depth' , 'AWC')])) <=.5

#use all stand structure vars
cormat[1:5,1:5] <- TRUE

#use only one side of matrix (avoid duplicated models)
cormat[!lower.tri(cormat)] <- NA

# build models with all combinations of abiotic variables, with a maximum of 4 variables per model
# this is very quick (a few seconds)

plot_dredge <- dredge(plot_climate_global, 
                      fixed = c('avg_height','avg_Diam', 'avg_ENN',
                                'BA.Plot', 'Pct_pimo'),
                      subset = cormat, 
                      m.lim = c(5, 8), 
                      trace = 2)


saveRDS(plot_dredge, "./model output/plot_dredge_outputs.rds")

best_plot_call <- as.character(unlist(attr((plot_dredge)[1], "model.calls")))
plot_model_dredge <- eval(parse(text = best_plot_call))
summary(plot_model_dredge)
r.squaredGLMM(plot_model_dredge)
summary(plot_model_dredge)
AICc(plot_model_dredge)

plot_dredge[1:10] #top 10 models
write.csv(plot_dredge[1:10], "./model output/plot_level_mod_selection_table.csv")


#validation and fitting interaction terms

plot(resid(plot_model_dredge)~predict(plot_model_dredge)) #looks good!
abline(h=0)

summary(plot_model_dredge)
vif.mer(plot_model_dredge)
r.squaredGLMM(plot_model_dredge)

plot_model_int_depth <- update(plot_model_dredge, . ~ . + BA.Plot:Avg_depth)

saveRDS(plot_model_dredge, "./model output/plot_model.rds")

plot_model_int_cwd <- update(plot_model_dredge, . ~ . + BA.Plot:cwd_normal_cum)
plot_model_int_tmax <- update(plot_model_dredge, . ~ . + BA.Plot:avg_summer_tmax)

AICc(plot_model_dredge)
AICc(plot_model_int_depth)
AICc(plot_model_int_cwd)
AICc(plot_model_int_tmax)

plot(allEffects(plot_model_int_depth))
plot(allEffects(plot_model_int_cwd))
plot(allEffects(plot_model_int_tmax))

summary(plot_model_int_depth)
summary(plot_model_int_cwd)
summary(plot_model_int_tmax)

#------------------------------------------------------------------
########## Bootstrapping
library("boot")

#function to use for bootstrapping, returns the fixed effects from each model run
bs <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- lmer(formula, data=d)
  return(fixef(fit)) 
} 


# bootstrapping with 5000 replications
#this only takes a minute or so

#calls bs function with the given formula, returns all the fixed effects
bootstrap_plot_level <- boot(data= plot_data_scaled, 
                             statistic=bs, 
                             R=5000, 
                             formula=plot_model_dredge@call,
                             parallel = "no")


saveRDS(bootstrap_plot_level, "./model output/bootstrap_plot.rds")

bootstrap_plot_level
# summary(bootstrap_plot_level)
# plot(bootstrap_plot_level, index=1) # intercept
# plot(bootstrap_plot_level, index=2) #
# plot(bootstrap_plot_level, index=3) #

# get 95% confidence intervals 

#initialize the data frame to receive the CI inforation
ci_plot <- cbind(attr(bootstrap_plot_level$t0, "names"), summary(bootstrap_plot_level))
ci_plot$lo <- NA
ci_plot$hi <- NA

#calculate CI for each parameter
for( i in 1:nrow(ci_plot)){
  ci_boot <- boot.ci(bootstrap_plot_level, type = "bca", index = i)
  #extract low and hi, add to ci_plot df
  ci_plot$lo[i] <- ci_boot$bca[1,4]
  ci_plot$hi[i] <- ci_boot$bca[1,5]
  ci_plot$mean[i] <- mean(bootstrap_plot_level$t[, i])
}

saveRDS(ci_plot, "./model output/ci_plot.rds")
write.csv(ci_plot, "./model output/ci_plot.csv")

