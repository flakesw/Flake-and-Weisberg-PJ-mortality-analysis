#### Analysis of field data
library(plyr)
library(reshape2)
library(effects)
library(car)
library(Hmisc)
library(arm)
library(vegan)
library(lme4)
library(MuMIn)
library(coefplot2)
library(boot)
source("calculate_awc.R")
options(na.action = na.fail)


vif.mer <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)] }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}


alldata_scaled <- read.csv("pj_resample_plot_level_data_scaled.csv")

alldata_unscaled <- read.csv("pj_resample_plot_level_data_scaled.csv")

#---------------------------------------------------------------------------------
####################### Plot analysis
###### Select climate variables
######################################

#global model with all abiotic variables
plot_climate_global <- lmer(Avg_delta_pdc ~ cwd_normal_cum + fdsi_anom + 
                              cwd_peak_anom + tmean + ppt + 
                              vpd_normal_annual_max + vpd_peak_anom + avg_summer_tmax + 
                              peak_tmax + min_ann_ppt +  
                              Pndjfm + Avg_depth + AWC + (1|Cluster), data = alldata, na.action = na.fail)


#correlation matrix to select only variables pairwise correlated with r <=0.5
cormat <- abs(cor(alldata[, c('cwd_normal_cum', 'fdsi_anom',
                              'cwd_peak_anom' , 'tmean', 'ppt',
                              'vpd_normal_annual_max' , 'vpd_peak_anom' , 'avg_summer_tmax', 'peak_tmax', 
                              'min_ann_ppt',  'Pndjfm' , 'Avg_depth' , 'AWC')])) <=.5
cormat[!lower.tri(cormat)] <- NA

# build models with all combinations of abiotic variables, with a maximum of 4 variables per model
# this take a while!

plot_climate_set3 <- dredge(plot_climate_global, subset = cormat, m.lim = c(0, 4))
(plot_climate_set3)[c(1:10)]

## top 10 models selected
mod1 <- lmer(formula = Weighted_pdc ~ Avg_depth + avg_summer_tmax + cwd_normal_cum + 
               vpd_peak_anom + (1 | Cluster), data = alldata)
mod2 <- lmer(formula = Weighted_pdc ~ Avg_depth + avg_summer_tmax + 
               cwd_normal_cum + (1 | Cluster), data = alldata)
mod3 <- lmer(formula = Weighted_pdc ~ Avg_depth + + avg_summer_tmax + AWC + cwd_normal_cum  
             + (1 | Cluster), data = alldata)
mod4 <- lmer(formula = Weighted_pdc ~ Avg_depth + cwd_normal_cum + 
               vpd_normal_annual_max + vpd_peak_anom + (1 | Cluster), data = alldata)
mod5 <- lmer(formula = Weighted_pdc ~ Avg_depth + cwd_normal_cum + 
               peak_tmax + vpd_normal_annual_max  + (1 | Cluster), data = alldata)
mod6 <- lmer(formula = Weighted_pdc ~ Avg_depth + cwd_normal_cum + 
               vpd_normal_annual_max + (1 | Cluster), data = alldata, na.action = na.fail)
mod7 <- lmer(formula = Weighted_pdc ~ Avg_depth + AWC + cwd_normal_cum + 
               vpd_normal_annual_max + (1 | Cluster), data = alldata)
mod8 <- lmer(formula = Weighted_pdc ~ Avg_depth + cwd_normal_cum + peak_tmax + vpd_peak_anom +
               (1 | Cluster), data = alldata)
mod9 <- lmer(formula = Weighted_pdc ~ Avg_depth + cwd_normal_cum + vpd_peak_anom +
               (1 | Cluster), data = alldata)
mod10 <- lmer(formula = Weighted_pdc ~ Avg_depth + cwd_normal_cum + peak_tmax +  
                (1 | Cluster), data = alldata)
max(vif.mer(mod10))


alldata$Cluster <- (as.character((alldata$Cluster)))
#### analysis with stand structure vars

#add in the stand structure variables to abiotic model
plot_model <- lmer(Weighted_pdc ~ Avg_depth + avg_summer_tmax + cwd_normal_cum + avg_height + 
                     vpd_peak_anom + avg_Diam + avg_ENN + avg_BA_4m  + Pct_pimo + Pct_pimo*vpd_peak_anom + (1|Cluster), data = alldata)

# summary(plot_model)
# AICc(plot_model)
# 
# plot(resid(plot_model)~predict(plot_model))
# abline(h=0)
# 
# summary(plot_model)
# vif.mer(plot_model)
# r.squaredGLMM(plot_model)
# 
# plot(allEffects(plot_model_fix))
# plot(Effect(focal.predictors = c("avg_BA_4m", "vpd_peak_anom"), mod = plot_model_fix, 
#             ylim = c(-40, 0), partial.residuals = TRUE))

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
#this takes a while

#calls bs function with the given formula, returns all the fixed effects
bootstrap_plot_level <- boot(data=alldata, 
                           statistic=bs, 
                           R=5000, 
                           formula=Weighted_pdc ~ Avg_depth + avg_BA_4m + avg_summer_tmax + cwd_normal_cum + avg_height +
                             vpd_peak_anom + avg_Diam + avg_ENN   + Pct_pimo + avg_BA_4m * vpd_peak_anom + (1|Cluster))


# 
# bootstrap_plot_level
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


#--------------------------------------------------
## Effects plot
#------------------------------------------------

#generates figure 4a

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

climvars <- all_clim[all_clim$site %in% plot_vars$site, ]

# for unscaling variables -- these things come from a full model in the model set, you'll have to change
# the number or the path depending upon what kind of model object and the order of your models, etc.
# You can also just get them from the raw data.
avg_BA_mean <- mean(alldata_unscaled$avg_BA)
avg_BA_sd <- sd(alldata_unscaled$avg_BA)

avg_ENN_mean <-  mean(alldata_unscaled$avg_ENN)
avg_ENN_sd <- sd(alldata_unscaled$avg_ENN)

avg_BA_4m_mean <-  mean(alldata_unscaled$avg_BA_4m)
avg_BA_4m_sd <- sd(alldata_unscaled$avg_BA_4m)

Pct_pimo_mean <-  mean(plot_vars_back$Pct_pimo)
Pct_pimo_sd <- sd(plot_vars_back$Pct_pimo)

Avg_depth_mean <-  mean(plot_vars_back$Avg_depth)
Avg_depth_sd <- sd(plot_vars_back$Avg_depth)

avg_summer_tmax_mean <-  mean(climvars$avg_summer_tmax)
avg_summer_tmax_sd <- sd(climvars$avg_summer_tmax)

cwd_normal_cum_mean <-  mean(climvars$cwd_normal_cum)
cwd_normal_cum_sd <- sd(climvars$cwd_normal_cum)

AWC_mean <-  mean(soils_back$AWC)
AWC_sd <- sd(soils_back$AWC)

vpd_peak_anom_mean <-  mean(climvars$vpd_peak_anom)
vpd_peak_anom_sd <- sd(climvars$vpd_peak_anom)

SDI.Plot_mean <-  mean(plot_vars_back$SDI.Plot)
SDI.Plot_sd <- sd(plot_vars_back$SDI.Plot)

################################
## Calcuate partial residuals, effects, for avg_BA
################################

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
  
  # Calculate partial residuals for weighted pdc and its interaction with soil depth. This is a mess, 
  # but it's just residuals + beta(variable)*x
  
  # part_resids <- resid + model@frame[, "avg_BA"] * fixef(model)["avg_BA"]+ ranefs + fixef(model)["(Intercept)"] 
  
  
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


# create plots of log-partial-residuals y-axis versus natural x-scale



tiff(filename="figure4a.tif", 
     type="cairo",
     units="in", 
     width = 4, 
     height = 2, 
     pointsize=15, 
     res=600,
     compression = "lzw")



layout(matrix(c(1,2,3), nrow=1, ncol=3, byrow = TRUE))
par(mar = c(1,1,0,0), oma = c(3,3.7,1,1), family = "serif", bty = 'n', cex = 0.5)

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
          col = addTrans("blue",30), border = NA)
  
  axis(1, at = c(0, .1, .2, .3)) #draw the x-axis
  
  #draw the y-axis just for the leftmost panel
  if(i == 1){axis(2)}
  
  #label the levels of the interacting variable
  text(x = 0.2, y = -0, labels = paste0(values[i]*100, "% Peak VPD \n Anomaly = ", 
                                        round(((unname(quantile(alldata$vpd_peak_anom, values[i])) * vpd_peak_anom_sd + vpd_peak_anom_mean)), 3)),
       cex = 0.8)
  
}


#label the x- and y-axes
mtext(side = 1, text = expression(paste("Avg 4m BA (", m^2, ")")), cex= 0.7, outer = TRUE, line = 1.5)
mtext(side = 2, text = "Predicted Change in \nPlot Weighted Canopy", cex = 0.7, outer = TRUE, line = 1.1)

dev.off()


