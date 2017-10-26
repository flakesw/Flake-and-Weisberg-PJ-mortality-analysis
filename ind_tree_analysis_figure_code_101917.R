# Create figures for Flake and Weisberg PJ Mortality
library("effects")
library("lme4")

source("pj_mortality_functions_090617.R")

##############################################
# Figure SXXX
# Does PDC predict mortality?

# #only climate vars
p_mort_clim <- glmer(Died ~
                       Avg_depth + cwd_normal_cum + fdsi_anom + (1| site), nAGQ = 1, family = "binomial",
                     data = all_tree_vars[all_tree_vars$Spp == "PIMO", ])
summary(p_mort_clim)
plot(allEffects(p_mort_clim))
coefplot2::coefplot2(p_mort_clim)
AICc(p_mort_clim)
vif.mer(p_mort_clim)
r.squaredGLMM(p_mort_clim)

#only stand structure vars
p_mort_stand <- glmer(Died ~ ENN_dist + Neighbor_larger + BA_4m +
                        (1|Cluster/site), nAGQ = 1, family = "binomial",
                      data = all_tree_vars[all_tree_vars$Spp == "PIMO", ])
summary(p_mort_stand)
plot(allEffects(p_mort_stand))
coefplot2::coefplot2(p_mort_stand)
AICc(p_mort_stand)
vif.mer(p_mort_stand)
r.squaredGLMM(p_mort_stand)



#only tree-level vars
p_mort_tree <- glmer(Died ~ Height + Diam +
                       (1|Cluster/site), nAGQ = 1, family = "binomial",
                     data = all_tree_vars[all_tree_vars$Spp == "PIMO", ])
summary(p_mort_tree)
plot(allEffects(p_mort_tree))
coefplot2::coefplot2(p_mort_tree)
AICc(p_mort_tree)
vif.mer(p_mort_tree)
r.squaredGLMM(p_mort_tree)



#intercept-only model
p_mort_int <- glm(Died ~ 1, family = "binomial",
                  data = all_tree_vars[all_tree_vars$Spp == "PIMO", ])

# random effects only model
p_mort_int_mix <- glmer(Died ~ 1 + (1| site), nAGQ = 1, family = "binomial",
                        data = all_tree_vars[all_tree_vars$Spp == "PIMO", ])

## Testing effect of different forms of PDC05 (percent dead crown in 2005)

pdc_mort <- update(p_mort, . ~ . + PDC05)
relgrad <- with(pdc_mort@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
plot(allEffects(pdc_mort))

pdc_mort_poly <- update(p_mort, . ~ . + poly(PDC05, 2))
relgrad <- with(pdc_mort_poly@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
plot(allEffects(pdc_mort_poly))

pdc_mort_exp <- update(p_mort, . ~ . + exp(PDC05))
relgrad <- with(pdc_mort_exp@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))

#these  models aren't great but close enough to be useful, I think


## Model comparisons
AICc(p_mort)
AICc(p_mort_clim)
AICc(p_mort_stand)
AICc(p_mort_tree)
AICc(p_mort_int)
AICc(p_mort_int_mix)

AICc(pdc_mort)
AICc(pdc_mort_poly)
AICc(pdc_mort_exp)

r.squaredGLMM(p_mort)
r.squaredGLMM(p_mort_clim)
r.squaredGLMM(p_mort_stand)
r.squaredGLMM(p_mort_tree)
pR2(p_mort_int)
r.squaredGLMM(p_mort_int_mix)


r.squaredGLMM(pdc_mort)
r.squaredGLMM(pdc_mort_poly)
r.squaredGLMM(pdc_mort_exp)

# hists and confusion matrices; not finished
hist(fitted(p_mort)[all_tree_vars$Spp == "PIMO" & all_tree_vars$Died == "Y"])
hist(fitted(p_mort)[all_tree_vars$Spp == "PIMO" & all_tree_vars$Died == "N"])
hist(fitted(p_mort)[which(all_tree_vars$Spp == "PIMO" & all_tree_vars$Died == "Y")])
hist(fitted(p_mort)[which(all_tree_vars$Spp == "PIMO" & all_tree_vars$Died == "N")])

######################################################
# Figure SXXX
# Plot the ROCs, adding the AUC for each
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
tiff(filename="./plots/classification_accuracy.tif", 
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
#### Generate figure 4
###############################################################
bootstrap_pmort <- readRDS("./model output/bootstrap_pmort.rds")

pdc_mort_poly <- update(p_mort, . ~ . + poly(PDC05, 2))

####partial effect plot of 2005 PDC on mortality risk 2005-2015

minx <- -mean(all_tree_vars_orig$PDC05)/sd(all_tree_vars_orig$PDC05)
maxx <- (95-mean(all_tree_vars_orig$PDC05))/sd(all_tree_vars_orig$PDC05)



#PDC05 effect from model
mort_eff <- Effect(pdc_mort_poly, focal.predictors = c("PDC05"), xlevels = list(PDC05 = seq(minx,maxx,length.out = 100)))

#prediction and 1-se envelope
x <- unlist(mort_eff$x)*sd(all_tree_vars_orig$PDC05) + mean(all_tree_vars_orig$PDC05) #unscale x variable
up <- exp(unlist(mort_eff$fit) + (mort_eff$se))
low <- exp(unlist(mort_eff$fit) - (mort_eff$se))
y <- exp(unlist(mort_eff$fit))

tiff(filename="./plots/Figure_4_pdc_mort_effects.tif", 
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
     ylim = c(0, 0.3),
     xlim = c(0, 100),
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
axis(2, cex.axis = 0.7, at = c(0, 0.1, 0.2, .3))
mtext(text = "Percent dead canopy in 2005", side = 1, cex = 0.8, line = 2.3)
mtext(text = "p(mortality), 2005-2015", side = 2, cex = 0.8, line = 2.3)


dev.off()


#------------------------------------------------------------
## Plots for bootstraps of variable values
#-------------------------------------------------
#Generate Figure 3

bootstrap_pimo_int <- readRDS("./model output/bootstrap_pimo_int.rds")
bootstrap_juos <- readRDS("./model output/bootstrap_juos.rds")
bootstrap_pmort <- readRDS("./model output/bootstrap_pmort.rds")

#order for variables to be plotted
var_names <- character()
var_names[1] <- "Intercept"
var_names[2] <- "Neigh. Rat. (log)"
var_names[3] <- "ENN Dist (log)"
var_names[4] <- "4m BA (log)"
var_names[5] <- "Height"
var_names[6] <- "Diam (log)"
var_names[7] <- "4m BA : CWD"
var_names[8] <- "ENN Dist : Neigh. Rat."

attr(bootstrap_pimo_int$t0, "names")
attr(bootstrap_juos$t0, "names")
attr(bootstrap_pmort$t0, "names")

#manually select the order of variables to line up with the variable names
order_p <- c(1, 9, 7, 5, 8, 6, 10)
order_j <- c(1, 9, 7, 5, 8, 6, 10)
order_pmort <- c(1, 9, 7, 5, 8, 6)

opar <- par(no.readonly = TRUE)

par(opar)

tiff(filename="./plots/Figure_3_coefficient_density_plots.tif", 
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


#loop over each variable
for (i in 1:length(var_names)){
  
  #pull out density plots for selected variables, controlled by the order_x lists defined above
  if(i <= 6){
    d_p <- density(bootstrap_pimo_int$t[ , order_p[i]])
    d_j <- density(bootstrap_juos$t[, order_j[i]])
    d_pm <- density(bootstrap_pmort$t[, order_pmort[i]])}
  
  #changing assignment for the interaction terms
  if(i >= 7){
    d_p <- density(bootstrap_pimo_int$t[ , order_p[7]])
    d_j <- density(bootstrap_juos$t[, order_j[7]])
  }
  
  if(i == 1){ #different axes for intercept
    
    #plot pinyon and juniper dieback density plots
    plot(NA,
         xlim = c(min(c(min(d_p$x), min(d_j$x))), max(c(max(d_p$x), max(d_j$x)))),
         ylim = c(0,max(c(max(d_p$y), max(d_j$y)))),
         xaxt = 'n', 
         yaxt = 'n',
         bty = 'n',
         xlab = "",
         ylab = "")
    lines(d_p, col = "grey50") #pinyon in darker grey
    polygon(d_p, col = addTrans("grey50", 90))
    lines(d_j, col = "grey90", new = FALSE) #juniper in lighter grey
    polygon(d_j, col = addTrans("grey90", 60))
    axis(1)
    
    #label the row with the i'th value of the var_names vector
    mtext(var_names[i],side=3,line=-1.8, 
          at=par("usr")[1] - 3,
          cex=.75,
          adj = 1) #in this case, Intercept
    
    mtext("Canopy dieback", side = 3, cex = 0.9) #title for left column
    
    
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
    
    mtext("Pinyon mortality", side = 3, cex = 0.9) #label for right column
    
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
    
    mtext(var_names[i],side=3,line=-1.5, 
          at=par("usr")[1] - .7,
          cex=.75,
          adj = 1) # label row with variable name
    
    
    
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
    
    abline(v = 0, lwd = 2, lty = 2)
    
    if (i == 6){axis(1, at = c(-1, 0, 1, 2))}
  }
  
  if(i == 7){ #interaction term for pinyon
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
          at=par("usr")[1] - .7,
          cex=.75,
          adj = 1)
    
  }
  
  if(i == 8){ #interaction term for juniper dieback
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
          at=par("usr")[1] - .7,
          cex=.75,
          adj = 1)
    
    axis(1)
  }
  
}




mtext(text = "Bootstrap coefficient estimates (scaled)", side = 1,
      line = 2.3, at = 0, cex = 0.9)

dev.off()


#####################################################################
### Generate figure 2
### multipanel figure
#####################################################################
## this is the ugliest, WETtest piece of code you've ever seen but it should work.
## run this whole chunk until dev.off() way down there about 350 lines
pardefault <- par(no.readonly = T)


# for unscaling variables -- these things come from a full model in the model set, you'll have to change
# the number or the path depending upon what kind of model object and the order of your models, etc.
# You can also just get them from the raw data.
library(lme4)


all_tree_vars <- read.csv("./Clean data/all_tree_vars_scaled_and_transformed.csv")
all_tree_vars_trans <- read.csv("./Clean data/all_tree_vars_unscaled_trans.csv")
all_tree_vars_unscaled <- read.csv( "./Clean data/all_tree_vars_unscaled_untrans.csv")

alldata <- read.csv("./Clean data/plot_level_vars_scaled.csv")

alldata_unscaled <- read.csv("./Clean data/plot_level_vars.csv")


BA_4m_mean <- mean((all_tree_vars_trans[all_tree_vars_trans$Spp == "PIMO", "BA_4m_log"]))
BA_4m_sd <- sd((all_tree_vars_trans[all_tree_vars_trans$Spp == "PIMO", "BA_4m_log"]))

ENN_dist_mean <- mean((all_tree_vars_trans[all_tree_vars_trans$Spp == "JUOS", "ENN_dist_log"]))
ENN_dist_sd <- sd((all_tree_vars_trans[all_tree_vars_trans$Spp == "JUOS", "ENN_dist_log"]))

Neighbor_larger_mean <- mean((all_tree_vars_trans[all_tree_vars_trans$Spp == "JUOS", "Neighbor_larger_log"]))
Neighbor_larger_sd <- sd((all_tree_vars_trans[all_tree_vars_trans$Spp == "JUOS", "Neighbor_larger_log"]))

cwd_normal_cum_mean <- mean(all_tree_vars_unscaled[all_tree_vars_unscaled$Spp == "PIMO", "cwd_normal_cum"])
cwd_normal_cum_sd<- sd(all_tree_vars_unscaled[all_tree_vars_unscaled$Spp == "PIMO", "cwd_normal_cum"])



#open device

tiff(filename="./plots/Figure 2 interaction plots.tif",
     type="cairo",
     units="in",
     width = 6,
     height = 6,
     pointsize=15,
     res=600,
     compression = "lzw")



layout(matrix(c(1,2,3,4,5,6), nrow=2, ncol=3, byrow = TRUE))
par(mar = c(3,.5,1.2,0), oma = c(2,4,1,1), family = "serif", bty = 'n', 
    cex = 0.7, cex.lab = 0.7, cex.axis = 0.7)

#############################################################
### Pimo part -- figure 4a
#############################################################

model <- readRDS("./model output/pimo_model.rds")
# 
# num_tree <- (nrow(all_tree_vars[all_tree_vars$Spp == "PIMO", ]))
# ranefs <- vector(mode = "numeric", length = num_tree) 

###BA_4meter
calculate_preds_BA_4m <- function(quant){
  #set all terms = 0, which is their mean value because all variables are standardized
  C <- data.frame(matrix(ncol = length(names(coef(model)$Cluster)), nrow = 1000))
  C <- as.data.frame(apply(C, c(1,2), FUN = function(x){return(0)}))
  names(C) <- names(coef(model)$Cluster)
  C$"(Intercept)" <- 1
  C$cwd_normal_cum <- unname(quantile(all_tree_vars$cwd_normal_cum, quant))
  #min value calculated to correspond to a real BA of -.001
  C$BA_4m_log <- seq(-1.068, max(subset(all_tree_vars, select = "BA_4m_log")), length.out = 1000)
  C$"cwd_normal_cum:BA_4m_log" = C$cwd_normal_cum * C$BA_4m_log
  
  se <- vector(mode='numeric', length=nrow(C))
  
  for (i in 1:nrow(C)){
    se[i] <- as.numeric(calc.se(C[i, ]))
  }
  
  
  preds <- data.frame(
    predictions = as.matrix(C) %*% as.matrix(fixef(model)),
    lo = as.matrix(C) %*% as.matrix(fixef(model)) - 1.96*se,
    hi = as.matrix(C) %*% as.matrix(fixef(model)) + 1.96*se,
    # BA_4m = model@frame$`BA_4m`,
    BA_4m_sim = C$BA_4m_log,
    BA_4m_sim_natural_scale = (exp((C$BA_4m_log * BA_4m_sd) + BA_4m_mean)/10000) # in m2
  )
  
  return(preds)
  
}

calculate_partial_residuals_BA_4m <- function(quant){
  model_frame <- model@frame
  resids <- resid(model)
  mod_matrix <- model.matrix(model)
  fixed <- fixef(model)
  
  
  if(quant == 0.1){
    cwd_low <-  unname(quantile(all_tree_vars$cwd_normal_cum, 0))
    cwd_high <-  unname(quantile(all_tree_vars$cwd_normal_cum, .30))
  }
  if(quant == 0.5){
    cwd_low <-  unname(quantile(all_tree_vars$cwd_normal_cum, 0.3))
    cwd_high <-  unname(quantile(all_tree_vars$cwd_normal_cum, .7))
  }
  if(quant == 0.9){
    cwd_low <-  unname(quantile(all_tree_vars$cwd_normal_cum, 0.7))
    cwd_high <-  unname(quantile(all_tree_vars$cwd_normal_cum, 1))
  }
  
  select_data <- which(model_frame$cwd_normal_cum > cwd_low & model_frame$cwd_normal_cum < cwd_high)
  
  model_frame <- as.data.frame(model_frame[select_data, ])
  resids <- resids[select_data]
  mod_matrix <- mod_matrix[select_data, ]
  ranefs <- vector(mode = "numeric", length = length(select_data))
  
  for (i in (1:length(select_data))){
    j <- select_data[i]
    ranefs[i] <- ranef(model)$site[rownames(ranef(model)$site) == paste(model@frame$site[j], model@frame$Cluster[j], sep=":"), ]
  }
  
  B <- as.matrix(fixed[c("BA_4m_log")])
  X <- as.matrix(mod_matrix[, c("BA_4m_log")])
  
  part_resids <- resids + t(B) %*% t(X) + fixed["(Intercept)"] + ranefs 
  BA_4m <- mod_matrix[, "BA_4m_log"]
  BA_4m_back_trans <- (exp((BA_4m * BA_4m_sd) + BA_4m_mean)/10000)
  return(cbind(t(part_resids), BA_4m_back_trans))
  
}

for(i in 1:3){
  values = c(0.1,0.5,0.9)
  # Use the function to get all the stuff we need (predictions, partial residuals, and SEs)
  pred <- calculate_preds_BA_4m(values[i])
  #part_resids <- calculate_partial_residuals_BA_4m(values[i])
  
  
  plot(NA,
       ylim = c(-50, 10),
       xlim = c(0.001, 3),
       xlab = "",
       ylab = "",
       log = "x",
       cex.lab = 1.5,
       xaxt = "n",
       yaxt = "n"
  )
  
  #points(part_resids[, 1] ~ part_resids[, 2], pch = 21, col = "grey70", cex = 0.5)
  
  lines(pred$predictions ~ pred$BA_4m_sim_natural_scale, lwd = 2)
  lines((pred$lo) ~ pred$BA_4m_sim_natural_scale, lwd = 1.3) 
  lines((pred$hi) ~ pred$BA_4m_sim_natural_scale, lwd = 1.3)
  polygon(c(pred$BA_4m_sim_natural_scale, rev(pred$BA_4m_sim_natural_scale)), c((pred$hi), rev((pred$lo))),
          col = addTrans("grey",30), border = NA) #fills in the area between high and low confidence intervals
  
  
  
  text(x = 0.05, y = 5, labels = paste0(values[i]*100, "% CWD = \n", 
                                         round(((unname(quantile(all_tree_vars$cwd_normal_cum, values[i])) * cwd_normal_cum_sd + cwd_normal_cum_mean)), 0), " mm"),
       cex = 0.7)
  
  axis(1, at=c(0.001, 0.01, .1, 1))
  
  if(i == 1){
    mtext(side = 3, at = 0.003, text= "(a)")
    axis(2)
    mtext(side = 2, text = "Predicted Change in Canopy (%)", cex = 0.7, outer = FALSE, line = 1.7)}
  
  #label x axis
  if(i == 2){mtext(side = 1, text = expression(paste("BA, 4m buffer (", m^2, ")")), cex= 0.7, outer = FALSE, line = 2)}
  
}


#######################################################
### juos part -- figure 4b
#######################################################


model <- readRDS("./model output/juos_model.rds")

###BA_4meter
calculate_preds_Neighbor_larger <- function(quant){
  #set all terms = 0, which is their mean value because all variables are standardized
  C <- data.frame(matrix(ncol = length(names(coef(model)$Cluster)), nrow = 1000))
  C <- as.data.frame(apply(C, c(1,2), FUN = function(x){return(0)}))
  names(C) <- names(coef(model)$Cluster)
  C$"(Intercept)" <- 1
  C$ENN_dist_log <- unname(quantile(all_tree_vars$ENN_dist_log, quant))
  #min value calculated to correspond to a real BA of -.001
  C$Neighbor_larger_log <- minmax.sequence(all_tree_vars, "Neighbor_larger_log", 1000)
  C$`ENN_dist_log:Neighbor_larger_log` = C$ENN_dist_log * C$Neighbor_larger_log
  
  se <- vector(mode='numeric', length=nrow(C))
  
  for (i in 1:nrow(C)){
    se[i] <- as.numeric(calc.se(C[i, ]))
  }
  
    preds <- data.frame(
    predictions = as.matrix(C) %*% as.matrix(fixef(model)),
    lo = as.matrix(C) %*% as.matrix(fixef(model)) - 1.96*se,
    hi = as.matrix(C) %*% as.matrix(fixef(model)) + 1.96*se,
    # BA_4m = model@frame$`BA_4m`,
    Neighbor_larger_sim = C$Neighbor_larger_log,
    Neighbor_larger_natural_scale = (exp((C$Neighbor_larger_log * Neighbor_larger_sd) + Neighbor_larger_mean))
  )
  
  return(preds)
  
}

calculate_partial_residuals_NL <- function(quant){
  model_frame <- model@frame
  resids <- resid(model)
  mod_matrix <- model.matrix(model)
  fixed <- fixef(model)
  
  
  if(quant == 0.1){
    enn_low <-  unname(quantile(all_tree_vars$ENN_dist_log, 0))
    enn_high <-  unname(quantile(all_tree_vars$ENN_dist_log, .30))
  }
  if(quant == 0.5){
    enn_low <-  unname(quantile(all_tree_vars$ENN_dist_log, 0.3))
    enn_high <-  unname(quantile(all_tree_vars$ENN_dist_log, .7))
  }
  if(quant == 0.9){
    enn_low <-  unname(quantile(all_tree_vars$ENN_dist_log, 0.7))
    enn_high <-  unname(quantile(all_tree_vars$ENN_dist_log, 1))
  }
  
  select_data <- which(model_frame$ENN_dist_log > enn_low & model_frame$ENN_dist_log < enn_high)
  
  model_frame <- as.data.frame(model_frame[select_data, ])
  resids <- resids[select_data]
  mod_matrix <- mod_matrix[select_data, ]
  ranefs <- vector(mode = "numeric", length = length(select_data))
  
  for (i in (1:length(select_data))){
    j <- select_data[i]
    ranefs[i] <- ranef(model)$site[rownames(ranef(model)$site) == paste(model@frame$site[j], model@frame$Cluster[j], sep=":"), ]
  }
  
  B <- as.matrix(fixed[c("Neighbor_larger_log")])
  X <- as.matrix(mod_matrix[, c("Neighbor_larger_log")])
  
  part_resids <- resids + t(B) %*% t(X) + fixed["(Intercept)"] + ranefs 
  NL <- mod_matrix[, "Neighbor_larger_log"]
  NL_back_trans <- (exp((NL * Neighbor_larger_sd) + Neighbor_larger_mean))
  return(cbind(t(part_resids), NL_back_trans))
  
}

i <- 1
for(i in 1:3){
  values = c(0.1,0.5,0.9)
  # Use the function to get all the stuff we need (predictions, partial residuals, and SEs)
  pred <- calculate_preds_Neighbor_larger(values[i])
  
  #part_resids <- calculate_partial_residuals_NL(values[i])
  
  plot(NA,
       ylim = c(-25, 10),
       xlim = c(min(pred$Neighbor_larger_natural_scale), max(pred$Neighbor_larger_natural_scale)),
       xlab = "",
       ylab = "",
       log = "x",
       cex.lab = 1.5,
       xaxt = "n",
       yaxt = "n"
  )
  #points(part_resids[, 1] ~ part_resids[, 2], pch = 21, col = "grey70", cex = 0.5)
  
  lines(pred$predictions ~ I(exp((pred$Neighbor_larger_sim * Neighbor_larger_sd) + Neighbor_larger_mean)), lwd = 2)
  lines((pred$lo) ~ I(exp((pred$Neighbor_larger_sim * Neighbor_larger_sd) + Neighbor_larger_mean)), lwd = 1.3) 
  lines((pred$hi) ~ I(exp((pred$Neighbor_larger_sim * Neighbor_larger_sd) + Neighbor_larger_mean)), lwd = 1.3)
  polygon(c(I(exp((pred$Neighbor_larger_sim * Neighbor_larger_sd) + Neighbor_larger_mean)), rev(I(exp((pred$Neighbor_larger_sim * Neighbor_larger_sd) + Neighbor_larger_mean)))), c((pred$hi), rev((pred$lo))),
          col = addTrans("grey",30), border = NA) #fills in the area between high and low confidence intervals
  
  axis(1, at=c(0.001,.1, 10), labels = c(.001, .1, 10))
  
  text(x = 0.9, y = 5, labels = paste0(values[i]*100, "% ENN Distance = \n", 
                                        round(exp((unname(quantile(all_tree_vars$ENN_dist, values[i])) * ENN_dist_sd + ENN_dist_mean)), 1), " m"),
       cex = 0.7)
  
  
  if(i == 1){
    mtext(side = 3, at = 0.003, text= "(b)")
    axis(2)
    mtext(side = 2, text = "Predicted Change in Canopy (%)", cex = 0.7, outer = FALSE, line = 1.7)}
  
  #label x axis
  if(i == 2){mtext(side = 1, text = "Neighbor Size Ratio (log)", cex= 0.7, outer = FALSE, line = 2)}
  
}

dev.off()

par <- pardefault


#-------------------------------------------------------------------------------------------
# Cwd effects plot
#-------------------------------------------------------------------------------------------

tiff(filename="./plots/Figure X plot-level cwd effect plot.tif",
     type="cairo",
     units="in",
     width = 4,
     height = 3,
     pointsize=15,
     res=600,
     compression = "lzw")

par(mar = c(3,.5,1.2,0), oma = c(2,4,1,1), family = "serif", bty = 'n', 
    cex = 0.7, cex.lab = 0.7, cex.axis = 0.7)


### set some parameters to be use throughout
ymax <- 20
ymin <- -80

# model used for effects plot
model <- readRDS("./model output/plot_model.rds")

# ranefs <- vector(mode = "numeric", length = 98) #initialize a vector for the ranefs
#
# #pull out the random intercept for each plot
# for (i in (1:98)){
#   ranefs[i] <- ranef(model)$Cluster[rownames(ranef(model)$Cluster) == model@frame$Cluster[i], ]
# }

#set all terms = 0, which is their mean value because all variables are standardized
C <- data.frame(matrix(ncol = length(names(coef(model)$Cluster)), nrow = 1000))
C <- as.data.frame(apply(C, c(1,2), FUN = function(x){return(0)}))
names(C) <- names(coef(model)$Cluster)
C$"(Intercept)" <- 1
C$cwd_normal_cum <- minmax.sequence(alldata, "cwd_normal_cum", nrow(C))

se <- vector(mode='numeric', length=nrow(C))

for (i in 1:nrow(C)){
  se[i] <- as.numeric(calc.se(C[i, ]))
}

pred <- data.frame(
  predictions = as.matrix(C) %*% as.matrix(fixef(model)),
  lo = as.matrix(C) %*% as.matrix(fixef(model)) - 1.96*se,
  hi = as.matrix(C) %*% as.matrix(fixef(model)) + 1.96*se,
  cwd_normal_cum_sim = C$cwd_normal_cum
)
cwd_normal_cum_orig = ((model@frame$cwd_normal_cum * cwd_normal_cum_sd) + cwd_normal_cum_mean)
part_resids = model@frame$Avg_delta_pdc
#set up the empty plot window
plot(NA,
     ylim = c(ymin, ymax),
     xlim = c((min(cwd_normal_cum_orig)), max(cwd_normal_cum_orig)),
     xlab = "",
     ylab = "",
     # log = "x",
     cex.lab = 1,
     xaxt = "n",
     yaxt = "n"
)

points(part_resids ~ cwd_normal_cum_orig, pch = 20, cex = 1, col = "grey70")
#draw prediction line and confidence envelope of predictions
lines(pred$predictions ~ I(((pred$cwd_normal_cum_sim * cwd_normal_cum_sd) + cwd_normal_cum_mean)), lwd = 2)
lines((pred$lo) ~ I(((pred$cwd_normal_cum_sim * cwd_normal_cum_sd) + cwd_normal_cum_mean)), lwd = 1.3)
lines((pred$hi) ~ I(((pred$cwd_normal_cum_sim * cwd_normal_cum_sd) + cwd_normal_cum_mean)), lwd = 1.3)

#fill in the confidence envelope
polygon(c(I((pred$cwd_normal_cum_sim * cwd_normal_cum_sd) + cwd_normal_cum_mean), rev(I((pred$cwd_normal_cum_sim * cwd_normal_cum_sd) + cwd_normal_cum_mean))), c((pred$hi), rev((pred$lo))),
        col = addTrans("grey",30), border = NA)

#label the levels of the interacting variable

axis(1, cex.axis = 1) #draw the x-axis
axis(2, cex.axis = 1)

mtext(side = 2, text = "Mean change in canopy (%)", cex = 0.7, outer = FALSE, line = 2.5)
mtext(side = 1, text = "Normal annual\nclimatic water deficit (mm)", cex = 0.7, outer = FALSE, line = 2.8)

dev.off()
