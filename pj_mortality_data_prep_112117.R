# PJ mortality data prep
# Author: Sam Flake
#         sflake@gmail.com or swflake@ncsu.edu
# Description: This script does some data wrangling to turn the raw individual files
#   data into usable combined data format. It also does some aggregation to generate
#   plot-level variables.


library(plyr)

source("pj_mortality_functions_090617.R")

setwd("C:/Users/Sam/Google Drive/Projects/MS Thesis/Flake and Weisberg PJ mortality analysis")

#read data
all_clim <- read.csv("./Raw data/ALL_climate_variables.csv") #climate from PRISM and derived variables from T. Dilts
all_clim$peak_tmax <- all_clim$peak_tmax / all_clim$avg_summer_tmax # calculate relativized peak t-max
fdsi <- read.csv("./Raw data/fdsi_annual_timeseries.csv") #fdsi and some other variables derived from PRISM
map_mat <- read.csv("./Raw data/map_mat_normals.csv") #MAP and MAT from PRISM
winterppt <- read.csv("./Raw data/ppt_winter_normals.csv") #winter PPT from PRISM
names(winterppt)[2] <- "site"
plot_vars <- read.csv("./Raw data/all_vars_EXPORT.csv") #plot-level structure variables
tree_vars <- read.csv("./Raw data/all_trees_with_delta_and_ENN_041916.csv") #tree data with delta PDC and ENN distances calculated

soil_in <- read.csv("./Raw data/soils_missing_imputed_012016.csv") #raw soil data

alldata <- join(all_clim, map_mat, by = c("site"), type = "inner")
alldata <- join(alldata, winterppt[, 2:3], by = c("site"), type = "inner")

tree_vars$Diam <- 2 * sqrt(tree_vars$BA_cm / 3.14159) # calculate effective diam from cross-sectional area

########################################################################
########## Get plot-level vars
########################################################################
soils <- calculate.awc(soil_in)
soils_back <- soils

#get means of some tree-level variables

names(plot_vars)[2] <- "site"
names(soils)[1] <- "site"
names(tree_vars)[3] <- "site"

options(na.action = na.omit)
i <- 1
for (i in 1:nrow(plot_vars)){
  plot_vars$avg_height[i] <- mean(tree_vars[tree_vars$site == as.character(plot_vars$site[i]), "Height"], na.omit = TRUE)
  plot_vars$avg_Diam[i] <- mean(tree_vars[tree_vars$site == as.character(plot_vars$site[i]), "Diam"], na.omit = TRUE)
  plot_vars$avg_ENN[i] <- mean(tree_vars[tree_vars$site == as.character(plot_vars$site[i]), "ENN_dist"], na.omit = TRUE)
  plot_vars$avg_BA[i] <- mean(as.numeric(tree_vars[tree_vars$site == as.character(plot_vars$site[i]), "BA_cm"]), na.rm = TRUE)
  plot_vars$avg_BA_2m[i] <- mean(as.numeric(tree_vars[tree_vars$site == as.character(plot_vars$site[i]), "BA_2m"]), na.rm = TRUE)
  plot_vars$avg_BA_4m[i] <- mean(as.numeric(tree_vars[tree_vars$site == as.character(plot_vars$site[i]), "BA_4m"]), na.rm = TRUE)
  plot_vars$avg_BA_6m[i] <- mean(as.numeric(tree_vars[tree_vars$site == as.character(plot_vars$site[i]), "BA_6m"]), na.rm = TRUE)
  plot_vars$fdsi_anom[i] <- min(fdsi[fdsi$year %in% c(2005:2015) & fdsi$plot == as.character(plot_vars$site[i]), "fdsi"])
  plot_vars$BA.Plot[i] <- sum(tree_vars[tree_vars$site == as.character(plot_vars$site[i]), "BA_cm"], na.omit = TRUE) / 1e4
  plot_vars$Dens.Plot[i] <- nrow(tree_vars[tree_vars$site == as.character(plot_vars$site[i]), ])
}

plot_vars_back <- plot_vars


alldata <- join(plot_vars, alldata,  by = "site", type = "inner")
alldata <- join(alldata, soils[, c(1,7)], by = "site", type = "inner")

write.csv(alldata, "./Clean data/plot_level_variables_all_columns.csv")

alldata_reduced <- alldata[, -c(1, 5, 9:12, 13, 14, 17:22, 24:26, 30, 31, 32, 34, 35, 37, 39, 43, 44)]

write.csv(alldata_reduced, "./Clean data/plot_level_vars.csv")

plot_data_scaled <- alldata_reduced

#scale all the plot-level data, except the response variables
plot_data_scaled[, -c(2:5)][, sapply(plot_data_scaled[, -c(2:5)], is.numeric)] <- #scale everything
  scale(alldata_reduced[, -c(2:5)][, sapply(alldata_reduced[, -c(2:5)], is.numeric)])

write.csv(plot_data_scaled, "./Clean data/plot_level_vars_scaled.csv")

#---------------------------------------------------------------------------
###################### Tree-level data prep

## Data prep
all_tree_vars <- join(tree_vars, alldata_reduced, by = "site", type = "inner") #create tree dataframe

all_tree_vars <- all_tree_vars[all_tree_vars$LiveDead == "L", ] #only trees alive in 2005
all_tree_vars$site <- as.factor(all_tree_vars$site)

# Figure out which trees died from 2005-2015 (were live in 2005 and dead in 2015)
# LiveDead is 2005; Live is 2015. Confusing I know
all_tree_vars$Died <- "N"

for(i in 1:nrow(all_tree_vars)){
  if (all_tree_vars$Live[i] == "N" & all_tree_vars$LiveDead[i] == "L"){
    all_tree_vars$Died[i] <- "Y"
  }
}

all_tree_vars$Died <- as.factor(all_tree_vars$Died)


# Assign mortality agent presence/absence based on mortality agent lisings in the 
# raw data
all_tree_vars$ips <- ifelse(all_tree_vars$MortalityCause1 =="IPS" | all_tree_vars$MortalityCause2 =="IPS" |all_tree_vars$MortalityCause3 =="IPS" |
                              all_tree_vars$MortalityCause4 =="IPS" |all_tree_vars$MortalityCause5 =="IPS", "Y", "N")
all_tree_vars[is.na(all_tree_vars$ips), "ips"] <- "N"
all_tree_vars$rtb <- ifelse(all_tree_vars$MortalityCause1 =="RTB" | all_tree_vars$MortalityCause2 =="RTB" |all_tree_vars$MortalityCause3 =="RTB" |
                              all_tree_vars$MortalityCause4 =="RTB" |all_tree_vars$MortalityCause5 =="RTB", "Y", "N")
all_tree_vars[is.na(all_tree_vars$rtb), "rtb"] <- "N"
all_tree_vars$ptb <- ifelse(all_tree_vars$MortalityCause1 =="PTB" | all_tree_vars$MortalityCause2 =="PTB" |all_tree_vars$MortalityCause3 =="PTB" |
                              all_tree_vars$MortalityCause4 =="PTB" |all_tree_vars$MortalityCause5 =="PTB", "Y", "N")
all_tree_vars[is.na(all_tree_vars$ptb), "ptb"] <- "N"
all_tree_vars$pns <- ifelse(all_tree_vars$MortalityCause1 =="PNS" | all_tree_vars$MortalityCause2 =="PNS" |all_tree_vars$MortalityCause3 =="PNS" |
                              all_tree_vars$MortalityCause4 =="PNS" |all_tree_vars$MortalityCause5 =="PNS", "Y", "N")
all_tree_vars[is.na(all_tree_vars$pns), "pns"] <- "N"
all_tree_vars$psf <- ifelse(all_tree_vars$MortalityCause1 =="PSF" | all_tree_vars$MortalityCause2 =="PSF" |all_tree_vars$MortalityCause3 =="PSF" |
                              all_tree_vars$MortalityCause4 =="PSF" |all_tree_vars$MortalityCause5 =="PSF", "Y", "N")
all_tree_vars[is.na(all_tree_vars$psf), "psf"] <- "N"
all_tree_vars$ptm <- ifelse(all_tree_vars$MortalityCause1 =="PTM" | all_tree_vars$MortalityCause2 =="PTM" |all_tree_vars$MortalityCause3 =="PTM" |
                              all_tree_vars$MortalityCause4 =="PTM" |all_tree_vars$MortalityCause5 =="PTM", "Y", "N")
all_tree_vars[is.na(all_tree_vars$ptm), "ptm"] <- "N"
all_tree_vars$pmb <- ifelse(all_tree_vars$MortalityCause1 =="PMB" | all_tree_vars$MortalityCause2 =="PMB" |all_tree_vars$MortalityCause3 =="PMB" |
                              all_tree_vars$MortalityCause4 =="PMB" |all_tree_vars$MortalityCause5 =="PMB", "Y", "N")
all_tree_vars[is.na(all_tree_vars$pmb), "pmb"] <- "N"
all_tree_vars$mt <- ifelse(all_tree_vars$MortalityCause1 =="MT" | all_tree_vars$MortalityCause2 =="MT" |all_tree_vars$MortalityCause3 =="MT" |
                             all_tree_vars$MortalityCause4 =="MT" |all_tree_vars$MortalityCause5 =="MT" | all_tree_vars$DMR>0, "Y", "N")
all_tree_vars[is.na(all_tree_vars$mt), "mt"] <- "N" 

all_tree_vars$bb <- ifelse(all_tree_vars$MortalityCause1 =="BB" | all_tree_vars$MortalityCause2 =="BB" |all_tree_vars$MortalityCause3 =="BB" |
                             all_tree_vars$MortalityCause4 =="BB" |all_tree_vars$MortalityCause5 =="BB", "Y", "N")
all_tree_vars[is.na(all_tree_vars$bb), "bb"] <- "N" 


# save a backup before transforming
all_tree_vars_orig <- all_tree_vars

write.csv(all_tree_vars_orig, "./Clean data/all_tree_vars_unscaled_untrans_all_columns.csv")

all_tree_vars <- all_tree_vars[, -c(1, 5, 6, 10, 11, 12, 14:21, 23:37, 56:59, 62,63,64)]

write.csv(all_tree_vars, "./Clean data/all_tree_vars_unscaled_untrans.csv")

#log-transfor some variables
all_tree_vars$Neighbor_larger_log <- log(all_tree_vars$Neighbor_larger)
all_tree_vars$Diam_log <- log(all_tree_vars$Diam)
all_tree_vars$ENN_dist_log <- log(all_tree_vars$ENN_dist)
all_tree_vars$BA_2m_log <- log(all_tree_vars$BA_2m + .01)
all_tree_vars$BA_4m_log <- log(all_tree_vars$BA_4m + .01)
all_tree_vars$BA_6m_log <- log(all_tree_vars$BA_6m + .01)

write.csv(all_tree_vars, "./Clean data/all_tree_vars_unscaled_trans.csv")

#scale all the numeric variables except the response variable
all_tree_vars[, -54][, sapply(all_tree_vars[, -54], is.numeric)] <- #scale everything
  scale(all_tree_vars[, -54][, sapply(all_tree_vars[, -54], is.numeric)])

write.csv(all_tree_vars, "./Clean data/all_tree_vars_scaled_and_transformed.csv")
