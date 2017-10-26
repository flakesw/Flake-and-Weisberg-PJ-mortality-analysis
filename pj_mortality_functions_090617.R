## PJ Mortality functions
# Author: Sam Flake
#         sflake@gmail.com or swflake@ncsu.edu
#
# These functions are needed for "ind_tree_analysis" and for plotting figures


addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  # from Sasha Epskamp, https://github.com/SachaEpskamp/qgraph/blob/master/R/addTrans.R
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}

#calculates VIF for mixed models
vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  # from AU Franks, https://github.com/aufrank/R-hacks/blob/master/mer-utils.R
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

# Get soils data
# From Jabloun and Sahli 2006 (Model 4)

# water content at field capacity
calculate.fc <- function(Sa, Si, Cl){
  fc <- 118.932*(Cl) + 119.0866*(Si) + 119.1104*(Sa) +
    162.31731/(Cl) - 46.21921/(Si) - 5.12991/(Sa) +
    18.1733 * log(Cl) + 0.0013*(Cl)*(Si) + 0.0022*(Si)*(Sa) -
    11939.3493
  return(fc)
}

# water content at permament wilt point
calculate.pwp <- function(Sa, Si, Cl){
  pwp <- -1.5722 * (Si) - 0.5423*(Sa) - 0.0072*((Cl)^2)+ 0.0072*((Si)^2) -
    0.0059*((Sa)^2) + 160.14591/(Cl) + 6.60011/(Sa) +
    0.0022*(Cl)*(Si) - 0.0039*(Cl)*(Sa) + 92.3851
  return(pwp)
}


# AWC as difference between water capacity at field capacity 
# and permanent wilt point
calculate.awc <- function(soil_in){
Sa <- soil_in$Percent.Sand
Si <- soil_in$Percent.Silt
Cl <- soil_in$Percent.Clay

fc <- calculate.fc(Sa, Si, Cl)
pwp <- calculate.pwp(Sa, Si, Cl)

awc <- fc - pwp
soil_data <- data.frame(Plot = as.character(soil_in[, 1]), 
                           Pct_Sa = Sa, 
                           Pct_Si = Si, 
                           Pct_Cl = Cl, 
                           Field_cap = fc,
                           Perm_wilt = pwp, 
                           AWC = awc)
						   
return(soil_data)
}


minmax.sequence <- function(data, variable, length){
  
  seq(min(subset(data, select = variable)), 
    max(subset(data, select = variable)), 
    length.out=length)
  
}

# function to calculate standard errors for a given vcov matrix and row of x values, to get
# point estimates of prediction error for a given prediction. cf. source code for effects::effect.default
# This gets used for the confidence interval of the effects plot

calc.se <- function(x){
  sqrt(as.matrix(x) %*% vcov(model) %*% t(x))
} 
