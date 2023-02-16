####################################################################
####################################################################
##    
##    Check and plot covariates (Temp and Depth) effects from VAST models
##    Daniel Vilas (danielvilasgonzalez@gmail.com/dvilasg@uw.edu)
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c("effects",'ggplot2')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#install coldpool to extract SBT for the EBS
if (!('VAST' %in% installed.packages())) {
  devtools::install_github("james-thorson/VAST@main", INSTALL_opts="--no-staged-install")};library(VAST)

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#setwd
#out_dir<-'E:/UW/Adapting Monitoring to a Changing Seascape/'
out_dir<-'C:/Users/danie/Desktop/'
setwd(out_dir)

#list of sp
#splist<-list.dirs('./slope shelf EBS NBS VAST/',full.names = FALSE,recursive = FALSE)
#splist<-sort(splist[-1])

#loop over species to fit models
#for (sp in sp.list) {
  
  #sp<-'Gadus macrocephalus'
  sp<-'Gadus macrocephalus 1'
  
  #list of models
  #models<-list.dirs(paste0('./slope shelf EBS NBS VAST/',sp),full.names = FALSE,recursive = FALSE)
  models<-list.dirs(paste0('./',sp),full.names = FALSE,recursive = FALSE)
  
#####################
# Effects package
#####################
  
  for (m in models) {
    
    m<-models[1]
    load(paste0('./',sp,'/',m,'/fit.RData'))
    
    # Must add data-frames to global environment (hope to fix in future)
    covariate_data_full = fit$effects$covariate_data_full
    summary(covariate_data_full)
    catchability_data_full = fit$effects$catchability_data_full
    
    # Plot 1st linear predictor, but could use `transformation` to apply link function
    pred = Effect.fit_model( fit,
                             focal.predictors = c("ScaleLogDepth"),
                             which_formula = "X1",
                             xlevels = 100,
                             transformation = list(link=identity, inverse=exp) )
    plot(pred)
    
    

######################################
# For S.Brodie
######################################

load('C:/Users/danie/Desktop/Gadus macrocephalus 1/depth/fit.RData')

#Get parameters from VAST model object
Gammas <- fit$ParHat[c("gamma1_cp","gamma2_cp")]
#Gammas_1 <- Gammas[[1]] [1,1,1:2]
#gamma1_ctp - effect of density covariates that change over time on encounter prob
#gamma2_ctp - effect of density covariates that change over time on positive catches

df<-fit$covariate_data


Gammas_1 <- Gammas[[1]] #[i,1,1:2]
Gammas_2 <- Gammas[[2]] #[i,1,1:2]
Xlim = range(df$CPE)
Xcenter = mean(df$CPE)
X = seq(Xlim[1], Xlim[2], length=1e4)
Y1 = (X-Xcenter)*Gammas_1[1]#  + (X-Xcenter)^2*Gammas_1[1]
plot( x=X, y=exp(Y1-max(Y1)), col="black", type="l", lty="solid", lwd=2, xlab="Depth",ylab="Response", main=paste0("Occurrence "))
Y2 = (X-Xcenter)*Gammas_2[1] + (X-Xcenter)*Gammas_2[2]
plot( x=X, y=exp(Y2-max(Y2)), col=c("black","grey"), type="l", lty="solid", lwd=2, xlab="Depth",ylab="Response", main=paste0("Abundance ",spp.names[i]))


#####################
# Effects package
#####################

library(effects)  # Used to visualize covariate effects

# Must add data-frames to global environment (hope to fix in future)
covariate_data_full = fit$effects$covariate_data_full
catchability_data_full = fit$effects$catchability_data_full

# Plot 1st linear predictor, but could use `transformation` to apply link function
pred = Effect.fit_model( fit,
                         focal.predictors = c("CPE"),
                         which_formula = "X1",
                         xlevels = 100,
                         transformation = list(link=identity, inverse=identity) )
plot(pred)

#####################
# pdp package
#####################

library(pdp)

# Make function to interface with pdp
pred.fun = function( object, newdata ){
  predict( x=object,
           Lat_i = object$data_frame$Lat_i,
           Lon_i = object$data_frame$Lon_i,
           t_i = object$data_frame$t_i,
           a_i = object$data_frame$a_i,
           what = "P1_iz",
           new_covariate_data = newdata,
           do_checks = FALSE )
}

# Run partial
Partial = partial( object = fit,
                   pred.var = "CPE",
                   pred.fun = pred.fun,
                   train = fit$covariate_data )

# Make plot using ggplot2
library(ggplot2)
autoplot(Partial)
