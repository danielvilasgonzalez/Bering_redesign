####################################################################
####################################################################
##    
##    simulate datafrom fit object, simulate data and get design sample index
##    get bias
##    danielvilasgonzalez@gmail.com/dvilasg@uw.edu
##
####################################################################
####################################################################


#clear all objects
rm(list = ls(all.names = TRUE)) 

#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c('ggplot2','units')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#install coldpool to extract SBT for the EBS
if (!('VAST' %in% installed.packages())) {
  devtools::install_github("james-thorson/VAST@main", INSTALL_opts="--no-staged-install")};library(VAST)

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#setwd
out_dir<-'C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/'
setwd(out_dir)

#list of sp
spp<-list.dirs('./data processed/species/',full.names = FALSE,recursive = FALSE)

#selected species
spp<-c('Limanda aspera',
       'Gadus chalcogrammus',
       'Gadus macrocephalus',
       'Atheresthes stomias',
       'Reinhardtius hippoglossoides',
       'Lepidopsetta polyxystra',
       'Hippoglossoides elassodon',
       'Pleuronectes quadrituberculatus',
       'Hippoglossoides robustus',
       'Boreogadus saida',
       'Eleginus gracilis',
       'Anoplopoma fimbria',
       'Chionoecetes opilio',
       'Paralithodes platypus',
       'Paralithodes camtschaticus')

#loop over spp
for (sp in spp) {
  
  sp<-"Gadus macrocephalus"
  
  #load fit file
  load(paste0('./shelf EBS NBS VAST/',sp,'/fit.RData')) #fit
  
  #reload model
  fit_new<-reload_model(x = fit)
    
  ##################################################
  ####   Simulate 1000 iterations of data
  ##################################################
  sim_data <- array(data = NA, dim = c(nrow(fit$spatial_list$n_g),
                                       length(unique(fit$year_labels)),
                                       1000))
  
  for (isim in 1:1000) {
    
    isim<-1
    
    Sim1 <- FishStatsUtils::simulate_data(fit = fit, 
                                          type = 1, 
                                          random_seed = isim)
    
    sim_data[, , isim] <- matrix(data = Sim1$b_i[pred_TF == 1] * 0.001, 
                                 nrow = nrow(grid_goa), 
                                 ncol = length(unique(data$YEAR)))
    if(isim%%100 == 0) print(paste("Done with", ispp, "Iteration", isim))
  }
  