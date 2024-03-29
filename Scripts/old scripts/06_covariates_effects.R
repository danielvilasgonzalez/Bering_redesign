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
pack_cran<-c("effects",'ggplot2','splines','cowplot')

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
splist<-list.dirs('./data processed/',full.names = FALSE,recursive = FALSE)

#folder region
#fol_region<-'slope EBS VAST'

#####################
# Effects package
#####################

#formulas
fs<-c('X1','X2')

#loop over species to fit models
#for (sp in sp.list) {
  
  sp<-'Gadus macrocephalus'
  
  #check effects
  #load(paste(out_dir,fol_region,sp,'diagnostics.RData',sep='/'))
  
  #load data_geostat to unscale
  data_geostat<-readRDS(paste(out_dir,fol_region,sp,'/data_geostat_temp.rds',sep='/'))
  
  #list of models
  models<-list.dirs(paste(out_dir,fol_region,sp,sep = '/'),full.names = FALSE,recursive = FALSE)
  models<-models[models!='null']

#pdf file for each sp  
  pdf(file = paste(out_dir,fol_region,sp,paste0(sp,'_effects.pdf'),sep='/'),   # The directory you want to save the file in
      width = 10, # The width of the plot in inches
      height = 10,
      onefile = TRUE)
  
  #grid list
  grid_list<-list()

  for (m in models) {
    m<-models[1]
    
    #plot list
    plot_list<-list()
    
    #load fit file
    #load(paste(out_dir,fol_region,sp,m,'/fit.RData',sep = '/'))
    load('C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/slope EBS VAST/Gadus macrocephalus/depth3d/fit.RData')
    load('C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/shelf EBS NBS VAST/Gadus macrocephalus/temp3d/b2_19822022fit.RData')
    
    
    if (grepl('depth',m)) {
      
      #loop over formulas
      for (f in fs) {
      
        #formula
        f<-'X2'
        
        #for title purposes
        f1<-ifelse(f=='X1','presence','positive')
        
        # Must add data-frames to global environment (hope to fix in future)
        covariate_data_full = fit$effects$covariate_data_full
        catchability_data_full = fit$effects$catchability_data_full
        
        #calculate effect of covariate
        pred = Effect.fit_model( fit,
                                 focal.predictors = c("BotTemp"),
                                 which_formula = f,
                                 xlevels = 100,
                                 transformation = list(link=identity, inverse=indentity) )
        
        #plot response
        p<-plot(pred,
                main=paste('SBT effect',f1),
                axes=list(grid=TRUE, 
                          x=list(rug=TRUE,
                                 BotTemp=list(lab="SBT (°C)")),
                          y=list(rug=TRUE,
                                 lab="Catch (kg)")))
        
        #add plot to list
        plot_list[[paste0(m,'depth_',f)]]<-p
        
        #p<-NULL
        
      }
    }
      
    if (grepl('temp',m)) {
        
        #loop over formulas
        for (f in fs) {
          
          #formula
          #f<-'X2'
          
          #for title purposes
          f1<-ifelse(f=='X1','presence','positive')
          
          #calculate effect of covariate
          pred = Effect.fit_model( fit,
                                   focal.predictors = c("ScaleBotTemp"),
                                   which_formula = f,
                                   xlevels = 100,
                                   transformation = list(link=identity, inverse=indentity) )
          
          #plot response
          p<-plot(pred,
                  main=paste(m,'model -',f1),
                  axes=list(grid=TRUE, 
                            x=list(rug=FALSE,
                                   ScaleBotTemp=list(lab="Temp (scale °C)")),
                            y=list(rug=FALSE,
                                   lab="Catch (kg)")))
          
          #add plot to list
          plot_list[[paste0(m,'temp_',f)]]<-p
          
          #p<-NULL
          
          }
      }
        
        
    #if (grepl('depth',m) & grepl("temp",m)) {
      #multiplot
      plot_row<-plot_grid(plotlist = plot_list,nrow = 2,ncol=2)
    #}   else {
      #multiplot
    #  plot_row<-plot_grid(plotlist = plot_list,nrow = 2,ncol=2)
    #} 
  
    # now add the title
    title <- ggdraw()+ 
             draw_label(sp,
                        x = 0,
                        size = 16)+
             theme(plot.margin = margin(0, 0, 0, 7))#;plot(title)
    
    #print plot
    
      gr<-plot_grid(title,
                    plot_row,
                    ncol = 1,
                    rel_heights = c(0.1, 1))
    
    grid_list[[m]]<-gr
    
    print(gr)
    #gr<-NULL
    #close pdf
    #dev.off()
    gc()  
}
  
  dev.off()
  #pdf file for each sp  
  # pdf(file = paste0(getwd(),'/slope shelf EBS NBS VAST/',sp,'/',sp,'_effects.pdf'),   # The directory you want to save the file in
  #     width = 10, # The width of the plot in inches
  #     height = 10,
  #     onefile = TRUE)
  # plot_list
  # dev.off()

#}  
  
    #try to use transform argument for plot effects
    # #You can use the attributes to unscale:
    # scaled<-scale(covariate_data_full$LogDepth)
    # unscaled<-scaled * attr(scaled, 'scaled:scale') + attr(scaled, 'scaled:center')
    # expunscaled<-exp(unscaled)
    # 
    # y<-pred$data$ScaleLogDepth
    # exp_unscale<-function(x) as.vector(exp(x * attr(scaled, 'scaled:scale') + attr(scaled, 'scaled:center')))
    # exp_unscale(scaled)
    # log_scale<-function(x) as.vector(scale(log(x)))
library(marginaleffects)
    
  # Plot 1st linear predictor, but could use `transformation` to apply link function
  quant = function(x) seq(min(x),max(x),length=21)
  #because of the scale thing is not exactly a column and it comes with additional attributes
  df<-data.frame(fit$covariate_data$ScaleLogDepth)
  rownames(df)<-rownames(fit$covariate_data)
  colnames(df)<-'ScaleLogDepth'
  newdata = datagrid( newdata=df[,'ScaleLogDepth',drop=FALSE], ScaleLogDepth = quant )
  pred = predictions( fit, newdata=newdata, covariate="X1",conf_level = 0.95)
 
  library(ggplot2)
  library(gridExtra)
  ggplot( pred, aes(ScaleLogDepth, predicted)) +
    geom_line( aes(y=predicted), color="blue", linewidth=1 ) +
    geom_ribbon( aes( x=ScaleLogDepth, ymin=conf.low, ymax=conf.high), fill=rgb(0,0,1,0.2) ) +
    facet_wrap(vars(category), scales="free", ncol=2) +
    labs(y="Predicted response")
  
  
######################################
# For S.Brodie
######################################

load('./slope shelf EBS NBS VAST/Gadus macrocephalus//depth3d/fit.RData')

#Get parameters from VAST model object
Gammas <- fit$ParHat[c("gamma1_cp","gamma2_cp")]
#Gammas_1 <- Gammas[[1]] [1,1,1:2]
#gamma1_ctp - effect of density covariates that change over time on encounter prob
#gamma2_ctp - effect of density covariates that change over time on positive catches

df<-fit$covariate_data


Gammas_1 <- Gammas[[1]] #[i,1,1:2]
Gammas_2 <- Gammas[[2]] #[i,1,1:2]
Xlim = range(data_geostat$ScaleLogDepth)
Xcenter = mean(data_geostat$ScaleLogDepth)
X = seq(Xlim[1], Xlim[2], length=1e4)
Y1 = (X-Xcenter)*Gammas_1[1]  + ((X-Xcenter)^2)*Gammas_1[2] #+ ((X-Xcenter)^3)*Gammas_1[3]
plot( x=X, y=exp(Y1-max(Y1)), col="black", type="l", lty="solid", lwd=2, xlab="Depth",ylab="Response", main=paste0("Occurrence "))
Y2 = (X-Xcenter)*Gammas_2[1] + (X-Xcenter)*Gammas_2[2]
plot( x=X, y=exp(Y2-max(Y2)), col=c("black","grey"), type="l", lty="solid", lwd=2, xlab="Depth",ylab="Response", main=paste0("Abundance ",spp.names[i]))


