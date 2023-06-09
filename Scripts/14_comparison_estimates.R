####################################################################
####################################################################
##    
##    extract estimates (index and CV)
##    compute RRMSE and pbias
##    plot
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

###################################
# grid to get index from simulated densities at cell
###################################

#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
grid$cell<-1:nrow(grid)

###################################
# Sampling designs (from 11_sampling_strata_optimization.R)
###################################

#sampling designs
samp_df<-expand.grid(strat_var=c('Depth_varTemp','varTemp','Depth'),
                     target_var=c('sumDensity'), #,'sqsumDensity'
                     n_samples=c(520), #c(300,500) 520 (EBS+NBS+CRAB);26 (CRAB); 350 (EBS-CRAB); 494 (NBS-CRAB)
                     n_strata=c(15),
                     stringsAsFactors = FALSE) #c(5,10,15)

#add scenario number and baseline designs
samp_df$samp_scn<-paste0(paste0('scn',1:nrow(samp_df)))
samp_df<-rbind(samp_df,c('baseline','current',520,15,'scnbase'),
               c('baseline w/o corner','current',494,15,'scnbase_bis'))

###################################
# SBT scenarios
###################################

#load SBT scenarios table
load('./tables/SBT_projection.RData')#df_sbt

#name SBT scenarios
df_sbt$sbt<-paste0('SBT',df_sbt$sbt_n)
df_sbt$sbt2<-paste0(df_sbt$sbt,'_',df_sbt$Scenario)

###################################
###################################
# HISTORICAL EVALUATION
###################################
###################################

#loop over spp
#for (sp in spp) {
  
  sp<-"Gadus macrocephalus"
  
  #create folder to store results
  dir.create(paste0('./output/species/',sp,'/historical evaluation/'))
  
  #load indices for each sbt and scn
  load(file = paste0("./output/species/",sp,'/historical design-based indices/indices.RData')) #sp_index
  
  ############################
  # INDEX
  ############################
  
  #get indices values
  ind<-index_hist['index',,,,]
  
  #array to dataframe
  ind1<-as.data.frame.table(ind)
  names(ind1)<-c('year','approach','sim','scn','index')
  ind1$year<-gsub('y','',ind1$year)
  
  #aggregate df to get mean, q95 and q5 for each group  
  df<-aggregate(index ~ year + scn + approach,ind1,FUN = function(x) c(mean = mean(x), q95 = quantile(x,probs=0.95) , q5 = quantile(x,probs=0.05)) )
  colnames(df$index)<-c('mean','q95','q5')

  #sort factors for plotting purposes
  df$scn<-factor(df$scn,
                 levels = c('scnbase','scnbase_bis','scn3','scn2','scn1'))
    
  #save plot
  p<-
   ggplot()+
     #geom_line(data=df,aes(x=year,y=index[,'mean'],color=scn,group=interaction(scn,approach),linetype=approach))+
     geom_ribbon(data=df,aes(x=year,ymax=index[,'q95'],ymin=index[,'q5'],group=interaction(scn,approach),fill=scn),alpha=0.1)+
     geom_line(data=df,aes(x=year,y=index[,'mean'],color=scn,group=interaction(scn,approach),linetype=approach),linewidth=0.7)+
     #geom_ribbon(data=subset(df,scenario=='scnbase'),aes(x=year,ymax=index[,'q95'],ymin=index[,'q5'],group=scenario,fill=scenario),alpha=0.3)+
     #geom_line(data=subset(df,scenario=='scnbase'),aes(x=year,y=index[,'mean'],color=scenario,group=scenario,alpha=scenario))+
     #geom_line(data=subset(df,scn=='scnbase'),aes(x=year,y=index[,'mean'],group=scn),color='black',alpha=0.8)+
     #geom_ribbon(data=subset(df,scn=='scnbase'),aes(x=year,ymax=index[,'q95'],ymin=index[,'q5'],group=scn),fill='black',alpha=0.5)+
     labs(y='t',x='')+
     scale_fill_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                       labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
     scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                        labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
     scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
                        labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
     theme_bw()+
     scale_linetype_manual(values = c('current'='solid',
                                      'buffer'='dashed',
                                      'random'='dotted'))+
     scale_y_continuous(expand = c(0,0))+
     theme(axis.text.x = element_text(angle=90,vjust=0.5),panel.grid.minor = element_line(linetype=2,color='grey'))+ 
     expand_limits(y = 0)

  #save plot
  ragg::agg_png(paste0('./figures/species/',sp,'/hist_indices.png'), width = 10, height = 4, units = "in", res = 300)
  p
  dev.off()
      
  #save df results
  save(df, file = paste0("./output/species/",sp,'/historical evaluation/hist_indices.RData'))
    
  ############################
  # CV
  ############################
  
  #get CV   
  ind<-index_hist['cv',,,,]
  
  #array to dataframe
  ind1<-as.data.frame.table(ind)
  names(ind1)<-c('year','approach','sim','scn','index')
  ind1$year<-gsub('y','',ind1$year)
  
  #aggregate df to get mean, q95 and q5 for each group  
  df<-aggregate(index ~ year + scn+ approach,ind1,FUN = function(x) c(mean = mean(x), q95 = quantile(x,probs=0.95) , q5 = quantile(x,probs=0.05)) )
  colnames(df$index)<-c('mean','q95','q5')
  
  #sort for plotting purposes  
  df$scn<-factor(df$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))
  
  #plot  
  p<-
   ggplot()+
     geom_line(data=df,aes(x=year,y=index[,'mean'],color=scn,group=interaction(scn,approach),linetype=approach))+
     #geom_line(data=subset(df,scn=='scnbase'),aes(x=year,y=index[,'mean'],group=1),color='black',linewidth=1.2)+
     #geom_line(data=subset(df,scn=='scnbase9'),aes(x=year,y=index[,'mean'],group=1),color='black',linetype='dashed')+
     #geom_line(data=subset(df,scn=='scnbase25'),aes(x=year,y=index[,'mean'],group=1),color='black',linetype='dotted')+
     labs(y='CV',x='')+
     # scale_fill_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='black','scnbase_bis'='#e15759'),
     #                   labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
     scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                        labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
     # scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
     #                    labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
     theme_bw()+
     scale_linetype_manual(values = c('current'='solid',
                                      'buffer'='dashed',
                                      'random'='dotted'))+
     scale_y_continuous(expand = c(0,0),limits=c(0,0.11))+
     theme(axis.text.x = element_text(angle=90,vjust=0.5),panel.grid.minor = element_line(linetype=2,color='grey'))+ 
     expand_limits(y = 0)
    
  #save plot
  ragg::agg_png(paste0('./figures/species/',sp,'/hist_indices_cv.png'), width = 10, height = 4, units = "in", res = 300)
  p
  dev.off()
    
  #merge results to sampling design table
  df1<-merge(df,samp_df,by.x='scn',by.y='samp_scn',all.x=TRUE)
    
  #sort and corrections for plotting purposes
  df1$scn<-factor(df1$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))
  df1$strat_var<-gsub('_','\n',df1$strat_var)
  df1$strat_var<-factor(df1$strat_var,levels=c('baseline','baseline\nbis','Depth','varTemp','Depth\nvarTemp'))
  df1$approach<-factor(df1$approach,levels=c('current','buffer','random'))
    
  #plot
  p<-
   ggplot()+
     geom_boxplot(data=df1,aes(x=strat_var,y=index[,'mean'],fill=scn,group =interaction(scn,approach),linetype=approach),alpha=0.8)+
     labs(y='CV',x='')+
     theme_bw()+
     theme(panel.grid.minor = element_line(linetype=2,color='grey'))+
     expand_limits(y = 0)+
     scale_fill_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                       labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
     scale_x_discrete(#expand=c(0.1,0.01),
                      labels=c('baseline','baseline w/o corner','depth','var temp','depth + var temp'))+
     scale_linetype_manual(values = c('current'='solid',
                                      'buffer'='dashed',
                                      'random'='dotted'))+
     scale_y_continuous(expand = c(0,0),limits=c(0,0.12))
    
  #save plot
  ragg::agg_png(paste0('./figures/species/',sp,'/hist_indices_cv_box.png'), width = 7, height = 4, units = "in", res = 300)
  p
  dev.off()
    
  #save object
  save(df1, file = paste0("./output/species/",sp,'/historical evaluation/hist_cv.RData'))
    
#}

  ############################
  # RRMSE - PBIAS
  ############################
      
  #to store results  
  df<-data.frame(matrix(NA,nrow = 0,ncol=5))
  names(df)<-c('sim','scn','approach','rrmse','pbias')
  
  #get simulations
  ld<-list.dirs(paste0('./output/species/',sp,'/simulated historical data/'),full.names = FALSE,recursive = FALSE)
  
  #load design-based index
  load(paste0("./output/species/",sp,'/historical design-based indices/indices.RData'))  #sim_data
  
  #loop over simulations (n_sim=100) to calculate RMSE and pbias
  for (sim in ld) {
    
    #sim<-ld[1]
    
    #print scenario to check progress
    cat(paste(" #############  ", sim, " #############\n"))
    
    #load index
    load(paste0("./output/species/",sp,'/simulated historical data/',sim,'/sim_dens.RData'))  #sim_data
    
    #get true index
    ind_true<-(colSums(data.frame(sim_dens) * t(grid$Area_in_survey_km2))/1000)
    
    #get historical index
    ind_sim<-index_hist['index',,,sim,]
    ind_sim1<-as.data.frame.table(ind_sim)
    names(ind_sim1)<-c('year','approach','scn','index')
    ind_sim1$year<-gsub('y','y',ind_sim1$year)
    
    #loop over sampling designs
    for (scn in unique(ind_sim1$scn)) {
      
    #scn<-unique(ind_sim1$scn)[1]
      
      #loop over approaches to allocate samples
      for (apr in c('buffer','current','random')) {
        
      #ss<-'current'
      
      #subset
      ind_sim2<-ind_sim1[which(ind_sim1$approach==apr & ind_sim1$scn==scn),'index']
        
      #calculate RRMSE
      rmse_i<-sqrt(mean((ind_true - ind_sim2)^2))/mean(ind_true)
      
      #calculate PBIAS
      pbias_i<-mean((ind_true - ind_sim2) / abs(ind_true))
      
      #append results
      df<-rbind(df,data.frame(sim=sim,scn=scn,approach=apr,rrmse=rmse_i,pbias=pbias_i))
      
      }
    }
  }
  
  #save object
  save(df, file = paste0("./output/species/",sp,'/historical evaluation/rrmse_pbias.RData'))
  #load(paste0("./output/species/",sp,'/historical evaluation/rrmse_pbias.RData'))
  
  #merge to get sampling scenarios data
  df1<-merge(df,samp_df,by.x='scn',by.y='samp_scn',all.x=TRUE)
  df1$scn<-factor(df1$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))
  df1$approach<-factor(df1$approach,levels=c('current','buffer','random'))
  
  #plot
  p<-
   ggplot()+
     geom_boxplot(data=df1,aes(x=scn,y=rrmse,fill=scn,group =interaction(scn,approach),linetype=approach),alpha=0.8)+
     labs(y='RRMSE',x='')+
     theme_bw()+
     theme(panel.grid.minor = element_line(linetype=2,color='grey'))+
     expand_limits(y = 0)+
     scale_linetype_manual(values = c('current'='solid',
                                      'buffer'='dashed',
                                      'random'='dotted'))+
     scale_fill_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                       labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
     scale_x_discrete(labels=c('baseline','baseline w/o corner','depth','var temp','depth + var temp'))+
     scale_y_continuous(expand = c(0,0),limits=c(0,0.17))
  
  #save plot
  ragg::agg_png(paste0('./figures/species/',sp,'/hist_indices_rmse_box.png'), width = 7, height = 4, units = "in", res = 300)
  p
  dev.off()

###################################
###################################
# PROJECTED EVALUATION
###################################
###################################

#loop over spp
#for (sp in spp) {
  
  sp<-"Gadus macrocephalus"
  
  #create folder to store results
  dir.create(paste0('./output/species/',sp,'/projected evaluation/'))
  
  #load indices for each sbt and scn
  load(file = paste0("./output/species/",sp,'/projected design-based indices/indices.RData')) #sp_index
  
  ############################
  # INDEX
  ############################
  
  #get index
  ind<-index_proj['index',,,,,]/0.1 #missing 10
  ind1<-as.data.frame.table(ind)
  names(ind1)<-c('year','approach','sim','scn','sbt','index')
  ind1$year<-gsub('y','',ind1$year)
  
  #aggregate values by groups
  df<-aggregate(index ~ year + scn +sbt+ approach ,ind1,FUN = function(x) c(mean = mean(x), q95 = quantile(x,probs=0.95) , q5 = quantile(x,probs=0.05)) )
  colnames(df$index)<-c('mean','q95','q5')

  #name SBT scenario
  df_sbt$sbt<-paste0('SBT',df_sbt$sbt_n)
  df_sbt$sbt2<-paste0(df_sbt$sbt,'_',df_sbt$Scenario)

  #get sbt names and sort for plotting purposes  
  df1<-merge(df,df_sbt,by='sbt',all.x=TRUE)
  df1$sbt2<-factor(df1$sbt2,levels = unique(df1$sbt2)[c(1,5:12,2:4)])
  df1$Scenario<-factor(df1$Scenario,levels = unique(df1$Scenario)[c(1,5:12,2:4)])
  
  #if removing cold scenarios
  df2<-df1[which(!grepl('cold',df1$sbt2)),]
  #df2$Scenario<-factor(df2$Scenario,levels = unique(df2$Scenario)[c(1,5:12,2:4)])
  
  #plot
  p<-
   ggplot()+
     geom_ribbon(data=df2,aes(x=year,ymax=index[,'q95'],ymin=index[,'q5'],group=interaction(scn,approach),fill=scn),alpha=0.1)+
     geom_line(data=df2,aes(x=year,y=index[,'mean'],color=scn,group=interaction(scn,approach),linetype=approach),linewidth=0.7)+
     labs(y='t',x='')+
     theme_bw()+
     theme(axis.text.x = element_text(angle=90,vjust=0.5),
           panel.spacing = unit(0, "lines"),
           strip.background = element_rect('white'),
           axis.line = element_line(colour = "grey"),
           panel.grid.major.y =element_line(colour = "grey"),
           strip.text = element_text(size=11),
           panel.grid.minor = element_line(linetype=2,color='grey'))+
     scale_fill_manual(breaks = c('scnbase','scnbase_bis','scn3','scn2','scn1'),
                       labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
                       values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                       name='sampling design')+ #labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
     scale_color_manual(breaks = c('scnbase','scnbase_bis','scn3','scn2','scn1'),
                        labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
                        values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                        name='sampling design')+ #labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
     scale_linetype_manual(values = c('current'='solid',
                                      'buffer'='dashed',
                                      'random'='dotted'))+
     expand_limits(y = 0)+
     scale_y_continuous(expand = c(0.01,0.01))+
     scale_x_discrete(expand = c(0.05,0.05))+
     facet_wrap(~ Scenario,nrow=2,ncol=4)
  
  #save plot
  ragg::agg_png(paste0('./figures/species/',sp,'/proj_indices.png'), width = 12, height = 5, units = "in", res = 300)
  p
  dev.off()
  
  #save df object
  save(df1, file = paste0("./output/species/",sp,'/projected evaluation/indices.RData'))
  
  ############################
  # CV
  ############################
  
  #get CV
  ind<-index_proj['cv',,,,,]
  ind1<-as.data.frame.table(ind)
  names(ind1)<-c('year','approach','sim','scn','sbt','index')
  ind1$year<-gsub('y','',ind1$year)
  
  #aggregate by groups
  df<-aggregate(index ~ year + scn + sbt + approach,ind1,FUN = function(x) c(mean = mean(x), q95 = quantile(x,probs=0.95) , q5 = quantile(x,probs=0.05)) )
  colnames(df$index)<-c('mean','q95','q5')
  
  #sbt name and sort for plotting purposes
  df1<-merge(df,df_sbt,by='sbt',all.x=TRUE)
  df1$sbt2<-factor(df1$sbt2,levels = unique(df1$sbt2)[c(1,5:12,2:4)])
  df1$Scenario<-factor(df1$Scenario,levels = unique(df1$Scenario)[c(1,5:12,2:4)])
  df$scn<-factor(df$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))
  
  #removing cold scenarios
  df2<-df1[which(!grepl('cold',df1$sbt2)),]
  df2$scn<-factor(df2$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))
  
  #plot
  p<-
   ggplot()+
     geom_line(data=df2,aes(x=year,y=index[,'mean'],color=scn,group =interaction(scn,approach),linetype=approach))+
     labs(y='CV',x='')+
     theme_bw()+
     theme(axis.text.x = element_text(angle=90,vjust=0.5),
           panel.spacing = unit(0, "lines"),
           strip.background = element_rect('white'),
           axis.line = element_line(colour = "grey"),
           panel.grid.major.y =element_line(colour = "grey"),
           strip.text = element_text(size=11),
           panel.grid.minor = element_line(linetype=2,color='grey'))+
     scale_color_manual(breaks = c('scnbase','scnbase_bis','scn3','scn2','scn1'),
                        labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
                        values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                        name='sampling design')+ #labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
     scale_linetype_manual(values = c('current'='solid',
                                      'buffer'='dashed',
                                      'random'='dotted'))+
     expand_limits(y = 0)+
     scale_y_continuous(expand = c(0,0.01))+
     scale_x_discrete(expand = c(0.05,0.05))+
     facet_wrap(~ Scenario,nrow=2,ncol=4)
  
  #save plot
  ragg::agg_png(paste0('./figures/species/',sp,'/proj_indices_cv.png'), width = 12, height = 5, units = "in", res = 300)
  p
  dev.off()
  
  #plot
  p<-
   ggplot()+
     geom_boxplot(data=df2,aes(x=scn,y=index[,'mean'],fill=scn,group =interaction(scn,approach),linetype=approach),alpha=0.8)+
     labs(y='CV',x='')+
     theme_bw()+
     theme(panel.grid.minor = element_line(linetype=2,color='grey'))+#,
     scale_fill_manual(breaks = c('scnbase','scnbase_bis','scn3','scn2','scn1'),
                       labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
                       values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                       name='sampling design')+ #labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
     scale_linetype_manual(values = c('current'='solid',
                                      'buffer'='dashed',
                                      'random'='dotted'))+
     expand_limits(y = 0)+
     scale_y_continuous(expand = c(0.001,0.001))+
     scale_x_discrete(expand = c(0.1,0.1),breaks = c('scnbase','scnbase_bis','scn3','scn2','scn1'),
                      labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'))#+
     
  
    #plot with SBT facets
    ggplot()+
    geom_boxplot(data=df2,aes(x=scn,y=index[,'mean'],fill=scn,group =interaction(scn,approach),linetype=approach),alpha=0.8)+
    labs(y='CV',x='')+
    theme_bw()+
      theme(axis.text.x = element_text(angle=90,vjust=0.5),
            panel.spacing = unit(0, "lines"),
            strip.background = element_rect('white'),
            axis.line = element_line(colour = "grey"),
            panel.grid.major.y =element_line(colour = "grey"),
            strip.text = element_text(size=11),
            panel.grid.minor = element_line(linetype=2,color='grey'))+
    scale_fill_manual(breaks = c('scnbase','scnbase_bis','scn3','scn2','scn1'),
                      labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
                      values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                      name='sampling design')+ #labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
    scale_linetype_manual(values = c('current'='solid',
                                     'buffer'='dashed',
                                     'random'='dotted'))+
    expand_limits(y = 0)+
    scale_y_continuous(expand = c(0.001,0.001))+
    scale_x_discrete(expand = c(0.1,0.1),breaks = c('scnbase','scnbase_bis','scn3','scn2','scn1'),
                     labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'))+
    facet_wrap(~ Scenario,nrow=2,ncol=4)
  
  #save plot
  ragg::agg_png(paste0('./figures/species/',sp,'/proj_indices_cv_box.png'), width = 7, height = 4, units = "in", res = 300)
  p
  dev.off()
  
  #save object
  save(df1, file = paste0("./output/species/",sp,'/projected evaluation/cv.RData'))


 ############################
 # RRMSE PBIAS
 ############################
  
 #to store results  
 df<-data.frame(matrix(NA,nrow = 0,ncol=6))
 names(df)<-c('sim','scn','approach','sbt','rrmse','pbias')

 #get simulations
 ld<-list.dirs(paste0('./output/species/',sp,'/simulated projected data/'),full.names = FALSE,recursive = FALSE)

 #loop over simulations
 for (sim in ld) {
  
    #sim<-ld[1]
    
    #print scenario to check progress
    cat(paste(" #############  ", sim, " #############\n"))
    
    for (sbt in paste0('SBT',1:12)) {
    
      
    #sbt<-'SBT2'  
    
    #load index
    load(paste0("./output/species/",sp,'/simulated projected data/',sim,'/sim_data_',sbt,'.RData'))  #sim_data
    
    #get true index
    ind_true<-sim_data$sim_ind/1000
    
    #load design-based index
    load(paste0("./output/species/",sp,'/projected design-based indices/indices.RData'))  #sim_data
    #dimnames(index_hist)
    ind_sim<-index_proj['index',,,sim,,]
    ind_sim1<-as.data.frame.table(ind_sim)
    names(ind_sim1)<-c('year','approach','scn','sbt','index')
    ind_sim1$year<-gsub('y','y',ind_sim1$year)
    
    #loop over sampling scenarios
    for (scn in unique(ind_sim1$scn)) {
      
      for (apr in c('current','buffer','random')) {
        
        #scn<-unique(ind_sim1$scn)[1]
        #ss<-'current'
        
        ind_sim2<-ind_sim1[which(ind_sim1$sbt==sbt & ind_sim1$approach==apr & ind_sim1$scn==scn),'index']
          
        #calculate RRMSE
        rmse_i<-sqrt(mean((ind_true - ind_sim2)^2))/mean(ind_true)
          
        #calculate PBIAS
        pbias_i<-mean((ind_true - ind_sim2) / abs(ind_true))
          
        #append results
        df<-rbind(df,data.frame(sim=sim,scn=scn,approach=apr,sbt=sbt,rrmse=rmse_i,pbias=pbias_i))
      }  
     }
   }
 }


 #save object
 save(df, file = paste0("./output/species/",sp,'/projected evaluation/rrmse_pbias.RData'))
 
 #merge to get names and sort for plotting pruposes
 df1<-merge(df,samp_df,by.x='scn',by.y='samp_scn',all.x=TRUE)
 df1$scn<-factor(df1$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))
 df1$approach<-factor(df1$approach,levels=c('current','buffer','random'))
 df2<-df1[which(df1$sbt %in% paste0('SBT',c(1,3:5,7,9,11,12))),]
  
  
 #plot
 p<-
  ggplot()+
   geom_boxplot(data=df2,aes(x=scn,y=rrmse,fill=scn,group =interaction(scn,approach),linetype=approach),alpha=0.8)+
   labs(y='RRMSE',x='')+
   theme_bw()+
   theme(panel.grid.minor = element_line(linetype=2,color='grey'))+#,
   expand_limits(y = 0)+
   scale_linetype_manual(values = c('current'='solid',
                                    'buffer'='dashed',
                                    'random'='dotted'))+
   scale_fill_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                     labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
   scale_x_discrete(labels=c('baseline','baseline w/o corner','depth','var temp','depth + var temp'))+
   scale_y_continuous(expand = c(0,0),limits=c(0,0.31))#+

 #save plot
 ragg::agg_png(paste0('./figures/species/',sp,'/proj_indices_rmse_box.png'), width = 7, height = 4, units = "in", res = 300)
 p
 dev.off()

#}

      