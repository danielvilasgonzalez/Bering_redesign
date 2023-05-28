####################################################################
####################################################################
##    
##    compare sampling design index estimates with complete index
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
# SCENARIOS from 11_sampling_strata_optimization.R
###################################

#sampling scenarios
samp_df<-expand.grid(strat_var=c('Lat_varTemp','Lat_meanTempF','Depth_meanTempF','Depth_varTemp','meanTempF_varTemp','meanTempF','varTemp','Depth'),
                     target_var=c('sumDensity'), #,'sqsumDensity'
                     n_samples=c(350,500), #c(300,500) 520 (EBS+NBS+CRAB);26 (CRAB); 350 (EBS-CRAB); 494 (NBS-CRAB)
                     n_strata=c(5,10,15)) #c(5,10,15)

#add scenario number
samp_df$samp_scn<-paste0(paste0('scn',1:nrow(samp_df)))

###################################
# SBT projections
###################################

#save SBT table
load('./tables/SBT_projection.RData')#df_sbt

#name scenario
df_sbt$sbt<-paste0('SBT',df_sbt$sbt_n)
df_sbt$sbt2<-paste0(df_sbt$sbt,'_',df_sbt$Scenario)

###################################
# HISTORICAL EVALUATION
###################################

#loop over spp
#for (sp in spp) {
  
  sp<-"Gadus macrocephalus"
  
  #create folder to store results
  dir.create(paste0('./output/species/',sp,'/historical evaluation/'))
  
  #load indices for each sbt and scn
  load(file = paste0("./output/species/",sp,'/historical design-based indices/indices.RData')) #sp_index
  
  #for (scn in dimnames(index_hist)[[5]]) {
    
    #scn<-dimnames(index_hist)[[5]][1]
  
  ############################
  # INDEX
  ############################
  
    ind<-index_hist[,'index',,,]
    ind1<-as.data.frame.table(ind)
    names(ind1)<-c('rep','year','sim','scn','index')
    ind1$year<-gsub('y','',ind1$year)
    
    df<-aggregate(index ~ year + scn,ind1,FUN = function(x) c(mean = mean(x), q95 = quantile(x,probs=0.95) , q5 = quantile(x,probs=0.05)) )
    colnames(df$index)<-c('mean','q95','q5')
    
    # df$ind<-df$index[,'mean']
    # df$q5<-df$index[,'q5']
    # df$q95<-df$index[,'q95']
    
    #check index by scn
    ggplot()+
      geom_line(data=subset(df,scn!='scnbase'),aes(x=year,y=index[,'mean'],color=scn,group =scn),alpha=0.5)+
      geom_line(data=subset(df,scn=='scnbase'),aes(x=year,y=index[,'mean'],group=1),color='black')+
      #geom_line(data=subset(df,scn=='scnbase9'),aes(x=year,y=index[,'mean'],group=1),color='black',linetype='dashed')+
      #geom_line(data=subset(df,scn=='scnbase25'),aes(x=year,y=index[,'mean'],group=1),color='black',linetype='dotted')+
      labs(y='t',x='year')+
      theme_bw()+
      theme(axis.text.x = element_text(angle=90),legend.position = 'none')+ 
      expand_limits(y = 0)
    
     ggplot()+
       geom_ribbon(data=df,aes(x=year,ymax=index[,'q95'],ymin=index[,'q5'],group=scn,fill=scn),alpha=0.5)+
       geom_line(data=df,aes(x=year,y=index[,'mean'],color=scn,group=scn))+
       labs(y='t',x='year')+
       theme_bw()+
       theme(axis.text.x = element_text(angle=90))+ 
       expand_limits(y = 0)+
       facet_wrap(~ scn)

    #save object
    save(df, file = paste0("./output/species/",sp,'/historical evaluation/indices.RData'))
    
    ############################
    # CV
    ############################
     
    ind<-index_hist[,'cv',,,]
    ind1<-as.data.frame.table(ind)
    names(ind1)<-c('rep','year','sim','scn','index')
    ind1$year<-gsub('y','',ind1$year)
    
    df<-aggregate(index ~ year + scn,ind1,FUN = function(x) c(mean = mean(x), q95 = quantile(x,probs=0.95) , q5 = quantile(x,probs=0.05)) )
    colnames(df$index)<-c('mean','q95','q5')
    
    # df$ind<-df$index[,'mean']
    # df$q5<-df$index[,'q5']
    # df$q95<-df$index[,'q95']
    
    df$scn<-factor(df$scn,levels=c('scnbase','scnbase9','scnbase25',paste0('scn',1:48)))
    
    #check cv by scn
    ggplot()+
      geom_line(data=subset(df,scn!='scnbase'),aes(x=year,y=index[,'mean'],color=scn,group =scn),alpha=0.5)+
      geom_line(data=subset(df,scn=='scnbase'),aes(x=year,y=index[,'mean'],group=1),color='black')+
      #geom_line(data=subset(df,scn=='scnbase9'),aes(x=year,y=index[,'mean'],group=1),color='black',linetype='dashed')+
      #geom_line(data=subset(df,scn=='scnbase25'),aes(x=year,y=index[,'mean'],group=1),color='black',linetype='dotted')+
      labs(y='CV',x='year')+
      theme_bw()+
      theme(axis.text.x = element_text(angle=90),legend.position = 'none')+ 
      expand_limits(y = 0)
    
    samp_df$strat_var<-as.character(samp_df$strat_var)
    samp_df1<-rbind(c('baseline',NA,520,15,'scnbase'),c('baseline',NA,511,15,'scnbase9'),c('baseline',NA,495,15,'scnbase25'),samp_df)
    
    df1<-merge(df,samp_df1,by.x='scn',by.y='samp_scn',all.x=TRUE)
    df1$scn<-factor(df1$scn,levels=c('scnbase','scnbase9','scnbase25',paste0('scn',1:48)))
    df1$n_samples<-factor(df1$n_samples,levels=c(520,511,495,350,500))
    df1$strat_var<-gsub('_','\n',df1$strat_var)
      
    ggplot()+
      geom_boxplot(data=df1,aes(x=scn,y=index[,'mean'],fill=scn,group =scn),alpha=0.5)+
      #geom_boxplot(data=subset(df,scn=='scnbase'),aes(x=scn,y=index[,'mean'],group=1),color='black')+
      #geom_line(data=subset(df,scn=='scnbase9'),aes(x=year,y=index[,'mean'],group=1),color='black',linetype='dashed')+
      #geom_line(data=subset(df,scn=='scnbase25'),aes(x=year,y=index[,'mean'],group=1),color='black',linetype='dotted')+
      labs(y='CV',x='sampling designs')+
      theme_bw()+
      theme(axis.text.x = element_text(angle=90),legend.position = 'none',
            panel.spacing = unit(0, "lines"),
            strip.background = element_blank(),
            axis.line = element_line(colour = "grey"),
            panel.grid.major.y =element_line(colour = "grey"),
            strip.text = element_text(size=11),
            strip.placement = "outside",
            #axis.text.x = element_text(angle = 90, hjust = 1),
            panel.background = element_rect(fill = 'white', colour = 'white'))+ 
      expand_limits(y = 0)+
      scale_y_continuous(expand = c(0,0))+
      facet_wrap(vars(n_samples,strat_var), strip.position = "top", scales = "free_x",  nrow=1)
      #coord_cartesian(ylim = c(0,0.1),clip = 'off') + 
      #annotate(geom = "text", x = c('scnbase9','scn15','scn45'), y = -0.035, label = c('baseline','350','500'), size = 4) +
      #scale_x_discrete(labels=c("0.5" = "Dose 0.5", "1" = "Dose 1", "2" = "Dose 2"))
    
    
    #save object
    save(df1, file = paste0("./output/species/",sp,'/historical evaluation/cv.RData'))
    
#}

############################
# RRMSE PBIAS
############################
    
#to store results  
df<-data.frame(matrix(NA,nrow = 0,ncol=5))
names(df)<-c('sim','scn','rep','rrmse','pbias')

#get simulations
ld<-list.dirs(paste0('./output/species/',sp,'/simulated historical data/'),full.names = FALSE,recursive = FALSE)

#loop over simulations
for (sim in ld) {
  
  #sim<-ld[1]
  
  #print scenario to check progress
  cat(paste(" #############  ", sim, " #############\n"))
  
  #load index
  load(paste0("./output/species/",sp,'/simulated historical data/',sim,'/sim_data.RData'))  #sim_data
  
  #get true index
  ind_true<-sim_data$sim_ind[,,1]/1000
  
  #load design-based index
  load(paste0("./output/species/",sp,'/historical design-based indices/indices.RData'))  #sim_data
  #dimnames(index_hist)
  ind_sim<-index_hist[,'index',,sim,]
  ind_sim1<-as.data.frame.table(ind_sim)
  names(ind_sim1)<-c('rep','year','scn','index')
  ind_sim1$year<-gsub('y','y',ind_sim1$year)
  ind_sim1$rep<-as.integer(ind_sim1$rep)
  
  #loop over sampling scenarios
  for (scn in unique(ind_sim1$scn)) {
    
    #scn<-unique(ind_sim1$scn)[1]
    
    #when base scenario there is no replicates
    if (grepl('base',scn)) {
      repls<-1
    } else {
      repls<-1:100
    }  
    
    #loop over replicates
    for (rep in repls) {
      
    #rep<-1
    
    ind_sim2<-ind_sim1[which(ind_sim1$rep==rep & ind_sim1$scn==scn),'index']
      
    #calculate RRMSE
    rmse_i<-sqrt(mean((ind_true - ind_sim2)^2))/mean(ind_true)
    
    #calculate PBIAS
    pbias_i<-mean((ind_true - ind_sim2) / abs(ind_true))
    
    #append results
    df<-rbind(df,data.frame(sim=sim,scn=scn,rep=rep,rrmse=rmse_i,pbias=pbias_i))
    
    }
  }
}

#save object
save(df, file = paste0("./output/species/",sp,'/historical evaluation/rrmse_pbias.RData'))
load(paste0("./output/species/",sp,'/historical evaluation/rrmse_pbias.RData'))

df1<-merge(df,samp_df1,by.x='scn',by.y='samp_scn',all.x=TRUE)
df1$scn<-factor(df1$scn,levels=c('scnbase','scnbase9','scnbase25',paste0('scn',1:48)))
df1$n_samples<-factor(df1$n_samples,levels=c(520,511,495,350,500))
df1$strat_var<-gsub('_','\n',df1$strat_var)

#plot
ggplot()+
  geom_boxplot(data=df1,aes(x=scn,y=rrmse,fill=scn,group =scn),alpha=0.5)+
  #geom_boxplot(data=subset(df,scn=='scnbase'),aes(x=scn,y=index[,'mean'],group=1),color='black')+
  #geom_line(data=subset(df,scn=='scnbase9'),aes(x=year,y=index[,'mean'],group=1),color='black',linetype='dashed')+
  #geom_line(data=subset(df,scn=='scnbase25'),aes(x=year,y=index[,'mean'],group=1),color='black',linetype='dotted')+
  labs(y='RRMSE',x='sampling designs')+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90),legend.position = 'none',
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        panel.grid.major.y =element_line(colour = "grey"),
        strip.text = element_text(size=11),
        strip.placement = "outside",
        #axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = 'white', colour = 'white'))+ 
  expand_limits(y = 0)+
  scale_y_continuous(expand = c(0,0))+
  facet_wrap(vars(n_samples,strat_var), strip.position = "top", scales = "free_x",  nrow=1)
#coord_cartesian(ylim = c(0,0.1),clip = 'off') + 
#annotate(geom = "text", x = c('scnbase9','scn15','scn45'), y = -0.035, label = c('baseline','350','500'), size = 4) +
#scale_x_discrete(labels=c("0.5" = "Dose 0.5", "1" = "Dose 1", "2" = "Dose 2"))


###################################
# PROJECTED EVALUATION
###################################

#loop over spp
#for (sp in spp) {
  
  sp<-"Gadus macrocephalus"
  
  #create folder to store results
  dir.create(paste0('./output/species/',sp,'/projected evaluation/'))
  
  #load indices for each sbt and scn
  load(file = paste0("./output/species/",sp,'/projected design-based indices/indices.RData')) #sp_index
  
  #for (scn in dimnames(index_hist)[[5]]) {
  
  #scn<-dimnames(index_hist)[[5]][1]
  
  ############################
  # INDEX
  ############################
  
  ind<-index_proj[,'index',,,,]
  ind1<-as.data.frame.table(ind)
  names(ind1)<-c('rep','year','sim','scn','sbt','index')
  ind1$year<-gsub('y','',ind1$year)
  
  df<-aggregate(index ~ year + scn +sbt ,ind1,FUN = function(x) c(mean = mean(x), q95 = quantile(x,probs=0.95) , q5 = quantile(x,probs=0.05)) )
  colnames(df$index)<-c('mean','q95','q5')
  
  # df$ind<-df$index[,'mean']
  # df$q5<-df$index[,'q5']
  # df$q95<-df$index[,'q95']
  
  #save SBT table
  load('./tables/SBT_projection.RData')#df_sbt
  
  #name scenario
  df_sbt$sbt<-paste0('SBT',df_sbt$sbt_n)
  df_sbt$sbt2<-paste0(df_sbt$sbt,'_',df_sbt$Scenario)
  df_sbt<-df_sbt[,c('sbt','sbt2')]
  
  df1<-merge(df,df_sbt,by='sbt',all.x=TRUE)
  df1$sbt2<-factor(df1$sbt2,levels = unique(df1$sbt2)[c(1,5:12,2:4)])
  
  #check index by scn
  ggplot()+
    geom_line(data=subset(df1,scn!='scnbase'),aes(x=year,y=index[,'mean'],color=scn,group =scn),alpha=0.5)+
    geom_line(data=subset(df1,scn=='scnbase'),aes(x=year,y=index[,'mean'],group=1),color='black')+
    #geom_line(data=subset(df,scn=='scnbase9'),aes(x=year,y=index[,'mean'],group=1),color='black',linetype='dashed')+
    #geom_line(data=subset(df,scn=='scnbase25'),aes(x=year,y=index[,'mean'],group=1),color='black',linetype='dotted')+
    labs(y='t',x='year')+
    theme_bw()+
    theme(axis.text.x = element_text(angle=90),legend.position = 'none')+ 
    expand_limits(y = 0)+
    facet_wrap(~ sbt2)
  
  
  # ggplot()+
  #   geom_ribbon(data=df,aes(x=year,ymax=index[,'q95'],ymin=index[,'q5'],group=scn,fill=scn),alpha=0.5)+
  #   geom_line(data=df,aes(x=year,y=index[,'mean'],color=scn,group=scn))+
  #   labs(y='t',x='year')+
  #   theme_bw()+
  #   theme(axis.text.x = element_text(angle=90))+ 
  #   expand_limits(y = 0)+
  #   facet_wrap(~ scn)
  
  #save object
  save(df1, file = paste0("./output/species/",sp,'/projected evaluation/indices.RData'))
  
  ############################
  # CV
  ############################
  
  ind<-index_proj[,'cv',,,,]
  ind1<-as.data.frame.table(ind)
  names(ind1)<-c('rep','year','sim','scn','sbt','index')
  ind1$year<-gsub('y','',ind1$year)
 
  df<-aggregate(index ~ year + scn + sbt,ind1,FUN = function(x) c(mean = mean(x), q95 = quantile(x,probs=0.95) , q5 = quantile(x,probs=0.05)) )
  colnames(df$index)<-c('mean','q95','q5')
  
  df1<-merge(df,df_sbt,by='sbt',all.x=TRUE)
  df1$sbt2<-factor(df1$sbt2,levels = unique(df1$sbt2)[c(1,5:12,2:4)])
  
  # df$ind<-df$index[,'mean']
  # df$q5<-df$index[,'q5']
  # df$q95<-df$index[,'q95']
  
  df$scn<-factor(df$scn,levels=c('scnbase','scnbase9','scnbase25',paste0('scn',1:48)))
  
  #check cv by scn
  ggplot()+
    geom_line(data=subset(df1,scn!='scnbase'),aes(x=year,y=index[,'mean'],color=scn,group =scn),alpha=0.5)+
    geom_line(data=subset(df1,scn=='scnbase'),aes(x=year,y=index[,'mean'],group=1),color='black')+
    #geom_line(data=subset(df,scn=='scnbase9'),aes(x=year,y=index[,'mean'],group=1),color='black',linetype='dashed')+
    #geom_line(data=subset(df,scn=='scnbase25'),aes(x=year,y=index[,'mean'],group=1),color='black',linetype='dotted')+
    labs(y='CV',x='year')+
    theme_bw()+
    theme(axis.text.x = element_text(angle=90),legend.position = 'none')+ 
    expand_limits(y = 0)+
    facet_wrap(~sbt2)
  
  samp_df$strat_var<-as.character(samp_df$strat_var)
  samp_df1<-rbind(c('baseline',NA,520,15,'scnbase'),c('baseline',NA,511,15,'scnbase9'),c('baseline',NA,495,15,'scnbase25'),samp_df)
  
  df1<-merge(df,samp_df1,by.x='scn',by.y='samp_scn',all.x=TRUE)
  df1$scn<-factor(df1$scn,levels=c('scnbase','scnbase9','scnbase25',paste0('scn',1:48)))
  df1$n_samples<-factor(df1$n_samples,levels=c(520,511,495,350,500))
  df1$strat_var<-gsub('_','\n',df1$strat_var)
  
  ggplot()+
    geom_boxplot(data=df1,aes(x=scn,y=index[,'mean'],fill=scn,group =scn),alpha=0.5)+
    #geom_boxplot(data=subset(df,scn=='scnbase'),aes(x=scn,y=index[,'mean'],group=1),color='black')+
    #geom_line(data=subset(df,scn=='scnbase9'),aes(x=year,y=index[,'mean'],group=1),color='black',linetype='dashed')+
    #geom_line(data=subset(df,scn=='scnbase25'),aes(x=year,y=index[,'mean'],group=1),color='black',linetype='dotted')+
    labs(y='CV',x='sampling designs')+
    theme_bw()+
    theme(axis.text.x = element_text(angle=90),legend.position = 'none',
          panel.spacing = unit(0, "lines"),
          strip.background = element_blank(),
          axis.line = element_line(colour = "grey"),
          panel.grid.major.y =element_line(colour = "grey"),
          strip.text = element_text(size=11),
          strip.placement = "outside",
          #axis.text.x = element_text(angle = 90, hjust = 1),
          panel.background = element_rect(fill = 'white', colour = 'white'))+ 
    expand_limits(y = 0)+
    scale_y_continuous(expand = c(0,0))+
    facet_wrap(vars(sbt), strip.position = "top", scales = "free_x",  nrow=2)
  #coord_cartesian(ylim = c(0,0.1),clip = 'off') + 
  #annotate(geom = "text", x = c('scnbase9','scn15','scn45'), y = -0.035, label = c('baseline','350','500'), size = 4) +
  #scale_x_discrete(labels=c("0.5" = "Dose 0.5", "1" = "Dose 1", "2" = "Dose 2"))
  
  
  #save object
  save(df1, file = paste0("./output/species/",sp,'/projected evaluation/cv.RData'))
  
#}

############################
# RRMSE PBIAS
############################
  
#to store results  
df<-data.frame(matrix(NA,nrow = 0,ncol=5))
names(df)<-c('sim','scn','rep','rrmse','pbias')

#get simulations
ld<-list.dirs(paste0('./output/species/',sp,'/simulated projected data/'),full.names = FALSE,recursive = FALSE)

#loop over simulations
for (sim in ld) {
  
  #sim<-ld[1]
  
  #print scenario to check progress
  cat(paste(" #############  ", sim, " #############\n"))
  
  #load index
  load(paste0("./output/species/",sp,'/projected historical data/',sim,'/sim_data.RData'))  #sim_data
  
  #get true index
  ind_true<-sim_data$sim_ind[,,3]
  
  #load design-based index
  load(paste0("./output/species/",sp,'/projected design-based indices/indices.RData'))  #sim_data
  #dimnames(index_hist)
  ind_sim<-index_hist[,'index',,sim,]
  ind_sim1<-as.data.frame.table(ind_sim)
  names(ind_sim1)<-c('rep','year','scn','index')
  ind_sim1$year<-gsub('y','y',ind_sim1$year)
  ind_sim1$rep<-as.integer(ind_sim1$rep)
  
  #loop over sampling scenarios
  for (scn in unique(ind_sim1$scn)) {
    
    #scn<-unique(ind_sim1$scn)[49]
    
    #when base scenario there is no replicates
    if (grepl('base',scn)) {
      repls<-1
    } else {
      repls<-1:100
    }  
    
    #loop over replicates
    for (rep in repls) {
      
      #rep<-1
      
      ind_sim2<-ind_sim1[which(ind_sim1$rep==rep & ind_sim1$scn==scn),'index']
      
      #calculate RRMSE
      rmse_i<-sqrt(mean((ind_true - ind_sim2)^2))/mean(ind_true)
      
      #calculate PBIAS
      pbias_i<-mean((ind_true - ind_sim2) / abs(ind_true))
      
      #append results
      df<-rbind(df,data.frame(sim=sim,scn=scn,rep=rep,rrmse=rmse_i,pbias=pbias_i))
      
    }
  }
}

#save object
save(df, file = paste0("./output/species/",sp,'/projected evaluation/rrmse_pbias.RData'))



##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################

  
  
  
  
  
  #load projections
  load( file = paste0("./output/species/",sp,'/fit_projection.RData')) #pr_list
  
  #convert array to df
  df_cv<-as.data.frame.table(sp_index[,'cv',,,])
  names(df_cv)<-c('iter','year','sbt','scn','cv')
  
  #merge to get name sampling scenario
  df_cv1<-merge(df_cv,samp_df,by.x='scn',by.y='samp_scn',all.x = TRUE)
  df_cv1$scn2<-paste0(df_cv1$scn,' - ',df_cv1$strat_var)
  df_cv1$sbt<-factor(df_cv1$sbt,
                     levels = paste0('SBT',1:12))
  
  #merge to get SBT scenario
  df_cv2<-merge(df_cv1,df_sbt[,c("sbt","sbt2")],by = 'sbt')
  
  #sort
  df_cv2$sbt2<-factor(df_cv2$sbt2,levels=unique(df_cv2$sbt2)[c(1,5:12,2:4)])
  
  #plot CV distribution by year/sbt/scn
  ggplot()+
    geom_boxplot(data=df_cv2,aes(x=year,y=cv,color=scn2),position = position_dodge(width = 0.9))+
    #geom_point(data=df_cv,aes(x=year,y=cv,color=scn),position = position_dodge(width = 0.9))+
    facet_wrap(~sbt2)+
    theme_bw()+
    labs(color='sampling designs')
  
  #aggregate annual average
  df_cv3<-aggregate(df_cv2$cv,by=list(year=df_cv2$year,
                                      scn2=df_cv2$scn2,
                                      scn=df_cv2$scn,
                                      sbt2=df_cv2$sbt2,
                                      sbt=df_cv2$sbt),
                    FUN=mean)
  
  #sort to plot sbt correctly
  df_cv3$sbt2<-factor(df_cv3$sbt2,
                      levels = unique(df_cv3$sbt2))
  
  #plot average CV by year/sbt/scn
  ggplot()+
    #geom_linerange(data=diff_index,aes(x=factor(year),ymin=0,ymax=rmse,color=scn),position = position_dodge(width = 0.5))+
    geom_line(data=df_cv3,aes(x=factor(year),y=x,color=scn2,group=scn2))+ #,position = position_dodge(width = 0.5)
    geom_point(data=df_cv3,aes(x=factor(year),y=x,fill=scn2),position = position_dodge(width = 0.5),shape=21,color='black')+
    #geom_hline(yintercept=0,linetype='dashed')+
    #geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn,shape=dummy),color='black',alpha=0.5,position = position_dodge(width = 0.9))+
    facet_wrap(~sbt2)+
    theme_bw()+
    theme()+
    labs(color='sampling design',fill='sampling design',shape='model-based',x='year',y='cv')#+  
  #scale_y_continuous(limits = c(0,25000))
  
  #aggregate to get average by scn/sbt
  df_cv4<-aggregate(df_cv3$x,by=list(
    scn2=df_cv3$scn2,
    scn=df_cv3$scn,
    sbt2=df_cv3$sbt2,
    sbt=df_cv3$sbt),
    FUN=mean)
  
  
  #plot average CV by sbt/scn
  ggplot()+
    #geom_linerange(data=diff_index,aes(x=factor(year),ymin=0,ymax=rmse,color=scn),position = position_dodge(width = 0.5))+
    #geom_line(data=diff_index5,aes(x=scn,y=rmse,color=scn2,group=scn2))+ #,position = position_dodge(width = 0.5)
    geom_point(data=df_cv4,aes(x=scn,y=x,fill=scn2),position = position_dodge(width = 0.5),shape=21,color='black',size=3)+
    #geom_hline(yintercept=0,linetype='dashed')+
    #geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn,shape=dummy),color='black',alpha=0.5,position = position_dodge(width = 0.9))+
    facet_wrap(~sbt2)+
    theme_bw()+
    theme()+
    labs(color='sampling design',fill='sampling design',shape='model-based',x='sampling scn',y='cv')#+  
  #scale_y_continuous(limits = c(0,25000))
  
  #plot average CV by sbt/scn
  ggplot()+
    #geom_linerange(data=diff_index,aes(x=factor(year),ymin=0,ymax=rmse,color=scn),position = position_dodge(width = 0.5))+
    #geom_line(data=diff_index5,aes(x=scn,y=rmse,color=scn2,group=scn2))+ #,position = position_dodge(width = 0.5)
    geom_boxplot(data=df_cv4,aes(x=reorder(scn, x, FUN = mean),y=x,color=scn2),position = position_dodge(width = 0.5),shape=21,size=1,alpha=0.8)+
    #geom_hline(yintercept=0,linetype='dashed')+
    #geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn,shape=dummy),color='black',alpha=0.5,position = position_dodge(width = 0.9))+
    #facet_wrap(~sbt2)+
    theme_bw()+
    theme()+
    labs(color='sampling design',fill='sampling design',x='sampling scn',y='cv')#+  
  #scale_y_continuous(limits = c(0,25000))
  
  
  #to store index
  df_index<-data.frame(matrix(NA,nrow = 0,ncol = 4))
  colnames(df_index)<-c('sbt','scn','year','index')
  diff_index<-df_true<-df_index
  
  #loop over sbt
  for (sbt in names(pr_list)) {
    
    #sbt<-names(pr_list)[1]
    
    #print scenario to check progress
    cat(paste(" #############     PROJECTING    #############\n",
              " #############   Species", sp, match(sp,spp), 'out of',length(spp),  "  #############\n",
              #" #############  Sampling ", samp, " #############\n",
              " #############  ", sbt, " #############\n"))
    
    #get projection
    pr<-pr_list[[sbt]]
    
    #loop over scn 
    for (scn in dimnames(sp_index)[[5]]) {
      
      #scn<-dimnames(sp_index)[[5]][1]
      
      #loop over years
      for (year in paste0('y',2023:2027)) {
        
        #year<-'y2023'
        
        #get index from projection
        df_true<-rbind(df_true,
                       data.frame(sbt=sbt,
                                  scn=scn,
                                  year=1982:2027,
                                  index=drop_units(pr$Index_ctl[1,,1]),row.names = NULL))
        
        #get index from stratified random sampling
        df_index<-rbind(df_index,
                        data.frame(sbt=sbt,
                                   scn=scn,
                                   year=as.integer(gsub('y','',year)),
                                   index=sp_index[,'index',year,sbt,scn]))
        
      }
    }
  }
  
  #check indices for each SBT and scenario
  df_true1<-subset(df_true,year %in% 2023:2027)
  df_true1<-df_true1[!duplicated(df_true1),]
  df_true1$dummy<-''
  
  
  #to sort facets
  df_true1$sbt<-factor(df_true1$sbt,
                       levels = paste0('SBT',1:12))
  df_index$sbt<-factor(df_index$sbt,
                       levels = paste0('SBT',1:12))
  
  #df_true1$sbt<-gsub('scn','',df_true1$sbt)
  #df_index$sbt<-gsub('scn','',df_index$sbt)
  
  df_index<-merge(df_index,df_sbt,by='sbt',all.x=TRUE)
  df_true1<-merge(df_true1,df_sbt,by='sbt',all.x=TRUE)
  
  
  #sort
  df_index$sbt2<-factor(df_index$sbt2,levels=unique(df_index$sbt2)[c(1,5:12,2:4)])  
  df_true1$sbt2<-factor(df_true1$sbt2,levels=unique(df_true1$sbt2)[c(1,5:12,2:4)])
  
  df_index<-merge(df_index,samp_df,by.x='scn',by.y='samp_scn',all.x=TRUE)
  df_true1<-merge(df_true1,samp_df,by.x='scn',by.y='samp_scn',all.x=TRUE)
  df_index$scn_strata<-paste(df_index$scn,df_index$strat_var)
  df_true1$scn_strata<-paste(df_true1$scn,df_true1$strat_var)
  
  #plot true and sampling index for checking
  ggplot()+
    geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn_strata,shape=dummy),color='black',position = position_dodge(width = 0.9),shape=8)+
    geom_boxplot(data=df_index,aes(x=factor(year),y=index,color=scn_strata),position = position_dodge(width = 0.9),alpha=0.1)+
    facet_wrap(~sbt2)+
    theme_bw()+
    theme()+
    labs(color='sampling design',shape='obs',x='year',y='index (t)')#+  
  #scale_y_continuous(limits = c(0,25000))
  
  #aggregate to get mean index by year/scn/sbt
  df_index1<-aggregate(df_index$index,by=list(year=df_index$year,
                                              scn=df_index$scn,
                                              sbt=df_index$sbt,
                                              sbt2=df_index$sbt2),
                       FUN=mean)
  
  #loop over sbt scn and year
  for (sbt in unique(df_index1$sbt)) {
    for (scn in unique(df_index1$scn)) {
      for (y in unique(df_index1$year)) {
        
        #sbt<-unique(df_index1$sbt)[1];scn<-unique(df_index1$scn)[1];y<-unique(df_index1$year)[1]
        
        #true index  
        true_i<-df_true1[which(df_true1$sbt==sbt & df_true1$scn == scn & df_true1$year==y),'index']/1000
        
        #sampling index
        samp_i<-df_index1[which(df_index1$sbt==sbt & df_index1$scn == scn & df_index1$year==y),'x']
        
        #true index for each 100 iterations
        true_ii<-rep(true_i,times=100)
        
        #subset index for year/sbt/scn
        samp_ii<-df_index[which(df_index$sbt==sbt & df_index$scn == scn & df_index$year==y),'index']
        
        #calculate RRMSE
        rmse_i<-sqrt(mean((true_ii - samp_ii)^2))/mean(true_ii)
        
        #calculate PBIAS
        pbias_i<-mean((true_ii - samp_ii) / abs(true_ii))
        
        #store results
        diff_index<-rbind(diff_index,
                          data.frame(sbt=paste0(sbt),
                                     scn=scn,
                                     year=y,
                                     diff=true_i-samp_i,
                                     rmse=rmse_i,
                                     pbias=pbias_i))  
      }
    }
  }
  
  #merge to get sampling scenario names
  diff_index2<-merge(diff_index,samp_df,by.x='scn',by.y='samp_scn',all.x = TRUE)
  
  ###################################
  # PLOTS
  ###################################
  
  #to sort facets
  diff_index$sbt<-factor(diff_index$sbt,
                         levels = paste0('SBT',1:12))
  
  #plot differences between true and sampling index
  ggplot()+
    geom_linerange(data=diff_index,aes(x=factor(year),ymin=0,ymax=diff,color=scn),position = position_dodge(width = 0.5))+
    geom_point(data=diff_index,aes(x=factor(year),y=diff,fill=scn),position = position_dodge(width = 0.5),shape=21,color='black')+
    geom_hline(yintercept=0,linetype='dashed')+
    #geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn,shape=dummy),color='black',alpha=0.5,position = position_dodge(width = 0.9))+
    facet_wrap(~sbt)+
    theme_bw()+
    theme()+
    labs(color='sampling design',fill='sampling design',shape='model-based',x='year',y='mean difference model-design index (t)')#+  
  #scale_y_continuous(limits = c(0,25000))
  
  #get sampling scenario names
  diff_index2$scn2<-paste0(diff_index2$scn,' - ',diff_index2$strat_var)
  
  #merge to get sbt projection names
  diff_index22<-merge(diff_index2,df_sbt[,c("sbt","sbt2")],by = 'sbt')
  
  #sort sbt facet
  diff_index22$sbt2<-factor(diff_index22$sbt2,
                            levels = unique(diff_index22$sbt2)[c(1,5:12,2:4)])
  
  #plot RRMSE by year/sbt/scn
  ggplot()+
    #geom_linerange(data=diff_index,aes(x=factor(year),ymin=0,ymax=rmse,color=scn),position = position_dodge(width = 0.5))+
    geom_line(data=diff_index22,aes(x=factor(year),y=rmse,color=scn2,group=scn2))+ #,position = position_dodge(width = 0.5)
    geom_point(data=diff_index22,aes(x=factor(year),y=rmse,fill=scn2),position = position_dodge(width = 0.5),shape=21,color='black')+
    #geom_hline(yintercept=0,linetype='dashed')+
    #geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn,shape=dummy),color='black',alpha=0.5,position = position_dodge(width = 0.9))+
    facet_wrap(~sbt2)+
    theme_bw()+
    theme()+
    labs(color='sampling design',fill='sampling design',shape='model-based',x='year',y='rmse of index')#+  
  #scale_y_continuous(limits = c(0,25000))
  
  #plot PBIAS by year/sbt/scn
  ggplot()+
    #geom_linerange(data=diff_index,aes(x=factor(year),ymin=0,ymax=rmse,color=scn),position = position_dodge(width = 0.5))+
    geom_line(data=diff_index22,aes(x=factor(year),y=pbias,color=scn2,group=scn2))+ #,position = position_dodge(width = 0.5)
    geom_point(data=diff_index22,aes(x=factor(year),y=pbias,fill=scn2),position = position_dodge(width = 0.5),shape=21,color='black')+
    geom_hline(yintercept=0,linetype='dashed')+
    #geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn,shape=dummy),color='black',alpha=0.5,position = position_dodge(width = 0.9))+
    facet_wrap(~sbt2)+
    theme_bw()+
    theme()+
    labs(color='sampling design',fill='sampling design',shape='model-based',x='year',y='pbias of index')#+  
  #scale_y_continuous(limits = c(0,25000))
  
  #to store index aggregated years
  diff_index3<-data.frame(matrix(NA,nrow = 0,ncol = 3))
  colnames(diff_index3)<-c('sbt','scn','index')
  
  #loop over sbt/scn to get aggregated estimates
  for (sbt in unique(df_index1$sbt)) {
    for (scn in unique(df_index1$scn)) {
      
      #sbt<-unique(df_index1$sbt)[1];scn<-unique(df_index1$scn)[1];#y<-unique(df_index1$year)[1]
      
      #true values
      true_i<-df_true1[which(df_true1$sbt==sbt & df_true1$scn == scn),]#/1000
      true_i<-true_i[order(true_i$year),'index']/1000
      true_ii<-rep(true_i,each=100)
      
      #sampling values
      samp_i<-df_index[which(df_index$sbt==sbt & df_index$scn == scn),]
      samp_i<-samp_i[order(samp_i$year),'index']
      
      #RRMSE
      rmse_i<-sqrt(mean((true_ii - samp_i)^2))/mean(true_ii)
      
      #PBIAS
      pbias_i<-mean((true_ii - samp_i) / abs(true_ii))
      
      #store results
      diff_index3<-rbind(diff_index3,
                         data.frame(sbt=paste0(sbt),
                                    scn=scn,
                                    #year=y,
                                    #diff=true_i-samp_i,
                                    rmse=rmse_i,
                                    pbias=pbias_i))  
      
    }
  }
  
  #merge for getting sampling scenarios and sbt projections names
  diff_index4<-merge(diff_index3,samp_df,by.x='scn',by.y='samp_scn',all.x = TRUE)
  diff_index4$scn2<-paste0(diff_index4$scn,' - ',diff_index4$strat_var)
  diff_index4$sbt<-factor(diff_index4$sbt,
                          levels = paste0('SBT',1:12))
  diff_index5<-merge(diff_index4,df_sbt[,c("sbt","sbt2")],by = 'sbt')
  diff_index5$sbt2<-factor(diff_index5$sbt2,
                           levels = unique(diff_index22$sbt2)[c(1,5:12,2:4)])
  
  #plot RRMSE annual averages by scn/sbt
  ggplot()+
    #geom_linerange(data=diff_index,aes(x=factor(year),ymin=0,ymax=rmse,color=scn),position = position_dodge(width = 0.5))+
    #geom_line(data=diff_index5,aes(x=scn,y=rmse,color=scn2,group=scn2))+ #,position = position_dodge(width = 0.5)
    geom_point(data=diff_index5,aes(x=scn,y=rmse,fill=scn2),position = position_dodge(width = 0.5),shape=21,color='black',size=3)+
    #geom_hline(yintercept=0,linetype='dashed')+
    #geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn,shape=dummy),color='black',alpha=0.5,position = position_dodge(width = 0.9))+
    facet_wrap(~sbt2)+
    theme_bw()+
    theme()+
    labs(color='sampling design',fill='sampling design',shape='model-based',x='sampling scn',y='rmse of index')#+  
  #scale_y_continuous(limits = c(0,25000))
  
  
  ggplot()+
    #geom_linerange(data=diff_index,aes(x=factor(year),ymin=0,ymax=rmse,color=scn),position = position_dodge(width = 0.5))+
    #geom_line(data=diff_index5,aes(x=scn,y=rmse,color=scn2,group=scn2))+ #,position = position_dodge(width = 0.5)
    geom_boxplot(data=diff_index5,aes(x=reorder(scn,rmse,FUN=mean),y=rmse,color=scn2),position = position_dodge(width = 0.5),size=1)+
    #geom_hline(yintercept=0,linetype='dashed')+
    #geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn,shape=dummy),color='black',alpha=0.5,position = position_dodge(width = 0.9))+
    #facet_wrap(~sbt2)+
    theme_bw()+
    theme()+
    labs(color='sampling design',fill='sampling design',shape='model-based',x='sampling scn',y='rmse of index')#+  
  #scale_y_continuous(limits = c(0,25000))
  
  
  
  #plot PBIAS annual averages by scn/sbt
  ggplot()+
    #geom_linerange(data=diff_index,aes(x=factor(year),ymin=0,ymax=rmse,color=scn),position = position_dodge(width = 0.5))+
    #geom_line(data=diff_index5,aes(x=scn,y=rmse,color=scn2,group=scn2))+ #,position = position_dodge(width = 0.5)
    geom_point(data=diff_index5,aes(x=scn,y=pbias,fill=scn2),position = position_dodge(width = 0.5),shape=21,color='black',size=3)+
    geom_hline(yintercept=0,linetype='dashed')+
    #geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn,shape=dummy),color='black',alpha=0.5,position = position_dodge(width = 0.9))+
    facet_wrap(~sbt2)+
    theme_bw()+
    theme()+
    labs(color='sampling design',fill='sampling design',shape='model-based',x='sampling scn',y='pbias')#+  
  #scale_y_continuous(limits = c(0,25000))
  
}


###################################
# PROJECTED EVALUATION
###################################

#loop over spp
for (sp in spp) {
  
  sp<-"Gadus macrocephalus"
  
  #load indices for each sbt and scn
  load(file = paste0("./output/species/",sp,'/historical design-based indices/indices.RData')) #sp_index
  
  #load projections
  load( file = paste0("./output/species/",sp,'/fit_projection.RData')) #pr_list
  
  #convert array to df
  df_cv<-as.data.frame.table(sp_index[,'cv',,,])
  names(df_cv)<-c('iter','year','sbt','scn','cv')
  
  #merge to get name sampling scenario
  df_cv1<-merge(df_cv,samp_df,by.x='scn',by.y='samp_scn',all.x = TRUE)
  df_cv1$scn2<-paste0(df_cv1$scn,' - ',df_cv1$strat_var)
  df_cv1$sbt<-factor(df_cv1$sbt,
                          levels = paste0('SBT',1:12))
  
  #merge to get SBT scenario
  df_cv2<-merge(df_cv1,df_sbt[,c("sbt","sbt2")],by = 'sbt')
  
  #sort
  df_cv2$sbt2<-factor(df_cv2$sbt2,levels=unique(df_cv2$sbt2)[c(1,5:12,2:4)])
  
  #plot CV distribution by year/sbt/scn
  ggplot()+
    geom_boxplot(data=df_cv2,aes(x=year,y=cv,color=scn2),position = position_dodge(width = 0.9))+
    #geom_point(data=df_cv,aes(x=year,y=cv,color=scn),position = position_dodge(width = 0.9))+
    facet_wrap(~sbt2)+
    theme_bw()+
    labs(color='sampling designs')
  
  #aggregate annual average
  df_cv3<-aggregate(df_cv2$cv,by=list(year=df_cv2$year,
                                           scn2=df_cv2$scn2,
                                            scn=df_cv2$scn,
                                           sbt2=df_cv2$sbt2,
                                      sbt=df_cv2$sbt),
                    FUN=mean)
  
  #sort to plot sbt correctly
  df_cv3$sbt2<-factor(df_cv3$sbt2,
                       levels = unique(df_cv3$sbt2))
  
  #plot average CV by year/sbt/scn
  ggplot()+
    #geom_linerange(data=diff_index,aes(x=factor(year),ymin=0,ymax=rmse,color=scn),position = position_dodge(width = 0.5))+
    geom_line(data=df_cv3,aes(x=factor(year),y=x,color=scn2,group=scn2))+ #,position = position_dodge(width = 0.5)
    geom_point(data=df_cv3,aes(x=factor(year),y=x,fill=scn2),position = position_dodge(width = 0.5),shape=21,color='black')+
    #geom_hline(yintercept=0,linetype='dashed')+
    #geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn,shape=dummy),color='black',alpha=0.5,position = position_dodge(width = 0.9))+
    facet_wrap(~sbt2)+
    theme_bw()+
    theme()+
    labs(color='sampling design',fill='sampling design',shape='model-based',x='year',y='cv')#+  
  #scale_y_continuous(limits = c(0,25000))
  
  #aggregate to get average by scn/sbt
  df_cv4<-aggregate(df_cv3$x,by=list(
                                      scn2=df_cv3$scn2,
                                      scn=df_cv3$scn,
                                      sbt2=df_cv3$sbt2,
                                      sbt=df_cv3$sbt),
                    FUN=mean)
  
  
  #plot average CV by sbt/scn
  ggplot()+
    #geom_linerange(data=diff_index,aes(x=factor(year),ymin=0,ymax=rmse,color=scn),position = position_dodge(width = 0.5))+
    #geom_line(data=diff_index5,aes(x=scn,y=rmse,color=scn2,group=scn2))+ #,position = position_dodge(width = 0.5)
    geom_point(data=df_cv4,aes(x=scn,y=x,fill=scn2),position = position_dodge(width = 0.5),shape=21,color='black',size=3)+
    #geom_hline(yintercept=0,linetype='dashed')+
    #geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn,shape=dummy),color='black',alpha=0.5,position = position_dodge(width = 0.9))+
    facet_wrap(~sbt2)+
    theme_bw()+
    theme()+
    labs(color='sampling design',fill='sampling design',shape='model-based',x='sampling scn',y='cv')#+  
  #scale_y_continuous(limits = c(0,25000))
  
  #plot average CV by sbt/scn
  ggplot()+
    #geom_linerange(data=diff_index,aes(x=factor(year),ymin=0,ymax=rmse,color=scn),position = position_dodge(width = 0.5))+
    #geom_line(data=diff_index5,aes(x=scn,y=rmse,color=scn2,group=scn2))+ #,position = position_dodge(width = 0.5)
    geom_boxplot(data=df_cv4,aes(x=reorder(scn, x, FUN = mean),y=x,color=scn2),position = position_dodge(width = 0.5),shape=21,size=1,alpha=0.8)+
    #geom_hline(yintercept=0,linetype='dashed')+
    #geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn,shape=dummy),color='black',alpha=0.5,position = position_dodge(width = 0.9))+
    #facet_wrap(~sbt2)+
    theme_bw()+
    theme()+
    labs(color='sampling design',fill='sampling design',x='sampling scn',y='cv')#+  
  #scale_y_continuous(limits = c(0,25000))
  
  
  #to store index
  df_index<-data.frame(matrix(NA,nrow = 0,ncol = 4))
  colnames(df_index)<-c('sbt','scn','year','index')
  diff_index<-df_true<-df_index
  
  #loop over sbt
  for (sbt in names(pr_list)) {
    
    #sbt<-names(pr_list)[1]
    
    #print scenario to check progress
    cat(paste(" #############     PROJECTING    #############\n",
              " #############   Species", sp, match(sp,spp), 'out of',length(spp),  "  #############\n",
              #" #############  Sampling ", samp, " #############\n",
              " #############  ", sbt, " #############\n"))
  
    #get projection
    pr<-pr_list[[sbt]]
    
    #loop over scn 
    for (scn in dimnames(sp_index)[[5]]) {
      
      #scn<-dimnames(sp_index)[[5]][1]
      
      #loop over years
      for (year in paste0('y',2023:2027)) {
      
        #year<-'y2023'
        
        #get index from projection
        df_true<-rbind(df_true,
                       data.frame(sbt=sbt,
                                  scn=scn,
                                  year=1982:2027,
                                  index=drop_units(pr$Index_ctl[1,,1]),row.names = NULL))
        
        #get index from stratified random sampling
        df_index<-rbind(df_index,
                        data.frame(sbt=sbt,
                                   scn=scn,
                                   year=as.integer(gsub('y','',year)),
                                   index=sp_index[,'index',year,sbt,scn]))

      }
    }
  }
  
  #check indices for each SBT and scenario
  df_true1<-subset(df_true,year %in% 2023:2027)
  df_true1<-df_true1[!duplicated(df_true1),]
  df_true1$dummy<-''
  
  
  #to sort facets
  df_true1$sbt<-factor(df_true1$sbt,
                       levels = paste0('SBT',1:12))
  df_index$sbt<-factor(df_index$sbt,
                       levels = paste0('SBT',1:12))
  
  #df_true1$sbt<-gsub('scn','',df_true1$sbt)
  #df_index$sbt<-gsub('scn','',df_index$sbt)
  
  df_index<-merge(df_index,df_sbt,by='sbt',all.x=TRUE)
  df_true1<-merge(df_true1,df_sbt,by='sbt',all.x=TRUE)
  
  
  #sort
  df_index$sbt2<-factor(df_index$sbt2,levels=unique(df_index$sbt2)[c(1,5:12,2:4)])  
  df_true1$sbt2<-factor(df_true1$sbt2,levels=unique(df_true1$sbt2)[c(1,5:12,2:4)])
  
  df_index<-merge(df_index,samp_df,by.x='scn',by.y='samp_scn',all.x=TRUE)
  df_true1<-merge(df_true1,samp_df,by.x='scn',by.y='samp_scn',all.x=TRUE)
  df_index$scn_strata<-paste(df_index$scn,df_index$strat_var)
  df_true1$scn_strata<-paste(df_true1$scn,df_true1$strat_var)
  
  #plot true and sampling index for checking
  ggplot()+
    geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn_strata,shape=dummy),color='black',position = position_dodge(width = 0.9),shape=8)+
    geom_boxplot(data=df_index,aes(x=factor(year),y=index,color=scn_strata),position = position_dodge(width = 0.9),alpha=0.1)+
    facet_wrap(~sbt2)+
    theme_bw()+
    theme()+
    labs(color='sampling design',shape='obs',x='year',y='index (t)')#+  
    #scale_y_continuous(limits = c(0,25000))
  
  #aggregate to get mean index by year/scn/sbt
  df_index1<-aggregate(df_index$index,by=list(year=df_index$year,
                                   scn=df_index$scn,
                                   sbt=df_index$sbt,
                                   sbt2=df_index$sbt2),
                       FUN=mean)
  
  #loop over sbt scn and year
  for (sbt in unique(df_index1$sbt)) {
    for (scn in unique(df_index1$scn)) {
      for (y in unique(df_index1$year)) {
        
      #sbt<-unique(df_index1$sbt)[1];scn<-unique(df_index1$scn)[1];y<-unique(df_index1$year)[1]
        
      #true index  
      true_i<-df_true1[which(df_true1$sbt==sbt & df_true1$scn == scn & df_true1$year==y),'index']/1000
      
      #sampling index
      samp_i<-df_index1[which(df_index1$sbt==sbt & df_index1$scn == scn & df_index1$year==y),'x']
      
      #true index for each 100 iterations
      true_ii<-rep(true_i,times=100)
      
      #subset index for year/sbt/scn
      samp_ii<-df_index[which(df_index$sbt==sbt & df_index$scn == scn & df_index$year==y),'index']
      
      #calculate RRMSE
      rmse_i<-sqrt(mean((true_ii - samp_ii)^2))/mean(true_ii)
      
      #calculate PBIAS
      pbias_i<-mean((true_ii - samp_ii) / abs(true_ii))
      
      #store results
      diff_index<-rbind(diff_index,
                        data.frame(sbt=paste0(sbt),
                                   scn=scn,
                                   year=y,
                                   diff=true_i-samp_i,
                                   rmse=rmse_i,
                                   pbias=pbias_i))  
      }
    }
  }
  
  #merge to get sampling scenario names
  diff_index2<-merge(diff_index,samp_df,by.x='scn',by.y='samp_scn',all.x = TRUE)
  
  ###################################
  # PLOTS
  ###################################
  
  #to sort facets
  diff_index$sbt<-factor(diff_index$sbt,
                                      levels = paste0('SBT',1:12))
  
  #plot differences between true and sampling index
  ggplot()+
    geom_linerange(data=diff_index,aes(x=factor(year),ymin=0,ymax=diff,color=scn),position = position_dodge(width = 0.5))+
    geom_point(data=diff_index,aes(x=factor(year),y=diff,fill=scn),position = position_dodge(width = 0.5),shape=21,color='black')+
    geom_hline(yintercept=0,linetype='dashed')+
    #geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn,shape=dummy),color='black',alpha=0.5,position = position_dodge(width = 0.9))+
    facet_wrap(~sbt)+
    theme_bw()+
    theme()+
    labs(color='sampling design',fill='sampling design',shape='model-based',x='year',y='mean difference model-design index (t)')#+  
    #scale_y_continuous(limits = c(0,25000))
  
  #get sampling scenario names
  diff_index2$scn2<-paste0(diff_index2$scn,' - ',diff_index2$strat_var)
  
  #merge to get sbt projection names
  diff_index22<-merge(diff_index2,df_sbt[,c("sbt","sbt2")],by = 'sbt')
  
  #sort sbt facet
  diff_index22$sbt2<-factor(diff_index22$sbt2,
                            levels = unique(diff_index22$sbt2)[c(1,5:12,2:4)])
  
  #plot RRMSE by year/sbt/scn
  ggplot()+
    #geom_linerange(data=diff_index,aes(x=factor(year),ymin=0,ymax=rmse,color=scn),position = position_dodge(width = 0.5))+
    geom_line(data=diff_index22,aes(x=factor(year),y=rmse,color=scn2,group=scn2))+ #,position = position_dodge(width = 0.5)
    geom_point(data=diff_index22,aes(x=factor(year),y=rmse,fill=scn2),position = position_dodge(width = 0.5),shape=21,color='black')+
    #geom_hline(yintercept=0,linetype='dashed')+
    #geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn,shape=dummy),color='black',alpha=0.5,position = position_dodge(width = 0.9))+
    facet_wrap(~sbt2)+
    theme_bw()+
    theme()+
    labs(color='sampling design',fill='sampling design',shape='model-based',x='year',y='rmse of index')#+  
  #scale_y_continuous(limits = c(0,25000))
  
  #plot PBIAS by year/sbt/scn
  ggplot()+
    #geom_linerange(data=diff_index,aes(x=factor(year),ymin=0,ymax=rmse,color=scn),position = position_dodge(width = 0.5))+
    geom_line(data=diff_index22,aes(x=factor(year),y=pbias,color=scn2,group=scn2))+ #,position = position_dodge(width = 0.5)
    geom_point(data=diff_index22,aes(x=factor(year),y=pbias,fill=scn2),position = position_dodge(width = 0.5),shape=21,color='black')+
    geom_hline(yintercept=0,linetype='dashed')+
    #geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn,shape=dummy),color='black',alpha=0.5,position = position_dodge(width = 0.9))+
    facet_wrap(~sbt2)+
    theme_bw()+
    theme()+
    labs(color='sampling design',fill='sampling design',shape='model-based',x='year',y='pbias of index')#+  
  #scale_y_continuous(limits = c(0,25000))
  
  #to store index aggregated years
  diff_index3<-data.frame(matrix(NA,nrow = 0,ncol = 3))
  colnames(diff_index3)<-c('sbt','scn','index')
  
  #loop over sbt/scn to get aggregated estimates
  for (sbt in unique(df_index1$sbt)) {
    for (scn in unique(df_index1$scn)) {
            
        #sbt<-unique(df_index1$sbt)[1];scn<-unique(df_index1$scn)[1];#y<-unique(df_index1$year)[1]
        
        #true values
        true_i<-df_true1[which(df_true1$sbt==sbt & df_true1$scn == scn),]#/1000
        true_i<-true_i[order(true_i$year),'index']/1000
        true_ii<-rep(true_i,each=100)
        
        #sampling values
        samp_i<-df_index[which(df_index$sbt==sbt & df_index$scn == scn),]
        samp_i<-samp_i[order(samp_i$year),'index']
        
        #RRMSE
        rmse_i<-sqrt(mean((true_ii - samp_i)^2))/mean(true_ii)
        
        #PBIAS
        pbias_i<-mean((true_ii - samp_i) / abs(true_ii))
        
        #store results
        diff_index3<-rbind(diff_index3,
                          data.frame(sbt=paste0(sbt),
                                     scn=scn,
                                     #year=y,
                                     #diff=true_i-samp_i,
                                     rmse=rmse_i,
                                     pbias=pbias_i))  
        
    }
  }
  
  #merge for getting sampling scenarios and sbt projections names
  diff_index4<-merge(diff_index3,samp_df,by.x='scn',by.y='samp_scn',all.x = TRUE)
  diff_index4$scn2<-paste0(diff_index4$scn,' - ',diff_index4$strat_var)
  diff_index4$sbt<-factor(diff_index4$sbt,
                         levels = paste0('SBT',1:12))
  diff_index5<-merge(diff_index4,df_sbt[,c("sbt","sbt2")],by = 'sbt')
  diff_index5$sbt2<-factor(diff_index5$sbt2,
                            levels = unique(diff_index22$sbt2)[c(1,5:12,2:4)])
  
  #plot RRMSE annual averages by scn/sbt
  ggplot()+
    #geom_linerange(data=diff_index,aes(x=factor(year),ymin=0,ymax=rmse,color=scn),position = position_dodge(width = 0.5))+
    #geom_line(data=diff_index5,aes(x=scn,y=rmse,color=scn2,group=scn2))+ #,position = position_dodge(width = 0.5)
    geom_point(data=diff_index5,aes(x=scn,y=rmse,fill=scn2),position = position_dodge(width = 0.5),shape=21,color='black',size=3)+
    #geom_hline(yintercept=0,linetype='dashed')+
    #geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn,shape=dummy),color='black',alpha=0.5,position = position_dodge(width = 0.9))+
    facet_wrap(~sbt2)+
    theme_bw()+
    theme()+
    labs(color='sampling design',fill='sampling design',shape='model-based',x='sampling scn',y='rmse of index')#+  
  #scale_y_continuous(limits = c(0,25000))
  
  
  ggplot()+
    #geom_linerange(data=diff_index,aes(x=factor(year),ymin=0,ymax=rmse,color=scn),position = position_dodge(width = 0.5))+
    #geom_line(data=diff_index5,aes(x=scn,y=rmse,color=scn2,group=scn2))+ #,position = position_dodge(width = 0.5)
    geom_boxplot(data=diff_index5,aes(x=reorder(scn,rmse,FUN=mean),y=rmse,color=scn2),position = position_dodge(width = 0.5),size=1)+
    #geom_hline(yintercept=0,linetype='dashed')+
    #geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn,shape=dummy),color='black',alpha=0.5,position = position_dodge(width = 0.9))+
    #facet_wrap(~sbt2)+
    theme_bw()+
    theme()+
    labs(color='sampling design',fill='sampling design',shape='model-based',x='sampling scn',y='rmse of index')#+  
  #scale_y_continuous(limits = c(0,25000))
  
  
  
  #plot PBIAS annual averages by scn/sbt
  ggplot()+
    #geom_linerange(data=diff_index,aes(x=factor(year),ymin=0,ymax=rmse,color=scn),position = position_dodge(width = 0.5))+
    #geom_line(data=diff_index5,aes(x=scn,y=rmse,color=scn2,group=scn2))+ #,position = position_dodge(width = 0.5)
    geom_point(data=diff_index5,aes(x=scn,y=pbias,fill=scn2),position = position_dodge(width = 0.5),shape=21,color='black',size=3)+
    geom_hline(yintercept=0,linetype='dashed')+
    #geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn,shape=dummy),color='black',alpha=0.5,position = position_dodge(width = 0.9))+
    facet_wrap(~sbt2)+
    theme_bw()+
    theme()+
    labs(color='sampling design',fill='sampling design',shape='model-based',x='sampling scn',y='pbias')#+  
  #scale_y_continuous(limits = c(0,25000))
  
}
      