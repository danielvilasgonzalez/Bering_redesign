####################################################################
####################################################################
##    
##    get estimates for historical and projection simulations
##    extract estimates (index and CV)
##    compute RRMSE
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

#install VAST in case it is not
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

###################################
# grid to get index from simulated densities at cell (available at the github VAST)
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
                     target_var=c('sumDensity'), 
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
df_sbt$old_sbt<-paste0('SBT',row.names(df_sbt))

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
  load(file = paste0("./output/species/",sp,'/historical design-based indices/indices.RData')) #index_hist
  
  ############################
  # INDEX
  ############################
  
  #load fit object to get the true index
  load(paste0('./shelf EBS NBS VAST/',sp,'./fit.RData'))
  
  #true index from VAST (unit t)
  index_true<-fit$Report$Index_ctl[1,,1] #index_true=CPUE_index$true_index[,1]
  
  #time series of index (kg) over replicates (0001:0100), approaches (current or random), and sampling scn (5: scn1,2,3,baseline,baseline_bis)
  ind<-index_hist['index',,,,]
  
  #array to dataframe
  ind1<-as.data.frame.table(ind)
  names(ind1)<-c('year','approach','sim','scn','index')
  ind1$year<-gsub('y','',ind1$year)
  
  #aggregate df to get mean, q95 and q5 for each group (year, sampling scenario and approach)
  df<-aggregate(index ~ year + scn + approach,ind1,FUN = function(x) c(mean = mean(x), q95 = quantile(x,probs=0.95) , q5 = quantile(x,probs=0.05)) )
  colnames(df$index)<-c('mean','q95','q5')

  #sort factors for plotting purposes
  df$scn<-factor(df$scn,
                 levels = c('scnbase','scnbase_bis','scn3','scn2','scn1'))
    
  #get the true index from VAST for plot
  index_true_plot=data.frame('year'=unique(df$year),
                             'index_true'=drop_units(index_true))
  
  #plot
  p<-
   ggplot()+
     geom_ribbon(data=df,aes(x=year,ymax=index[,'q95']/1000,ymin=index[,'q5']/1000,group=interaction(scn,approach),fill=scn),alpha=0.1)+
     geom_line(data=df,aes(x=year,y=index[,'mean']/1000,color=scn,group=interaction(scn,approach),linetype=approach),linewidth=0.7,alpha=0.8)+
     geom_line(data=index_true_plot,aes(x=year,y=index_true/1000,group=1),color='black',linewidth=0.7)+
     labs(y='t',x='')+
     scale_fill_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                       labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
     scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                        labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
     scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
                        labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
     theme_bw()+
     scale_linetype_manual(values = c('current'='solid',
                                      'random'='dashed'),
                           label=c('systematic','random'),
                           name='sample allocation')+
     scale_y_continuous(expand = c(0,0),
                        limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
     theme(axis.text.x = element_text(angle=90,vjust=0.5),panel.grid.minor = element_line(linetype=2,color='grey'))+ 
     expand_limits(y = 0)

  #save plot
  ragg::agg_png(paste0('./figures/species/',sp,'/hist_indices.png'), width = 10, height = 4, units = "in", res = 300)
  p
  dev.off()
      
  #save df results
  save(df, file = paste0("./output/species/",sp,'/historical evaluation/hist_indices.RData'))
    
  ############################
  # CV (coefficient of variance)
  ############################
  
  #get CV for each replicate, approach and sampling scenario
  #CV are calculated on the script #13 as, 
  ### CV <- sqrt(STRS_var) / STRS_mean 
  ### STRS_var<-sum(strs_var, by=year) ### strs_var<-var*area²/samples ### var<-var(CPUE)
  ### STRS_mean<-sum(index_strata, by=year) ### index_strata<-mean_strata*area ### mean_strata<-mean(CPUE)
  ind<-index_hist['cv',,,,]
  
  #array to dataframe
  ind1<-as.data.frame.table(ind)
  names(ind1)<-c('year','approach','sim','scn','index')
  ind1$year<-gsub('y','',ind1$year)
  
  #aggregate df to get mean, q95 and q5 for each group (year, sampling scenario and approach)
  df<-aggregate(index ~ year + scn+ approach,ind1,FUN = function(x) c(mean = mean(x), q95 = quantile(x,probs=0.95) , q5 = quantile(x,probs=0.05)) )
  colnames(df$index)<-c('mean','q95','q5')
  
  #sort for plotting purposes  
  df$scn<-factor(df$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))
  
  #plot  
  p<-
   ggplot()+
     geom_line(data=df,aes(x=year,y=index[,'mean'],color=scn,group=interaction(scn,approach),linetype=approach),linewidth=0.7,alpha=0.8)+
     labs(y='mean CV',x='')+
     scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                        labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
     theme_bw()+
    scale_linetype_manual(values = c('current'='solid',
                                     'random'='dashed'),
                          label=c('systematic','random'),
                          name='sample allocation')+
    
     scale_y_continuous(expand = c(0,0),
                        limits =  c(0, max(df$index[,'mean']) + mean(df$index[,'mean'])/10))+
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
  df1$strat_var<-factor(df1$strat_var,levels=c('baseline','baseline w/o corner','Depth','varTemp','Depth\nvarTemp'))
  df1$approach<-factor(df1$approach,levels=c('current','random'))
    
  #plot
  p<-
   ggplot()+
     geom_boxplot(data=df1,aes(x=strat_var,y=index[,'mean'],fill=scn,group =interaction(scn,approach),linetype=approach),alpha=1)+
     labs(y='CV',x='')+
     theme_bw()+
     theme(panel.grid.minor = element_line(linetype=2,color='grey'),axis.text.x.bottom = element_blank())+
     expand_limits(y = 0)+
      scale_fill_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                        breaks = c('scnbase','scnbase_bis',paste0('scn',3:1)),
                        labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
      scale_x_discrete(#expand=c(0.1,0.01),
                       labels=c('baseline','baseline w/o corner','depth','var temp','depth + var temp'))+
    scale_linetype_manual(values = c('current'='solid',
                                     'random'='dashed'),
                          label=c('systematic','random'),
                          name='sample allocation')+
    scale_y_continuous(expand = c(0,0),
                       limits =  c(0, max(df$index[,'mean']) + mean(df$index[,'mean'])/10))
    
  #save plot
  ragg::agg_png(paste0('./figures/species/',sp,'/hist_indices_cv_box.png'), width = 7, height = 4, units = "in", res = 300)
  p
  dev.off()
    
  #save object
  save(df1, file = paste0("./output/species/",sp,'/historical evaluation/hist_cv.RData'))
    
#}

  ############################
  # RRMSE (relative root mean square error) 
  ############################
  
  #true values CPUE and index from OM
  load(file=paste0('./output/species/',sp,'/optimization data/OM_CPUE_index.RData')) #CPUE_index
  summary(CPUE_index$CPUE)
  summary(CPUE_index$true_index)
  #sd(CPUE_index$EBS_NBS)
  
  #get simulations (as folders, from 0001 to 0100)
  ld<-list.dirs(paste0('./output/species/',sp,'/simulated historical data/'),full.names = FALSE,recursive = FALSE)
  
  #load design-based index time series for each replicate, sampling scenario and approach 
  load(paste0("./output/species/",sp,'/historical design-based indices/indices.RData'))  #index_hist

  #get estimated index SD for each survey across years, sampling scenario and approach
  index_hist1<-index_hist['index',,,,]
  index_sd<-apply(index_hist1,c(1,2,4),sd)
  #index_sd1<-apply(sqrt(index_hist['STRS_var',,,,]),c(1,2,4),mean) #

  #compute the true CV time series for each sampling scenario and approach as EstIndexSD/TrueIndex
  cv_true<-sweep(index_sd,MARGIN =1,STATS = index_true,FUN = '/') 
  
  #get the estimated CV time series for each replicate, sampling scenario and approach
  cv_sim<-index_hist['cv',,,,]
  
  #RRMSE array to store
  rrmse<-array(NA,
               dim=list(41,2,5),
               dimnames=list(paste0('y',1982:2022),dimnames(index_hist)[[3]],dimnames(index_hist)[[5]]))
  
  #loop over sampling scenarios and approaches
  for (scn in dimnames(index_hist)[[5]]) {
    for (apr in dimnames(index_hist)[[3]]) {

     #scn<-'scn1';apr<-'current'  
      
     #get estimated CV time series for each replica  
     cv_sim1<-cv_sim[,apr,,scn]
     #get true CV time series
     cv_true1<- cv_true[,apr,scn]
     
     #for (y in paste0('y',1982:2022)) {
       
     
     #y<-'y1982'
     
     #hist(cv_sim1[y,],nclass=20)
     #abline(v=cv_true1[y])
     #}
     
     #get RRMSE relative to the estimated CV as
     ### for each simulated CV replica difference true CV and squared, 
     rrmse[,apr,scn]<-sqrt(apply(sweep(cv_sim1,MARGIN = 1,STATS = cv_true1,'-')^2,
                  MARGIN = 1,'mean'))/apply(cv_sim1,MARGIN = 1,mean)
       
    }}
  
   #convert array to df
   rrmse<-as.data.frame.table(rrmse)
   names(rrmse)<-c('year','approach','scn','rrmse')
   rrmse$year<-gsub('y','',rrmse$year)
   
   #merge to get sampling scenarios data
   df1<-merge(rrmse,samp_df,by.x='scn',by.y='samp_scn',all.x=TRUE)
   df1$scn<-factor(df1$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))
   df1$approach<-factor(df1$approach,levels=c('current','random'))
   
   #plot boxplot
   p<-
     ggplot()+
     geom_boxplot(data=df1,aes(x=scn,y=rrmse,fill=scn,group =interaction(scn,approach),linetype=approach),alpha=1)+
     labs(y='RRMSE of CV',x='')+
     theme_bw()+
     theme(panel.grid.minor = element_line(linetype=2,color='grey'),axis.text.x.bottom = element_blank())+
     expand_limits(y = 0)+
     #scale_y_continuous(expand = c(0,0.001))+
     scale_fill_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                       labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
     scale_x_discrete(labels=c('baseline','baseline w/o corner','depth','var temp','depth + var temp'))+
       scale_linetype_manual(values = c('current'='solid',
                                        'random'='dashed'),
                             label=c('systematic','random'),
                             name='sample allocation')+
       scale_y_continuous(expand = c(0,0),
                          limits =  c(0, max(df1$rrmse) + mean(df1$rrmse)/10))
     
     #save plot
     ragg::agg_png(paste0('./figures/species/',sp,'/hist_rrmse_cv_box.png'), width = 7, height = 4, units = "in", res = 300)
     p
     dev.off()
   
   #plot time series
   p<-
     ggplot()+
     geom_line(data=df1,aes(x=year,y=rrmse,color=scn,group =interaction(scn,approach),linetype=approach),linewidth=0.7,alpha=0.8)+
     labs(y='RRMSE of CV',x='')+
     theme_bw()+
     theme(panel.grid.minor = element_line(linetype=2,color='grey'))+
       theme(axis.text.x = element_text(angle=90,vjust=0.5),panel.grid.minor = element_line(linetype=2,color='grey'))+ 
     expand_limits(y = 0)+
     scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                       labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
     scale_linetype_manual(values = c('current'='solid',
                                      'random'='dashed'),
                           label=c('systematic','random'),
                           name='sample allocation')+
     scale_y_continuous(expand = c(0,0),
                        limits =  c(0, max(df1$rrmse) + mean(df1$rrmse)/10))

   #save plot
   ragg::agg_png(paste0('./figures/species/',sp,'/hist_rrmse_cv.png'), width = 10, height = 4, units = "in", res = 300)
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
  load(file = paste0("./output/species/",sp,'/projected design-based indices/indices.RData')) #index_proj
  
  ############################
  # INDEX
  ############################
  
  #get index time series by SeaBottomTemp, sampling scenario, replicates and approaches
  ind<-index_proj['index',,,,,]
  ind1<-as.data.frame.table(ind)
  names(ind1)<-c('year','approach','sim','scn','sbt','index')
  ind1$year<-gsub('y','',ind1$year)
  
  
  
  #test plot to get projected design-based indices
  ##the color pattern reflects that the major difference is due to replicates
  ggplot()+
    geom_line(data=ind1,aes(x=year,y=index/1000,group=interaction(scn,approach,sim),color=scn),alpha=0.5)+
    facet_wrap(~sbt)
  
  #aggregate values by groups
  df<-aggregate(index ~ year + scn +sbt+ approach ,ind1,FUN = function(x) c(mean = mean(x), q95 = quantile(x,probs=0.95) , q5 = quantile(x,probs=0.05)) )
  colnames(df$index)<-c('mean','q95','q5')

  #name SBT scenario
  df_sbt$sbt<-paste0('SBT',df_sbt$sbt_n)
  df_sbt$sbt2<-paste0(df_sbt$sbt,'_',df_sbt$Scenario)

  #get sbt names and sort for plotting purposes  
  #df1<-merge(df,df_sbt,by='sbt',all.x=TRUE)
  df1<-merge(df,df_sbt,by.x='sbt',by.y='old_sbt',all.x=TRUE)
  df1$sbt2<-factor(df1$sbt2,levels = unique(df1$sbt2)[c(1,5:12,2:4)])
  df1$Scenario<-factor(df1$Scenario,levels = unique(df1$Scenario)[c(1,5:12,2:4)])
  
  #if removing cold scenarios
  df2<-df1[which(!grepl('cold',df1$sbt2)),]
  df2 <- df2[complete.cases(df2$sbt2), ]
  
  #plot
  p<-
   ggplot()+
     geom_ribbon(data=df2,aes(x=year,ymax=index[,'q95']/1000,ymin=index[,'q5']/1000,group=interaction(scn,approach),fill=scn),alpha=0.1)+
     geom_line(data=df2,aes(x=year,y=index[,'mean']/1000,color=scn,group=interaction(scn,approach),linetype=approach),linewidth=0.7,alpha=0.8)+
     geom_line()+
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
     expand_limits(y = 0)+
    scale_linetype_manual(values = c('current'='solid',
                                     'random'='dashed'),
                          label=c('systematic','random'),
                          name='sample allocation')+
    scale_y_continuous(expand = c(0,0),
                       limits =  c(0, max(df2$index[,'mean']/1000) + mean(df2$index[,'mean'])/1000))+
     scale_x_discrete(expand = c(0.05,0.05))+
     facet_wrap(~ Scenario,nrow=2,ncol=4)
  
  #save plot
  ragg::agg_png(paste0('./figures/species/',sp,'/proj_indices.png'), width = 12, height = 5, units = "in", res = 300)
  p
  dev.off()
  
  #save df object
  save(df1, file = paste0("./output/species/",sp,'/projected evaluation/indices.RData'))
  
  ############################
  # CV (coefficient variacion)
  ############################
  
  #get estimated CV 
  ind<-index_proj['cv',,,,,]
  ind1<-as.data.frame.table(ind)
  names(ind1)<-c('year','approach','sim','scn','sbt','index')
  ind1$year<-gsub('y','',ind1$year)
  
  #aggregate (mean and 95%CI) CV time series by groups (SeaBottomTemp, sampling scenario and approach) 
  df<-aggregate(index ~ year + scn + sbt + approach,ind1,FUN = function(x) c(mean = mean(x), q95 = quantile(x,probs=0.95) , q5 = quantile(x,probs=0.05)) )
  colnames(df$index)<-c('mean','q95','q5')
  
  #sbt name and sort for plotting purposes
  df1<-merge(df,df_sbt,by='sbt',all.x=TRUE)
  df1$sbt2<-factor(df1$sbt2,levels = unique(df1$sbt2)[c(1,5:12,2:4)])
  df1$Scenario<-factor(df1$Scenario,levels = unique(df1$Scenario)[c(1,5:12,2:4)])
  df$scn<-factor(df$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))
  
  #removing cold scenarios
  df2<-df1[which(!grepl('cold',df1$sbt2)),]
  df2 <- df2[complete.cases(df2), ]
  
  df2$scn<-factor(df2$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))
  
  #plot
  p<-
   ggplot()+
     geom_line(data=df2,aes(x=year,y=index[,'mean'],color=scn,group =interaction(scn,approach),linetype=approach),linewidth=0.7,alpha=0.8)+
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
                                     'random'='dashed'),
                          label=c('systematic','random'),
                          name='sample allocation')+
    scale_y_continuous(expand = c(0,0),
                       limits =  c(0, max(df2$index[,'mean']) + mean(df2$index[,'mean'])/10))+
     expand_limits(y = 0)+
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
     theme(panel.grid.minor = element_line(linetype=2,color='grey'),axis.text.x.bottom = element_blank())+#,
     scale_fill_manual(breaks = c('scnbase','scnbase_bis','scn3','scn2','scn1'),
                       labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
                       values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                       name='sampling design')+ #labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
    scale_linetype_manual(values = c('current'='solid',
                                     'random'='dashed'),
                          label=c('systematic','random'),
                          name='sample allocation')+
    scale_y_continuous(expand = c(0,0),
                       limits =  c(0, max(df2$index[,'mean']) + mean(df2$index[,'mean'])/10))+
    expand_limits(y = 0)+
    scale_x_discrete(expand = c(0.1,0.1),breaks = c('scnbase','scnbase_bis','scn3','scn2','scn1'),
                      labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'))#+
     
    #plot with SBT facets
    ggplot()+
    geom_boxplot(data=df2,aes(x=scn,y=index[,'mean'],fill=scn,group =interaction(scn,approach),linetype=approach),alpha=0.8)+
    labs(y='CV',x='')+
    theme_bw()+
      theme(axis.text.x = element_blank(),#element_text(angle=90,vjust=0.5),
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
                                       'random'='dashed'),
                            label=c('systematic','random'),
                            name='sample allocation')+
      scale_y_continuous(expand = c(0,0),
                         limits =  c(0, max(df2$index[,'mean']) + mean(df2$index[,'mean'])/10))+
    expand_limits(y = 0)+
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
  # RRMSE (relative root mean square error) 
  ############################
  
  #get replicates folders (0001 to 0100)
  ld<-list.dirs(paste0('./output/species/',sp,'/simulated projected data/'),full.names = FALSE,recursive = FALSE)

  #array true index 
  index_true<-array(NA,
                    dim=list(5,length(ld),length(dimnames(index_proj)[[6]])),
                    dimnames = list(gsub('y','',dimnames(index_proj)[[2]]),ld,paste0('SBT',1:12)))
  #array true cv
  cv_true<-array(NA,
                    dim=list(5,length(dimnames(index_proj)[[3]]),length(dimnames(index_proj)[[5]]),length(dimnames(index_proj)[[6]])),
                    dimnames = list(gsub('y','',dimnames(index_proj)[[2]]),dimnames(index_proj)[[3]],dimnames(index_proj)[[5]],paste0('SBT',1:12)))
  
  #for each projection and simulation, get the true index
  for (isim in ld) {
      for (sbt in dimnames(index_proj)[[6]]) {
  
        #isim<-ld[1]
        #sbt<-'SBT1'
        
        #load output densities and index from project_model() 
        load(paste0("./output/species/",sp,'/simulated projected data/',isim,'/sim_data_',sbt,'.RData')) 
        
        #get tru index
        true_ind<-sim_data$sim_ind
        #true_ind1<-colSums(sim_data$sim_dens*grid$Area_in_survey_km2)
        
        #store tru index for each projection and simulations
        index_true[,isim,sbt]<-true_ind
      }}
        
  #get the mean true index time series for each projection as the mean over simulations
  index_true1<-apply(index_true,c(1,3),mean)    
  sd_true1<-apply(index_true,c(1,3),sd)   
  
  #get estimated index SD for each survey across years, sampling scenario and approach
  index_proj1<-index_proj['index',,,,,]
  index_sd<-apply(index_proj1,c(1,2,4,5),sd)
  index_hist1<-index_hist['index',,,,]
  index_sd1<-apply(index_hist1,c(1,2,4),sd)
  
  for (scn in dimnames(index_proj)[[5]]) {
    for (apr in dimnames(index_proj)[[3]]) {
      
      hist<-index_sd1[,apr,scn]
      proj<-index_sd[,apr,scn,]
      
      
    }
  }
  
  
  #array true cv
  cv_truee<-array(NA,
                 dim=list(5,length(ld),length(dimnames(index_proj)[[3]]),length(dimnames(index_proj)[[5]]),length(dimnames(index_proj)[[6]])),
                 dimnames = list(gsub('y','',dimnames(index_proj)[[2]]),ld,dimnames(index_proj)[[3]],dimnames(index_proj)[[5]],paste0('SBT',1:12)))
  
  
  
  #loop over sampling scenarios
  for (sbt in dimnames(index_proj)[[6]]) {   
      for (scn in dimnames(index_proj)[[5]]) {
          for (apr in dimnames(index_proj)[[3]]) {
            for (isim in ld) {
            #sd<-get sd for each scenario over replicas
            #sd<-apply(sqrt(index_proj['STRS_var',,apr,,scn,sbt]),1,mean)
            #sd1<-apply(index_proj['index',,apr,,scn,sbt],1,sd)
            #sd2<-sd_true1[,sbt]
             
            sd<-index_sd[,apr,scn,sbt]
            
            #get true CV for each scenario
            cv_truee[,isim,apr,scn,sbt]<-sd/index_true[,isim,sbt]
            }
      }
    }
  }
  
  
  #loop over sampling scenarios
  for (sbt in dimnames(index_proj)[[6]]) {   
    for (scn in dimnames(index_proj)[[5]]) {
      for (apr in dimnames(index_proj)[[3]]) {
        
        #sd<-get sd for each scenario over replicas
        #sd<-apply(sqrt(index_proj['STRS_var',,apr,,scn,sbt]),1,mean)
        #sd1<-apply(index_proj['index',,apr,,scn,sbt],1,sd)
        #sd2<-sd_true1[,sbt]
        
        sd<-index_sd[,apr,scn,sbt]
        
        #get true CV for each scenario
        cv_true[,apr,scn,sbt]<-sd/index_true[,isim,sbt]
      }
    }
  }
  
  
  
  mean(cv_true) #too high
  hist(cv_true)
  #get the estimated CV for each survey
  cv_sim<-index_proj['cv',,,,,]
  mean(cv_sim) #too high  
  
  #RRMSE and MAE array to store
  rrmse<-array(NA,
               dim=list(5,2,5,12),
               dimnames=list(paste0('y',2023:2027),dimnames(index_hist)[[3]],dimnames(index_hist)[[5]],paste0('SBT',1:12)))
  mae<-rrmse
  
  #loop over surveys
  for (scn in dimnames(index_proj)[[5]]) {
    for (apr in dimnames(index_proj)[[3]]) {
      for (sbt in dimnames(index_proj)[[6]]) {
      
      #scn<-dimnames(index_hist)[[5]][1];apr<-dimnames(index_hist)[[3]][1];sbt<-'SBT1'

      #get estimated CV  
      cv_sim1<-cv_sim[,apr,,scn,sbt]
      
      #get true CV
      cv_true1<- cv_true[,apr,scn,sbt]
      cv_truee1<- cv_truee[,,apr,scn,sbt]
      
      cv_sim1-cv_truee1
      
      sqrt(apply((cv_sim1-cv_truee1)^2,
                 MARGIN = 1,mean))/apply(cv_sim1,MARGIN = 1,mean)
      
      #get RRMSE relative to the estimated CV
      rrmse[,apr,scn,sbt]<-sqrt(apply(sweep(cv_sim1,MARGIN = 2,STATS = cv_true1,'-')^2,
                             MARGIN = 1,mean))/apply(cv_sim1,MARGIN = 1,mean)
      #MAE for testing
      mae[,apr,scn,sbt]<-apply(abs(sweep(cv_sim1,MARGIN = 2,STATS = cv_true1,'-')),
                                 MARGIN = 1,mean)/apply(cv_sim1,MARGIN = 1,mean)
    }}}

  mae<-as.data.frame.table(mae)  
  summary(mae)

  rrmse<-as.data.frame.table(rrmse)
  summary(rrmse)
  names(rrmse)<-c('year','approach','scn','sbt','rrmse')
  rrmse$year<-gsub('y','',rrmse$year)

  #sbt name and sort for plotting purposes
  df1<-merge(rrmse,df_sbt,by='sbt',all.x=TRUE)
  df1$sbt2<-factor(df1$sbt2,levels = unique(df1$sbt2)[c(1,5:12,2:4)])
  df1$Scenario<-factor(df1$Scenario,levels = unique(df1$Scenario)[c(1,5:12,2:4)])
  df$scn<-factor(df$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))
  
  #removing cold scenarios
  df2<-df1[which(!grepl('cold',df1$sbt2)),]
  df2 <- df2[complete.cases(df2), ]
  df2$scn<-factor(df2$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))
  
  #plot boxplot distribution
  p<-
  ggplot()+
    geom_boxplot(data=df2,aes(x=scn,y=rrmse,fill=scn,group =interaction(scn,approach),linetype=approach),alpha=1)+
    labs(y='RRMSE of CV',x='')+
    theme_bw()+
    theme_bw()+
    theme(axis.text.x = element_text(),axis.text.x.bottom = element_blank(),
          panel.spacing = unit(0, "lines"),
          strip.background = element_rect('white'),
          axis.line = element_line(colour = "grey"),
          panel.grid.major.y =element_line(colour = "grey"),
          strip.text = element_text(size=11),
          panel.grid.minor = element_line(linetype=2,color='grey'))+
    theme(panel.grid.minor = element_line(linetype=2,color='grey'))+
    expand_limits(y = 0)+
    scale_linetype_manual(values = c('current'='solid',
                                     'random'='dashed'),
                          label=c('systematic','random'),
                          name='sample allocation')+
    scale_y_continuous(expand = c(0,0),
                       limits =  c(0, max(df2$rrmse) + mean(df2$rrmse)/10))+
    scale_fill_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                      labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
    scale_x_discrete(expand = c(0.1,0.1))+
    #scale_x_discrete(labels=gsub(' ','\n',c('baseline','baseline w/o corner','depth','var temp','depth + var temp')))+
    facet_wrap(~ sbt2,nrow=2,ncol=4) #, scales="free"
  
  #save plot
  ragg::agg_png(paste0('./figures/species/',sp,'/proj_indices_rmse_box.png'), width = 12, height = 4, units = "in", res = 300)
  p
  dev.off()
  
  #plot time series
  p<-
  ggplot()+
    geom_line(data=df2,aes(x=year,y=rrmse,color=scn,group =interaction(scn,approach),linetype=approach),linewidth=0.7,alpha=0.8)+
    labs(y='RRMSE of CV',x='')+
    theme_bw()+
    theme(axis.text.x = element_text(angle=90,vjust=0.5),
          panel.spacing = unit(0, "lines"),
          strip.background = element_rect('white'),
          axis.line = element_line(colour = "grey"),
          panel.grid.major.y =element_line(colour = "grey"),
          strip.text = element_text(size=11),
          panel.grid.minor = element_line(linetype=2,color='grey'))+
    theme(panel.grid.minor = element_line(linetype=2,color='grey'))+
    expand_limits(y = 0)+
    scale_linetype_manual(values = c('current'='solid',
                                     'random'='dashed'),
                          label=c('systematic','random'),
                          name='sample allocation')+
    scale_y_continuous(expand = c(0,0),
                       limits =  c(0, max(df2$rrmse) + mean(df2$rrmse)/10))+
    scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                      labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
    scale_x_discrete(expand = c(0.05,0.05))+
    #scale_x_discrete(labels=c('baseline','baseline w/o corner','depth','var temp','depth + var temp'))+
    facet_wrap(~ sbt2,nrow=2,ncol=4) #, scales="free"
  
  #save plot
  ragg::agg_png(paste0('./figures/species/',sp,'/proj_indices_rmse.png'), width = 12, height = 4, units = "in", res = 300)
  p
  dev.off()
  
  