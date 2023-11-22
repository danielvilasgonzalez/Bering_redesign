####################################################################
####################################################################
##    
##    compare and plot design estimates vs true estimates
##    
##    danielvilasgonzalez@gmail.com/dvilasg@uw.edu
##
####################################################################
####################################################################

######################
# SETTINGS
######################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c('raster','units','ggplot2')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#setwd
out_dir<-'C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/'
out_dir<-'/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/'
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
       'Paralithodes camtschaticus',
       'Chionoecetes bairdi')

#remove Anoploma and Reinhardtius because habitat preference reasons
spp<-setdiff(spp, c('Anoplopoma fimbria','Reinhardtius hippoglossoides'))

spp1<-c('Yellowfin sole',
       'Alaska pollock',
       'Pacific cod',
       'Arrowtooth flounder',
       #'Greenland turbot',
       'Northern rock sole',
       'Flathead sole',
       'Alaska plaice',
       'Bering flounder',
       'Arctic cod',
       'Saffon cod',
       #'Sablefish',
       'Snow crab',
       'Blue king crab',
       'Red king crab',
       'Tanner crab')
 
df_spp<-data.frame('spp'=spp,
                   'common'=spp1) 
  
#years
yrs<-c(1982:2019,2021,2022)
n_yrs<-length(yrs)

#how manyt projected years we want
n_proj<-5

#project_yrs
project_yrs<-((yrs[length(yrs)])+1):(yrs[length(yrs)]+n_proj)

#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
grid$cell<-1:nrow(grid)

#load baseline strata data
load('./output/baseline_strata.RData') #baseline_strata

#add percent of total area per strata
baseline_strata$strata_areas$pct<-baseline_strata$strata_areas$Area_in_survey_km2/sum(baseline_strata$strata_areas$Area_in_survey_km2)

###################################
# Sampling designs (from script #11) 
###################################

#sampling scenarios
samp_df<-expand.grid(strat_var=c('Depth_varTemp','varTemp','Depth'),
                     target_var=c('sumDensity'), #,'sqsumDensity'
                     n_samples=c(520), #c(300,500) 520 (EBS+NBS+CRAB);26 (CRAB); 350 (EBS-CRAB); 494 (NBS-CRAB)
                     n_strata=c(15),
                     stringsAsFactors = FALSE) #c(5,10,15)

#add scenario number
samp_df$samp_scn<-paste0(paste0('scn',1:nrow(samp_df)))
samp_df<-rbind(samp_df,c('baseline','systematic',520,15,'scnbase'),
               c('baseline w/o corner','systematic',494,15,'scnbase_bis'))

###################################
# SBT scenarios
###################################

#load SBT scenarios table
load('./tables/SBT_projection.RData')#df_sbt

#name SBT scenarios
df_sbt$sbt<-paste0('SBT',df_sbt$sbt_n)
df_sbt$sbt2<-paste0(df_sbt$sbt,'_',df_sbt$Scenario)

#number of historical simulations and projected simulations
n_sim<- 100

######################
# HISTORICAL INDEX
######################

#INDEX TRUE (MODEL$BASED)
#load(file = paste0("./output/species/dens_index_hist_OM.RData"))  #dens_index_hist_OM, 

#INDEX SIMULATED (DESIGN$BASED)
load(file = paste0("./output/index_hist.RData" )) #index_hist
ind<-index_hist[,'STRS_mean',,,,]

#array to dataframe
ind1<-as.data.frame.table(ind)
names(ind1)<-c('spp','year','approach','sim','scn','index')
ind1$year<-gsub('y','',ind1$year)

#aggregate df to get mean, q95 and q5 for each group (year, sampling scenario and approach)
df<-aggregate(index ~ spp + year + scn + approach,ind1,FUN = function(x) c(mean = mean(x), q95 = quantile(x,probs=0.95) , q5 = quantile(x,probs=0.05)) )
colnames(df$index)<-c('mean','q95','q5')

#sort factors for plotting purposes
df$scn<-factor(df$scn,
               levels = c('scnbase','scnbase_bis','scn3','scn2','scn1'))
df<-merge(df,df_spp,by='spp')
df$year<-as.integer(df$year)
df$com_sci<-paste0(df$common,'\n(',df$spp,')')

#INDEX TRUE (MODEL$BASED)
load(file = paste0("./output/species/dens_index_hist_OM.RData"))  #dens_index_hist_OM, 

true_ind<-data.frame(matrix(NA,nrow = length(yrs),ncol = length(spp)))
rownames(true_ind)<-yrs
colnames(true_ind)<-spp

for (sp in spp) {
  #sp<-spp[1]
  true_ind[,sp]<-drop_units(dens_index_hist_OM[[sp]]$index[,as.character(yrs),1])
}

true_ind$year<-as.character(yrs)
true_ind1<-reshape2::melt(true_ind,id.vars='year')
names(true_ind1)<-c('year','spp','value')
true_ind1<-merge(true_ind1,df_spp,by='spp')
true_ind1$year<-as.integer(true_ind1$year)
true_ind1$dummy<-'true index'
true_ind1$com_sci<-paste0(true_ind1$common,'\n(',true_ind1$spp,')')

#plot
#p<-
  ggplot()+
  geom_ribbon(data=df,aes(x=year,ymax=index[,'q95']/1000,ymin=index[,'q5']/1000,group=interaction(scn,approach,com_sci),fill=scn),alpha=0.1)+
  geom_line(data=df,aes(x=year,y=index[,'mean']/1000,color=scn,group=interaction(scn,approach),linetype=approach),linewidth=0.7,alpha=0.8)+
  geom_point(data=true_ind1,aes(x=year,y=value/1000,group=com_sci,shape=dummy),fill='black',color='black',alpha=0.7,size=1.5)+
  labs(y='t',x='')+
  scale_fill_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                    labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
  scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                     labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
  scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
                     labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
  theme_bw()+
    scale_linetype_manual(values = c('sys'='solid',
                                     'rand'='dashed',
                                     'sb'='dotted'),
                          label=c('systematic','random','spatially-balanced'),
                          name='station allocation')+
  scale_shape_manual(values=c('true index'=4),name='')+
  scale_x_continuous(expand=c(0,0),
                     breaks = c(1985,1990,1995,2000,2005,2010,2015,2020),
                     minor_breaks = setdiff(1982:2022,c(1982,1985,1990,1995,2000,2005,2010,2015,2020,2022)))+
  scale_y_continuous(expand = c(0,0),limits = c(0,NA),labels = scales::comma)+
  #                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
  theme(panel.grid.minor = element_line(linetype=2,color='grey90'),strip.background = element_rect(fill='white'),
        legend.position=c(.92,.30),legend.key.size = unit(20, 'points'),legend.text = element_text(size=10), #legend.position=c(.85,.19)
        legend.title = element_text(size=14),strip.text = element_text(size=12))+ #axis.text.x = element_text(angle=90,vjust=0.5),
  expand_limits(y = 0)+
  facet_wrap(~com_sci,scales='free',dir='v',nrow = 3)
  #pacific cod #facet_wrap(~com_sci,scales='free',dir='v',nrow = 5)

  
  ragg::agg_png(paste0('./figures/ms_hist_indices2.png'), width = 12, height = 15, units = "in", res = 300)
  p
  dev.off()
  
  
######################
# HISTORICAL CV
######################

#get CV for each replicate, approach and sampling scenario
#CV are calculated on the script #13 as, 
### CV <- sqrt(STRS_var) / STRS_mean 
### STRS_var<-sum(strs_var, by=year) ### strs_var<-var*area²/samples ### var<-var(CPUE)
### STRS_mean<-sum(index_strata, by=year) ### index_strata<-mean_strata*area ### mean_strata<-mean(CPUE)
ind<-index_hist[,'CV_sim',,,,]

#array to dataframe
ind1<-as.data.frame.table(ind)
names(ind1)<-c('spp','year','approach','sim','scn','index')
ind1$year<-gsub('y','',ind1$year)

#aggregate df to get mean, q95 and q5 for each group (year, sampling scenario and approach)
df<-aggregate(index ~ spp + year + scn + approach,ind1,FUN = function(x) c(mean = mean(x), q95 = quantile(x,probs=0.95) , q5 = quantile(x,probs=0.05)) )
colnames(df$index)<-c('mean','q95','q5')

#sort factors for plotting purposes
df$scn<-factor(df$scn,
               levels = c('scnbase','scnbase_bis','scn3','scn2','scn1'))
df<-merge(df,df_spp,by='spp')


#plot  
#p<-
df$year<-as.numeric(df$year)
df$com_sci<-paste0(df$common,'\n(',df$spp,')')

ggplot()+
  geom_line(data=df,aes(x=year,y=index[,'mean'],color=scn,group=interaction(scn,approach),linetype=approach),linewidth=0.7,alpha=0.8)+
  labs(y='CV',x='')+
  scale_fill_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                    labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
  scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                     labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
  scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
                     labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
  theme_bw()+
  scale_x_continuous(expand=c(0,0),
                     breaks = c(1985,1990,1995,2000,2005,2010,2015,2020),
                     minor_breaks = setdiff(1982:2022,c(1982,1985,1990,1995,2000,2005,2010,2015,2020,2022)))+
   scale_linetype_manual(values = c('sys'='solid',
                                    'rand'='dashed',
                                    'sb'='dotted'),
                         label=c('systematic','random','spatially-balanced'),
                         name='station allocation')+
  theme(panel.grid.minor = element_line(linetype=2,color='grey90'),strip.background = element_rect(fill='white'),
        legend.position=c(.92,.30),legend.key.size = unit(20, 'points'),legend.text = element_text(size=10), #legend.position=c(.85,.19)
        legend.title = element_text(size=14),strip.text = element_text(size=12))+ #axis.text.x = element_text(angle=90,vjust=0.5),
  expand_limits(y = 0)+
  facet_wrap(~com_sci,scales='free',dir='v',nrow = 3)
  # scale_y_continuous(expand = c(0,0),
  #                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
  #theme(panel.grid.minor = element_line(linetype=2,color='grey'))+ #axis.text.x = element_text(angle=90,vjust=0.5),
  #expand_limits(y = 0)+
  #facet_wrap(~common,scales='free')
  

#save plot
# ragg::agg_png(paste0('./figures/species/',sp,'/hist_indices_cv.png'), width = 10, height = 4, units = "in", res = 300)
# p
# dev.off()

#merge results to sampling design table
df1<-merge(df,samp_df,by.x='scn',by.y='samp_scn',all.x=TRUE)

#sort and corrections for plotting purposes
df1$scn<-factor(df1$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))
df1$strat_var<-gsub('_','\n',df1$strat_var)
df1$strat_var<-factor(df1$strat_var,levels=c('baseline','baseline w/o corner','Depth','varTemp','Depth\nvarTemp'))
df1$approach<-factor(df1$approach,levels=c('sys','rand','sb'))
df1$com_sci<-paste0(df1$common,'\n(',df1$spp,')')

#plot
#p<-
  ggplot()+
  geom_boxplot(data=df1,aes(x=strat_var,y=index[,'mean'],fill=scn,group =interaction(scn,approach),linetype=approach),alpha=1)+
  labs(y='CV',x='')+
    scale_fill_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                      labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
    scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                       labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
    scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
                       labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
    theme_bw()+
    scale_linetype_manual(values = c('sys'='solid',
                                     'rand'='dashed',
                                     'sb'='dotted'),
                          label=c('systematic','random','spatially-balanced'),
                          name='station allocation')+
    scale_shape_manual(values=c('true index'=21),name='')+
    scale_y_continuous(expand = c(0.01,0),limits = c(0,NA))+
    #                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
    theme(panel.grid.minor = element_line(linetype=2,color='grey90'),strip.background = element_rect(fill='white'),
          legend.position=c(.90,.30),legend.key.size = unit(20, 'points'),legend.text = element_text(size=10),
          legend.title = element_text(size=14),strip.text = element_text(size=15),
          axis.text.x = element_blank())+ #axis.text.x = element_text(angle=90,vjust=0.5),
    #expand_limits(y = 0)+
    facet_wrap(~com_sci,scales='free_y',dir='v',nrow = 3)
  
                     #limits =  c(0, max(df$index[,'mean']) + mean(df$index[,'mean'])/10))

#save plot
ragg::agg_png(paste0('./figures/ms_hist_indices_cv_box.png'), width = 12, height = 12, units = "in", res = 300)
p
dev.off()

######################
# RRMSE of CV
######################

#get the estimated CV time series for each replicate, sampling scenario and approach
cv_sim<-index_hist[,'CV_sim',,,,]

#get estimated index SD for each survey across years, sampling scenario and approach
index_hist1<-index_hist[,'STRS_mean',,,,]
index_sd<-
  apply(index_hist1,c(1,2,3,5),sd)
#index_sd1<-apply(sqrt(index_hist['STRS_var',,,,]),c(1,2,4),mean) #

cv_true<-list()

for (sp in spp) {
  #p<-spp[1]
  true_ind2<-true_ind[,sp]
  index_sd1<-index_sd[sp,,,]
  cv_true[[sp]]<-sweep(index_sd1,MARGIN =1,STATS = true_ind2,FUN = '/') 
}

#RRMSE array to store
rrmse<-array(NA,
             dim=list(length(spp),length(yrs),dim(index_hist)[[4]],dim(index_hist)[[6]]),
             dimnames=list(spp,paste0('y',yrs),dimnames(index_hist)[[4]],dimnames(index_hist)[[6]]))

for (sp in spp) {
#loop over species, sampling designs and approaches
for (scn in dimnames(index_hist)[[6]]) {
  for (apr in dimnames(index_hist)[[4]]) {
    
    
    #sp<-spp[1];scn<-'scn1';apr<-'sys'
    
    #scn<-'scn1';apr<-'current'  
    
    #get estimated CV time series for each replica  
    cv_sim1<-cv_sim[sp,,apr,,scn]
    #get true CV time series
    cv_true1<-cv_true[[sp]]
    cv_true2<- cv_true1[,apr,scn]
    
    #for (y in paste0('y',1982:2022)) {
    
    
    #y<-'y1982'
    
    #hist(cv_sim1[y,],nclass=20)
    #abline(v=cv_true1[y])
    #}
    
    #get RRMSE relative to the estimated CV as
    ### for each simulated CV replica difference true CV and squared, 
    rrmse[sp,,apr,scn]<-sqrt(apply(sweep(cv_sim1,MARGIN = 1,STATS = cv_true2,'-')^2,
                                MARGIN = 1,'mean'))/apply(cv_sim1,MARGIN = 1,mean)
    
  }}}

#convert array to df
rrmse<-as.data.frame.table(rrmse)
names(rrmse)<-c('spp','year','approach','scn','rrmse')
rrmse$year<-gsub('y','',rrmse$year)

#merge to get sampling scenarios data
df1<-merge(rrmse,samp_df,by.x='scn',by.y='samp_scn',all.x=TRUE)
df1$scn<-factor(df1$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))
df1$approach<-factor(df1$approach,levels=c('sys','rand','sb'))
df1<-merge(df1,df_spp,by.x='spp',by.y='spp')
df1$com_sci<-paste0(df1$common,'\n(',df1$spp,')')

#plot boxplot
#p<-
  ggplot()+
  geom_boxplot(data=df1,aes(x=scn,y=rrmse,fill=scn,group =interaction(scn,approach,spp),linetype=approach),alpha=1)+
  labs(y='RRMSE of CV',x='')+
    scale_fill_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                      labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
    scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                       labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
   scale_linetype_manual(values = c('sys'='solid',
                                     'rand'='dashed',
                                     'sb'='dotted'),
                          label=c('systematic','random','spatially-balanced'),
                          name='station allocation')+
    scale_shape_manual(values=c('true index'=21),name='')+
    theme_bw()+
    scale_y_continuous(expand = c(0.01,0),limits = c(0,0.5))+
    #                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
    theme(panel.grid.minor = element_line(linetype=2,color='grey90'),strip.background = element_rect(fill='white'),
          legend.position=c(.90,.30),legend.key.size = unit(20, 'points'),legend.text = element_text(size=10),
          legend.title = element_text(size=14),strip.text = element_text(size=15),
          axis.text.x = element_blank())+ #axis.text.x = element_text(angle=90,vjust=0.5),
    #expand_limits(y = 0)+
    facet_wrap(~com_sci,scales='free_y',dir='v',nrow = 3)
  
  
  #save plot210notio
  ragg::agg_png(paste0('./figures/ms_hist_indices_rrmse_box.png'), width = 12, height = 12, units = "in", res = 300)
  p
  dev.off()
  
######################
# PROJECTED INDEX
######################

#INDEX TRUE (MODEL$BASED)
load(file = paste0("./output/species/dens_index_proj_OM.RData"))  #dens_index_hist_OM, 

#INDEX SIMULATED (DESIGN$BASED)
load(file = paste0("./output/species/projected design estimates.RData" )) #index_hist
ind<-index_proj[,'STRS_mean',,,,,,]

#array to dataframe
ind1<-as.data.frame.table(ind)
names(ind1)<-c('spp','year','approach','sim','sbt','scn','index')
ind1$year<-gsub('y','',ind1$year)

#aggregate df to get mean, q95 and q5 for each group (year, sampling scenario and approach)
df<-aggregate(index ~ spp + year + scn + approach + sbt,ind1,FUN = function(x) c(mean = mean(x), q95 = quantile(x,probs=0.95) , q5 = quantile(x,probs=0.05)) )
colnames(df$index)<-c('mean','q95','q5')

#sort factors for plotting purposes
df$scn<-factor(df$scn,
               levels = c('scnbase','scnbase_bis','scn3','scn2','scn1'))
df<-merge(df,df_spp,by='spp')

#INDEX TRUE (MODEL$BASED)
load(file = paste0("./output/species/dens_index_proj_OM.RData"))  #dens_index_hist_OM, 

# true_ind<-data.frame(matrix(NA,nrow = length(yrs),ncol = length(spp)))
# rownames(true_ind)<-yrs
# colnames(true_ind)<-spp



########## calculate tru index based on true densities on projection

dens_index_proj_OM<-list()

for (sp in spp) {
  
  for (sbt in df_sbt$sbt_n) {
    
    cat(paste(" #############   Species", sp, match(sp,spp), 'out of',length(spp),  "  #############\n",
              " #############  SBT", sbt, "of",8, " #############\n"))
    
    #sp<-spp[1];sbt<-df_sbt$sbt_n[1]
    
    load(paste0("./output/species/",sp,'/simulated projected data/fit_projection_SBT',sbt,'.RData'))
    
    #store index and dens
    index<-pm$Index_ctl[,,1]
    dens<-pm$D_gct[,1,]
    
    dens_index_proj_OM[[paste0(sp,"_SBT",sbt)]]<-list('index'=index,'dens'=dens)
  }
}

save(dens_index_proj_OM, file = paste0("./output/species/dens_index_proj_OM.RData")) 
load(paste0("./output/species/dens_index_proj_OM.RData"))

#remove fit
#rm(fit)
#dyn.unload('C:/Program Files/R/R-4.2.2/library/VAST/executables/VAST_v13_1_0_TMBad.dll')

true_ind<-array(NA,
                dim=c(n_proj,length(spp),8),
                dimnames = list(2023:2027,spp,df_sbt$sbt_n))

for (sp in spp) {
  
  for (sbt in df_sbt$sbt_n) {
    
    #dens_index_proj_OM
    
    x<-drop_units(dens_index_proj_OM[[paste0(sp,'_SBT',sbt)]]$index)
    
    #sp<-spp[1]
    true_ind[,sp,sbt]<-x[(length(x)-n_proj+1):length(x)]
    
    
  }

}

true_ind1<-as.data.frame.table(true_ind)
names(true_ind1)<-c('year','spp','sbt','value')
true_ind1<-merge(true_ind1,df_spp,by='spp')
true_ind1$year<-as.integer(as.character(true_ind1$year))
true_ind1$dummy<-'true index'
true_ind1$com_sci<-paste0(true_ind1$common,'\n(',true_ind1$spp,')')
true_ind1$sbt<-paste0('SBT',true_ind1$sbt)

# true_ind$year<-as.character(yrs)
# true_ind1<-reshape2::melt(true_ind,id.vars='year')
# names(true_ind1)<-c('year','spp','value')
# true_ind1<-merge(true_ind1,df_spp,by='spp')


#INDEX SIMULATED (DESIGN$BASED)
load(file = paste0("./output/species/projected design estimates.RData" )) #index_proj
ind<-index_proj[,'STRS_mean',,,,,,]

#array to dataframe
ind1<-as.data.frame.table(ind)
names(ind1)<-c('spp','year','approach','sim','sbt','scn','index')
ind1$year<-gsub('y','',ind1$year)

#aggregate df to get mean, q95 and q5 for each group (year, sampling scenario and approach)
df<-aggregate(index ~ spp + year + scn + approach + sbt,ind1,FUN = function(x) c(mean = mean(x), q95 = quantile(x,probs=0.95) , q5 = quantile(x,probs=0.05)) )
colnames(df$index)<-c('mean','q95','q5')

#sort factors for plotting purposes
df$scn<-factor(df$scn,
               levels = c('scnbase','scnbase_bis','scn3','scn2','scn1'))
df<-merge(df,df_spp,by='spp')
df$year<-as.integer(df$year)
df$com_sci<-paste0(df$common,'\n(',df$spp,')')

#for (sp in spp) {
  
sp<-'Gadus macrocephalus'
df1<-subset(df,spp==sp)
true_ind11<-subset(true_ind1,spp==sp)

sel_sbt<-paste0('SBT',c(3,6,8))

true_ind111<-subset(true_ind11,sbt %in% sel_sbt)
df11<-subset(df1,sbt %in% sel_sbt)

#plot
p<-
  ggplot()+
  geom_ribbon(data=df11,aes(x=year,ymax=index[,'q95']/1000,ymin=index[,'q5']/1000,group=interaction(scn,approach,com_sci,sbt),fill=scn),alpha=0.1)+
  geom_line(data=df11,aes(x=year,y=index[,'mean']/1000,color=scn,group=interaction(scn,approach),linetype=approach),linewidth=0.7,alpha=0.8)+
  geom_point(data=true_ind111,aes(x=year,y=value/1000,group=interaction(com_sci,sbt),shape=dummy),fill='black',color='black',alpha=0.7,size=1.5)+
  labs(y='t',x='')+
  scale_fill_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                    labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
  scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                     labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
  scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
                     labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
  theme_bw()+
  scale_linetype_manual(values = c('systematic'='solid',
                                   'random'='dashed'),
                        label=c('systematic','random'),
                        name='sample allocation')+
  scale_shape_manual(values=c('true index'=4),name='')+
  scale_x_continuous(expand=c(0,0))+#,
                     # breaks = c(1985,1990,1995,2000,2005,2010,2015,2020),
                     # minor_breaks = setdiff(1982:2022,c(1982,1985,1990,1995,2000,2005,2010,2015,2020,2022)))+
  scale_y_continuous(expand = c(0,0),limits = c(0,NA),labels = scales::comma)+
  #                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
  theme(panel.grid.minor = element_line(linetype=2,color='grey90'),strip.background = element_rect(fill='white'),plot.margin =margin(1,15,1,1),
        panel.grid.minor.x = element_blank(),#aspect.ratio = 1,
        legend.position = 'none',#legend.position=c(.85,.19),
        legend.key.size = unit(20, 'points'),legend.text = element_text(size=10),
        legend.title = element_text(size=14),strip.text = element_text(size=12))+ #axis.text.x = element_text(angle=90,vjust=0.5),
  expand_limits(y = 0)+
  facet_wrap(~com_sci+sbt,scales='free',dir='v',nrow = 5)

#save plot
ragg::agg_png(paste0('./figures/example_proj_indices2.png'), width = 4, height = 8, units = "in", res = 300)
p
dev.off()



  ######################
  # PROJECTED CV
  ######################
  
  #get CV for each replicate, approach and sampling scenario
  #CV are calculated on the script #13 as, 
  ### CV <- sqrt(STRS_var) / STRS_mean 
  ### STRS_var<-sum(strs_var, by=year) ### strs_var<-var*area²/samples ### var<-var(CPUE)
  ### STRS_mean<-sum(index_strata, by=year) ### index_strata<-mean_strata*area ### mean_strata<-mean(CPUE)
  ind<-index_proj[,'CV_sim',,,,,,]
  
  #array to dataframe
  ind1<-as.data.frame.table(ind)
  names(ind1)<-c('spp','year','approach','sim','sbt','scn','index')
  ind1$year<-gsub('y','',ind1$year)
  
  #aggregate df to get mean, q95 and q5 for each group (year, sampling scenario and approach)
  df<-aggregate(index ~ spp + year + scn + approach + sbt,ind1,FUN = function(x) c(mean = mean(x), q95 = quantile(x,probs=0.95) , q5 = quantile(x,probs=0.05)) )
  colnames(df$index)<-c('mean','q95','q5')
  
  #sort factors for plotting purposes
  df$scn<-factor(df$scn,
                 levels = c('scnbase','scnbase_bis','scn3','scn2','scn1'))
  df<-merge(df,df_spp,by='spp')
  
  
  #plot  
  #p<-
  
  ggplot()+
    geom_line(data=df,aes(x=year,y=index[,'mean'],color=scn,group=interaction(scn,approach,sbt),linetype=approach),linewidth=0.7,alpha=0.8)+
    labs(y='CV',x='')+
    scale_fill_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                      labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
    scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                       labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
    scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
                       labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
    theme_bw()+
    scale_linetype_manual(values = c('systematic'='solid',
                                     'random'='dashed'),
                          label=c('systematic','random'),
                          name='sample allocation')+
    # scale_y_continuous(expand = c(0,0),
    #                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
    theme(axis.text.x = element_text(angle=90,vjust=0.5),panel.grid.minor = element_line(linetype=2,color='grey'))+ 
    expand_limits(y = 0)+
    facet_wrap(~common,scales='free')
  
  
    #save plot
  # ragg::agg_png(paste0('./figures/species/',sp,'/hist_indices_cv.png'), width = 10, height = 4, units = "in", res = 300)
  # p
  # dev.off()
  
  #merge results to sampling design table
  df1<-merge(df,samp_df,by.x='scn',by.y='samp_scn',all.x=TRUE)
  
  #sort and corrections for plotting purposes
  df1$scn<-factor(df1$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))
  df1$strat_var<-gsub('_','\n',df1$strat_var)
  df1$strat_var<-factor(df1$strat_var,levels=c('baseline','baseline w/o corner','Depth','varTemp','Depth\nvarTemp'))
  df1$approach<-factor(df1$approach,levels=c('systematic','random'))
  df1$com_sci<-paste0(df1$common,'\n(',df1$spp,')')
  
  #plot
  #p<-
  df2<-subset(df1,spp==sp & sbt %in% paste0('SBT',c(3,6,8)))
  
  p<-
    ggplot()+
    geom_boxplot(data=df2,aes(x=strat_var,y=index[,'mean'],fill=scn,group =interaction(scn,approach),linetype=approach),alpha=1)+
    labs(y='CV',x='')+
    scale_fill_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                      labels = c('baseline','baseline w/o corner
                                 ','depth','var temp','depth + var temp'),name='sampling design')+
    scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                       labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
    scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
                       labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
    theme_bw()+
    scale_linetype_manual(values = c('systematic'='solid',
                                     'random'='dashed'),
                          label=c('systematic','random'),
                          name='sample allocation')+
    scale_shape_manual(values=c('true index'=21),name='')+
    scale_y_continuous(expand = c(0.01,0),limits = c(0,NA))+
    #                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
    theme(panel.grid.minor = element_line(linetype=2,color='grey90'),strip.background = element_rect(fill='white'),
         legend.key.size = unit(20, 'points'),legend.text = element_text(size=10),
         legend.position = 'none',
          legend.title = element_text(size=14),strip.text = element_text(size=15), #legend.position='none',#c(.85,.19),
          axis.text.x = element_blank())+ #axis.text.x = element_text(angle=90,vjust=0.5),
    #expand_limits(y = 0)+
    facet_wrap(~com_sci+sbt,dir='v',nrow = 3)
  
  #limits =  c(0, max(df$index[,'mean']) + mean(df$index[,'mean'])/10))
  
  #save plot
  ragg::agg_png(paste0('./figures/example_proj_indices_cv_box.png'), width = 4, height = 8, units = "in", res = 300)
  p
  dev.off()
  
  
  
  
  ######################
  # RRMSE of CV
  ######################
  
  #get the estimated CV time series for each replicate, sampling scenario and approach
  cv_sim<-index_proj[,'CV_sim',,,,,,]
  
  #get estimated index SD for each survey across years, sampling scenario and approach
  index_proj1<-index_proj[,'STRS_mean',,,,,,]
  index_sd<-
    apply(index_proj1,c(1,2,3,5,6),sd)
  #index_sd1<-apply(sqrt(index_hist['STRS_var',,,,]),c(1,2,4),mean) #
  
  cv_true<-array(NA,
                 dim=c(length(spp),length(dimnames(cv_sim)[[5]])),
                 dimnames = list(spp,dimnames(cv_sim)[[5]]))
  
  #RRMSE array to store
  cv_true<-array(NA,
               dim=list(length(spp),length(project_yrs),length(dimnames(cv_sim)[[5]]),2,5),
               dimnames=list(spp,paste0('y',project_yrs),dimnames(cv_sim)[[5]],dimnames(index_proj)[[4]],dimnames(index_proj)[[8]]))
  
  for (sp in spp) {
    for (sbt in dimnames(cv_sim)[[5]]) {
      
    #p<-spp[1]
      sbt1<-gsub('SBT','',sbt)
    true_ind2<-true_ind[,sp,sbt1]
    index_sd1<-index_sd[sp,,,sbt,]
    cv_true[sp,,as.character(sbt),,]<-sweep(index_sd1,MARGIN =1,STATS = true_ind2,FUN = '/') 
    }
  }
  
  #RRMSE array to store
  rrmse<-array(NA,
               dim=list(length(spp),length(project_yrs),length(dimnames(cv_sim)[[5]]),2,5),
               dimnames=list(spp,paste0('y',project_yrs),dimnames(cv_sim)[[5]],dimnames(index_proj)[[4]],dimnames(index_proj)[[8]]))

    for (sp in spp) {
      for (sbt in dimnames(cv_sim)[[5]]) {
    #loop over species, sampling designs and approaches
    for (scn in dimnames(index_proj)[[8]]) {
      for (apr in dimnames(index_proj)[[4]]) {
        
        
        #sp<-spp[1];scn<-'scn1';apr<-'systematic'
        
        #scn<-'scn1';apr<-'current'  
        
        #get estimated CV time series for each replica  
        cv_sim1<-cv_sim[sp,,apr,,sbt,scn]
        #get true CV time series
        cv_true1<-cv_true[sp,,as.character(sbt),apr,scn]
        cv_true2<- cv_true1
        
        #for (y in paste0('y',1982:2022)) {
        
        
        #y<-'y1982'
        
        #hist(cv_sim1[y,],nclass=20)
        #abline(v=cv_true1[y])
        #}
        
        #get RRMSE relative to the estimated CV as
        ### for each simulated CV replica difference true CV and squared, 
        rrmse[sp,,sbt,apr,scn]<-sqrt(apply(sweep(cv_sim1,MARGIN = 1,STATS = cv_true2,'-')^2,
                                       MARGIN = 1,'mean'))/apply(cv_sim1,MARGIN = 1,mean)
        
      }}}}
  
  #convert array to df
  rrmse<-as.data.frame.table(rrmse)
  names(rrmse)<-c('spp','year','sbt','approach','scn','rrmse')
  rrmse$year<-gsub('y','',rrmse$year)
  
  #merge to get sampling scenarios data
  df1<-merge(rrmse,samp_df,by.x='scn',by.y='samp_scn',all.x=TRUE)
  df1$scn<-factor(df1$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))
  df1$approach<-factor(df1$approach,levels=c('systematic','random'))
  df1<-merge(df1,df_spp,by.x='spp',by.y='spp')
  df1$com_sci<-paste0(df1$common,'\n(',df1$spp,')')
  
  df2<-subset(df1,spp=='Gadus macrocephalus' & sbt %in% sel_sbt)
  
  
  
  #plot boxplot
  p<-
    ggplot()+
    geom_boxplot(data=df2,aes(x=scn,y=rrmse,fill=scn,group =interaction(scn,approach,spp,sbt),linetype=approach),alpha=1)+
    labs(y='RRMSE of CV',x='')+
    scale_fill_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
                      labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
    # scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
    #                    labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
    scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
                       labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
    theme_bw()+
    scale_linetype_manual(values = c('systematic'='solid',
                                     'random'='dashed'),
                          label=c('systematic','random'),
                          name='sample allocation')+
    scale_shape_manual(values=c('true index'=21),name='')+
    scale_y_continuous(expand = c(0.01,0),limits = c(0,NA))+
    #                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
    theme(panel.grid.minor = element_line(linetype=2,color='grey90'),strip.background = element_rect(fill='white'),
          legend.position = 'none',#legend.position=c(.85,.19),
          legend.key.size = unit(20, 'points'),legend.text = element_text(size=10),
          legend.title = element_text(size=14),strip.text = element_text(size=15),
          axis.text.x = element_blank())+ #axis.text.x = element_text(angle=90,vjust=0.5),
    #expand_limits(y = 0)+
    facet_wrap(~com_sci+sbt,scales='free_y',dir='v',nrow = 5)
  
  
  #save plot
  ragg::agg_png(paste0('./figures/example_proj_indices_rrmse_box.png'), width = 4, height = 8, units = "in", res = 300)
  p
  dev.off()
  
  
  
