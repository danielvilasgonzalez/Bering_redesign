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
                     n_samples=c(350), #c(300,500) 520 (EBS+NBS+CRAB);26 (CRAB); 350 (EBS-CRAB); 494 (NBS-CRAB)
                     n_strata=c(10)) #c(5,10,15)

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
# loop
###################################

#loop over spp
for (sp in spp) {
  
  sp<-"Gadus macrocephalus"
  
  #load indices for each sbt and scn
  load(file = paste0("./output/species/",sp,'/indices.RData')) #sp_index
  
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
  
  #plot CV distribution by year/sbt/scn
  ggplot()+
    geom_boxplot(data=df_cv2,aes(x=year,y=cv,color=scn2),position = position_dodge(width = 0.9))+
    #geom_point(data=df_cv,aes(x=year,y=cv,color=scn),position = position_dodge(width = 0.9))+
    facet_wrap(~sbt2)+
    theme_bw()+
    labs(color='sampling scn')
  
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
    labs(color='design-based',fill='design-based',shape='model-based',x='year',y='cv')#+  
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
    labs(color='design-based',fill='design-based',shape='model-based',x='sampling scn',y='cv')#+  
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
  
  #plot true and sampling index for checking
  ggplot()+
    geom_boxplot(data=df_index,aes(x=factor(year),y=index,color=scn),position = position_dodge(width = 0.9))+
    geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn,shape=dummy),color='black',alpha=0.5,position = position_dodge(width = 0.9))+
    facet_wrap(~sbt)+
    theme_bw()+
    theme()+
    labs(color='design-based',shape='model-based',x='year',y='index (t)')#+  
    #scale_y_continuous(limits = c(0,25000))
  
  #aggregate to get mean index by year/scn/sbt
  df_index1<-aggregate(df_index$index,by=list(year=df_index$year,
                                   scn=df_index$scn,
                                   sbt=df_index$sbt),
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
    labs(color='design-based',fill='design-based',shape='model-based',x='year',y='mean difference model-design index (t)')#+  
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
    labs(color='design-based',fill='design-based',shape='model-based',x='year',y='rmse of index')#+  
  #scale_y_continuous(limits = c(0,25000))
  
  #plot PBIAS by year/sbt/scn
  ggplot()+
    #geom_linerange(data=diff_index,aes(x=factor(year),ymin=0,ymax=rmse,color=scn),position = position_dodge(width = 0.5))+
    geom_line(data=diff_index22,aes(x=factor(year),y=pbias,color=scn2,group=scn2))+ #,position = position_dodge(width = 0.5)
    geom_point(data=diff_index22,aes(x=factor(year),y=pbias,fill=scn2),position = position_dodge(width = 0.5),shape=21,color='black')+
    #geom_hline(yintercept=0,linetype='dashed')+
    #geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn,shape=dummy),color='black',alpha=0.5,position = position_dodge(width = 0.9))+
    facet_wrap(~sbt2)+
    theme_bw()+
    theme()+
    labs(color='design-based',fill='design-based',shape='model-based',x='year',y='pbias of index')#+  
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
    labs(color='design-based',fill='design-based',shape='model-based',x='sampling scn',y='rmse of index')#+  
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
    labs(color='design-based',fill='design-based',shape='model-based',x='sampling scn',y='pbias')#+  
  #scale_y_continuous(limits = c(0,25000))
  
}
      