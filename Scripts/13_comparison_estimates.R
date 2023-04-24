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

#loop over spp
for (sp in spp) {
  
  sp<-"Gadus macrocephalus"
  
  #load indices for each sbt and scn
  load(file = paste0("./output/species/",sp,'/indices.RData')) #sp_index
  
  #load projections
  load( file = paste0("./output/species/",sp,'/fit_projection.RData')) #pr_list
  
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
  
  
  df_true1$sbt<-factor(df_true1$sbt,
                       levels = paste0('SBT',1:12))
  df_index$sbt<-factor(df_index$sbt,
                       levels = paste0('SBT',1:12))
  
  #df_true1$sbt<-gsub('scn','',df_true1$sbt)
  #df_index$sbt<-gsub('scn','',df_index$sbt)
  
  ggplot()+
    geom_boxplot(data=df_index,aes(x=factor(year),y=index,color=scn),position = position_dodge(width = 0.9))+
    geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn,shape=dummy),color='black',alpha=0.5,position = position_dodge(width = 0.9))+
    facet_wrap(~sbt)+
    theme_bw()+
    theme()+
    labs(color='design-based',shape='model-based',x='year',y='index (t)')#+  
    #scale_y_continuous(limits = c(0,25000))
  
  df_index1<-aggregate(df_index$index,by=list(year=df_index$year,
                                   scn=df_index$scn,
                                   sbt=df_index$sbt),
                       FUN=mean)
  
  for (sbt in unique(df_index1$sbt)) {
    for (scn in unique(df_index1$scn)) {
      for (y in unique(df_index1$year)) {
        
      #sbt<-unique(df_index1$sbt)[1];scn<-unique(df_index1$scn)[1];y<-unique(df_index1$year)[1]
        
      true_i<-df_true1[which(df_true1$sbt==sbt & df_true1$scn == scn & df_true1$year==y),'index']/1000
      samp_i<-df_index1[which(df_index1$sbt==sbt & df_index1$scn == scn & df_index1$year==y),'x']
      
      
      diff_index<-rbind(diff_index,
                        data.frame(sbt=paste0(sbt),
                                   scn=scn,
                                   year=y,
                                   index=true_i-samp_i))  
      
      }
    }
  }
  
  
  diff_index$sbt<-factor(diff_index$sbt,
                       levels = paste0('SBT',1:12))
  
  ggplot()+
    geom_linerange(data=diff_index,aes(x=factor(year),ymin=0,ymax=index,color=scn),position = position_dodge(width = 0.5))+
    geom_point(data=diff_index,aes(x=factor(year),y=index,fill=scn),position = position_dodge(width = 0.5),shape=21,color='black')+
    geom_hline(yintercept=0,linetype='dashed')+
    #geom_point(data=df_true1,aes(x=factor(year),y=index/1000,group=scn,shape=dummy),color='black',alpha=0.5,position = position_dodge(width = 0.9))+
    facet_wrap(~sbt)+
    theme_bw()+
    theme()+
    labs(color='design-based',fill='design-based',shape='model-based',x='year',y='mean difference model-design index (t)')#+  
    #scale_y_continuous(limits = c(0,25000))
  
  
  
  
  
}
      