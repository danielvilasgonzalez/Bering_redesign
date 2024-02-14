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
pack_cran<-c('raster','units','ggplot2','data.table','sf','reshape2','data.table')

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
crabs<-c('BB_RKC',"PBL_BKC" ,"PBL_RKC" ,"STM_BKC" ,"SNW_CRB" ,"TNR_CRB")

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

#common species name
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
       'Snow crab_EBSNBS',
       'Blue king crab_EBSNBS',
       'Red king crab_EBSNBS',
       'Tanner crab_EBSNBS')
 
df_spp<-data.frame('spp'=spp,
                   'common'=spp1) 
df_sppcrab<-data.frame('spp'=crabs,
           'common'=c('Bristol Bay\nred king crab','Pribilof Islands\nred king crab','Pribilof Islands\nblue king crab',
                      'St. Matthew Island\nblue king crab','Snow crab','Tanner crab'))  
df_spp1<-rbind(df_spp,df_sppcrab)
df_spp1<-df_spp1[-c(11:14),]
df_spp1<-df_spp1[order(df_spp1$common),]
df_spp1$label<-letters[1:nrow(df_spp1)]

#years
#yrs<-c(1982:2019,2021,2022)
yrs<-c(1982:2022)
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
grid2<-grid

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
samp_df<-rbind(samp_df,c('existing','systematic',520,15,'scnbase'),
               c('existing w/o corner' ,'systematic',494,15,'scnbase_bis'))

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

####################################
# CRABS
####################################

# Set the path to the geodatabase file
gdb_path <- "./shapefiles/CrabStrataShapefiles_GAPgrid.gdb/"

# List the layers/tables in the geodatabase file
gdb_layers <- st_layers(gdb_path)

# Select the specific table you want to read
selected_layer <- gdb_layers$name[1]

# Read the selected table from the geodatabase
#"BBRKC_strata"        "Pribilof_BKC_strata" "Pribilof_RKC_strata" "StMatt_BKC_strata"   "Norton_RKC_Strata"  
#[6] "EBS_CO_CB_strata"   
gdb_table1 <- st_read(dsn = gdb_path, layer = gdb_layers$name[1])
gdb_table2 <- st_read(dsn = gdb_path, layer = gdb_layers$name[2])
gdb_table3 <- st_read(dsn = gdb_path, layer = gdb_layers$name[3])
gdb_table4 <- st_read(dsn = gdb_path, layer = gdb_layers$name[4])
gdb_table5 <- st_read(dsn = gdb_path, layer = gdb_layers$name[5])
gdb_table6 <- st_read(dsn = gdb_path, layer = gdb_layers$name[6])

crabs<-c('BB_RKC','PBL_BKC','PBL_RKC','STM_BKC','SNW_CRB','TNR_CRB')
names(area)<-crabs
crabs_spp<-c('Paralithodes camtschaticus','Paralithodes platypus','Paralithodes camtschaticus','Paralithodes platypus','Chionoecetes opilio','Chionoecetes bairdi')

areacrab<-list()

for (i in 1:length(gdb_layers$name)) {
  
  #i<-3
  
  #load grid of NBS and EBS
  load('./extrapolation grids/northern_bering_sea_grid.rda')
  load('./extrapolation grids/eastern_bering_sea_grid.rda')
  grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
  grid$cell<-1:nrow(grid)
  #grid2<-grid
  coordinates(grid)<-~Lon + Lat
  proj4string(grid) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
  pts<-spTransform(grid,CRSobj = crs(gdb_table1))
  #pts<-as.data.frame(pts)
  
  #st to sp
  gdb_tablea<-as(get(paste0('gdb_table',i)),'Spatial')
  #plot(pts)
  #plot(gdb_tablea)
  
  if (length(gdb_tablea$Shape_Area)!=1) {
    iarea<-sum(gdb_tablea$Shape_Area)/1000000
  } else{
    iarea<-gdb_tablea$Shape_Area/1000000
  }
  
  areacrab[[i]]<-iarea
  
  
  #plot(gdb_tablea[2])
  
  #points over polygon
  xx<-over(pts,gdb_tablea)
  
  # Check which points fall within the sf object using st_within
  points_within_polygon <- as.data.frame(pts)[!is.na(xx), 'cell']
  
  #name<-
  grid2$newcolumn<-FALSE
  grid2$newcolumn[points_within_polygon]<-TRUE
  names(grid2)[ncol(grid2)]<-gdb_layers$name[i]
}

#correct areas
areacrab[[5]]<-areacrab[[6]]
names(areacrab)<-crabs

# Convert the list to a dataframe with one column
dfarea <- data.frame(area = unlist(areacrab))

# Convert the column to factors with column names as levels
dfarea$sp <- rownames(dfarea)

#get cells in each region
BB_RKC_cells<-grid2[which(grid2$BBRKC_strata==TRUE),'cell']
PBL_KC_cells<-grid2[which(grid2$Pribilof_BKC_strata==TRUE),'cell']
STM_BKC_cells<-grid2[which(grid2$StMatt_BKC_strata==TRUE),'cell']
EBS_C_cells<-grid2[which(grid2$EBS_CO_CB_strata==TRUE),'cell']

######################
# HISTORICAL INDEX
######################

#dataframe to store estimated indices
ind2<-data.frame(matrix(NA,nrow = 0,ncol = 7))
names(ind2)<-c('spp','year','approach','sur','scn','index','sim')

#list of files (100 files, one for each simulated data)
files<-list.files('./output/ms_sim_survey_hist/',pattern = 'index_hist',recursive = TRUE,full.names = TRUE)
files<-files[!grepl('*/index_hist_crab.RData$',files)]
  
  
#loop over simulated data - files  
for (sim in 1:100) {
    
  #sim<-1
  
  #print
  cat(paste0('##### ',' sim', sim))
  
  #load file  
  load(files[(sim*2)-1])
  index_hist_crab<-index_hist
  load(files[sim*2])
  
  ind<-index_hist[,'STRS_mean',,,,]
  ind_crab<-index_hist_crab[,'STRS_mean',,,,]
  
  #array to dataframe
  ind<-as.data.frame.table(ind)
  ind_crab<-as.data.frame.table(ind_crab)
  
  ind1<-rbind(ind,ind_crab) 
  names(ind1)<-c('spp','year','approach','sur','scn','index')
  ind1$year<-gsub('y','',ind1$year)
  ind1$sim<-sim
    
  #append
  ind2<-rbind(ind2,ind1)
    
  #remove object
  rm(index_hist);rm(index_hist_crab)
}

#save simulated index
save(ind2,file = './output/estimated_index_hist.RData') #ind2
#load('./output/estimated_index_hist.RData')

#aggregate df to get mean, q95 and q5 for each group (sp, year, sampling scenario and approach)
df<-aggregate(index ~ spp + year + scn + approach,ind2,FUN = function(x) c(mean = mean(x), q95 = quantile(x,probs=0.95) , q5 = quantile(x,probs=0.05)) )
colnames(df$index)<-c('mean','q95','q5')

#subset by species (removing crab species to avoid confusion with crab stocks)
df<-df[which(df$spp %in% df_spp1$spp),]

#sort factors for plotting purposes
df$scn<-factor(df$scn,
               levels = c('scnbase','scnbase_bis','scn3','scn2','scn1'))

#merge df and sp df
df<-merge(df,df_spp1,by='spp')
df$year<-as.integer(df$year)

#load true index and density of histoorical
load(file = paste0("./output/species/dens_index_hist_OM.RData"))  #dens_index_hist_OM, 
names(dens_index_hist_OM)

#loop over crab stocks to store true density and index (spatially clipping)
for (c in crabs) {
  
  #c<-crabs[4]
  
  if (grepl('PBL',c)) {
    cells<-PBL_KC_cells
  } else if (grepl('CRB',c)) {
    cells<-EBS_C_cells
  } else if (grepl('BB',c)) {
    cells<-BB_RKC_cells
  } else if (grepl('STM',c)) {
    cells<-STM_BKC_cells
  }
  
  #c<-crabs[1]
  sp<-crabs_spp[match(c,crabs)]
  dens_spp<-dens_index_hist_OM[[sp]]$dens[cells,,]
  dim(dens_spp)
  dens_spp1<-sweep(dens_spp,1,STATS = grid2[which(grid2$cell %in% cells),'Area_in_survey_km2'],FUN = '*')
  dens_index_hist_OM[[c]]<-list('dens'=dens_spp1,
                                'index'=colSums(dens_spp1))
  
}

#save true dens and index of historical including crab stocks
save(dens_index_hist_OM,file = paste0("./output/species/dens_index_hist_OM.RData"))  #dens_index_hist_OM, 
#load(file = paste0("./output/species/dens_index_hist_OM.RData"))  #dens_index_hist_OM, 


dens_index_hist_OM$TNR_CRB$index
dens_index_hist_OM$SNW_CRB$index


#df to store results
true_ind<-data.frame(matrix(NA,nrow = length(yrs),ncol = length(c(spp,crabs))))
rownames(true_ind)<-yrs
colnames(true_ind)<-c(spp,crabs)

#loop over species and crab stock to extract the true index
for (sp in c(spp,crabs)) {
  
  #sp<-'SNW_CRB'
  
  if (sp %in% crabs) {
    true_ind[,sp]<-dens_index_hist_OM[[sp]]$index
  } else {
    true_ind[,sp]<-drop_units(dens_index_hist_OM[[sp]]$index[,as.character(yrs),1])
  }
  
}

#arrange true index data
true_ind$year<-as.character(yrs)
true_ind<-subset(true_ind,year!='2020')
true_ind1<-reshape2::melt(true_ind,id.vars='year')
names(true_ind1)<-c('year','spp','value')
unique(true_ind1$spp)
true_ind1<-true_ind1[which(true_ind1$spp %in% df_spp1$spp),]
true_ind1<-merge(true_ind1,df_spp1,by='spp')
true_ind1$year<-as.integer(true_ind1$year)
true_ind1$dummy<-'true index'

#save true ind
save(true_ind,file = paste0("./output/true_ind_hist.RData"))  
load(file = paste0("./output/true_ind_hist.RData"))  #true_ind

#to adjust y axis limits
df$value<-df$index[,'q95']/1000
y_scale<-aggregate(value ~ common, df,max)
y_scale$scale<-y_scale$value+y_scale$value*0.2
y_scale$text<-y_scale$value+y_scale$value*0.15
y_scale$apr<-'sys'
y_scale$year<-2010
y_scale$scn<-'scn1'

#fxn to turn axis into scientific
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

#sort approach (station allocation)
df$approach <- factor(df$approach, levels = c("sys", "sb", "rand"))

head(df)

unique(df$spp)
df[which(df$spp=='SNW_CRB'),]
#df[which(df$spp=='SNW_CRB'),'value']<-df[which(df$spp=='SNW_CRB'),'value']*1000
unique(true_ind1$spp)
#true_ind1[which(true_ind1$spp=='SNW_CRB'),'value']<-true_ind1[which(true_ind1$spp=='SNW_CRB'),'value']*100
#true_ind1[which(true_ind1$spp=='STM_BKC'),'value']<-true_ind1[which(true_ind1$spp=='STM_BKC'),'value']*100

#plot
p<-
  ggplot()+
  geom_ribbon(data=df,aes(x=year,ymax=index[,'q95']/1000,ymin=index[,'q5']/1000,group=interaction(scn,approach,common),fill=scn),alpha=0.1)+
  geom_line(data=df,aes(x=year,y=index[,'mean']/1000,color=scn,group=interaction(scn,approach,common),linetype=approach),linewidth=0.7,alpha=0.8)+
  geom_point(data=true_ind1,aes(x=year,y=value/1000,group=common,shape=dummy),fill='black',color='black',alpha=0.7,size=1.5)+
  labs(y='t',x='')+
    scale_fill_manual(values=c('scn1'='#4b7a99','scn2'='#679bc3','scn3'='#8db6c3','scnbase'='#474554','scnbase_bis'='#878787'),
                      labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
    scale_color_manual(values=c('scn1'='#4b7a99','scn2'='#679bc3','scn3'='#8db6c3','scnbase'='#474554','scnbase_bis'='#878787'),
                       labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
  scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
                     labels = c('existing','existing w/o corner' ,'depth','var temp','depth + var temp'),name='stratification')+
  theme_bw()+
    scale_linetype_manual(values = c('sys'='solid',
                                     'sb'='dashed',
                                     'rand'='dotted'),
                          label=c('systematic','balanced random','random'),
                          name='station allocation')+
  scale_shape_manual(values=c('true index'=4))+
    #coord_trans(y = "exp")+
  scale_x_continuous(expand=c(0,0),
                     breaks = c(1985,1990,1995,2000,2005,2010,2015,2020),
                     minor_breaks = setdiff(1982:2022,c(1982,1985,1990,1995,2000,2005,2010,2015,2020,2022)))+
  scale_y_continuous(expand = c(0,0),limits = c(0,NA),labels = scientific_10)+
  #                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
  theme(panel.grid.minor = element_line(linetype=2,color='grey90',),#strip.background = element_rect(fill='white'),
        legend.key.size = unit(12, 'points'),legend.direction = 'vertical',legend.text = element_text(size=9), #legend.position=c(.85,.19)
        legend.title = element_text(size=10),legend.spacing.x = unit(0.05, "cm"),legend.box.spacing =  unit(0.01, "cm"), #,strip.text = element_text(size=12)
        strip.background = element_blank(),legend.background = element_blank(),legend.box = 'vertical',legend.position = 'right',#legend.justification = 'right',legend.position='bottom',#legend.position=c(.84,.05),
        strip.text = element_blank(),axis.title.x = element_blank())+ #axis.text.x = element_text(angle=90,vjust=0.5),
  expand_limits(y = 0)+
  geom_text(data=y_scale,aes(label = common, y = text),x = Inf, vjust = 'inward', hjust = 1.1,size=4, lineheight = 0.8) + #,fontface='italic'
  geom_text(data=df,aes(label = label),x = 1984, y = Inf, vjust = 1.5,size=5) + #,fontface='italic'
  guides(fill=guide_legend(ncol=1,order = 1),color=guide_legend(ncol=1,order = 1),linetype=guide_legend(ncol=1,order = 2),shape='none')+
  #facet_wrap(~com_sci,scales='free',dir='v',nrow = 3)
  #pacific cod 
  geom_blank(data=y_scale,aes(x=year,y=scale,fill=scn,group =interaction(scn,apr)))+
  facet_wrap(~common,scales='free_y',dir='h',nrow = 5)

#save index plot
ragg::agg_png(paste0('./figures/ms_hist_indices_v2.png'), width = 14, height = 9, units = "in", res = 300)
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

#dataframe to store estimated indices
cv2<-data.frame(matrix(NA,nrow = 0,ncol = 7))
names(cv2)<-c('spp','year','approach','sur','scn','cv','sim')

#list of files (100 files, one for each simulated data)
files<-list.files('./output/ms_sim_survey_hist/',pattern = 'index_hist',recursive = TRUE,full.names = TRUE)
  
#loop over simulated data - files  
for (sim in 1:100) {
    
    #sim<-100
    
    #print
    cat(paste0('##### ',' sim', sim))
    
    #load file  
    load(files[(sim*2)-1])
    index_hist_crab<-index_hist
    load(files[sim*2])
    
    cv<-index_hist[,'CV_sim',,,,]
    cv_crab<-index_hist_crab[,'CV_sim',,,,]
    
    #array to dataframe
    cv<-as.data.frame.table(cv)
    cv_crab<-as.data.frame.table(cv_crab)
    
    cv1<-rbind(cv,cv_crab) 
    names(cv1)<-c('spp','year','approach','sur','scn','index')
    cv1$year<-gsub('y','',cv1$year)
    cv1$sim<-sim
    
    #append
    cv2<-rbind(cv2,cv1)
    
    #remove object
    rm(index_hist);rm(index_hist_crab)
}
  
#save cv sim data  
save(cv2,file = './output/estimated_cvsim_hist.RData')
#load(file = './output/estimated_cvsim_hist.RData')

#year to number
cv2$year<-as.numeric(as.character(cv2$year))

# #subset by species (removing crab species to avoid confusion with crab stocks)
# cv2<-cv2[which(cv2$spp %in% df_spp1$spp),]
# 
# #sort factors for plotting purposes
# cv2$scn<-factor(cv2$scn,
#                levels = c('scnbase','scnbase_bis','scn3','scn2','scn1'))
# 
# #merge df and sp df
# cv2<-merge(cv2,df_spp1,by='spp')

#plot over time
ggplot()+
  geom_line(data=cv2,aes(x=year,y=index,color=scn,group=interaction(scn,approach,sim),linetype=approach),linewidth=0.7,alpha=0.8)+
  labs(y='CV',x='')+
  scale_fill_manual(values=c('scn1'='#4b7a99','scn2'='#679bc3','scn3'='#8db6c3','scnbase'='#474554','scnbase_bis'='#878787'),
                    labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
  scale_color_manual(values=c('scn1'='#4b7a99','scn2'='#679bc3','scn3'='#8db6c3','scnbase'='#474554','scnbase_bis'='#878787'),
                     labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
  scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
                     labels = c('existing','existing w/o corner' ,'depth','var temp','depth + var temp'),name='stratification')+
  theme_bw()+
  scale_x_continuous(expand=c(0,0),
                     breaks = c(1985,1990,1995,2000,2005,2010,2015,2020),
                     minor_breaks = setdiff(1982:2022,c(1982,1985,1990,1995,2000,2005,2010,2015,2020,2022)))+
  scale_linetype_manual(values = c('sys'='solid',
                                   'sb'='dashed',
                                   'rand'='dotted'),
                        label=c('systematic','spatially-balanced','random'),
                        name='station allocation')+
  theme(panel.grid.minor = element_line(linetype=2,color='grey90'),strip.background = element_rect(fill='white'),
        legend.position=c(.92,.30),legend.key.size = unit(20, 'points'),legend.text = element_text(size=10), #legend.position=c(.85,.19)
        legend.title = element_text(size=14),strip.text = element_text(size=12))+ #axis.text.x = element_text(angle=90,vjust=0.5),
  expand_limits(y = 0)+
  facet_wrap(~spp,scales='free',dir='v',nrow = 3)

#aggregate df to get mean, q95 and q5 for each group (year, sampling scenario and approach)
df<-aggregate(index ~ spp + year + scn + approach,cv2,FUN = function(x) c(mean = mean(x), q95 = quantile(x,probs=0.95) , q5 = quantile(x,probs=0.05)) )
colnames(df$index)<-c('mean','q95','q5')

#sort factors for plotting purposes
df$scn<-factor(df$scn,
               levels = c('scnbase','scnbase_bis','scn3','scn2','scn1'))
df<-df[which(df$spp %in% df_spp1$spp),]
df<-merge(df,df_spp1,by='spp')
df$approach<-factor(df$approach,levels=c('sys','sb','rand'))

#plot  
#p<-
df$year<-as.numeric(df$year)
#df$com_sci<-paste0(df$common,'\n(',df$spp,')')

#save plot
# ragg::agg_png(paste0('./figures/species/',sp,'/hist_indices_cv.png'), width = 10, height = 4, units = "in", res = 300)
# p
# dev.off()

#merge results to sampling design table
df1<-merge(df,samp_df,by.x='scn',by.y='samp_scn',all.x=TRUE)

#sort and corrections for plotting purposes
df1$scn<-factor(df1$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))
#df1$strat_var<-gsub('_','\n',df1$strat_var)
#df1$strat_var<-factor(df1$strat_var,levels=c('existing','existing w/o corner' ,'Depth','varTemp','Depth\nvarTemp'))


#for geom_blank(0 and adjust scale)
df1$value<-df1$index[,'mean']
y_scale<-aggregate(value ~ common, df1,max)
y_scale$scale<-y_scale$value+y_scale$value*0.2
y_scale$text<-y_scale$value+y_scale$value*0.15
y_scale$apr<-'sys'
y_scale$scn<-'scn1'
y_scale$year<-2022

#plot
p<-
  ggplot()+
  geom_boxplot(data=df1,aes(x=scn,y=value,fill=scn,group =interaction(scn,approach,spp),linetype=approach),alpha=1,outlier.alpha = 0.3,outlier.size = 1.2,outlier.stroke = 0)+ #x=reorder(scn,value)
  labs(y='CV',x='')+
    stat_summary(data=df1,aes(x=scn,y=value,fill=scn,group =interaction(scn,approach,spp),linetype=approach),alpha=1,position = position_dodge(),geom = "crossbar", fun = "median", linetype = "solid", width = .7,linewidth=0.4)+
    scale_fill_manual(values=c('scn1'='#4b7a99','scn2'='#679bc3','scn3'='#8db6c3','scnbase'='#474554','scnbase_bis'='#878787'),
                   labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
    # scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
    #                    labels = c('existing','existing w/o corner' ,'depth','var temp','depth + var temp'),name='stratification')+
    #scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=,'scnbase_bis'=1),
    #                   labels = c('existing','existing w/o corner' ,'depth','var temp','depth + var temp'),name='stratification')+
    theme_bw()+ 
    facet_wrap(~common,scales='free_y',dir='h',nrow = 5)+
    scale_linetype_manual(values = c('sys'='solid',
                                     'sb'='dashed',
                                     'rand'='dotted'),
                          label=c('systematic','balanced random','random'),
                          name='station allocation')+
    scale_shape_manual(values=c('true index'=21),name='')+
    scale_y_continuous(labels=function(x) sprintf('%.2f',x),expand = c(NA,0.1),limits = c(0,NA))+ #expand = c(NA,0.1),limits = c(0,NA)
    #                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
    theme(panel.grid.minor = element_line(linetype=2,color='grey90',),#strip.background = element_rect(fill='white'),
          legend.key.size = unit(12, 'points'),legend.direction = 'vertical',legend.text = element_text(size=9),axis.text.x = element_blank(), #legend.position=c(.85,.19)
          legend.title = element_text(size=10),legend.spacing.x = unit(0.05, "cm"),legend.box.spacing =  unit(0.01, "cm"), #,strip.text = element_text(size=12)
          strip.background = element_blank(),legend.background = element_blank(),legend.box = 'vertical',legend.position = 'right',#legend.justification = 'right',legend.position='bottom',#legend.position=c(.84,.05),
          strip.text = element_blank(),axis.title.x = element_blank())+ #axis.text.x = element_text(angle=90,vjust=0.5),
    expand_limits(y = 0)+
    geom_text(data=y_scale,aes(label = common, y = text),x = Inf, vjust = 'inward', hjust = 1.1,size=4, lineheight = 0.8) + #,fontface='italic'
    geom_text(data=df1,aes(label = paste0(label,'      ')),x = 'scnbase', y = Inf, vjust = 1.5,size=5) + #,fontface='italic'
    #  geom_text(data=df1,aes(label = common),x = Inf, y = -Inf, hjust = 1.1, vjust = -0.8,size=4) + #,fontface='italic'
  guides(fill=guide_legend(ncol=1,order = 1),color=guide_legend(ncol=1,order = 1),linetype=guide_legend(ncol=1,order = 2),shape='none')+
  #facet_wrap(~com_sci,scales='free',dir='v',nrow = 3)
  geom_blank(data=y_scale,aes(x=scn,y=scale,fill=scn,group =interaction(scn,apr)))
  
#save plot
ragg::agg_png(paste0('./figures/ms_hist_indices_cv_box_v2.png'), width = 14, height = 8, units = "in", res = 300)
p
dev.off()

df<-df1


p<-
  ggplot()+
  geom_line(data=df,aes(x=year,y=index[,'mean'],color=scn,group=interaction(scn,approach),linetype=approach),linewidth=0.7,alpha=0.8)+
  labs(y='CV',x='')+
  scale_fill_manual(values=c('scnbase'='#474554','scnbase_bis'='#878787','scn3'='#8db6c3','scn2'='#679bc3','scn1'='#4b7a99'),
                    labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
  scale_color_manual(values=c('scnbase'='#474554','scnbase_bis'='#878787','scn3'='#8db6c3','scn2'='#679bc3','scn1'='#4b7a99'),
                     labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
  # scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
  #                    labels = c('existing','existing w/o corner' ,'depth','var temp','depth + var temp'),name='stratification')+
  theme_bw()+ 
  scale_x_continuous(expand=c(0,0),
                     breaks = c(1985,1990,1995,2000,2005,2010,2015,2020),
                     minor_breaks = setdiff(1982:2022,c(1982,1985,1990,1995,2000,2005,2010,2015,2020,2022)))+
  
  facet_wrap(~common,scales='free_y',dir='h',nrow = 5)+
   scale_linetype_manual(values = c('sys'='solid',
                                     'sb'='dashed',
                                     'rand'='dotted'),
                          label=c('systematic','balanced random','random'),
                          name='station allocation')+
  scale_shape_manual(values=c('true index'=21),name='')+
  scale_y_continuous(labels=function(x) sprintf('%.2f',x),expand = c(NA,0.1),limits = c(0,NA))+ #expand = c(NA,0.1),limits = c(0,NA)
  #                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
  theme(panel.grid.minor = element_line(linetype=2,color='grey90',),#strip.background = element_rect(fill='white'),
        legend.key.size = unit(12, 'points'),legend.direction = 'vertical',legend.text = element_text(size=9), #legend.position=c(.85,.19)
        legend.title = element_text(size=10),legend.spacing.x = unit(0.05, "cm"),legend.box.spacing =  unit(0.01, "cm"), #,strip.text = element_text(size=12)
        strip.background = element_blank(),legend.background = element_blank(),legend.box = 'vertical',legend.position = 'right',#legend.justification = 'right',legend.position='bottom',#legend.position=c(.84,.05),
        strip.text = element_blank(),axis.title.x = element_blank())+ #axis.text.x = element_text(angle=90,vjust=0.5),
  expand_limits(y = 0)+
  geom_text(data=y_scale,aes(label = common, y = text),x = Inf, vjust = 'inward', hjust = 1.1,size=4, lineheight = 0.8) + #,fontface='italic'
  geom_text(data=df,aes(label = label),x = 1984, y = Inf, vjust = 1.5,size=5) + #,fontface='italic'
  #  geom_text(data=df1,aes(label = common),x = Inf, y = -Inf, hjust = 1.1, vjust = -0.8,size=4) + #,fontface='italic'
  guides(fill=guide_legend(ncol=1,order = 1),color=guide_legend(ncol=1,order = 1),linetype=guide_legend(ncol=1,order = 2),shape='none')+
  geom_blank(data=y_scale,aes(x=year,y=scale,fill=scn,group =interaction(scn,apr)))

# scale_y_continuous(expand = c(0,0),
#                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
#theme(panel.grid.minor = element_line(linetype=2,color='grey'))+ #axis.text.x = element_text(angle=90,vjust=0.5),
#expand_limits(y = 0)+
#facet_wrap(~common,scales='free')

ragg::agg_png(paste0('./figures/ms_hist_cv_timeseries_v2.png'), width = 14, height = 9, units = "in", res = 300)
p
dev.off()


df1$scn<-factor(df1$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))


#plot
p<-
  ggplot()+
  geom_boxplot(data=df1,aes(x=scn,y=value,fill=scn,group =interaction(scn,approach),linetype=approach),alpha=1,outlier.alpha = 0.3,outlier.size = 1.2,outlier.stroke = 0)+ #x=reorder(scn,value)
  stat_summary(data=df1,aes(x=scn,y=value,fill=scn,group =interaction(scn,approach),linetype=approach),alpha=1,position = position_dodge(),geom = "crossbar", fun = "median", linetype = "solid", width = .7,linewidth=0.4)+
  labs(y='CV',x='')+
  scale_fill_manual(values=c('scnbase'='#474554','scnbase_bis'='#878787','scn3'='#8db6c3','scn2'='#679bc3','scn1'='#4b7a99'),
                    labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
  scale_color_manual(values=c('scnbase'='#474554','scnbase_bis'='#878787','scn3'='#8db6c3','scn2'='#679bc3','scn1'='#4b7a99'),
                     labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
  # scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
  #                    labels = c('existing','existing w/o corner' ,'depth','var temp','depth + var temp'),name='stratification')+
  theme_bw()+ 
  #facet_wrap(~common,scales='free_y',dir='h',nrow = 5)+
  scale_linetype_manual(values = c('sys'='solid',
                                   'sb'='dashed',
                                   'rand'='dotted'),
                        label=c('systematic','balanced random','random'),
                        name='station allocation')+
  scale_shape_manual(values=c('true index'=21),name='')+
  scale_y_continuous(labels=function(x) sprintf('%.2f',x),expand = c(NA,0.1),limits = c(0,0.9))+ #expand = c(NA,0.1),limits = c(0,NA)
  #                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
  theme(panel.grid.minor = element_line(linetype=2,color='grey90'),#strip.background = element_rect(fill='white'),
        legend.position=c(0.75,0.88),legend.key.size = unit(12, 'points'),legend.text = element_text(size=9), #legend.position=c(.85,.19)
        legend.title = element_text(size=10),legend.spacing.y = unit(0.05, "cm"),legend.box.spacing =  unit(0.01, "cm"), #,strip.text = element_text(size=12)
        strip.background = element_blank(),legend.box.background = element_rect(color='black'),legend.direction = 'vertical',legend.box = 'horizontal',legend.background = element_blank(),
        strip.text = element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank())+ #axis.text.x = element_text(angle=90,vjust=0.5),
  expand_limits(y = 0)+
  #geom_text(data=df1,aes(label = common),x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,size=4) + #,fontface='italic'
  #  geom_text(data=df1,aes(label = common),x = Inf, y = -Inf, hjust = 1.1, vjust = -0.8,size=4) + #,fontface='italic'
  guides(fill=guide_legend(ncol=1,order=1),color=guide_legend(ncol=1,order=1),linetype=guide_legend(ncol=1,order = 2))#+
  #facet_wrap(~com_sci,scales='free',dir='v',nrow = 3)
  #geom_blank(data=y_scale,aes(x=scn,y=scale,fill=scn,group =interaction(scn,apr)))

  
  ragg::agg_png(paste0('./figures/ms_hist_indices_cv_box_allspp_v2.png'), width = 7, height = 6, units = "in", res = 300)
  p
  dev.off()
  
  
######################
# RRMSE of CV
######################

#get the estimated CV time series for each replicate, sampling scenario and approach
load(file = './output/estimated_cvsim_hist.RData') #cv2
#cv2[approach=='sb' & scn=='scnbase' & spp=='Limanda aspera' & year=='2019' & sim=='1']

#year to character
cv2$year<-as.character(cv2$year)
summary(cv2)

#rename
names(cv2)[6]<-'cvsim'

#get estimated index SD for each survey across years, sampling scenario and approach
load('./output/estimated_index_hist.RData') #ind2
summary(ind2)
index_sd<-aggregate(index ~ spp + year + scn + approach + sim,ind2,FUN = function(x) c(sd = sd(x)))
summary(index_sd)
summary(index_sd[which(index_sd$approach=='sys' & index_sd$scn=='scnbase_bis'),])

#get true ind
load('./output/true_ind_hist.RData') #ind2

#true index reshape
true_ind2<-reshape2::melt(true_ind,id.vars='year')
names(true_ind2)[2]<-'spp'

#merge sd and true index
df<-merge(index_sd,true_ind2,by=c('year','spp'))
df$cvtrue<-df$index/df$value
names(df)[6:7]<-c('indsim','indtrue')
df[which(df$approach=='sys' & df$scn=='scnbase'),]

#year to character
df$year<-as.character(df$year)

#use datatable to fasten the process
setDT(df)
setDT(cv2)
#cv2<-data.table(cv2)

dim(cv2)
df1<-merge(cv2,df,by=c('spp','scn','approach','sim','year'),all.x=TRUE,allow.cartesian=TRUE)
#dim(df)
dim(df1)

df1$sqdiffcv<-(df1$cvsim-df1$cvtrue)^2
df1[which(df1$approach=='sys' & df1$scn=='scnbase'),]

df2<-df1[, .(mean_sqdiffcv = mean(sqdiffcv,na.rm=FALSE),mean_cvsim=mean(cvsim,na.rm=FALSE),mean_cvtrue=mean(cvtrue,na.rm=FALSE)), by = .(spp,scn,approach,sim,year)]
df2$sqrtmean_sqdiffcv<-sqrt(df2$mean_sqdiffcv)
df2$rrmse<-df2$sqrtmean_sqdiffcv/df2$mean_cvsim

#remove 3sd of rrmse
df2<-df2[which(df2$rrmse <= mean(df2$rrmse)+3*sd(df2$rrmse)),]
df3<-df2[which(df2$spp %in% df_spp1$spp),]
df3[which(df3$approach=='sys' & df3$scn=='scnbase'),]

#select species and sort 
df3<-merge(df3,df_spp1,by='spp')
df3$scn<-factor(df3$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))
df3$approach<-factor(df3$approach,levels=c('sys','sb','rand'))

save(df3,file = './output/rrmse_cv_hist.RData')

y_scale<-aggregate(rrmse ~ spp+common+label, df3,max)
y_scale$scale<-y_scale$rrmse+y_scale$rrmse*0.2
y_scale$text<-y_scale$rrmse+y_scale$rrmse*0.15
y_scale$apr<-'sys'
y_scale$scn<-'scn1'
y_scale$year<-2022

#load rrmse
load(file = './output/rrmse_cv_hist.RData') #df3,

#plot boxplot
#p<-
  ggplot()+
  geom_boxplot(data=df3,aes(x=scn,y=rrmse,fill=scn,group =interaction(scn,approach,spp),linetype=approach),alpha=1)+
    expand_limits(y = 0)+
    geom_hline(yintercept=0.3,col='red')+
    facet_wrap(~spp,scales='free_y',dir='h',nrow = 5)
  
  
  ggplot()+
    geom_line(data=df3,aes(x=year,y=rrmse,color=scn,group =interaction(scn,approach,spp,sim),linetype=approach),alpha=1)+
    scale_color_manual(values=c('scnbase'='#474554','scnbase_bis'='#878787','scn3'='#8db6c3','scn2'='#679bc3','scn1'='#4b7a99'),
                       labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
    expand_limits(y = 0)+
    geom_hline(yintercept=0.3,col='red')+
    facet_wrap(~spp,scales='free_y',dir='h',nrow = 5)
  
  df4<-subset(df3,spp='BB_RKC')
  
  
  ggplot()+
    geom_line(data=df4,aes(x=year,y=rrmse,color=scn,group =interaction(scn,approach,spp,sim),linetype=approach),alpha=1)+
    scale_color_manual(values=c('scnbase'='#474554','scnbase_bis'='#878787','scn3'='#8db6c3','scn2'='#679bc3','scn1'='#4b7a99'),
                       labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
    expand_limits(y = 0)+
    geom_hline(yintercept=0.3,col='red')+
    facet_wrap(~sim,scales='free_y',dir='h',nrow = 5)
  
  
  
  
  
  
  
  
  df3$common<-factor(df3$common,levels=c(df_spp1$common))
  y_scale$common<-factor(y_scale$common,levels=c(df_spp1$common))
  

  #plot
  p<-
    ggplot()+
    geom_boxplot(data=df3,aes(x=scn,y=rrmse,fill=scn,group =interaction(scn,approach,spp),linetype=approach),alpha=1,outlier.alpha = 0.3,outlier.size = 1.2,outlier.stroke = 0)+ #x=reorder(scn,value)
      stat_summary(data=df3,aes(x=scn,y=rrmse,fill=scn,group =interaction(scn,approach,spp),linetype=approach),alpha=1,position = position_dodge(),geom = "crossbar", fun = "median", linetype = "solid", width = .7,linewidth=0.3)+
    scale_fill_manual(values=c('scnbase'='#474554','scnbase_bis'='#878787','scn3'='#8db6c3','scn2'='#679bc3','scn1'='#4b7a99'),
                      labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
    scale_color_manual(values=c('scnbase'='#474554','scnbase_bis'='#878787','scn3'='#8db6c3','scn2'='#679bc3','scn1'='#4b7a99'),
                       labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
    # scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
    #                    labels = c('existing','existing w/o corner' ,'depth','var temp','depth + var temp'),name='stratification')+
    theme_bw()+ 
    scale_linetype_manual(values = c('sys'='solid',
                                     'sb'='dashed',
                                     'rand'='dotted'),
                          label=c('systematic','balanced random','random'),
                          name='station allocation')+
    scale_shape_manual(values=c('true index'=21),name='')+
      scale_y_continuous(labels=function(x) sprintf('%.2f',x),expand = c(NA,0.1),limits = c(0,NA))+ #expand = c(NA,0.1),limits = c(0,NA)
    #                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
      theme(panel.grid.minor = element_line(linetype=2,color='grey90',),#strip.background = element_rect(fill='white'),
            legend.key.size = unit(12, 'points'),legend.direction = 'vertical',legend.text = element_text(size=9),axis.text.x = element_blank(), #legend.position=c(.85,.19)
            legend.title = element_text(size=10),legend.spacing.x = unit(0.05, "cm"),legend.box.spacing =  unit(0.01, "cm"), #,strip.text = element_text(size=12)
            strip.background = element_blank(),legend.background = element_blank(),legend.box = 'vertical',legend.position = 'right',#legend.justification = 'right',legend.position='bottom',#legend.position=c(.84,.05),
            strip.text = element_blank(),axis.title.x = element_blank())+ #axis.text.x = element_text(angle=90,vjust=0.5),
      expand_limits(y = 0)+
      geom_text(data=y_scale,aes(label = common, y = text),x = Inf, vjust = 'inward', hjust = 1.1,size=4, lineheight = 0.8) + #,fontface='italic'
      geom_text(data=df3,aes(label = paste0(label,'      ')),x = 'scnbase', y = Inf, vjust = 1.5,size=5) + #,fontface='italic'
      geom_blank(data=y_scale,aes(x=scn,y=scale,fill=scn,group =interaction(scn,apr)))+
      facet_wrap(~common,scales='free_y',dir='h',nrow = 5)+
    expand_limits(y = 0)+
      labs(y='RRMSE of CV',x='')+
    guides(fill=guide_legend(ncol=1,order=1),color=guide_legend(ncol=1,order=1),linetype=guide_legend(ncol=1,order = 2))#+

  #geom_blank(data=y_scale,aes(x=scn,y=scale,fill=scn,group =interaction(scn,apr)))
  
  
  ragg::agg_png(paste0('./figures/ms_hist_rrmse_cv_box.png'), width = 14, height = 8, units = "in", res = 300)
  p
  dev.off()
  
  
  p<-
    ggplot()+
    geom_boxplot(data=df3,aes(x=scn,y=rrmse,fill=scn,group =interaction(scn,approach),linetype=approach),alpha=1,outlier.alpha = 0.3,outlier.size = 1.2,outlier.stroke = 0)+ #x=reorder(scn,value)
    stat_summary(data=df3,aes(x=scn,y=rrmse,fill=scn,group =interaction(scn,approach),linetype=approach),alpha=1,position = position_dodge(),geom = "crossbar", fun = "median", linetype = "solid", width = .7,linewidth=0.3)+
    labs(y='RRMSE of CV',x='')+
    scale_fill_manual(values=c('scnbase'='#474554','scnbase_bis'='#878787','scn3'='#8db6c3','scn2'='#679bc3','scn1'='#4b7a99'),
                      labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
    scale_color_manual(values=c('scnbase'='#474554','scnbase_bis'='#878787','scn3'='#8db6c3','scn2'='#679bc3','scn1'='#4b7a99'),
                       labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
    # scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
    #                    labels = c('existing','existing w/o corner' ,'depth','var temp','depth + var temp'),name='stratification')+
    theme_bw()+ 
    #facet_wrap(~common,scales='free_y',dir='h',nrow = 5)+
    scale_linetype_manual(values = c('sys'='solid',
                                     'sb'='dashed',
                                     'rand'='dotted'),
                          label=c('systematic','balanced random','random'),
                          name='station allocation')+
    scale_shape_manual(values=c('true index'=21),name='')+
    scale_y_continuous(labels=function(x) sprintf('%.2f',x),expand = c(NA,0.1),limits = c(0,1.7))+ #expand = c(NA,0.1),limits = c(0,NA)
    #                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
    theme(panel.grid.minor = element_line(linetype=2,color='grey90'),#strip.background = element_rect(fill='white'),
          legend.position=c(0.75,0.88),legend.key.size = unit(12, 'points'),legend.text = element_text(size=9), #legend.position=c(.85,.19)
          legend.title = element_text(size=10),legend.spacing.y = unit(0.05, "cm"),legend.box.spacing =  unit(0.01, "cm"), #,strip.text = element_text(size=12)
          strip.background = element_blank(),legend.box.background = element_rect(color='black'),legend.direction = 'vertical',legend.box = 'horizontal',legend.background = element_blank(),
          strip.text = element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank())+ #axis.text.x = element_text(angle=90,vjust=0.5),
    expand_limits(y = 0)+
    #geom_text(data=df1,aes(label = common),x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,size=4) + #,fontface='italic'
    #  geom_text(data=df1,aes(label = common),x = Inf, y = -Inf, hjust = 1.1, vjust = -0.8,size=4) + #,fontface='italic'
    guides(fill=guide_legend(ncol=1,order=1),color=guide_legend(ncol=1,order=1),linetype=guide_legend(ncol=1,order = 2))#+
  #facet_wrap(~com_sci,scales='free',dir='v',nrow = 3)
  #geom_blank(data=y_scale,aes(x=scn,y=scale,fill=scn,group =interaction(scn,apr)))
  
  
  #save plot210notio
  ragg::agg_png(paste0('./figures/ms_hist_indices_rrmse_box_allsp.png'), width = 9, height = 8, units = "in", res = 300)
  p
  dev.off()

  ######################
  # BIAS of CV
  ######################
  
  #get the estimated CV time series for each replicate, sampling scenario and approach
  load(file = './output/estimated_cvsim_hist.RData') #cv2
  #cv2[approach=='sb' & scn=='scnbase' & spp=='Limanda aspera' & year=='2019' & sim=='1']
  
  #year to character
  cv2$year<-as.character(cv2$year)
  summary(cv2)
  
  #rename
  names(cv2)[6]<-'cvsim'
  
  #get estimated index SD for each survey across years, sampling scenario and approach
  load('./output/estimated_index_hist.RData') #ind2
  summary(ind2)
  index_sd<-aggregate(index ~ spp + year + scn + approach + sim,ind2,FUN = function(x) c(sd = sd(x)))
  summary(index_sd)
  summary(index_sd[which(index_sd$approach=='sys' & index_sd$scn=='scnbase_bis'),])
  
  #get true ind
  load('./output/true_ind_hist.RData') #ind2
  
  #true index reshape
  true_ind2<-reshape2::melt(true_ind,id.vars='year')
  names(true_ind2)[2]<-'spp'
  
  #merge sd and true index
  df<-merge(index_sd,true_ind2,by=c('year','spp'))
  df$cvtrue<-df$index/df$value
  names(df)[6:7]<-c('indsim','indtrue')
  df[which(df$approach=='sys' & df$scn=='scnbase'),]
  
  #year to character
  df$year<-as.character(df$year)
  
  #use datatable to fasten the process
  data.table::setDT(df)
  data.table::setDT(cv2)
  #cv2<-data.table(cv2)
  
  dim(cv2)
  df1<-merge(cv2,df,by=c('spp','scn','approach','sim','year'),all.x=TRUE,allow.cartesian=TRUE)
  #dim(df)
  dim(df1)
  
  df1$sqdiffcv<-(df1$cvsim-df1$cvtrue)^2
  df1[which(df1$approach=='sys' & df1$scn=='scnbase'),]
  
  df2<-df1[, .(mean_sqdiffcv = mean(sqdiffcv,na.rm=FALSE),mean_cvsim=mean(cvsim,na.rm=FALSE),mean_cvtrue=mean(cvtrue,na.rm=FALSE)), by = .(spp,scn,approach,sim,year)]
  df2$sqrtmean_sqdiffcv<-sqrt(df2$mean_sqdiffcv)
  df2$rrmse<-df2$sqrtmean_sqdiffcv/df2$mean_cvsim
  
  #remove 3sd of rrmse
  df2<-df2[which(df2$rrmse <= mean(df2$rrmse)+3*sd(df2$rrmse)),]
  df3<-df2[which(df2$spp %in% df_spp1$spp),]
  
  #select species and sort 
  df3<-merge(df3,df_spp1,by='spp')
  df3$scn<-factor(df3$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))
  df3$approach<-factor(df3$approach,levels=c('sys','sb','rand'))
  
  save(df3,file = './output/rrmse_cv_hist.RData')
  
  y_scale<-aggregate(rrmse ~ spp+common+label, df3,max)
  y_scale$scale<-y_scale$rrmse+y_scale$rrmse*0.2
  y_scale$text<-y_scale$rrmse+y_scale$rrmse*0.15
  y_scale$apr<-'sys'
  y_scale$scn<-'scn1'
  y_scale$year<-2022
  
  
  
  
  
  #get the estimated CV time series for each replicate, sampling scenario and approach
  #load(file = './output/estimated_cvsim_hist.RData') #cv2
  #true cv
  load(file = './output/rrmse_cv_hist.RData') #df3
  
  df3$cvbias<-df3$mean_cvsim-df3$mean_cvtrue
  df3$rbias<-100*(df3$cvbias/df3$mean_cvtrue)

  y_scale<-aggregate(rbias ~ spp+common+label, df3,max)
  y_scale$scale<-y_scale$rbias+y_scale$rbias*0.2
  y_scale$text<-y_scale$rbias+y_scale$rbias*0.15
  y_scale$apr<-'sys'
  y_scale$scn<-'scn1'
  y_scale$year<-2022
  y_scale$common<-factor(y_scale$common,levels=c(df_spp1$common))
  
  #plot
  #p<-
    ggplot()+
    geom_boxplot(data=df3,aes(x=scn,y=rbias,fill=scn,group =interaction(scn,approach,spp),linetype=approach),alpha=1,outlier.alpha = 0.3,outlier.size = 1.2,outlier.stroke = 0)+ #x=reorder(scn,value)
    stat_summary(data=df3,aes(x=scn,y=rbias,fill=scn,group =interaction(scn,approach,spp),linetype=approach),alpha=1,position = position_dodge(),geom = "crossbar", fun = "median", linetype = "solid", width = .7,linewidth=0.3)+
    scale_fill_manual(values=c('scnbase'='#474554','scnbase_bis'='#878787','scn3'='#8db6c3','scn2'='#679bc3','scn1'='#4b7a99'),
                      labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
    scale_color_manual(values=c('scnbase'='#474554','scnbase_bis'='#878787','scn3'='#8db6c3','scn2'='#679bc3','scn1'='#4b7a99'),
                       labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
    # scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
    #                    labels = c('existing','existing w/o corner' ,'depth','var temp','depth + var temp'),name='stratification')+
    theme_bw()+ 
    scale_linetype_manual(values = c('sys'='solid',
                                     'sb'='dashed',
                                     'rand'='dotted'),
                          label=c('systematic','balanced random','random'),
                          name='station allocation')+
    scale_shape_manual(values=c('true index'=21),name='')+
    scale_y_continuous(labels=function(x) sprintf('%.2f',x),expand = c(NA,0.1),limits = c(0,NA))+ #expand = c(NA,0.1),limits = c(0,NA)
    #                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
    theme(panel.grid.minor = element_line(linetype=2,color='grey90',),#strip.background = element_rect(fill='white'),
          legend.key.size = unit(12, 'points'),legend.direction = 'vertical',legend.text = element_text(size=9),axis.text.x = element_blank(), #legend.position=c(.85,.19)
          legend.title = element_text(size=10),legend.spacing.x = unit(0.05, "cm"),legend.box.spacing =  unit(0.01, "cm"), #,strip.text = element_text(size=12)
          strip.background = element_blank(),legend.background = element_blank(),legend.box = 'vertical',legend.position = 'right',#legend.justification = 'right',legend.position='bottom',#legend.position=c(.84,.05),
          strip.text = element_blank(),axis.title.x = element_blank())+ #axis.text.x = element_text(angle=90,vjust=0.5),
    expand_limits(y = 0)+
    geom_text(data=y_scale,aes(label = common, y = text),x = Inf, vjust = 'inward', hjust = 1.1,size=4, lineheight = 0.8) + #,fontface='italic'
    geom_text(data=df3,aes(label = paste0(label,'      ')),x = 'scnbase', y = Inf, vjust = 1.5,size=5) + #,fontface='italic'
    geom_blank(data=y_scale,aes(x=scn,y=scale,fill=scn,group =interaction(scn,apr)))+
    facet_wrap(~common,scales='free_y',dir='h',nrow = 5)+
    expand_limits(y = 0)+
    labs(y='RBIAS of CV',x='')+
    guides(fill=guide_legend(ncol=1,order=1),color=guide_legend(ncol=1,order=1),linetype=guide_legend(ncol=1,order = 2))#+
  
  
  
   
######################
# PROJECTED INDEX
######################


  ############################
  #combine true indeces
  ############################
  #files true indices
  files<-list.files('./output/species/',pattern = ' ms_sim_proj_ind',full.names = TRUE)
  #files true densities from projected
  densfiles<-list.files('./output/species/',pattern = 'ms_sim_proj_dens',full.names = TRUE)
  
  #df to store values
  proj_ind2<-data.frame(matrix(NA,nrow = 0,ncol = 5))
  names(proj_ind2)<-c('sp','year','sim','index','sbt')
  
  #loop over files
  for (f in 1:length(files)) {
    
    #f<-1
    
    cat(paste0('##### ',' - SBT',f))
    
    #load files
    load(files[f])
    load(densfiles[f])
    
    #get true index and conver array to df
    head(proj_ind)
    proj_ind1<-as.data.frame.table(proj_ind)
    names(proj_ind1)<-c('sp','year','sim','index')
    proj_ind1$sbt<-paste0('SBT',f)
    
    #filter densities for crab species
    simdata1<-simdata[,crabs_spp,,]
    
    #loop over crab species
    for (c in crabs) {
      
      if (grepl('PBL',c)) {
        cells<-PBL_KC_cells
      } else if (grepl('CRB',c)) {
        cells<-EBS_C_cells
      } else if (grepl('BB',c)) {
        cells<-BB_RKC_cells
      } else if (grepl('STM',c)) {
        cells<-STM_BKC_cells
      }
      
      #filter by spp
      sp<-crabs_spp[match(c,crabs)]
      dens_spp<-
        simdata[cells,sp,,]
      dim(dens_spp)
      
      #area cell over density cell
      dens_spp<-sweep(dens_spp,1,STATS = grid2[which(grid2$cell %in% cells),'Area_in_survey_km2'],FUN = '*')
      
      #reshape df
      proj_ind1crab<-reshape2::melt(colSums(dens_spp))
      names(proj_ind1crab)<-c('year','sim','index')
      proj_ind1crab$sp<-c
      proj_ind1crab$sbt<-paste0('SBT',f)
      proj_ind1crab<-proj_ind1crab[,c('sp','year','sim','index','sbt')]
      
      #append
      proj_ind1<-rbind(proj_ind1,proj_ind1crab)
    }
    
    #append
    proj_ind2<-rbind(proj_ind2,proj_ind1)
    
    #plot indices
    print(
      ggplot()+
        geom_line(data=proj_ind1,aes(x=year,y=index,group=sim))+
        facet_wrap(~sp,scales = 'free_y',ncol = 3)+
        ggtitle(paste0('SBT',f)))
  }
  
  #save true indices with crab stocks
  save(proj_ind2,file = './output/true_ind_proj.RData') 
  
  ############################
  #combining indices
  ############################  
  ind2<-data.frame(matrix(NA,nrow = 0,ncol = 8))
  names(ind2)<-c('spp','year','approach','sur','scn','index','sbt','sim')
  
  for (sbt in 1:8) {
        
    #sbt<-1
    
     files<-list.files('./output/ms_sim_survey_proj/',pattern = paste0('SBT',sbt),recursive = TRUE,full.names = TRUE)

     for (sim in 1:100) {
       
       #sim<-1
       
       cat(paste0('##### ',' - SBT',sbt,' - sim', sim))
       
       load(files[sim*2])
       #dimnames(index_proj)
       ind<-index_proj[,'STRS_mean',,,,]
       load(files[(sim*2)-1])
       ind_crab<-index_proj[,'STRS_mean',,,,]
       
       #array to dataframe
       ind<-as.data.frame.table(ind)
       ind_crab<-as.data.frame.table(ind_crab)
       ind1<-rbind(ind,ind_crab)
       head(ind1) 
       names(ind1)<-c('spp','year','approach','sur','scn','index')
       ind1$year<-gsub('y','',ind1$year)
       ind1$sbt<-paste0('SBT',sbt)
       ind1$sim<-sim
       
       #append
       ind2<-rbind(ind2,ind1)
       
       #remove object
       rm(index_proj)
     }
  }
  
  save(ind2,file = './output/estimated_index_proj.RData')
  #load('./output/estimated_index_proj.RData')

  head(ind2)
  df_spp1
  #arctic cod
  ind2<-ind2[!(ind2$spp=='Boreogadus saida' & ind2$sbt=='SBT8' & ind2$sim==1),]
  #ind2<-subset(ind2,spp!='Boreogadus saida' & sbt!='SBT8' & sim!=1)
  df<-aggregate(index ~ spp + year + scn + approach + sbt , #+ sim
                ind2,
                FUN = function(x) c(mean = mean(x), q95 = quantile(x,probs=0.95) , q5 = quantile(x,probs=0.05)) )
  colnames(df$index)<-c('mean','q95','q5')
  
  #sort factors for plotting purposes
  df$scn<-factor(df$scn,
                 levels = c('scnbase','scnbase_bis','scn3','scn2','scn1'))
  df<-merge(df,df_spp1,by='spp')
  
  summary(df)
  df$year<-as.numeric(df$year)
  
  df1<-df[which(df$sbt=='SBT3'),]
  
  #to adjust y axis limits
  df1$value<-df1$index[,'q95']/1000
  y_scale<-aggregate(value ~ common, df1,max)
  y_scale$scale<-y_scale$value+y_scale$value*0.2
  y_scale$text<-y_scale$value+y_scale$value*0.15
  y_scale$apr<-'sys'
  y_scale$year<-2023
  y_scale$scn<-'scn1'
  
  #fxn to turn axis into scientific
  scientific_10 <- function(x) {
    parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
  }
  
  #sort approach (station allocation)
  df1$approach <- factor(df1$approach, levels = c("sys", "sb", "rand"))
  
  df1[which(df1$common=='Arctic cod'),]
  df[which(df$common=='Arctic cod'),]
  
  p<-
  ggplot()+
    geom_ribbon(data=df1,aes(x=year,ymax=index[,'q95']/1000,ymin=index[,'q5']/1000,group=interaction(scn,approach,common,sbt),fill=scn),alpha=0.1)+
    geom_line(data=df1,aes(x=year,y=index[,'mean']/1000,color=scn,group=interaction(scn,approach,common,sbt),linetype=approach),linewidth=0.7,alpha=0.8)+
    scale_fill_manual(values=c('scn1'='#4b7a99','scn2'='#679bc3','scn3'='#8db6c3','scnbase'='#474554','scnbase_bis'='#878787'),
                      labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
    scale_color_manual(values=c('scn1'='#4b7a99','scn2'='#679bc3','scn3'='#8db6c3','scnbase'='#474554','scnbase_bis'='#878787'),
                       labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
    scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
                       labels = c('existing','existing w/o corner' ,'depth','var temp','depth + var temp'),name='stratification')+
    theme_bw()+
    labs(y='t',x='')+
    scale_linetype_manual(values = c('sys'='solid',
                                     'sb'='dashed',
                                     'rand'='dotted'),
                          label=c('systematic','balanced random','random'),
                          name='station allocation')+
    scale_shape_manual(values=c('true index'=4))+
    #coord_trans(y = "exp")+
    scale_x_continuous(expand=c(0,0),
                       breaks = c(2023,2024,2025,2026,2027),
                       minor_breaks = c(2024,2026))+
    scale_y_continuous(expand = c(0,0),limits = c(0,NA),labels = scientific_10)+
    #                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
    theme(panel.grid.minor = element_line(linetype=2,color='grey90',),#strip.background = element_rect(fill='white'),
          legend.key.size = unit(12, 'points'),legend.direction = 'vertical',legend.text = element_text(size=9), #legend.position=c(.85,.19)
          legend.title = element_text(size=10),legend.spacing.x = unit(0.05, "cm"),legend.box.spacing =  unit(0.01, "cm"), #,strip.text = element_text(size=12)
          strip.background = element_blank(),legend.background = element_blank(),legend.box = 'vertical',legend.position = 'right',#legend.justification = 'right',legend.position='bottom',#legend.position=c(.84,.05),
          strip.text = element_blank())+ #axis.text.x = element_text(angle=90,vjust=0.5),
    expand_limits(y = 0)+
    geom_text(data=y_scale,aes(label = common, y = text),x = Inf, vjust = 'inward', hjust = 1.1,size=4, lineheight = 0.8) + #,fontface='italic'
    geom_text(data=df,aes(label = label),x = 2023.2, y = Inf, vjust = 1.5,size=5) + #,fontface='italic'
    guides(fill=guide_legend(ncol=1,order = 1),color=guide_legend(ncol=1,order = 1),linetype=guide_legend(ncol=1,order = 2),shape='none')+
    #facet_wrap(~com_sci,scales='free',dir='v',nrow = 3)
    #pacific cod 
    geom_blank(data=y_scale,aes(x=year,y=scale,fill=scn,group =interaction(scn,apr)))+
    facet_wrap(~common,scales='free_y',dir='h',nrow = 5)
  
#save plot
ragg::agg_png(paste0('./figures/SBT3_proj_indices.png'), width = 14, height = 8, units = "in", res = 300)
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
  
  ############################
  #combining CVs
  ############################  
  ind2<-data.frame(matrix(NA,nrow = 0,ncol = 8))
  names(ind2)<-c('spp','year','approach','sur','scn','cv','sbt','sim')

for (sbt in 1:8) {
  
  files<-list.files('./output/ms_sim_survey_proj/',pattern = paste0('SBT',sbt),recursive = TRUE,full.names = TRUE)
  
  for (sim in 1:100) {
    
    cat(paste0('##### ',' - SBT',sbt,' - sim', sim))
    
    load(files[sim*2])
    #dimnames(index_proj)
    ind<-index_proj[,'CV_sim',,,,]
    load(files[(sim*2)-1])
    ind_crab<-index_proj[,'CV_sim',,,,]
    
    #array to dataframe
    ind<-as.data.frame.table(ind)
    ind_crab<-as.data.frame.table(ind_crab)
    ind1<-rbind(ind,ind_crab)
    head(ind1) 
    names(ind1)<-c('spp','year','approach','sur','scn','cv')
    ind1$year<-gsub('y','',ind1$year)
    ind1$sbt<-paste0('SBT',sbt)
    ind1$sim<-sim
    
    #append
    ind2<-rbind(ind2,ind1)
    
    #remove object
    rm(index_proj)
  }
}

cvsim<-ind2
save(cvsim,file = './output/cvsim_proj.RData')
load('./output/cvsim_proj.RData') #cvsim


# Convert your data frame to a data.table
setDT(cvsim) #library(data.table)

# Calculate only the mean grouped by specified columns
result_mean <- cvsim[, .(mean = mean(cv)), by = .(spp, year, scn, approach, sbt)]
result_mean<-merge(result_mean,df_spp1,by='spp')

####################
# by sbt one sp
####################

# sp<-'Gadus macrocephalus'
# result_mean$approach<-factor(result_mean$approach,
#                               levels = c('sys','sb','rand'))
# result_mean$scn<-factor(result_mean$scn,
#                          levels = c('scnbase','scnbase_bis','scn3','scn2','scn1'))
# result_mean1<-subset(result_mean,spp==sp)
# result_mean2<-merge(data.frame(result_mean1),df_spp,by.x='spp',by.y='spp')
# result_mean2$label<-letters[as.numeric(gsub('SBT','',result_mean2$sbt))]
# # result_mean2$approach<-factor(result_mean2$approach,
# #                          levels = c('sys','sb','rand'))
# #   
# #p<-
# ggplot()+
#   geom_boxplot(data=data.frame(result_mean2),aes(x=scn,y=mean,fill=scn,group=interaction(scn,approach,sbt,spp),linetype=approach))+
#   stat_summary(data=data.frame(result_mean2),aes(x=scn,y=mean,fill=scn,group=interaction(scn,approach,sbt,spp),linetype=approach),position = position_dodge(),geom = "crossbar", fun = "median", linetype = "solid", width = .7,linewidth=0.4)+
#   labs(y='CV',x='')+
#   scale_fill_manual(values=c('scn1'='#4b7a99','scn2'='#679bc3','scn3'='#8db6c3','scnbase'='#474554','scnbase_bis'='#878787'),
#                     labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
#   scale_color_manual(values=c('scn1'='#4b7a99','scn2'='#679bc3','scn3'='#8db6c3','scnbase'='#474554','scnbase_bis'='#878787'),
#                      labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
#   theme_bw()+
#   #scale_x_continuous(expand=c(0,0),
#                      #breaks = c(1985,1990,1995,2000,2005,2010,2015,2020),
#   #)+ #minor_breaks = setdiff(1982:2022,c(1982,1985,1990,1995,2000,2005,2010,2015,2020,2022))
#   scale_linetype_manual(values = c('sys'='solid',
#                                    'sb'='dashed',
#                                    'rand'='dotted'),
#                         label=c('fixed','random-balanced','random'),
#                         name='station allocation')+
#   theme(panel.grid.minor = element_line(linetype=2,color='grey90'),#strip.background = element_rect(fill='white'),
#         legend.key.size = unit(12, 'points'),legend.text = element_text(size=9), #legend.position=c(.85,.19)
#         legend.title = element_text(size=10),legend.spacing.x = unit(0.05, "cm"),legend.box.spacing =  unit(0.01, "cm"), #,strip.text = element_text(size=12)
#         strip.background = element_blank(),legend.background = element_blank(),#legend.justification = 'right',legend.poswition='bottom',#
#         strip.text = element_blank(),axis.text.x = element_blank())+ #axis.text.x = element_text(angle=90,vjust=0.5),
#   expand_limits(y = 0)+
#   ggtitle(label = sp)+
#   scale_y_continuous(expand = c(NA,0.1),limits = c(0,max(result_mean2$mean)+max(result_mean2$mean)*0.1))+ #expand = c(NA,0.1),limits = c(0,NA)
#   geom_text(data=data.frame(result_mean2),aes(label = sbt),x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,size=4) + #,fontface='italic'
#   #geom_text(data=data.frame(result_mean2),aes(label = paste0(label,'      ')),x = 'scnbase', y = Inf,vjust = 1.3,size=5) + #,fontface='italic'
#   guides(fill=guide_legend(ncol=1,order = 1),color=guide_legend(ncol=1,order = 1),linetype=guide_legend(ncol=1,order = 2),shape='none')+
#   facet_wrap(~sbt,dir='h',nrow = 2)
# 
# #save plot210notio
# ragg::agg_png(paste0('./figures/sbt_ss_cv_proj_box.png'), width = 12, height =6, units = "in", res = 300)
# p
# dev.off()

#######################################
# by sbt all species combined
#######################################

#to adjust y axis limits
result_mean$mean<-result_mean$mean
y_scale<-aggregate(mean ~ sbt, result_mean,max)
y_scale$value<-y_scale$mean
y_scale$scale<-y_scale$value+y_scale$value*0.2
y_scale$text<-y_scale$value+y_scale$value*0.15
y_scale$apr<-'sys'
y_scale$year<-2023
y_scale$scn<-'scn1'
y_scale$label<-letters[1:nrow(y_scale)]

#sort approach (station allocation)
result_mean$approach <- factor(result_mean$approach, levels = c("sys", "sb", "rand"))
result_mean$scn<-factor(result_mean$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))

#plot
p<-
ggplot()+
  geom_boxplot(data=result_mean,aes(x=scn,y=mean,fill=scn,group =interaction(scn,approach,sbt),linetype=approach),alpha=1,outlier.alpha = 0.3,outlier.size = 1.2,outlier.stroke = 0)+ #x=reorder(scn,value)
  stat_summary(data=result_mean,aes(x=scn,y=mean,fill=scn,group =interaction(scn,approach,sbt),linetype=approach),alpha=1,position = position_dodge(),geom = "crossbar", fun = "median", linetype = "solid", width = .7,linewidth=0.3)+
    labs(y='CV',x='')+
  scale_fill_manual(values=c('scnbase'='#474554','scnbase_bis'='#878787','scn3'='#8db6c3','scn2'='#679bc3','scn1'='#4b7a99'),
                    labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
  scale_color_manual(values=c('scnbase'='#474554','scnbase_bis'='#878787','scn3'='#8db6c3','scn2'='#679bc3','scn1'='#4b7a99'),
                     labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
  # scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
  #                    labels = c('existing','existing w/o corner' ,'depth','var temp','depth + var temp'),name='stratification')+
  theme_bw()+ 
  #facet_wrap(~common,scales='free_y',dir='h',nrow = 5)+
  scale_linetype_manual(values = c('sys'='solid',
                                   'sb'='dashed',
                                   'rand'='dotted'),
                        label=c('systematic','balanced random','random'),
                        name='station allocation')+
  scale_shape_manual(values=c('true index'=21),name='')+
  scale_y_continuous(labels=function(x) sprintf('%.2f',x),expand = c(NA,0.1),limits = c(0,NA))+ #expand = c(NA,0.1),limits = c(0,NA)
  #                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
  theme(panel.grid.minor = element_line(linetype=2,color='grey90',),#strip.background = element_rect(fill='white'),
        legend.key.size = unit(12, 'points'),legend.direction = 'vertical',legend.text = element_text(size=9), #legend.position=c(.85,.19)
        legend.title = element_text(size=10),legend.spacing.x = unit(0.05, "cm"),legend.box.spacing =  unit(0.01, "cm"), #,strip.text = element_text(size=12)
        strip.background = element_blank(),legend.background = element_blank(),legend.box = 'vertical',legend.position = 'right',#legend.justification = 'right',legend.position='bottom',#legend.position=c(.84,.05),
        strip.text = element_blank(), axis.text.x = element_blank(),axis.title.x = element_blank())+ #axis.text.x = element_text(angle=90,vjust=0.5),
  expand_limits(y = 0)+
  #geom_text(data=y_scale,aes(label = common, y = text),x = Inf, vjust = 'inward', hjust = 1.1,size=4, lineheight = 0.8) + #,fontface='italic'
  #geom_text(data=y_scale,aes(label = label),x = 'scn1', y = Inf, vjust = 1.5,size=5) + #,fontface='italic'
  #geom_text(data=y_scale,aes(label = paste0(sbt,'      ')),x = 'scnbase', y = Inf,vjust = 1.3,size=5) + #,fontface='italic'
  
  geom_text(data=y_scale,aes(label = sbt),x = Inf, hjust = 1.1,lineheight = 0.8,y = Inf, vjust = 1.5,size=5) + #,fontface='italic'
  #geom_text(data=y_scale,aes(label = label),x = 'scnbase', y = Inf, vjust = 1.5,size=5) + #,fontface='italic'
  geom_blank(data=y_scale,aes(x=scn,y=scale,fill=scn,group =interaction(scn,apr,sbt)))+
  geom_text(data=y_scale,aes(label = paste0(label,'     ')),x = 'scnbase', y = Inf, vjust = 1.5,size=5) + #,fontface='italic'
  guides(fill=guide_legend(ncol=1,order = 1),color=guide_legend(ncol=1,order = 1),linetype=guide_legend(ncol=1,order = 2),shape='none')+
  facet_wrap(~sbt,dir='h',nrow = 2)


#save plot210notio
ragg::agg_png(paste0('./figures/sbt_ms_cv_proj_box.png'), width = 12, height = 6, units = "in", res = 300)
p
dev.off()

####################
# by spp one sbt
####################

# sbt<-'SBT8'
# result_mean1<-subset(result_mean,sbt==sbt)
# result_mean2<-data.frame(result_mean1)
# #result_mean2$label<-letters[as.numeric(as.character(result_mean2$spp)]
# 
# #for limits purposes geom_blank
# y_scale<-aggregate(mean ~ common, result_mean2,max)
# y_scale$scale<-y_scale$mean+y_scale$mean*0.1
# y_scale$apr<-'sys'
# y_scale$scn<-'scn1'
# y_scale$year<-2022
# 
# #plot 
# ggplot()+
#   geom_boxplot(data=data.frame(result_mean2),aes(x=scn,y=mean,fill=scn,group=interaction(scn,approach,common),linetype=approach))+
#   stat_summary(data=data.frame(result_mean2),aes(x=scn,y=mean,fill=scn,group=interaction(scn,approach)),position = position_dodge(),geom = "crossbar", fun = "median", linetype = "solid", width = .7,linewidth=0.4)+
#   labs(y='CV',x='')+
#   scale_fill_manual(values=c('scn1'='#2C8472','scn2'='#44B8BE','scn3'='#8BE2FA','scnbase'='#5C5354','scnbase_bis'='#A19999'),
#                     labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
#   scale_color_manual(values=c('scn1'='#2C8472','scn2'='#44B8BE','scn3'='#8BE2FA','scnbase'='#5C5354','scnbase_bis'='#A19999'),
#                      labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
#   theme_bw()+
#   #scale_x_continuous(expand=c(0,0),
#   #breaks = c(1985,1990,1995,2000,2005,2010,2015,2020),
#   #)+ #minor_breaks = setdiff(1982:2022,c(1982,1985,1990,1995,2000,2005,2010,2015,2020,2022))
#   scale_linetype_manual(values = c('sys'='solid',
#                                    'sb'='dashed',
#                                    'rand'='dotted'),
#                         label=c('systematic','random-balanced','random'),
#                         name='station allocation')+
#   theme(panel.grid.minor = element_line(linetype=2,color='grey90'),#strip.background = element_rect(fill='white'),
#         legend.key.size = unit(12, 'points'),legend.text = element_text(size=9), #legend.position=c(.85,.19)
#         legend.title = element_text(size=10),legend.spacing.x = unit(0.05, "cm"),legend.box.spacing =  unit(0.01, "cm"), #,strip.text = element_text(size=12)
#         strip.background = element_blank(),legend.background = element_blank(),#legend.justification = 'right',legend.position='bottom',#
#         strip.text = element_blank(),axis.text.x = element_blank(),legend.position=c(.84,.10),legend.box = 'horizontal')+ #axis.text.x = element_text(angle=90,vjust=0.5),
#   expand_limits(y = 0)+
#   ggtitle(label = sbt)+
#   scale_y_continuous(expand = c(NA,1),limits = c(0,NA))+ #expand = c(NA,0.1),limits = c(0,NA)
#   geom_text(data=data.frame(result_mean2),aes(label = common),x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,size=4) + #,fontface='italic'
#   geom_text(data=data.frame(result_mean2),aes(label = paste0(label,'      ')),x = 'scnbase', y = Inf,vjust = 1.3,size=5) + #,fontface='italic'
#   guides(fill=guide_legend(ncol=1,order = 1),color=guide_legend(ncol=1,order = 1),linetype=guide_legend(ncol=1,order = 2),shape='none')+
#   facet_wrap(~common,scales = 'free_y',dir='v',nrow = 5)+
#   geom_blank(data=y_scale,aes(x=scn,y=scale,fill=scn,group =interaction(scn,apr)))


####################
# by spp all sbt
####################

#for limits purposes geom_blank
y_scale<-aggregate(mean ~ common, result_mean,max)
y_scale$scale<-y_scale$mean+y_scale$mean*0.2
y_scale$apr<-'sys'
y_scale$scn<-'scn1'
y_scale$year<-2022
y_scale$label<-letters[1:nrow(y_scale)]
y_scale$text<-y_scale$mean+y_scale$mean*0.15

#plot 
p<-
ggplot()+
  geom_boxplot(data=result_mean,aes(x=scn,y=mean,fill=scn,group =interaction(scn,approach,common),linetype=approach),alpha=1,outlier.alpha = 0.3,outlier.size = 1.2,outlier.stroke = 0)+ #x=reorder(scn,value)
  stat_summary(data=result_mean,aes(x=scn,y=mean,fill=scn,group =interaction(scn,approach,common),linetype=approach),alpha=1,position = position_dodge(),geom = "crossbar", fun = "median", linetype = "solid", width = .7,linewidth=0.3)+
  labs(y='CV',x='')+
  stat_summary(data=result_mean,aes(x=scn,y=mean,fill=scn,group =interaction(scn,approach,spp),linetype=approach),alpha=1,position = position_dodge(),geom = "crossbar", fun = "median", linetype = "solid", width = .7,linewidth=0.4)+
  scale_fill_manual(values=c('scn1'='#4b7a99','scn2'='#679bc3','scn3'='#8db6c3','scnbase'='#474554','scnbase_bis'='#878787'),
                    labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
  # scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
  #                    labels = c('existing','existing w/o corner' ,'depth','var temp','depth + var temp'),name='stratification')+
  #scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=,'scnbase_bis'=1),
  #                   labels = c('existing','existing w/o corner' ,'depth','var temp','depth + var temp'),name='stratification')+
  theme_bw()+ 
  facet_wrap(~common,scales='free_y',dir='h',nrow = 5)+
  scale_linetype_manual(values = c('sys'='solid',
                                   'sb'='dashed',
                                   'rand'='dotted'),
                        label=c('systematic','balanced random','random'),
                        name='station allocation')+
  scale_shape_manual(values=c('true index'=21),name='')+
  scale_y_continuous(labels=function(x) sprintf('%.2f',x),expand = c(NA,0.1),limits = c(0,NA))+ #expand = c(NA,0.1),limits = c(0,NA)
  #                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
  theme(panel.grid.minor = element_line(linetype=2,color='grey90',),#strip.background = element_rect(fill='white'),
        legend.key.size = unit(12, 'points'),legend.direction = 'vertical',legend.text = element_text(size=9),axis.text.x = element_blank(), #legend.position=c(.85,.19)
        legend.title = element_text(size=10),legend.spacing.x = unit(0.05, "cm"),legend.box.spacing =  unit(0.01, "cm"), #,strip.text = element_text(size=12)
        strip.background = element_blank(),legend.background = element_blank(),legend.box = 'vertical',legend.position = 'right',#legend.justification = 'right',legend.position='bottom',#legend.position=c(.84,.05),
        strip.text = element_blank(),axis.title.x = element_blank())+ #axis.text.x = element_text(angle=90,vjust=0.5),
  expand_limits(y = 0)+
  geom_text(data=y_scale,aes(label = common, y = text),x = Inf, vjust = 'inward', hjust = 1.1,size=4, lineheight = 0.8) + #,fontface='italic'
  geom_text(data=y_scale,aes(label = paste0(label,'     '),y=text),x = 'scnbase', vjust = 'inward',size=5) + #,fontface='italic'
  guides(fill=guide_legend(ncol=1,order = 1),color=guide_legend(ncol=1,order = 1),linetype=guide_legend(ncol=1,order = 2),shape='none')+
  #facet_wrap(~com_sci,scales='free',dir='v',nrow = 3)
  #pacific cod 
  geom_blank(data=y_scale,aes(x=scn,y=scale,fill=scn,group =interaction(scn,apr,common)))+
  guides(fill=guide_legend(ncol=1,order = 1),color=guide_legend(ncol=1,order = 1),linetype=guide_legend(ncol=1,order = 2),shape='none')+
  #geom_blank(data=y_scale,aes(x=scn,y=scale,fill=scn,group =interaction(scn,apr)))+
  facet_wrap(~common,scales='free_y',dir='h',nrow = 5)
  
#save plot
ragg::agg_png(paste0('./figures/spp_ms_proj_cv_box.png'), width = 14, height = 8, units = "in", res = 300)
p
dev.off()


####################
# all
####################

#plot
p<-
ggplot()+
  geom_boxplot(data=result_mean,aes(x=scn,y=mean,fill=scn,group =interaction(scn,approach),linetype=approach),alpha=1,outlier.alpha = 0.3,outlier.size = 1.2,outlier.stroke = 0)+ #x=reorder(scn,value))+
  stat_summary(data=result_mean,aes(x=scn,y=mean,fill=scn,group =interaction(scn,approach),linetype=approach),alpha=1,position = position_dodge(),geom = "crossbar", fun = "median", linetype = "solid", width = .7,linewidth=0.3)+
  labs(y='CV',x='')+
  scale_fill_manual(values=c('scnbase'='#474554','scnbase_bis'='#878787','scn3'='#8db6c3','scn2'='#679bc3','scn1'='#4b7a99'),
                    labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
  scale_color_manual(values=c('scnbase'='#474554','scnbase_bis'='#878787','scn3'='#8db6c3','scn2'='#679bc3','scn1'='#4b7a99'),
                     labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
  # scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
  #                    labels = c('existing','existing w/o corner' ,'depth','var temp','depth + var temp'),name='stratification')+
  theme_bw()+ 
  #facet_wrap(~common,scales='free_y',dir='h',nrow = 5)+
  scale_linetype_manual(values = c('sys'='solid',
                                   'sb'='dashed',
                                   'rand'='dotted'),
                        label=c('systematic','balanced random','random'),
                        name='station allocation')+
  scale_shape_manual(values=c('true index'=21),name='')+
  scale_y_continuous(labels=function(x) sprintf('%.2f',x),expand = c(NA,0.1),limits = c(0,1))+ #expand = c(NA,0.1),limits = c(0,NA)
  #                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
  theme(panel.grid.minor = element_line(linetype=2,color='grey90'),#strip.background = element_rect(fill='white'),
        legend.position=c(0.75,0.88),legend.key.size = unit(12, 'points'),legend.text = element_text(size=9), #legend.position=c(.85,.19)
        legend.title = element_text(size=10),legend.spacing.y = unit(0.05, "cm"),legend.box.spacing =  unit(0.01, "cm"), #,strip.text = element_text(size=12)
        strip.background = element_blank(),legend.box.background = element_rect(color='black'),legend.direction = 'vertical',legend.box = 'horizontal',legend.background = element_blank(),
        strip.text = element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank())+ #axis.text.x = element_text(angle=90,vjust=0.5),
  expand_limits(y = 0)+
  #geom_text(data=df1,aes(label = common),x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,size=4) + #,fontface='italic'
  #  geom_text(data=df1,aes(label = common),x = Inf, y = -Inf, hjust = 1.1, vjust = -0.8,size=4) + #,fontface='italic'
  guides(fill=guide_legend(ncol=1,order=1),color=guide_legend(ncol=1,order=1),linetype=guide_legend(ncol=1,order = 2))#+
#facet_wrap(~com_sci,scales='free',dir='v',nrow = 3)
#geom_blank(data=y_scale,aes(x=scn,y=scale,fill=scn,group =interaction(scn,apr)))

#save plot210notio
ragg::agg_png(paste0('./figures/ms_proj_indices_cv_box_allsp.png'), width = 9, height = 8, units = "in", res = 300)
p
dev.off()

######################
# RRMSE of CV
######################
  
#get simulated CV
cvsim<-ind2
cvsim<-data.frame(cvsim)

#get true index per simulated data and projection scn
load(file = './output/true_ind_proj.RData')  #proj_ind2,
setDT(proj_ind2)
#estimated index
load('./output/estimated_index_proj.RData') #ind2
setDT(ind2)
#cv true
rrmse<-data.frame(matrix(NA,nrow = 0,ncol = 12))
names(rrmse)<-c("spp","scn" ,"approach" ,"sim","year" ,       
                  "mean_sqdiffcv","mean_cvsim","mean_cvtrue","sqrtmean_sqdiffcv","rrmse",            
                  "common","label"   )

for (sbtscn in 1:8) {
  
  for (si in 1:100) {
    
    cat(paste0('##### ',' - SBT',sbtscn,' - sim', si))
    
    #si<-1;sbtscn<-1
    
    #cv sim
    cvsim1<-cvsim[which(sbt==paste0('SBT',sbtscn)& sim==si) ,]
    #estimated index
    ind3<-ind2[which(sbt==paste0('SBT',sbtscn) & sim==si),]
    #true index
    proj_ind3<-proj_ind2[which(sbt==paste0('SBT',sbtscn) & sim==si),]
    names(proj_ind3)[c(1,4)]<-c('spp','trueind')
    #estimated sd of est ind
    index_sd<-aggregate(index ~ spp + year + scn + approach + sim,ind3,FUN = function(x) c(sd = sd(x)))
    names(index_sd)[ncol(index_sd)]<-'index_sd'
    #index_sd[which(index_sd$approach=='sys' & index_sd$scn=='scnbase'),]
    
    #true index reshape
    # proj_ind4<-reshape2::melt(data.frame(proj_ind3),id.vars='year')
    # names(true_ind2)[2]<-'spp'
    
    #merge sd and true index
    df<-merge(index_sd,proj_ind3,by=c('year','spp','sim'))
    df$cvtrue<-df$index_sd/df$trueind
    df[which(df$approach=='sys' & df$scn=='scnbase'),]
    
    #year to character
    df$year<-as.character(df$year)
    
    #use datatable to fasten the process
    setDT(df)
    #setDT(cvsim1)
    #cv2<-data.table(cv2)
    
    df1<-merge(cvsim1,df,by=c('spp','scn','approach','sim','year','sbt'),all.x=TRUE,allow.cartesian=TRUE)
    #dim(df)
    dim(df1)
    dim(cvsim1)
    
    df1$sqdiffcv<-(df1$cv-df1$cvtrue)^2
    df1[which(df1$approach=='sys' & df1$scn=='scnbase'),]
    
    df2<-df1[, .(mean_sqdiffcv = mean(sqdiffcv,na.rm=FALSE),mean_cvsim=mean(cv,na.rm=FALSE),mean_cvtrue=mean(cvtrue,na.rm=FALSE)), by = .(spp,scn,approach,sim,year)]
    df2$sqrtmean_sqdiffcv<-sqrt(df2$mean_sqdiffcv)
    df2$rrmse<-df2$sqrtmean_sqdiffcv/df2$mean_cvsim
    
    #remove 3sd of rrmse
    df2<-df2[which(df2$rrmse <= mean(df2$rrmse)+3*sd(df2$rrmse)),]
    df3<-df2[which(df2$spp %in% df_spp1$spp),]
    
    #select species and sort 
    df3<-merge(df3,df_spp1,by='spp')
    df3$scn<-factor(df3$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))
    df3$approach<-factor(df3$approach,levels=c('sys','sb','rand'))
    df3$sbt<-paste0('SBT',sbtscn)
    
    #append
    rrmse<-rbind(rrmse,df3)
    
    }
  }
    
save(rrmse,file = './output/rrmse_cv_proj.RData')
load(file = './output/rrmse_cv_proj.RData') #rrmse,

df<-rrmse  
df2<-rrmse


# # Function to remove values 3 SD away from the mean for each category
# remove_outliers <- function(x) {
#   mean_val <- mean(x)
#   sd_val <- sd(x)
#   x[!(x > mean_val + 3 * sd_val | x < mean_val - 3 * sd_val)]
# }
# 
# # Apply the function to each category
# df_cleaned <- do.call(rbind, by(df, df$sbt, FUN = function(sub_df) {
#   sub_df$value <- remove_outliers(sub_df$value)
#   return(sub_df)
# }))
# 
# # Reset row names if needed
# rownames(df_cleaned) <- NULL
# 
# print(df_cleaned)
# 
# #remove 3sd of rrmse
# df2<-df2[which(df2$rrmse <= mean(df2$rrmse)+3*sd(df2$rrmse)),]
df3<-df2[which(df2$spp %in% df_spp1$spp),]

#######################################
# by sbt all species combined
#######################################

#to adjust y axis limits
df3$mean<-df3$rrmse
y_scale<-aggregate(mean ~ sbt, df3,max)
y_scale$value<-y_scale$mean
y_scale$scale<-y_scale$value+y_scale$value*0.2
y_scale$text<-y_scale$value+y_scale$value*0.15
y_scale$apr<-'sys'
y_scale$year<-2023
y_scale$scn<-'scn1'
y_scale$label<-letters[1:nrow(y_scale)]

#sort approach (station allocation)
df3$approach <- factor(df3$approach, levels = c("sys", "sb", "rand"))
df3$scn<-factor(df3$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))

#plot
p<-
  ggplot()+
  geom_boxplot(data=df3,aes(x=scn,y=mean,fill=scn,group =interaction(scn,approach,sbt),linetype=approach),alpha=1,outlier.alpha = 0.3,outlier.size = 1.2,outlier.stroke = 0)+ #x=reorder(scn,value)
  stat_summary(data=df3,aes(x=scn,y=mean,fill=scn,group =interaction(scn,approach,sbt),linetype=approach),alpha=1,position = position_dodge(),geom = "crossbar", fun = "median", linetype = "solid", width = .7,linewidth=0.3)+
  labs(y='RRMSE of CV',x='')+
  scale_fill_manual(values=c('scnbase'='#474554','scnbase_bis'='#878787','scn3'='#8db6c3','scn2'='#679bc3','scn1'='#4b7a99'),
                    labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
  scale_color_manual(values=c('scnbase'='#474554','scnbase_bis'='#878787','scn3'='#8db6c3','scn2'='#679bc3','scn1'='#4b7a99'),
                     labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
  # scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
  #                    labels = c('existing','existing w/o corner' ,'depth','var temp','depth + var temp'),name='stratification')+
  theme_bw()+ 
  #facet_wrap(~common,scales='free_y',dir='h',nrow = 5)+
  scale_linetype_manual(values = c('sys'='solid',
                                   'sb'='dashed',
                                   'rand'='dotted'),
                        label=c('systematic','balanced random','random'),
                        name='station allocation')+
  scale_shape_manual(values=c('true index'=21),name='')+
  scale_y_continuous(labels=function(x) sprintf('%.2f',x),expand = c(NA,0.1),limits = c(0,NA))+ #expand = c(NA,0.1),limits = c(0,NA)
  #                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
  theme(panel.grid.minor = element_line(linetype=2,color='grey90',),#strip.background = element_rect(fill='white'),
        legend.key.size = unit(12, 'points'),legend.direction = 'vertical',legend.text = element_text(size=9), #legend.position=c(.85,.19)
        legend.title = element_text(size=10),legend.spacing.x = unit(0.05, "cm"),legend.box.spacing =  unit(0.01, "cm"), #,strip.text = element_text(size=12)
        strip.background = element_blank(),legend.background = element_blank(),legend.box = 'vertical',legend.position = 'right',#legend.justification = 'right',legend.position='bottom',#legend.position=c(.84,.05),
        strip.text = element_blank(), axis.text.x = element_blank(),axis.title.x = element_blank())+ #axis.text.x = element_text(angle=90,vjust=0.5),
  expand_limits(y = 0)+
  geom_text(data=y_scale,aes(label = sbt),x = Inf, hjust = 1.1,lineheight = 0.8,y = Inf, vjust = 1.5,size=5) + #,fontface='italic'
  #geom_blank(data=y_scale,aes(x=scn,y=scale,fill=scn,group =interaction(scn,apr,sbt)))+
  geom_text(data=y_scale,aes(label = paste0(label,'     ')),x = 'scnbase', y = Inf, vjust = 1.5,size=5) + #,fontface='italic'
  guides(fill=guide_legend(ncol=1,order = 1),color=guide_legend(ncol=1,order = 1),linetype=guide_legend(ncol=1,order = 2),shape='none')+
  facet_wrap(~sbt,dir='h',scales = 'free_y',nrow = 2)


#save plot210notio
ragg::agg_png(paste0('./figures/sbt_ms_rrmse_proj_box.png'), width = 12, height = 6, units = "in", res = 300)
p
dev.off()

####################
# by spp one sbt
####################

# sbt<-'SBT8'
# result_mean1<-subset(result_mean,sbt==sbt)
# result_mean2<-data.frame(result_mean1)
# #result_mean2$label<-letters[as.numeric(as.character(result_mean2$spp)]
# 
# #for limits purposes geom_blank
# y_scale<-aggregate(mean ~ common, result_mean2,max)
# y_scale$scale<-y_scale$mean+y_scale$mean*0.1
# y_scale$apr<-'sys'
# y_scale$scn<-'scn1'
# y_scale$year<-2022
# 
# #plot 
# ggplot()+
#   geom_boxplot(data=data.frame(result_mean2),aes(x=scn,y=mean,fill=scn,group=interaction(scn,approach,common),linetype=approach))+
#   stat_summary(data=data.frame(result_mean2),aes(x=scn,y=mean,fill=scn,group=interaction(scn,approach)),position = position_dodge(),geom = "crossbar", fun = "median", linetype = "solid", width = .7,linewidth=0.4)+
#   labs(y='CV',x='')+
#   scale_fill_manual(values=c('scn1'='#2C8472','scn2'='#44B8BE','scn3'='#8BE2FA','scnbase'='#5C5354','scnbase_bis'='#A19999'),
#                     labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
#   scale_color_manual(values=c('scn1'='#2C8472','scn2'='#44B8BE','scn3'='#8BE2FA','scnbase'='#5C5354','scnbase_bis'='#A19999'),
#                      labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
#   theme_bw()+
#   #scale_x_continuous(expand=c(0,0),
#   #breaks = c(1985,1990,1995,2000,2005,2010,2015,2020),
#   #)+ #minor_breaks = setdiff(1982:2022,c(1982,1985,1990,1995,2000,2005,2010,2015,2020,2022))
#   scale_linetype_manual(values = c('sys'='solid',
#                                    'sb'='dashed',
#                                    'rand'='dotted'),
#                         label=c('systematic','random-balanced','random'),
#                         name='station allocation')+
#   theme(panel.grid.minor = element_line(linetype=2,color='grey90'),#strip.background = element_rect(fill='white'),
#         legend.key.size = unit(12, 'points'),legend.text = element_text(size=9), #legend.position=c(.85,.19)
#         legend.title = element_text(size=10),legend.spacing.x = unit(0.05, "cm"),legend.box.spacing =  unit(0.01, "cm"), #,strip.text = element_text(size=12)
#         strip.background = element_blank(),legend.background = element_blank(),#legend.justification = 'right',legend.position='bottom',#
#         strip.text = element_blank(),axis.text.x = element_blank(),legend.position=c(.84,.10),legend.box = 'horizontal')+ #axis.text.x = element_text(angle=90,vjust=0.5),
#   expand_limits(y = 0)+
#   ggtitle(label = sbt)+
#   scale_y_continuous(expand = c(NA,1),limits = c(0,NA))+ #expand = c(NA,0.1),limits = c(0,NA)
#   geom_text(data=data.frame(result_mean2),aes(label = common),x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,size=4) + #,fontface='italic'
#   geom_text(data=data.frame(result_mean2),aes(label = paste0(label,'      ')),x = 'scnbase', y = Inf,vjust = 1.3,size=5) + #,fontface='italic'
#   guides(fill=guide_legend(ncol=1,order = 1),color=guide_legend(ncol=1,order = 1),linetype=guide_legend(ncol=1,order = 2),shape='none')+
#   facet_wrap(~common,scales = 'free_y',dir='v',nrow = 5)+
#   geom_blank(data=y_scale,aes(x=scn,y=scale,fill=scn,group =interaction(scn,apr)))


####################
# by spp all sbt
####################

#to adjust y axis limits
df3$mean<-df3$rrmse
y_scale<-aggregate(mean ~ common, df3,max)
y_scale$value<-y_scale$mean
y_scale$scale<-y_scale$value+y_scale$value*0.2
y_scale$text<-y_scale$value+y_scale$value*0.15
y_scale$apr<-'sys'
y_scale$year<-2023
y_scale$scn<-'scn1'
y_scale$label<-letters[1:nrow(y_scale)]

#sort approach (station allocation)
df3$approach <- factor(df3$approach, levels = c("sys", "sb", "rand"))
df3$scn<-factor(df3$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))

#plot
p<-
  ggplot()+
  geom_boxplot(data=df3,aes(x=scn,y=mean,fill=scn,group =interaction(scn,approach,common),linetype=approach),alpha=1,outlier.alpha = 0.3,outlier.size = 1.2,outlier.stroke = 0)+ #x=reorder(scn,value)
  stat_summary(data=df3,aes(x=scn,y=mean,fill=scn,group =interaction(scn,approach,common),linetype=approach),alpha=1,position = position_dodge(),geom = "crossbar", fun = "median", linetype = "solid", width = .7,linewidth=0.3)+
  labs(y='RRMSE of CV',x='')+
  scale_fill_manual(values=c('scnbase'='#474554','scnbase_bis'='#878787','scn3'='#8db6c3','scn2'='#679bc3','scn1'='#4b7a99'),
                    labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
  scale_color_manual(values=c('scnbase'='#474554','scnbase_bis'='#878787','scn3'='#8db6c3','scn2'='#679bc3','scn1'='#4b7a99'),
                     labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
  # scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
  #                    labels = c('existing','existing w/o corner' ,'depth','var temp','depth + var temp'),name='stratification')+
  theme_bw()+ 
  #facet_wrap(~common,scales='free_y',dir='h',nrow = 5)+
  scale_linetype_manual(values = c('sys'='solid',
                                   'sb'='dashed',
                                   'rand'='dotted'),
                        label=c('systematic','balanced random','random'),
                        name='station allocation')+
  scale_shape_manual(values=c('true index'=21),name='')+
  scale_y_continuous(labels=function(x) sprintf('%.2f',x),expand = c(NA,0.1),limits = c(0,NA))+ #expand = c(NA,0.1),limits = c(0,NA)
  #                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
  theme(panel.grid.minor = element_line(linetype=2,color='grey90',),#strip.background = element_rect(fill='white'),
        legend.key.size = unit(12, 'points'),legend.direction = 'vertical',legend.text = element_text(size=9), #legend.position=c(.85,.19)
        legend.title = element_text(size=10),legend.spacing.x = unit(0.05, "cm"),legend.box.spacing =  unit(0.01, "cm"), #,strip.text = element_text(size=12)
        strip.background = element_blank(),legend.background = element_blank(),legend.box = 'vertical',legend.position = 'right',#legend.justification = 'right',legend.position='bottom',#legend.position=c(.84,.05),
        strip.text = element_blank(), axis.text.x = element_blank(),axis.title.x = element_blank())+ #axis.text.x = element_text(angle=90,vjust=0.5),
  expand_limits(y = 0)+
  #geom_text(data=y_scale,aes(label = common),x = Inf, hjust = 1.1,lineheight = 0.8,y = Inf, vjust = 1.5,size=5) + #,fontface='italic'
  geom_blank(data=y_scale,aes(x=scn,y=scale,fill=scn,group =interaction(scn,apr,common)))+
  #geom_text(data=y_scale,aes(label = paste0(label,'     ')),x = 'scnbase', y = Inf, vjust = 1.5,size=5) + #,fontface='italic'
  guides(fill=guide_legend(ncol=1,order = 1),color=guide_legend(ncol=1,order = 1),linetype=guide_legend(ncol=1,order = 2),shape='none')+
  geom_text(data=y_scale,aes(label = common, y = text),x = Inf, vjust = 'inward', hjust = 1.1,size=4, lineheight = 0.8) + #,fontface='italic'
  geom_text(data=y_scale,aes(label = paste0(label,'     '),y=text),x = 'scnbase', vjust = 'inward',size=5) + #,fontface='italic'
  facet_wrap(~common,dir='h',scales = 'free_y',nrow = 4)


#save plot
ragg::agg_png(paste0('./figures/spp_ms_proj_rrmse_box.png'), width = 14, height = 8, units = "in", res = 300)
p
dev.off()


####################
# all
####################

#plot
p<-
  ggplot()+
  geom_boxplot(data=df3,aes(x=scn,y=mean,fill=scn,group =interaction(scn,approach),linetype=approach),alpha=1,outlier.alpha = 0.3,outlier.size = 1.2,outlier.stroke = 0)+ #x=reorder(scn,value))+
  stat_summary(data=df3,aes(x=scn,y=mean,fill=scn,group =interaction(scn,approach),linetype=approach),alpha=1,position = position_dodge(),geom = "crossbar", fun = "median", linetype = "solid", width = .7,linewidth=0.3)+
  scale_fill_manual(values=c('scnbase'='#474554','scnbase_bis'='#878787','scn3'='#8db6c3','scn2'='#679bc3','scn1'='#4b7a99'),
                    labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
  scale_color_manual(values=c('scnbase'='#474554','scnbase_bis'='#878787','scn3'='#8db6c3','scn2'='#679bc3','scn1'='#4b7a99'),
                     labels = c('existing','existing w/o corner' ,'opt depth','opt varSBT','opt depth + varSBT'),name='stratification')+
  labs(y='RRMSE of CV',x='')+
  # scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
  #                    labels = c('existing','existing w/o corner' ,'depth','var temp','depth + var temp'),name='stratification')+
  theme_bw()+ 
  #facet_wrap(~common,scales='free_y',dir='h',nrow = 5)+
  scale_linetype_manual(values = c('sys'='solid',
                                   'sb'='dashed',
                                   'rand'='dotted'),
                        label=c('systematic','balanced random','random'),
                        name='station allocation')+
  scale_shape_manual(values=c('true index'=21),name='')+
  scale_y_continuous(labels=function(x) sprintf('%.2f',x),expand = c(NA,0.1),limits = c(0,1.2))+ #expand = c(NA,0.1),limits = c(0,NA)
  #                    limits =  c(0, max(df$index[,'mean']/1000) + mean(df$index[,'mean'])/10000))+
  theme(panel.grid.minor = element_line(linetype=2,color='grey90'),#strip.background = element_rect(fill='white'),
        legend.position=c(0.75,0.88),legend.key.size = unit(12, 'points'),legend.text = element_text(size=9), #legend.position=c(.85,.19)
        legend.title = element_text(size=10),legend.spacing.y = unit(0.05, "cm"),legend.box.spacing =  unit(0.01, "cm"), #,strip.text = element_text(size=12)
        strip.background = element_blank(),legend.box.background = element_rect(color='black'),legend.direction = 'vertical',legend.box = 'horizontal',legend.background = element_blank(),
        strip.text = element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank())+ #axis.text.x = element_text(angle=90,vjust=0.5),
  expand_limits(y = 0)+
  #geom_text(data=df1,aes(label = common),x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,size=4) + #,fontface='italic'
  #  geom_text(data=df1,aes(label = common),x = Inf, y = -Inf, hjust = 1.1, vjust = -0.8,size=4) + #,fontface='italic'
  guides(fill=guide_legend(ncol=1,order=1),color=guide_legend(ncol=1,order=1),linetype=guide_legend(ncol=1,order = 2))#+

#save plot210notio
ragg::agg_png(paste0('./figures/ms_proj_indices_rrmse_box_allsp.png'), width = 9, height = 8, units = "in", res = 300)
p
dev.off()

    
    
    
    
   