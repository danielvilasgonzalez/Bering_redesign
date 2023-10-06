############################################################
#
#     CHECK MULTIPLE OPTIONS REPLICATES/SURVEY OPTIONS:
#       A-30 simulated data and for each 100 survey sim. n=3000 per year, sampling design and sbt scenario
#       B-30 simulated data and for each 1 survey sim. n=30 per year, sampling design and sbt scenario
#       C-1 simulated data (take two/three individually to compare among them) and for each 100 sim. n=100 per year, samp and sbt
#     SURVEY SIMULATION RUNS EVERY YEAR - sampling location are allowed to vary every year
#     CALCULATE CV AMONG REPLICATES, TRUE CV (for comparison purposes) and CV WITHIN SAMPLING (strata)
#     CHECK COEFFICIENT OF VARIANCE (CV) ACROSS PROJECTED SIMULATIONS
#
############################################################

##############################################
# SETTINGS
##############################################

#libraries from cran to call or install/load
pack_cran<-c('ggplot2','units','splines','raster','sp')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#install VAST if it is not
if (!('VAST' %in% installed.packages())) {
  devtools::install_github("james-thorson/VAST@main", INSTALL_opts="--no-staged-install")};library(VAST)

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#setwd
out_dir<-'/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/'
setwd(out_dir)

############################################
############################################
# CALCULATE CV AMONG SIMULATIONS FOR TRUE INDECES
############################################
############################################

#example species
sp<-'Gadus macrocephalus'

##############################################
# EXTRACT PROJECTED TRUE INDECES (30sim x 5SBT)
##############################################

#list of files of projected objects
#pm is a list of length()=number of simulations (30)
lf<-list.files(paste0('/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/output/species/',sp,'/simulated projected data/'),pattern = 'test_fit_projection_SBT',full.names = TRUE)

#number simulations and projected years
n_proj<-5
n_sim<-30

#store projected index (true index)
ind_proj<-array(NA,
                dim = c(n_sim,n_proj,length(lf)),
                dimnames=list(1:n_sim,2023:2027,1:5))

#loop over projected list projections, each pm represent n_sim in each SBT scenario (here, we tested 5 SBT scenarios)
for (p in lf) {
  
  #p<-lf[1] 
  
  cat(paste('########',match(p,lf),'\n'))
  
  #load object
  load(p) #pm
  
  #loop over simulations
  for (pi in 1:n_sim){
    
    #pi<-1
    
    #extract timse series of index for each sim and sbt
    ts<-drop_units(pm[[pi]]$Index_ctl[,42:46,1])
    ind_proj[pi,,match( p,lf)]<-ts
  
  }
}

#store projected index
save(ind_proj,file=paste0("/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/output/species/",sp,'/simulated projected data/test_fit_projection_index.RData'))
load(file=paste0("/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/output/species/",sp,'/simulated projected data/test_fit_projection_index.RData'))
proj_true<-ind_proj
#plot true index

x<-as.data.frame.table(ind_proj)
names(x)<-c('replicate','year','sbt','index')

#CV timeseries
ggplot()+
  geom_line(data=x,aes(x=year,y=index,color=sbt,group=interaction(replicate,sbt)))+
  labs(x='',y='true index')+
  #scale_y_continuous(limits = c(0,1.5),expand = c(0,0))+
  facet_wrap(~sbt,nrow=1)+
  theme_bw()

#######################
# GET CV AMONG PROJECTED TRUE INDECES
#######################

#store (obs - mean)^2
ind_proj1<-array(NA,
                dim = c(n_sim,n_proj,length(lf)),
                dimnames=list(1:n_sim,2023:2027,1:5))

#store variance as sum(ob-mean^2/n_proj-1)
var_proj<-array(NA,
                 dim = c(n_proj,length(lf)),
                 dimnames=list(2023:2027,1:5))

#CV timeseries
# ggplot()+
#   geom_line(data=cv1,aes(x=year,y=value,color=sbt))+
#   labs(x='future scn',y='CV (true index) among replicates')+
#   scale_y_continuous(limits = c(0,1.5),expand = c(0,0))+
#   theme_bw()


#VARIANCE AMONG PROJECTIONS 
# Calculate Mean: For each time point, calculate the mean (average) value across all the projections. This will give you a baseline reference point.
# Calculate Variance: For each time point, calculate the squared difference between each projection's value and the mean value calculated in step 3. Sum up these squared differences for all projections and divide by the number of projections minus 1 (to get an unbiased estimate) to calculate the variance at that time point.
# Variance = Σ(projection_value - mean_value)^2 / (number_of_projections - 1)

#mean for each year across all simulations for each future scenario
mean_sbt<-apply(ind_proj, c(2,3), mean)

#loop to get the variance for each SBT and year 
for (p in 1:n_proj) {
  for (y in as.character(2023:2027)) {
    
    #p<-1;y<-'2023'
    
    #mean across simulations per year and imulations
    imean<-mean_sbt[y,p]
    
    #loop over simulations
    for (pi in 1:n_sim) {
      
      ind_proj1[pi,y,p]<-(ind_proj[pi,y,p]-imean)^2
      
    }
    
    var_proj[y,p]<-sum(ind_proj1[,y,p])/(n_proj-1)
    
  }
}

ind_proj1


#get CV
cv<-sqrt(var_proj)/mean_sbt
cv1<-reshape2::melt(cv)
names(cv1)<-c('year','sbt','value')
cv1$year<-as.integer(cv1$year)
cv1$sbt<-as.character(cv1$sbt)

#CV timeseries
ggplot()+
  geom_line(data=cv1,aes(x=year,y=value,color=sbt))+
  labs(x='future scn',y='CV (true index) among replicates')+
  scale_y_continuous(limits = c(0,1.5),expand = c(0,0))+
  theme_bw()

#distribution of annual CVs
ggplot()+
  geom_boxplot(data=cv1,aes(x=sbt,y=value,fill=sbt))+
  labs(x='future scn',y='CV (true index) among replicates')+
  scale_y_continuous(limits = c(0,1.5),expand = c(0,0))+
  theme_bw()


############################################
############################################
# CALCULATE CV AMONG SIMULATIONS FOR SAMPLING DESIGN INDECES
############################################
############################################

#setwd
out_dir<-'/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/'
setwd(out_dir)

###################################
# SBT projections
###################################

#save SBT table
load('./tables/SBT_projection.RData')#df_sbt

#name scenario
df_sbt$sbt<-paste0('SBT',df_sbt$sbt_n)
df_sbt$sbt2<-paste0(df_sbt$sbt,'_',df_sbt$Scenario)

df_sbt1<-subset(df_sbt,sbt %in% paste0('SBT',c(1,2,3,5,7))) #status quo, warm moderate variation, warm very high variation, gradually warm, severe gradually warm

n_sim<-30

#yrs
yrs<-1982:2022

#how manyt projected years we want
n_proj<-5

#project_yrs
project_yrs<-((yrs[length(yrs)])+1):(yrs[length(yrs)]+n_proj)

###################################
# LOAD GRID EBS (remember to keep the same order as in fit_model if multiple grids)
###################################

#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
grid$cell<-1:nrow(grid)
#add col and row number
x1<-grid[,c('Lon','Lat','cell')]
names(x1)<-c('x','y','z')
coordinates(x1)=~x + y
crs(x1)<-c(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
x2<-spTransform(x1,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
x3<-data.frame(x2)
x3$x<-as.integer(x3$x)
x3$y<-as.integer(x3$y)
lon<-sort(unique(x3$x),decreasing = FALSE) #1556
lat<-sort(unique(x3$y),decreasing = TRUE) #1507
lons<-data.frame(x=lon,col=1:length(lon))
lats<-data.frame(y=lat,row=1:length(lat))
x4<-merge(x3,lons,by='x',all.x=TRUE)
x5<-merge(x4,lats,by='y',all.x=TRUE)
x5<-x5[,c('y','x','z','col','row')]
colnames(x5)<-c('Lat','Lon','cell','col','row')
x5<-x5[,c('Lat','Lon','cell','col','row')]

###################################
# YEARS AND BASELINE
###################################

#load grid
load('./extrapolation grids/lastversion_grid_EBS.RData')
yrs<-1982:2022
grid_ebs<-grid.ebs_year[which(grid.ebs_year$region != 'EBSslope' & grid.ebs_year$Year %in% yrs),]
dim(grid_ebs)

#load baseline strata and specify corner stations
load('./output/baseline_strata.RData')
baseline_strata$locations[grep('GF|HG|JI|IH|ON|QP|PO',baseline_strata$locations$stationid),]
baseline_strata$locations$corner<-ifelse(grepl('GF|HG|JI|IH|ON|QP|PO',baseline_strata$locations$stationid),'TRUE','FALSE')

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
samp_df<-rbind(samp_df,c('baseline','current',520,15,'scnbase'),
               c('baseline w/o corner','current',494,15,'scnbase_bis'))

######################
# SIMULATE SAMPLING SURVEY AND GET SAMPLES AND ALLOCATIONS FOR EACH SIMULATION (30), YEAR (5), SBT (5)
######################

# 30 simulated data and for each: 100 survey simulations, 5 year,and 10 sampling designs (2*5)
# This is the full combination, but when evaluating, evaluate A, B and C
######################################

#number of simulated data (for test 30)
n_sim<-30

#number of simulated surveys
n_sur<-100

#example sp
sp<-"Gadus macrocephalus"

#loop over sampling designs
for (samp in unique(samp_df$samp_scn)) { #sampling designs
  
  #samp<-unique(samp_df$samp_scn)[1]
  
  s<-match(samp,samp_df$samp_scn)
  
  if (grepl('base',samp)) { #if it contains base is baseline scenario
    
    strata<-as.data.frame(baseline_strata$cell_strata)[,c('cell','Lat','Lon','Stratum')]
    names(strata)[4]<-c('strata')
    D8<-strata
    
    #allocations<-
    if(samp_df[s,'samp_scn']=='scnbase'){
      allocations<-data.frame('strata'=baseline_strata$n_samples$stratum,
                              'n_samples'=baseline_strata$n_samples$scnbase)
    }else{
      allocations<-data.frame('strata'=baseline_strata$n_samples$stratum,
                              'n_samples'=baseline_strata$n_samples$scnbase_bis)
    }
    
  } else {
    
    #load results_optimization
    load(file=paste0('./output/ms_optim_allocations_',samp,'.RData')) #list = c('result_list','ss_sample_allocations','ms_sample_allocations','samples_strata','cv_temp')
    load(file=paste0('./output/species/multisp_optimization_static_data.RData')) #df
    df<-df[,c("Lat",'Lon','cell')]
    
    #cell - strata
    strata<-all$result_list$solution$indices #all$result_list$sol_by_cell
    colnames(strata)<-c('cell','strata')
    strata<-strata[order(strata$cell),]
    
    #add a strata value to each cell
    D8<-merge(df,strata,by='cell',all.x=TRUE)
    D8<-D8[,c("cell","Lat","Lon","strata")]
    names(D8)[4]<-'strata'
    D8$strata<-as.numeric(D8$strata)
    D8$strata<-ifelse(is.na(D8$strata),999,D8$strata)
    
    #strata - n_samples
    allocations<-data.frame('strata'=all$samples_strata$strata,
                            'n_samples'=all$samples_strata$n_samples)
    
  }
  
  #to store all points
  all_points<-array(NA,
                    dim = list(sum(allocations$n_samples),4,length(project_yrs),2,n_sim,n_sur),
                    dimnames = list(c(1:sum(allocations$n_samples)),c('Lon','Lat','cell','strata'),1:length(project_yrs),c('systematic','random'),1:n_sim,1:n_sur))
  
  #to store survey samples
  sim_survey <- array(data = NA, dim = c(length(dimnames(all_points)[[1]]),
                                         length(sp)+2,
                                         length(project_yrs), #add strata number
                                         #length(dimnames(all_points)[[3]]),
                                         2,
                                         nrow(df_sbt1),
                                         n_sim,
                                         n_sur),
                      dimnames = list(dimnames(all_points)[[1]],
                                      c(sp,'strata','cell'),
                                      c(project_yrs),
                                      #dimnames(all_points)[[3]],
                                      c('systematic','random'),
                                      df_sbt1$sbt,
                                      1:n_sim,
                                      1:n_sur))
  
  
  #loop over SBT scenarios
  for (sbt in paste0(df_sbt1$sbt_n)) {
    
    #sbt<-(df_sbt1$sbt_n)[1]
    
    #load projections
    load(paste0('/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/output/species/',sp,'/simulated projected data/test_fit_projection_SBT',sbt,'.Rdata'))
    
    #loop over each simulated data
    for (isim in 1:n_sim) {
 
      #print simulation to check progress
      cat(paste(" #############  projected - sampling design",samp ,'- SBT',sbt,'- simulation ', isim, "of",n_sim,"#############\n"))
      
      #get simulated data
      sim_proj_dens<-drop_units(pm[[isim]]$D_gc[,1,42:46])

      #loop over simulated surveys (100)
      for (isur in 1:n_sur) { #simulations

        #loop over years
        for (y in 1:length(project_yrs)) { #years
          
          #y=1
          #year
          yy<-c(project_yrs)[y]
          
          #to store points
          dfcurrent<-data.frame(matrix(NA,nrow=0,ncol=4))
          colnames(dfcurrent)<-c('Lon','Lat','cell','strata')
          dfrandom<-dfcurrent
          
          #for while purposes
          flag<-TRUE
          
          #random sample for each strata wit distance constrains
          for(n_istrata in 1:nrow(allocations)) {
            
            #n_istrata<-1
            
            istrata<-allocations[n_istrata,'strata']
            
            #subset cells for strata
            df<-subset(as.data.frame(D8),strata==istrata)
            df1<-subset(as.data.frame(D8),strata==istrata & cell %in% baseline_strata$locations$cell)
            n_i<-allocations[n_istrata,'n_samples']
            
            #while error, keep running
            while(flag){
              
              #keep running if error
              tryCatch({
                
                #print loop state
                #cat(paste(" -- strata", n_istrata, 'of',nrow(allocations),  " -- CURRENT  --\n"))
                
                ##############################
                # CURRENT STRATIFIED APPROACH
                ##############################
                
                if(samp_df[s,'samp_scn']=='scnbase'){
                  pointsc<-baseline_strata$locations[which(baseline_strata$locations$stratum==istrata),c('longitude','latitude','cell','stratum')]
                  names(pointsc)<-c('Lon' ,'Lat', 'cell','strata')
                }else if (samp_df[s,'samp_scn']=='scnbase_bis'){
                  pointsc<-baseline_strata$locations[which(baseline_strata$locations$corner==FALSE & baseline_strata$locations$stratum==istrata),c('longitude','latitude','cell','stratum')]
                  names(pointsc)<-c('Lon' ,'Lat',  'cell','strata')
                } else {
                  
                  #if more required samples than available
                  if (nrow(df1)<n_i) {
                    
                    #vectors to store results
                    #dropcell<-c()
                    #selcell<-c()
                    
                    #duplicate df strata
                    dff<-df
                    
                    #df removing available samples from current design
                    dff<-subset(dff, !(cell %in% df1$cell))
                    
                    #cells to complete the required cells
                    ii<-n_i-nrow(df1)
                    
                    #loop over the required samples using random
                    selcell<-sample(dff$cell,ii)
                    
                    
                    #subset points if equal to required number of samples, if not ERROR and rerun iteration
                    if (nrow(subset(df,cell %in% c(selcell,df1$cell))) != n_i) {flag<-TRUE; stop("-- different number of points on current approach",call. = FALSE)}else{
                      pointsc<-subset(df,cell %in% c(selcell,df1$cell))[,c('Lon','Lat','cell','strata')]
                      names(pointsc)<-c('Lon','Lat','cell','strata')}
                    
                    #else there are enough samples to get from the current sampling design    
                  } else {
                    ss<-sample(1:nrow(df1),size = n_i,replace = FALSE)
                    pointsc<-df1[ss,c('Lon','Lat','cell','strata')]}}
                
                
                ##############################
                # STRATIFIED RANDOM APPROACH
                ##############################
                
                #cat(paste(" -- strata", n_istrata, 'of',nrow(allocations),  " -- RANDOM  --\n"))
                
                #duplicate df strata
                dff<-df
                
                #random selection of samples
                selcell<-sample(dff$cell,n_i)
                pointsr<-subset(df,cell %in% selcell)[,c('Lon','Lat','cell','strata')]
                names(pointsr)<-c('Lon','Lat','cell','strata')
                
                #append data if buffer and current df have equal to required number of samples
                if (nrow(pointsc) == n_i & nrow(pointsr) == n_i) {
                  dfcurrent<-rbind(dfcurrent,pointsc)
                  #dfbuffer<-rbind(dfbuffer,pointsb)
                  dfrandom<-rbind(dfrandom,pointsr)
                  flag<-FALSE}
              },
              
              #error message
              error=function(e) {
                message(paste0("ERROR RETRYING",':\n'),e)})
              if (!flag) next
            }
            
            #to start while again
            flag<-TRUE
          }
          
            #append points
            all_points[,,y,'systematic',isim,isur]<-unlist(dfcurrent)
            #all_points[,,y,'buffer']<-unlist(dfbuffer)
            all_points[,,y,'random',isim,isur]<-unlist(dfrandom)
            
            #get points iterations
            pointsc<-data.frame(unlist(all_points[,,y,'systematic',isim,isur]))
            #pointsb<-data.frame(unlist(all_points[,,y,'buffer']))                          
            pointsr<-data.frame(unlist(all_points[,,y,'random',isim,isur]))   
          
    
            #dimnames(sim_proj_dens)
            sim_proj_dens1c<-data.frame(sim_proj_dens[pointsc$cell,y])
            names(sim_proj_dens1c)<-sp
            
            #dimnames(sim_proj_dens)
            sim_proj_dens1r<-data.frame(sim_proj_dens[pointsr$cell,y])
            names(sim_proj_dens1r)<-sp
            
            #append survey densities for each iteration and simulated data
            sim_survey[,,as.character(yy),'systematic',paste0('SBT',sbt),as.character(isim),as.character(isur)]<-
              unlist(data.frame(cbind(sim_proj_dens1c,'strata'=pointsc$strata,'cell'=pointsc$cell)))
            #sim_survey[,,'buffer']<-cbind(sim_dens[pointsb$cell,],pointsb$strata)
            sim_survey[,,as.character(yy),'random',paste0('SBT',sbt),as.character(isim),as.character(isur)]<-
              unlist(data.frame(cbind(sim_proj_dens1r,'strata'=pointsr$strata,'cell'=pointsr$cell)))
  
        }
      }
    }
  }
  #store results
  save(sim_survey, file = paste0('./output/species/sim_proj_survey_test_',samp,'.RData'))  
  save(all_points, file = paste0('./output/species/allocations_proj_survey_test_',samp,'.RData'))  
  
}

# load(file = paste0('./output/species/sim_proj_survey_test_',samp,'.RData'))
# 
# sbt<-'SBT1'
# isim<-1
# cbind(unlist(sim_survey[,,as.character(year),'random',sbt,isim,1]),unlist(sim_survey[,,as.character(year),'systematic',sbt,isim,1]))
# 


######################
# GET INDEX (STRS_mean) FROM SIMULATED SURVEY
######################

#INDEX TRUE (MODEL$BASED)
#load(file = paste0("./output/species/dens_index_project_OM.RData"))  #dens_index_hist_OM, 

#array to store
index_proj<-array(NA,
                  dim=list(1,3,n_proj,2,n_sim,n_sur,nrow(df_sbt1),nrow(samp_df)),
                  dimnames=list(sp,c('STRS_mean','STRS_var','CV_sim'),paste0('y',project_yrs),c('systematic','random'),1:n_sim,1:n_sur,df_sbt1$sbt,samp_df$samp_scn))


#loop over sampling designs
for (samp in samp_df$samp_scn) {
  
  #samp<-samp_df$samp_scn[5]
  
  s<-match(samp,samp_df$samp_scn)
  
  #print scenario to check progress
  cat(paste(" #############  sampling design -", samp, " #############\n"))
  
  #store results
  load(file = paste0('./output/species/sim_proj_survey_test_',samp,'.RData'))  #sim_survey
  
  #when base sampling other files
  if (grepl('base',samp)) {
    
    #conditions on baseline scenarios
    if (samp == 'scnbase') {
      baseline_strata$locations2<-baseline_strata$locations
    } else if (samp == 'scnbase_bis') {
      baseline_strata$locations2<-baseline_strata$locations[which(baseline_strata$locations$corner=='FALSE'),]
    } 
    
    #sort by strata
    baseline_strata$strata_areas<-baseline_strata$strata_areas[order(baseline_strata$strata_areas$X1),]
    baseline_strata$locations2<-baseline_strata$locations2[order(baseline_strata$locations2$stratum),]
    
    #area by strata
    strata_areas <- baseline_strata$strata_areas
    
    #strata data
    survey_detail <- data.frame("Stratum" = baseline_strata$strata_areas$X1, #strata
                                'Nh' = baseline_strata$strata_areas$pct*53464, #number of cells
                                "nh" = data.frame(table(baseline_strata$locations2$stratum))[,c('Freq')]) #number of sample allocations
    
    #weight of strata for each
    survey_detail$Wh <- survey_detail$Nh / sum(survey_detail$Nh)
    survey_detail$wh <- with(survey_detail, nh/Nh)
    
  } else {
    
    #load optimization results
    load(paste0("./output/ms_optim_allocations_",samp_df[s,'samp_scn'],".RData")) #all
    
    #area
    area_cell<-merge(all$result_list$solution$indices, grid, by.x='ID',by.y='cell')
    
    #area by strata
    strata_areas <- aggregate(Area_in_survey_km2 ~ X1, 
                              FUN = sum,
                              data = area_cell)
    
    #strata data
    survey_detail <- data.frame("Stratum" = all$samples_strata$strata, #strata
                                'Nh' = as.integer(table(all$result_list$solution$indices$X1)), #number of cells
                                "nh" = all$samples_strata$n_samples) #number of sample allocations
    
    #weight of strata for each
    survey_detail$Wh <- survey_detail$Nh / sum(survey_detail$Nh)
    survey_detail$wh <- with(survey_detail, nh/Nh)
  }  
  
  #loop over sbt scenarios
  for (sbt in df_sbt1$sbt) {
    
    #print scenario to check progress
    cat(paste(" #############  SBT projection -", sbt, " #############\n"))
    
    #sbt<-'SBT1'
    
    #loop over simulated data
    for (isim in 1:n_sim) {
      
      #isim<-1
      
      #loop over simulated survey
      for (isur in 1:n_sur) {
        
        #loop over approaches
        for (apr in dimnames(sim_survey)[[4]]) {
          
          #apr<-dimnames(sim_survey)[[4]][1]
          
          #loop over sbt scenarios
          for (year in project_yrs) {
            
            #year<-'2023'
            
            #get strata mean over year  
            dens<-unlist(sim_survey[,,as.character(year),apr,sbt,isim,isur])
            y<-data.frame(dens[,c(sp,'strata')],check.names = FALSE)
            yy<-reshape2::melt(y,id.vars=c('strata'))
            
            #mean, sum and var by strata and year (variable)
            yyy<-aggregate(x=yy$value,
                           by=list(strata=yy$strata,sp=yy$variable),
                           FUN = function(x) c('mean' = mean(x,na.rm=T), 'sum' = sum(x),'var' = var(x,na.rm=T) ))
            
            #create df
            zzz<-data.frame('strata'=yyy$strata,'sp'=yyy$sp,'mean'=yyy$x[,c('mean')],'var'=yyy$x[,c('var')]) #/length(yy$value)
            zzzz<-merge(zzz,strata_areas,by.x='strata',by.y='X1',all.x=TRUE)
            zzzz<-merge(zzzz,survey_detail,by.x='strata',by.y='Stratum',all.x=TRUE)
            
            #add index strata for sum to compute index (mean strata density * area of strata) kg!
            zzzz$index_strata<-zzzz$mean*zzzz$Area_in_survey_km2
            
            #add strata var 
            zzzz$strs_var<-zzzz$var*(zzzz$Area_in_survey_km2^2)/zzzz$nh #sum(survey_detail$Nh) 
            
            #sum of strata var and mean density across years (kg/km2)
            zzzz1 <- aggregate(zzzz[,c('strs_var','index_strata')], by= list(zzzz$sp),FUN = sum)
            
            #get CV across years
            zzzz1$cv<- sqrt(zzzz1$strs_var) / zzzz1$index_strata
            
            #mean CV - to check 
            mean(zzzz1$cv,na.rm=TRUE) 
            
            #get outputs
            STRS_mean <- zzzz1$index_strata #absolute strata year biomass
            STRS_var <- zzzz1$strs_var #strata variance 
            CV <- sqrt(STRS_var) / STRS_mean #simulated CV
            
            #append results
            index_proj[,'STRS_mean',paste0('y',year),apr,isim,isur,sbt,samp]<-STRS_mean
            index_proj[,'STRS_var',paste0('y',year),apr,isim,isur,sbt,samp]<-STRS_var
            index_proj[,'CV_sim',paste0('y',year),apr,isim,isur,sbt,samp]<-CV
          }
        }
      }  
    }
  }
}

#save object
save(index_proj, file = paste0("./output/species/projected design estimates_test.RData" ))
load( file = paste0("./output/species/projected design estimates_test.RData" ))

#comparison estimates proj and hist
#CV (STRS_var/index)
mean(index_proj[,'CV_sim',,,,,,],na.rm=TRUE)
#index
mean(index_proj[,'STRS_mean',,,,,,])

######################
# CALCULATE CV ACROSS SIMULATIONS FOR SAMPLING DESIGN INDECES 
######################

#get sampling design based index
index_proj1<-index_proj[,'STRS_mean',,,,,,]

#mean for each year across all simulations for each future scenario
mean_sbt<-
  apply(index_proj1, c(1,2,3,5,6), mean)

#reshape df
ind<-as.data.frame.table(index_proj1)
names(ind)<-c('year','apr','sim','sur','sbt','samp','value')
ind$year<-gsub('y','',ind$year)
ind$year<-as.integer(ind$year)
mean<-as.data.frame.table(mean_sbt)
names(mean)<-c('year','apr','sim','sbt','samp','value')
mean$year<-gsub('y','',mean$year)
mean$year<-as.integer(mean$year)

#indeces 
ggplot()+
  geom_line(data=ind,aes(x=year,y=value,color=samp,linetype=apr,group=interaction(samp,apr,sbt,sim,sur)),alpha=0.7)+
  geom_line(data=mean,aes(x=year,y=value,group=interaction(samp,apr,sbt,sim)),color='black',alpha=0.8)+
  labs(x='',y='design-based index')+
  scale_linetype(name='samples allocation')+
  #scale_y_continuous(limits = c(0,1.5),expand = c(0,0))+
  theme_bw()+
  theme(aspect.ratio = 1)+
  scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
                     name='sampling design')+
  facet_wrap(~sbt,nrow=1)


#get sampling design based index
index_proj1<-index_proj[,'CV_sim',,,,,,]

#reshape df
ind<-as.data.frame.table(index_proj1)
names(ind)<-c('year','apr','sim','sur','sbt','samp','value')
ind$year<-gsub('y','',ind$year)
ind$year<-as.integer(ind$year)
mean<-as.data.frame.table(mean_sbt)
names(mean)<-c('year','apr','sim','sbt','samp','value')
mean$year<-gsub('y','',mean$year)
mean$year<-as.integer(mean$year)

#CV sim
ggplot()+
  #geom_line(data=ind,aes(x=year,y=value,color=samp,linetype=apr,group=interaction(samp,apr,sbt,sim,sur)),alpha=0.7)+
  #geom_line(data=mean,aes(x=year,y=value,group=interaction(samp,apr,sbt,sim)),color='black',alpha=0.8)+
  geom_boxplot(data=ind,aes(x=samp,y=value,color=samp,linetype=apr,group=interaction(samp,apr,sbt)),alpha=0.7)+
  labs(x='',y='CV sim')+
  scale_linetype(name='samples allocation')+
  #scale_y_continuous(limits = c(0,1.5),expand = c(0,0))+
  theme_bw()+
  theme(aspect.ratio = 1)+
  scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
                     name='sampling design')+
  facet_wrap(~sbt,nrow=1)



index_proj1<-index_proj[,'STRS_mean',,,,,,]

#mean for each year across all simulations for each future scenario
mean_sbt<-
  apply(index_proj1, c(1,2,3,5,6), mean)


#to store results
index_proj2<-array(NA,
                   dim = c(n_proj,2,n_sim,n_sur,5,5),
                   dimnames = list(paste0('y',2023:2027),c('systematic','random'),1:n_sim,1:n_sur,df_sbt1$sbt,samp_df$samp_scn))

var_proj<-array(NA,
                  dim=list(n_proj,2,n_sim,nrow(df_sbt1),nrow(samp_df)),
                  dimnames=list(paste0('y',project_yrs),c('systematic','random'),1:n_sim,df_sbt1$sbt,samp_df$samp_scn))

true_cv<-sd_est_ind<-var_proj
#VARIANCE AMONG PROJECTIONS 
# Calculate Mean: For each time point, calculate the mean (average) value across all the projections. This will give you a baseline reference point.
# Calculate Variance: For each time point, calculate the squared difference between each projection's value and the mean value calculated in step 3. Sum up these squared differences for all projections and divide by the number of projections minus 1 (to get an unbiased estimate) to calculate the variance at that time point.
# Variance = Σ(projection_value - mean_value)^2 / (number_of_projections - 1)



#loop to get the variance for each year and 
for (samp in samp_df$samp_scn) {
  
 for (apr in c('random','systematic')) {
   
   for (sbt in df_sbt1$sbt) {
  
    for (y in paste0('y',2023:2027)) {
      
      for (isim in as.character(1:n_sim)) {
        
      #samp<-'scn1';y<-'y2023';sbt<-'SBT1';apr<-'random';isim<-1
      
      #mean across simulations
      imean<-
        mean_sbt[y,apr,isim,sbt,samp]
      
      ind_sim<-as.vector(index_proj[sp,'STRS_mean',y,apr,isim,,sbt,samp])
      
      #var_proj[y,apr,isim,sbt,samp]<-sum(index_proj2[y,apr,isim,,sbt,samp])/(n_sur-1)
      
      #sd_est_ind[y,apr,isim,sbt,samp]<-sd(ind_sim)
      
          
      
        for (isur in 1:n_sur) {
          
          index_proj2[y,apr,isim,isur,sbt,samp]<-(index_proj1[y,apr,isim,isur,sbt,samp]-imean)^2
          
        }
      var_proj[y,apr,isim,sbt,samp]<-sum(index_proj2[y,apr,isim,,sbt,samp])/(n_sur-1)
      
      #SD for true CV
      sd_est_ind[y,apr,isim,sbt,samp]<-sqrt(var_proj[y,apr,isim,sbt,samp])
      
      #true CV
      true_cv[,apr,isim,sbt,samp]<-sd_est_ind[,apr,isim,sbt,samp]/proj_true[isim,,match(sbt,df_sbt1$sbt)]
      
      
     }
    }
   }
  }
}

#for true CV
true_cv

######################
# PLOT CV ESTIMATES ACROSS SIMULATIONS FOR SAMPLING DESIGN INDICES
######################

#cv<-sqrt(var_proj)/mean_sbt
true_cv1<-reshape2::melt(true_cv)
names(true_cv1)<-c('year','apr','sim','sbt','samp','value')
true_cv1$year<-gsub('y','',true_cv1$year)
true_cv1$year<-as.integer(true_cv1$year)
#cv1$sbt<-as.character(cv1$sbt)

true_cv2<-subset(true_cv1,sbt=='SBT1')

#CV timeseries
ggplot()+
  geom_line(data=true_cv1,aes(x=year,y=value,color=samp,linetype=apr,group=interaction(apr,samp,sbt,sim)))+
  labs(x='',y='true CV')+
  #scale_y_continuous(limits = c(0,1.5),expand = c(0,0))+
  theme_bw()+
  scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
                     name='sampling design')+
  facet_wrap(~sbt,nrow=1)



########################
########################

dimnames(index_proj1)
dimnames(var_proj)

cv<-array(NA,dim = dim(index_proj1),dimnames = dimnames(index_proj1))

for (isur in 1:n_sur) {

  #isur<-1
  cv[,,,isur,,]<-sqrt(var_proj)/index_proj1[,,,isur,,]  
  
}

#cv<-sqrt(var_proj)/mean_sbt
cv1<-reshape2::melt(cv)
names(cv1)<-c('year','apr','sim','isur','sbt','samp','value')
cv1$year<-gsub('y','',cv1$year)
cv1$year<-as.integer(cv1$year)
#cv1$sbt<-as.character(cv1$sbt)

######################
# PLOT CV ESTIMATES ACROSS SIMULATIONS FOR SAMPLING DESIGN INDICES
######################

#CV timeseries
ggplot()+
  geom_line(data=cv1,aes(x=year,y=value,color=samp,linetype=apr,group=interaction(apr,samp,sbt,sim)))+
  labs(x='',y='CV estimated index among replicates')+
  #scale_y_continuous(limits = c(0,1.5),expand = c(0,0))+
  theme_bw()+
  scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
                    name='sampling design')+
  facet_wrap(~sbt+sim)

cv2<-subset(cv1,sbt=='SBT1')

cv2$samp<-factor(cv2$samp,levels = c('scnbase','scnbase_bis','scn3','scn2', 'scn1'))

#CV timeseries
ggplot()+
  geom_line(data=cv2,aes(x=year,y=value,color=samp,linetype=apr,group=interaction(apr,samp,sbt,sim,isur)))+
  labs(x='',y='CV estimated index among replicates')+
  #scale_y_continuous(limits = c(0,1.5),expand = c(0,0))+
  theme_bw()+
  scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
                    name='sampling design')+
  theme(axis.text.x = element_text(angle=90))+
  facet_grid(sbt~sim,margins = 'vs')



cv1$samp<-factor(cv1$samp,levels = c('scnbase','scnbase_bis','scn3','scn2', 'scn1'))

#distribution of annual CVs
ggplot()+
  geom_boxplot(data=cv1,aes(x=samp,y=value,fill=samp,linetype=apr))+
  labs(x='',y='CV estimated index among replicates')+
  #scale_y_continuous(limits = c(0,1.5),expand = c(0,0))+
  theme_bw()+
  scale_fill_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
                    name='sampling design')+
  # scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
  #                    labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
  # scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
  #                    labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
   theme_bw()+
  scale_linetype_manual(values = c('systematic'='solid',
                                   'random'='dashed'),
                        label=c('systematic','random'),
                        name='sample allocation')+
  theme(axis.text.x = element_blank())+
  facet_wrap(~sbt+sim)


cv2<-
  aggregate(value ~ year + sbt + apr + samp, cv1,FUN='mean')

#distribution of annual CVs
ggplot()+
  geom_boxplot(data=cv1,aes(x=samp,y=value,fill=samp,linetype=apr))+
  labs(x='',y='mean CV estimated index among replicates')+
  #scale_y_continuous(limits = c(0,1.5),expand = c(0,0))+
  theme_bw()+
  scale_fill_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
                    name='sampling design')+
  # scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
  #                    labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
  # scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
  #                    labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
  theme_bw()+
  scale_linetype_manual(values = c('systematic'='solid',
                                   'random'='dashed'),
                        label=c('systematic','random'),
                        name='sample allocation')+
  theme(axis.text.x = element_blank())+
  facet_wrap(~sbt,nrow=1)

#distribution of annual CVs
ggplot()+
  geom_boxplot(data=cv2,aes(x=samp,y=value,fill=samp,linetype=apr))+
  labs(x='',y='mean CV estimated index among replicates')+
  #scale_y_continuous(limits = c(0,1.5),expand = c(0,0))+
  theme_bw()+
  scale_fill_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
                    name='sampling design')+
  # scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
  #                    labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
  # scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
  #                    labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
  theme_bw()+
  scale_linetype_manual(values = c('systematic'='solid',
                                   'random'='dashed'),
                        label=c('systematic','random'),
                        name='sample allocation')+
  theme(axis.text.x = element_blank())#+
  #facet_wrap(~sbt,nrow=1)

######################
# PLOT CV ESTIMATES ACROSS STRATA FOR SAMPLING DESIGN INDICES
######################

#INDEX SIMULATED (DESIGN$BASED)
load(file = paste0("./output/species/projected design estimates_test.RData" )) #index_hist
ind<-index_proj[,'CV_sim',,,,,,]

#array to dataframe
ind1<-as.data.frame.table(ind)
names(ind1)<-c('year','approach','sim','sur','sbt','samp','index')
ind1$year<-gsub('y','',ind1$year)

#aggregate df to get mean, q95 and q5 for each group (year, sampling scenario and approach)
df<-ind1
#sort factors for plotting purposes
df$samp<-factor(df$samp,
               levels = c('scnbase','scnbase_bis','scn3','scn2','scn1'))
#df<-merge(df,df_spp,by='spp')


#merge results to sampling design table
#df1<-merge(df,samp_df,by.x='scn',by.y='samp_scn',all.x=TRUE)

#sort and corrections for plotting purposes
#df1$scn<-factor(df1$scn,levels=c('scnbase','scnbase_bis',paste0('scn',3:1)))
#df1$strat_var<-gsub('_','\n',df1$strat_var)
#df1$strat_var<-factor(df1$strat_var,levels=c('baseline','baseline w/o corner','Depth','varTemp','Depth\nvarTemp'))
#df1$approach<-factor(df1$approach,levels=c('systematic','random'))
#df1$com_sci<-paste0(df1$common,'\n(',df1$spp,')')

df1<-subset(df,sbt=='SBT1')

#sort y by group
df1 <- df1[order(df1$sim, df1$index), ]

#plot
#distribution of annual CVs
ggplot()+
  geom_boxplot(data=df1,aes(x=samp,y=index,fill=samp,linetype=approach))+
  labs(x='',y='CV within index')+
  #scale_y_continuous(limits = c(0,1.5),expand = c(0,0))+
  theme_bw()+
  scale_fill_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
                    name='sampling design')+
  # scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
  #                    labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
  # scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
  #                    labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
  theme_bw()+
  scale_y_continuous(expand = c(0,0))+
  scale_linetype_manual(values = c('systematic'='solid',
                                   'random'='dashed'),
                        label=c('systematic','random'),
                        name='sample allocation')+
  theme(axis.text.x = element_blank())+
  facet_grid(sbt~sim,scales='free_x')

ggplot()+
  geom_line(data=df1,aes(x=year,y=index,color=samp,linetype=approach,group=interaction(sim,sur,samp,approach)),alpha=0.5)+
  labs(x='',y='CV within index')+
  #scale_y_continuous(limits = c(0,1.5),expand = c(0,0))+
  theme_bw()+
  scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
                    name='sampling design')+
  # scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
  #                    labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
  # scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
  #                    labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
  scale_linetype_manual(values = c('systematic'='solid',
                                   'random'='dashed'),
                        label=c('systematic','random'),
                        name='sample allocation')+
  theme(axis.text.x = element_text(angle=90))+
  scale_y_continuous(expand = c(0,0))+
  facet_grid(sbt~sim,scales='free_x')


ggplot()+
  geom_boxplot(data=df1,aes(x=samp,y=index,fill=samp,linetype=approach))+
  labs(x='',y='CV within index')+
  #scale_y_continuous(limits = c(0,1.5),expand = c(0,0))+
  theme_bw()+
  scale_fill_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),
                    name='sampling design')+
  # scale_color_manual(values=c('scn1'='#4e79a7','scn2'='#59a14f','scn3'='#edc948','scnbase'='#79706E','scnbase_bis'='#e15759'),
  #                    labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
  # scale_alpha_manual(values = c('scn1'=1,'scn2'=1,'scn3'=1,'scnbase'=0.8,'scnbase_bis'=1),
  #                    labels = c('baseline','baseline w/o corner','depth','var temp','depth + var temp'),name='sampling design')+
  theme_bw()+
  scale_y_continuous(expand = c(0,0))+
  scale_linetype_manual(values = c('systematic'='solid',
                                   'random'='dashed'),
                        label=c('systematic','random'),
                        name='sample allocation')+
  theme(axis.text.x = element_blank())

#####################
#####################
# from survey simulation variance 



#######################
# RRMSE of CV
#######################
