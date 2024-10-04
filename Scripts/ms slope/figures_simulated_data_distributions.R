#####################################
# Settings
#####################################

#libraries from cran to call or install/load
pack_cran<-c('shadowtext',
             'cowplot',
             'ggspatial',
             'raster',
             'rasterVis',
             'rgeos',
             'scales',
             'rnaturalearth',
             'grid',
             'ggplot2',
             'lubridate',
             'ragg',
             'rgdal',
             'dplyr',
             'sp',
             'ggtext',
             'magick',
             'shadowtext',
             'dplyr',
             'googledrive',
             'sf',
             'knitr',
             'data.table',
             'ggshadow',
             'devtools',
             'reshape2')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman");library(pacman)} else {
    library(pacman)
  }

#install akgfmaps to extract shapefile of Alaska
if (!('akgfmaps' %in% installed.packages())) {
  install_github('afsc-gap-products/akgfmaps')};library(akgfmaps)

#install VAST
if (!('VAST' %in% installed.packages())) {
  install_github("james-thorson/VAST@main", INSTALL_opts="--no-staged-install")};library(VAST)

#load/install packages
p_load(pack_cran,character.only = TRUE)

# Set the root directory for knitting to a work directory
dir <- '/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/'

#setwd
setwd(dir)

#selected species
sel_sp<-c("Gadus chalcogrammus", #Alaskan pollock
          "Gadus macrocephalus", #pacific cod
          "Reinhardtius hippoglossoides", #Greenland turbot
          "Chionoecetes opilio") #snow crab

#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
load('./extrapolation grids/bering_sea_slope_grid.rda')
colnames(bering_sea_slope_grid)[4]<-'Stratum'
bering_sea_slope_grid$Stratum<-NA
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),
                          data.frame(eastern_bering_sea_grid,region='EBS'),
                          data.frame(bering_sea_slope_grid,region='SBS')))
grid$cell<-1:nrow(grid)
grid2<-grid
grid<-grid2[,c('Lat','Lon','region','cell')]

###########################
# TEMP AMD REGIME
###########################

#load file grid
load('./grid_EBS_NBS.RData') #grid.ebs_year

#average SBT by region and year
temp_region<-aggregate(Temp ~ Year + region,grid.ebs_year,FUN=mean)
names(temp_region)<-c('year','region','temp')

#get average SBT for all region (EBSshelf+NBS+EBSslope)
alltemp<-aggregate(temp ~ year, temp_region,FUN=mean)
alltemp$region<-'all'
alltemp<-alltemp[,c("year","region","temp")]

#rbind
temp_region<-rbind(temp_region,alltemp)

# Scaling values by group
temp_region1 <- data.frame(temp_region %>%
                             group_by(region) %>%
                             mutate(scaledtemp = scale(temp)) %>%
                             ungroup())

#classify by cold/warm/normal years
temp_region1$typeyear <- ifelse(temp_region1$scaledtemp > 1, 
                                "warm", ifelse(temp_region1$scaledtemp < -1, 
                                               "cold", "normal"))

#cold/warm years
cyrs<-unique(temp_region1[which(temp_region1$typeyear == c('cold') & 
                                  temp_region1$region=='all'),'year'])
wyrs<-unique(temp_region1[which(temp_region1$typeyear == c('warm') & 
                                  temp_region1$year %in% c(1982:2022) & 
                                  temp_region1$region=='all'),'year'])

# Define the years you want to bold
bold_years <- c(2002, 2004, 2008, 2010, 2012, 2016)

# Create custom labels function
custom_labels <- function(x) {
  sapply(x, function(year) {
    if (year %in% bold_years) {
      paste0("**", year, "**")
    } else {
      as.character(year)
    }
  })
}

#plot
p<-
  ggplot() +
    geom_rect(aes(xmin = 2001.5, xmax = 2005.5, ymin = -Inf, ymax = Inf, fill = "red"), 
              alpha = 0.5) +
    geom_rect(aes(xmin = 2013.5, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "red"), 
              alpha = 0.5) +
    geom_rect(aes(xmin = 2005.5, xmax = 2013.5, ymin = -Inf, ymax = Inf, fill = "blue"), 
              alpha = 0.5) +
    geom_hline(yintercept = mean(subset(temp_region1, region == 'all' & year %in% 1982:2022)[,'temp']), alpha = 0.5, linetype = 'dashed') +
    geom_line(data = subset(temp_region1, region == 'all' & year %in% 1982:2022),
              aes(x = year, y = temp, color = region)) +
    geom_point(data = subset(temp_region1, region == 'all' & year %in% 1982:2022),
               aes(x = year, y = temp, color = region)) +
    scale_fill_manual(values = c("red" = "#cc1d1f", "blue" = "#1675ac"),
                      labels = c("normal-cold", "warm"),
                      name = "SBT regime") +
    scale_color_manual(values = c('all' = 'black'),
                       breaks = c('all'),
                       labels = c('Bering\nSea'),name='',guide=FALSE) +
    theme_bw() +
    scale_x_continuous(breaks = c(1985,1990,1995,2000,bold_years,2020),
                       minor_breaks = c(1982:2022), 
                       labels = custom_labels, limits = c(1982, 2022),
                       expand = c(0.01,0)) +
    theme(axis.text.x = element_markdown(angle = 45,hjust=1),
          axis.text.y = element_text(),text = element_text(size = 12)) +
    labs(x = '', y = 'mean SBT')

#save env plot
agg_png(paste0('./figures slope/temperature_regime1.png'), width = 7, height = 4, units = "in", res = 300)
print(p)
dev.off()

# Define a custom color scale function
custom_colors <- colorRampPalette(c("#1675ac", "white", "#cc1d1f"))

ggplot() +
  # Adding vertical rects per year, filled by temperature
  geom_rect(data = subset(temp_region1, region == 'all' & year %in% 1982:2022),
            aes(xmin = year - 0.5, xmax = year + 0.5, ymin = -Inf, ymax = Inf, 
                fill = temp)) +
  
  # Adding horizontal dashed line for mean temperature
  geom_hline(yintercept = mean(subset(temp_region1, region == 'all' & year %in% 1982:2022)[,'temp']), 
             alpha = 0.5, linetype = 'dashed') +
  
  # Plotting the temperature time series for the 'all' region
  geom_line(data = subset(temp_region1, region == 'all' & year %in% 1982:2022),
            aes(x = year, y = temp, color = region)) +
  geom_point(data = subset(temp_region1, region == 'all' & year %in% 1982:2022),
             aes(x = year, y = temp, color = region)) +
  
  # Define the color gradient based on the temperature values
  scale_fill_gradientn(colors = custom_colors(20), 
                       name = "SBT (°C)") +
  
  # Manually setting the color of the line for 'all' region
  scale_color_manual(values = c('all' = 'black'), 
                     breaks = c('all'), 
                     labels = c('Bering\nSea'), 
                     name = '', guide = FALSE) +
  
  # Customizing x-axis with custom breaks and labels
  scale_x_continuous(breaks = c(1985, 1990, 1995, 2000, bold_years, 2020), 
                     minor_breaks = c(1982:2022), 
                     labels = custom_labels, 
                     limits = c(1982, 2022), 
                     expand = c(0.01, 0)) +
  
  # Customizing theme elements
  theme_bw() +
  theme(axis.text.x = element_markdown(angle = 45, hjust = 1), 
        axis.text.y = element_text(), 
        text = element_text(size = 12)) +
  
  # Adding axis labels
  labs(x = '', y = 'mean SBT')

###########################
# DEPTH DISTRIBUTION OF OCCURRENCE
###########################

#load file observed dataframe input 
load('./obs_df.RData') #obs_df

#survey as factors and rename
obs_df$survey_name<-as.factor(obs_df$survey_name)
levels(obs_df$survey_name)<-c('EBSshelf','EBSslope','NBS')

#only presences, EBSshelf+EBSslope and selected species
obs_df<-subset(obs_df,cpue_kgkm2!=0 &
                 survey_name %in% c('EBSshelf','EBSslope') &
                 species %in% sel_sp)

#filter by cold and warm years
#obs_df1<-obs_df[which(obs_df$year %in% c(cyrs,wyrs)),]
obs_df1<-obs_df[which(obs_df$year %in% bold_years),]

#year type
obs_df1$year_type<-ifelse(obs_df1$year %in% c(2002,2004,2016), "warm",'cold')

#filtered to <600m
filt_obs_df1<-subset(obs_df1,depth_m<=600)

#plot Kernel density estimation for depth distributions
p<-
  ggplot(filt_obs_df1, aes(x = depth_m, fill = year_type)) +
    geom_density(alpha = 0.5) +
    labs(x = 'depth (m)', y = 'density of occurrence')+ 
         #title = 'Empirical depth distribution of occurrence in Warm vs. Cold Years') +
    scale_fill_manual(values = c('cold' = '#1675ac',
                                 'warm' = '#cc1d1f'),
                      labels = c("normal-cold", "warm"),
                      name = 'SBT regime') +
    theme_bw()+
    scale_x_continuous(limits = c(0,600))+
    theme(strip.text = element_text(size=12),strip.background = element_blank(),
          text= element_text(size=12))+
    facet_wrap(~species)

#save env plot
agg_png(paste0('./figures slope/depth_distribution_occurrence.png'), width = 7, height = 4, units = "in", res = 300)
print(p)
dev.off()

#load files density slope and shelf data 
load('./output slope/species/ms_sim_dens_all.RData')  #sim_dens1
dimnames(sim_dens1)

dens<-sim_dens1[,sel_sp,as.character(bold_years),]
dimnames(dens)

dens1<-as.data.frame.table(dens)
names(dens1)<-c('cell','species','Year','sim','dens')

#year type
dens1$year_type<-ifelse(dens1$Year %in% wyrs, "warm",'cold')

#load env
load('./grid_EBS_NBS_envs.RData') #xx
envs<-subset(xx,Year=='2002')
grid1<-merge(grid,envs,by=c('Lat','Lon'))

dens2<-merge(dens1,grid1[,c('cell','depth_m')],by=c('cell'))
dim(dens2)
dim(dens1)

#filtered to <600m
filt_dens2<-subset(dens2,depth_m<=600 & dens!=0)

#plot Kernel density estimation for depth distributions
p<-
  ggplot(filt_dens2, aes(x = depth_m, fill = year_type)) +
  geom_density(alpha = 0.5) +
  labs(x = 'depth (m)', y = 'density of occurrence')+ 
  #title = 'Empirical depth distribution of occurrence in Warm vs. Cold Years') +
  scale_fill_manual(values = c('cold' = '#1675ac',
                               'warm' = '#cc1d1f'),
                    labels = c("normal-cold", "warm"),
                    name = 'SBT regime') +
  theme_bw()+
  scale_x_continuous(limits = c(0,600))+
  theme(strip.text = element_text(size=12),strip.background = element_blank(),
        text= element_text(size=12))+
  facet_wrap(~species)


###########################
# FRACTION ABUNDANCE
###########################

#load file index of selected species
load('./ind_sel.RData') #ind_all

#reshape index
ind_all1<-melt(ind_all,id.vars=c('year','species'))

#keep NBS and EBSshelf
ind_all2<-subset(ind_all1,variable %in% c('NBS','EBSshelf'))

#get fraction
ind_all3<-ind_all2 %>% 
  group_by(species, year) %>%
  mutate(frac = value / sum(value))

# Create a dummy data frame for legend rectangles
legend_rects <- data.frame(
  ymin = -Inf,
  ymax = Inf,
  xmin = c(-Inf, 2005.5,2013.5),
  xmax = c(2005.5, 2013.5,Inf),
  fill = c('red', 'blue','red'),
  label = c('Period 1', 'Period 2','Period 1')
)

#plot
p<-
ggplot() +
  geom_rect(data = legend_rects, aes(xmin = xmin, xmax = xmax, 
                                     ymin = ymin, ymax = ymax, fill = label), 
            show.legend = TRUE,alpha=0.8) +
  # annotate("rect", xmin = 2002, xmax = 2005.5, 
  #          ymin = -Inf, ymax = Inf, alpha = 0.5, fill = 'red') +
  # annotate("rect", xmin = 2013.5, xmax = 2022.5, 
  #          ymin = -Inf, ymax = Inf, alpha = 0.5, fill = 'red') +
  # annotate("rect", xmin = 2005.5, xmax = 2013.5, 
  #          ymin = -Inf, ymax = Inf, alpha = 0.5, fill = 'blue') +
  geom_bar(data = ind_all3, aes(x = year, y = frac, fill = variable), 
           stat = 'identity',alpha=0.8) +
  scale_fill_manual(values = c('EBSshelf' = '#046407', 
                               'NBS' = '#B4AF46',
                               'Period 1' = '#cc1d1f',
                               'Period 2' = '#1675ac'),
                    breaks = c('EBSshelf', 'NBS', 'Period 2','Period 1'),
                    labels = c('EBS', 'NBS', 'normal-cold', 'warm'),
                    name = 'region and\nSBT regime') +
  scale_x_continuous(breaks = seq(from = 1985, to = 2020, by = 5),limits = c(2001.5,2022.5),
                     expand = c(0.01,0)) +
  scale_y_continuous(expand = c(0.02,0))+
  labs(y = 'abundance proportion') +
  theme_bw() +
  theme(strip.text = element_text(size = 12),
        strip.background = element_blank(),
        text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.8),
        axis.title.x = element_blank()) +
  facet_wrap(~ species, scales = 'free_x', nrow = 2) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))

#save env plot
agg_png(paste0('./figures slope/abundance_fraction.png'), width = 7, height = 4, units = "in", res = 300)
print(p)
dev.off()

#load files density slope and shelf data 
load('./output slope/species/ms_sim_dens_all.RData')  #sim_dens1
dimnames(sim_dens1)

#subset years and species densities
dens<-sim_dens1[,sel_sp,as.character(2002:2022),]
dimnames(dens)

#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
#load('./extrapolation grids/bering_sea_slope_grid.rda')
#colnames(bering_sea_slope_grid)[4]<-'Stratum'
#bering_sea_slope_grid$Stratum<-NA
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),
                          data.frame(eastern_bering_sea_grid,region='EBS')))
grid$cell<-1:nrow(grid)
grid2<-grid
grid<-grid2[,c('Area_in_survey_km2','region','cell')]

#data input
bio_df<-
expand.grid('sp'=dimnames(dens)[[2]],
            'year'=dimnames(dens)[[3]],
            'sim'=dimnames(dens)[[4]])

bio_df$EBSbio<-NA
bio_df$EBSpct<-NA
bio_df$NBSbio<-NA
bio_df$NBSpct<-NA

for (r in 1:nrow(bio_df)) {
  
  #r<-1
  
  cat(paste('#',r,"-",nrow(bio_df),'\n'))
  
  s<-bio_df[r,'sp']
  y<-bio_df[r,'year']
  sim<-bio_df[r,'sim']
  region<-bio_df[r,'region']
  
  bio<-dens[1:53464,s,y,sim]*grid$Area_in_survey_km2
  
  bio_df[r,'EBSbio']<-sum(bio[15181:53464])
  bio_df[r,'NBSbio']<-sum(bio[1:15180])
  bio_df[r,'EBSpct']<-sum(bio[15181:53464])/(sum(bio[15181:53464])+sum(bio[1:15180]))
  bio_df[r,'NBSpct']<-sum(bio[1:15180])/(sum(bio[15181:53464])+sum(bio[1:15180]))

}

#reshape
bio_df1<-bio_df[,c('sp','year','sim','EBSpct','NBSpct')]
bio_df2<-melt(bio_df1)

#mean over region and simulated data
sim_indmean<-aggregate(value ~ year + variable + sp,bio_df2,FUN=mean)
sim_indsd<-aggregate(value ~ year + variable + sp,bio_df2,FUN=sd)

#rename
names(sim_indmean)[ncol(sim_indmean)]<-'mean'
sim_indmean1<-cbind(sim_indmean,'sd'=sim_indsd$value)
sim_indmean1$year<-as.numeric(sim_indmean1$year)+2001

#sort levels
sim_indmean1$sp<-factor(sim_indmean1$sp,
                                levels=c("Chionoecetes opilio" ,"Gadus chalcogrammus","Gadus macrocephalus","Reinhardtius hippoglossoides" ))

#plot
p<-
  ggplot() +
  geom_rect(data = legend_rects, aes(xmin = xmin, xmax = xmax, 
                                     ymin = ymin, ymax = ymax, fill = label), 
            show.legend = TRUE,alpha=0.8) +
  # annotate("rect", xmin = 2002, xmax = 2005.5, 
  #          ymin = -Inf, ymax = Inf, alpha = 0.5, fill = 'red') +
  # annotate("rect", xmin = 2013.5, xmax = 2022.5, 
  #          ymin = -Inf, ymax = Inf, alpha = 0.5, fill = 'red') +
  # annotate("rect", xmin = 2005.5, xmax = 2013.5, 
  #          ymin = -Inf, ymax = Inf, alpha = 0.5, fill = 'blue') +
  geom_bar(data = sim_indmean1, aes(x = year, y = mean, fill = variable), 
           stat = 'identity',alpha=0.8) +
  geom_errorbar(data = subset(sim_indmean1,variable=='NBSpct'), aes(x = year, ymin = mean-sd, ymax=mean+sd), 
            stat = 'identity',alpha=0.5) +
  scale_fill_manual(values = c('EBSpct' = '#046407', 
                               'NBSpct' = '#B4AF46',
                               'Period 1' = '#cc1d1f',
                               'Period 2' = '#1675ac'),
                    breaks = c('EBSpct', 'NBSpct', 'Period 2','Period 1'),
                    labels = c('EBS', 'NBS', 'normal-cold', 'warm'),
                    name = 'region and\nSBT regime') +
  scale_x_continuous(breaks = seq(from = 1985, to = 2020, by = 5),limits = c(2001.5,2022.5),
                     expand = c(0.01,0)) +
  scale_y_continuous(expand = c(0.02,0))+
  labs(y = 'abundance proportion') +
  theme_bw() +
  theme(strip.text = element_text(size = 12),
        strip.background = element_blank(),
        text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.8),
        axis.title.x = element_blank()) +
  facet_wrap(~ sp, scales = 'free_x', nrow = 2) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))

#save env plot
agg_png(paste0('./figures slope/abundance_fraction_sim.png'), width = 7, height = 4, units = "in", res = 300)
print(p)
dev.off()


###########################
# SIMULATED METRICS
###########################

#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
load('./extrapolation grids/bering_sea_slope_grid.rda')
colnames(bering_sea_slope_grid)[4]<-'Stratum'
bering_sea_slope_grid$Stratum<-NA
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),
                          data.frame(eastern_bering_sea_grid,region='EBS'),
                          data.frame(bering_sea_slope_grid,region='SBS')))
grid$cell<-1:nrow(grid)
grid2<-grid
grid<-grid2[,c('Lat','Lon','region','cell')]

#load files density slope and shelf data 
load('./output slope/species/ms_sim_dens_all.RData')  #sim_dens1
load('./grid_EBS_NBS_envs.RData') #xx

#df to store results
metrics_df<-data.frame(matrix(NA,nrow = 0,ncol=16))
names(metrics_df)<-c("Year","COG_lat","COG_lon","COG_depth","COG_temp",
                     "total_bio",'mean_dens','mean_temp','depthq10','depthq90',"q10","q20","q80",
                     "q90","species",'sim')

#shelf slope spp
#shelfslope_spp<-unique(dens_all$species)

#yrs with slope data()
yrs<-as.character(c(2002,2004,2008,2010,2012,2016))

#subset dens
dens_all<-sim_dens1[,sel_sp,yrs,]

#calculate center of gravity, percentiles and effective area occupied
#loop over common spp in shelf and slope
#for (s in sel_sp) { 

for (y in yrs) { 
    
  #y<-yrs[1]
  
  #dens
  dens<-dens_all[,,y,]  
  
  #grid with envs
  grid.ebs_year1<-subset(xx,Year %in% y)
  grids<-merge(grid,grid.ebs_year1,by=c('Lon','Lat'))
  grids<-grids[order(grids$cell, decreasing = FALSE), ]   

  
  for (sim in 1:100) {
    for (s in sel_sp) {
    
      #sim<-1;s<-sel_sp[1]
      
      cat(paste(y,s,'\n'))

      #merge dens and grid
      dens1<-cbind(dens[,s,sim],grids)
      names(dens1)[1]<-'dens'
      
      #calculate quantiles using tapply
      quantiles <- aggregate(dens ~ Year, data = dens1, 
                             FUN = quantile,c(0.10,0.20,0.80,0.90))
      quantiles<- data.frame(as.matrix(quantiles))
      names(quantiles)<-c('Year','q10','q20','q80','q90')
      
      #get quantiles
      p1 <- quantiles$q10
      p2 <- quantiles$q90
      
      #filter
      dens1$quantile <- ifelse(dens1$dens < p2 & dens1$dens > p1, 1, 0)
      
      # Sort the dataframe by DepthGEBCO
      all_bio4 <- dens1[order(dens1$DepthGEBCO), ]
      
      # Calculate cumulative biomass
      all_bio4$bio<-all_bio4$dens*all_bio4$Area_in_survey_km2
      all_bio4$cumulative_biomass <- cumsum(all_bio4$bio)
      
      # Normalize the cumulative biomass to get a cumulative distribution (0 to 1)
      all_bio4$cumulative_biomass_scale <- all_bio4$cumulative_biomass / 
        max(all_bio4$cumulative_biomass)
      
      # Find the depths corresponding to the 10th and 90th percentiles
      ps <- c(0.10, 0.90)
      depth_10th <- all_bio4$DepthGEBCO[which.min(
        abs(all_bio4$cumulative_biomass_scale - ps[1]))]
      depth_90th <- all_bio4$DepthGEBCO[which.min(
        abs(all_bio4$cumulative_biomass_scale - ps[2]))]
      
      #COG and total abundance
      centroids <- dens1 %>%
        group_by(Year) %>%
        summarise(COG_lat=sum(dens*Lat,na.rm=TRUE)/sum(dens,na.rm=TRUE),
                  COG_lon=sum(dens*Lon,na.rm=TRUE)/sum(dens,na.rm=TRUE),
                  COG_depth=sum(dens*DepthGEBCO,na.rm=TRUE)/sum(dens,na.rm=TRUE),
                  COG_temp=sum(dens*Temp,na.rm = TRUE)/sum(dens,na.rm=TRUE),
                  total_bio=sum(dens*Area_in_survey_km2,na.rm=TRUE),
                  mean_dens=mean(dens,na.rm=TRUE),
                  mean_temp=mean(Temp,na.rm=TRUE),
                  depthq10=depth_10th,
                  depthq90=depth_90th)
      
      #cbind data
      metrics<-cbind(centroids,quantiles[,c('q10','q20','q80','q90')],'species'=s,'sim'=sim)
      
      #append
      metrics_df<-rbind(metrics_df,metrics) 
    }
  }
}

#save table
save(metrics_df,file='./output slope/simulated_metrics.RData')
load(file='./output slope/simulated_metrics.RData') #metrics_df

# Define a custom color scale function
custom_colors <- colorRampPalette(c("#1675ac", "white", "#cc1d1f"))

#selected species
sel_sp<-c("Gadus chalcogrammus", #Alaskan pollock
          "Gadus macrocephalus", #pacific cod
          "Reinhardtius hippoglossoides", #Greenland turbot
          "Chionoecetes opilio") #snow crab

# Load required package
library(dplyr)

# Calculate mean and SD for COG_lat and COG_depth grouped by Year
summary_stats <- metrics_df %>%
  group_by(Year,species) %>%
  summarise(
    mean_COG_lat = mean(COG_lat, na.rm = TRUE),
    mean_temp=mean(mean_temp,na.rm=TRUE),
    sd_COG_lat = sd(COG_lat, na.rm = TRUE),
    mean_COG_depth = mean(COG_depth, na.rm = TRUE),
    sd_COG_depth = sd(COG_depth, na.rm = TRUE),
    mean_depth10 = mean(depthq10, na.rm = TRUE),
    mean_depth90 = mean(depthq90, na.rm = TRUE),
  )

# View the summary
print(summary_stats)

###########################
# COG PLOT
###########################

#plot
p<-
ggplot(data = summary_stats) +
  #geom_errorbar(aes(xmin = mean_COG_depth - sd_COG_depth, xmax = mean_COG_depth + sd_COG_depth, 
  #                   y = mean_COG_lat, color = mean_temp)) +
  # Add vertical error bars for COG_lat (latitudinal uncertainty)
  #geom_errorbar(aes(ymin = mean_COG_lat - sd_COG_lat, ymax = mean_COG_lat + sd_COG_lat, 
  #                  x = mean_COG_depth, color = mean_temp)) +
  geom_shadowtext(aes(x = mean_COG_depth, y = mean_COG_lat, color = mean_temp, label = Year),
                  fontface = 'bold', bg.r = 0.05) +
  theme_bw() +
  labs(x = 'bathymetrical COG (m)', y = 'latitudinal COG') +
  # Add horizontal error bars for COG_depth (bathymetrical uncertainty)
  
  theme(aspect.ratio = 1,
        text = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) +
  scale_color_gradientn(colors = custom_colors(20), name = 'SBT (°C)',
                        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  scale_x_continuous(expand = expansion(mult = 0.1)) +  # Add buffer on the x-axis
  scale_y_continuous(expand = expansion(mult = 0.1)) +  # Add buffer on the y-axis
  facet_wrap(~species, scales = 'free')

#save env plot
agg_png(paste0('./figures slope/cog.png'), width = 6.5, height = 6.5, units = "in", res = 300)
print(p)
dev.off()


#plot using SBT color
p<-
ggplot(data = subset(metrics_df, species %in% sel_sp)) +
  #geom_shadowtext(aes(x = COG_depth, y = COG_lat, color = mean_temp, label = Year),
  #                fontface = 'bold', bg.r = 0.05) +
  geom_point(aes(x = COG_depth, y = COG_lat, fill = mean_temp),shape=21,color=rgb(0, 0, 0, alpha = 0.2)) +
  theme_bw() +
  labs(x = 'bathymetrical COG (m)', y = 'latitudinal COG (°)') +
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) +
  scale_fill_gradientn(colors = custom_colors(100), name = 'SBT (°C)',
                        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  scale_x_continuous(expand = expansion(mult = 0.1)) +  # Add buffer on the x-axis
  scale_y_continuous(expand = expansion(mult = 0.1)) +  # Add buffer on the y-axis
  facet_wrap(~species, scales = 'free')

unique(metrics_df$Year)


agg_png(paste0('./figures slope/cog_sim.png'), width = 6.5, height = 6.5, units = "in", res = 300)
print(p)
dev.off()

ggplot(data = subset(metrics_df, species %in% sel_sp)) +
  #geom_point(aes(x = COG_depth, y = COG_lat, color = mean_temp)) +
  stat_density_2d(aes(x = COG_depth, y = COG_lat,group=Year,color=mean_temp), contour = TRUE) + 
  theme_bw() +
  labs(x = 'bathymetrical COG (m)', y = 'latitudinal COG (°)') +
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) +
  scale_fill_gradientn(colors = custom_colors(100), name = 'SBT (°C)', 
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  scale_color_gradientn(colors = custom_colors(100), name = 'SBT (°C)', 
                        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  facet_wrap(~species, scales = 'free')

#scale to keep one legend
metrics_df1 <- metrics_df %>%
  group_by(species) %>%
  mutate(bio_scaled = scale(total_bio))

#plot using scaled bio color
ggplot(data=subset(metrics_df1,species %in% sel_sp))+
  geom_shadowtext(aes(x=COG_depth,y=COG_lat,color=bio_scaled,label=Year),
                  fontface='bold',bg.r = 0.05)+
  theme_bw()+
  labs(x='bathymetrical COG (m)',y='latitudinal COG')+
  theme(aspect.ratio = 1,strip.background = element_blank(),strip.text = element_text(size=12))+
  scale_color_viridis_c(name = 'relative\ntotal biomass',guide = guide_colorbar(frame.colour = "black", 
                                                                                ticks.colour = "black"))+
  facet_wrap(~species,scales='free')

#define a data frame for the rectangles
rect_data <- data.frame(xmin = c(2002, 2005.5, 2013.5),
                        xmax = c(2005.5, 2013.5, 2016),
                        fill = factor(c("warm", "cold", "warm")))

###########################
# EFFECTIVE AREA OCCUPIED
###########################

library(ggplot2)

# Create a data frame for the rectangles
rect_data <- data.frame(
  xmin = c(-Inf, 2013.5, 2005.5),
  xmax = c(2005.5, Inf, 2013.5),
  ymin = rep(-Inf, 3),
  ymax = rep(Inf, 3),
  fill = c("warm", "warm", "normal-cold")  # Labels for the legend
)

p <- 
  ggplot(data = metrics_df) + 
  # Use geom_rect to add rectangles with fill aesthetics
  geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), alpha = 0.5) +
  #geom_line(aes(x = Year, y = total_bio / mean_dens / 1000, group = sim), color = 'black', alpha = 0.3) +  # Convert to thousands
  #geom_point(aes(x = Year, y = total_bio / mean_dens / 1000, group = sim), color = 'black', alpha = 0.3) +  # Convert to thousands
  geom_boxplot(aes(x = Year, y = total_bio / mean_dens / 1000, group = Year), color = 'black', alpha = 0.3) +  # Convert to thousands
    
  labs(x = '', y = 'effective area occupied (thousands km²)') +  # Adjust y-axis label
  theme_bw() + 
  scale_x_continuous(breaks = bold_years,minor_breaks = c(2002:2016), expand = c(0, 0.5)) + 
  scale_y_continuous(labels = comma) +  # Add commas to y-axis numbers
  theme(aspect.ratio = 1, 
        strip.background = element_blank(), 
        strip.text = element_text(size = 12), 
        axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.8),
        text = element_text(size = 12)) +   
  facet_wrap(~species, nrow = 2, scales = 'free') + 
  scale_fill_manual(values = c("warm" = "#cc1d1f", "normal-cold" = "#1675ac"), 
                    name = "SBT regime") + 
  theme(legend.position = 'right')  # Position the legend on the right


#save env plot
agg_png(paste0('./figures slope/effective_area_sim.png'), width = 7.5, height = 7, units = "in", res = 300)
print(p)
dev.off()

###########################
# INTERDECILE DEPTH
###########################

#plot depth range
p <- 
ggplot(data = metrics_df) +
  #geom_point(aes(x = COG_depth, y = COG_lat, fill = mean_temp),shape=21,color=rgb(0, 0, 0, alpha = 0.2)) +
  
  geom_point(aes(x = depthq90-depthq10, y = species, fill = mean_temp),alpha=0.5,  # Correct placement of aes()
  #geom_point(aes(x = mean_depth90-mean_depth10, y = species, fill = mean_temp),  # Correct placement of aes()
                    shape = 21, size = 3,color=rgb(0, 0, 0, alpha = 0.2)) + #,alpha=0.3
  facet_wrap(~species,scales = 'free')+
  scale_fill_gradientn(colors = custom_colors(20),name = 'SBT (°C)',
                       guide = guide_colorbar(frame.colour = "black", 
                                              ticks.colour = "black"))+
  theme_bw()+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), 
    axis.title.y = element_blank(),
    strip.background = element_blank(),text = element_text(size=12),
    strip.text = element_text(size=12))+
  labs(x = 'depth at 90th percentile - depth at 10th percentile (m)')

#save env plot
agg_png(paste0('./figures slope/interdecile_depth_sim.png'), width = 7, height = 3, units = "in", res = 300)
print(p)
dev.off()
