
library(VAST)

load('/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/slope EBS VAST/Gadus macrocephalus/fit.RData')


fit$Report$Index_ctl
dim(fit$Report$Index_ctl)
dim(fit$Report$Index_gctl)
sum(fit$Report$Index_gctl[,,14,])

simulate_data()

#######PALETTE

# Sunrise Palette with 6 colors
sunrise_palette <- c("#FFD3B6", "#FFB88C", "#FF9A66", "#FF7B43", "#FF572C", "#FF3C1E")

# Twilight Palette with 6 colors
twilight_palette <- c("#0B0033", "#1C3A5A", "#37677A", "#5A8EAB", "#85B4D8", "#B2E1E8")




#####################################
# Polygon regions shapefiles (EBS, NBS and slope)
#####################################


#name shapefiles 
shfiles<-c('EBSshelfThorson','NBSThorson','EBSslopeThorson')

#loop over shapefiles
for (i in shfiles) {
  
  #shapefile EBS
  sh<-rgdal::readOGR(dsn='./shapefiles/',layer = i)
  
  if (i=='EBSslopeThorson') {
    
    #reproject shapefile
    proj4string(sh) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
    sh<-spTransform(sh,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))
    
  }
  
  #shapefile name
  shname<-paste0(gsub('Thorson','',i),'_sh')
  
  #assign shapefiles
  assign(shname,sh)
  
}

#merge shapefiles
bs_sh1<-raster::union(EBSshelf_sh,NBS_sh)
bs_sh<-raster::union(bs_sh1,EBSslope_sh)


#########
# correlation between density and temp across species

load(paste0('/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/output/multisp_optimization_static_data_ebsnbs_slope.RData'))

head(df)
names_spdens<-c(names(df[grepl('_sumDensity*$',names(df))]),'sumDensity')

cor_df<-data.frame(matrix(NA,0,3))
names(cor_df)<-c('env','sp','cor')


for (s in 1:length(names_spdens)) {
  
  #s<-1
  
  cor_df<-rbind(cor_df,c('meanTemp',names_spdens[s],round(cor(df$meanTemp,df[,names_spdens[s]]),digits = 3)))
  cor_df<-rbind(cor_df,c('varTemp',names_spdens[s],round(cor(df$varTemp,df[,names_spdens[s]]),digits = 3)))
  cor_df<-rbind(cor_df,c('depth',names_spdens[s],round(cor(df$Depth,df[,names_spdens[s]]),digits = 3)))
  
}

names(cor_df)<-c('env','sp','cor')



#########
# calculate the change on center of gravity over years and species  
setwd('/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/')

#load slope grid
load('./extrapolation grids/bering_sea_slope_grid.rda')
dim(bering_sea_slope_grid)
names(bering_sea_slope_grid)[4]<-'Stratum'
bering_sea_slope_grid$Stratum<-999
#gridslope<-data.frame(bering_sea_slope_grid,region='SLP')

#load EBS+NBS grid
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),
                          data.frame(eastern_bering_sea_grid,region='EBS'),
                          data.frame(bering_sea_slope_grid,region='SLP')))
grid$cell<-1:nrow(grid)
grid$cell<-as.numeric(grid$cell)

load('./data processed/grid_EBS_NBS.RData')
grid.ebs_year1<-unique(grid.ebs_year[,c("DepthGEBCO","region","Area_in_survey_km2","Lat","Lon")])
aggregate(Lat+Lon~ region, grid.ebs_year1,FUN=length)
grid.ebs_year2<-subset(grid.ebs_year1,region %in% c('EBSshelf','NBS'))

#####################
# with simdata
#####################

conv_tab<-read.csv('./tables/slope_ebsnbs_convspp.csv')
sp_shelfslope<-conv_tab[which(conv_tab$slope=='There is no evidence that the model is not converged' & conv_tab$EBS_NBS=='There is no evidence that the model is not converged'),'spp']

load('/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/output/species/ms_sim_dens.RData') #index_hist
sim_dens1_shelf<-sim_dens1
load('/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/output/species/ms_sim_dens_slope.RData') #index_hist
sim_dens1_slope<-sim_dens1
dim(sim_dens1_shelf);dim(sim_dens1_slope)

dimnames(sim_dens1)

cog_df<-data.frame(matrix(NA,nrow = 0,ncol=6))
names(cog_df)<-c('sp','sim','y','lat','lon', 'region')


for (s in sp_shelfslope) {
  
  #s<-sp_shelfslope[1];si<-dimnames(sim_dens1_shelf)[[4]][1];y<-dimnames(sim_dens1_slope)[[3]][1]
  
  for (si in dimnames(sim_dens1_shelf)[[4]]) {
  
    cat(paste('#################### ',s," - ",si,'\n'))
    
    for (y in dimnames(sim_dens1_slope)[[3]]) {
  
      
    sim_dens2<-c(sim_dens1_shelf[,s,y,si],sim_dens1_slope[,s,y,si])
    
      for (r in c('EBSshelf','NBS','EBSslope','all')) {
      
        if (r=='EBSshelf') {
          c<-which(grid.ebs_year1$region=='EBSshelf')
        } else if (r=='NBS') {
          c<-which(grid.ebs_year1$region=='NBS')
        } else if (r=="EBSslope") {
          c<-which(grid.ebs_year1$region=='EBSslope')
        } else if (r=='all') {
          c<-which(!is.na(grid.ebs_year1$region))
        }
  
        igrid<-grid.ebs_year1[c,]
        sim_dens3<-sim_dens2[c]
        
        cog_df<-rbind(cog_df,
                      data.frame('sp'=s,
                                 'sim'=si,
                                 'year'=y,
                                 'depth'=sum(sim_dens3*igrid$DepthGEBCO)/sum(sim_dens3),
                                 'lat'=sum(sim_dens3*igrid$Lat)/sum(sim_dens3),
                                 'lon'=sum(sim_dens3*igrid$Lon)/sum(sim_dens3),
                                 'region'=r))
        
      }
    }
  }
}

save(cog_df,file='./output/cog_spp.RData')
#load(file='./output/cog_spp.RData') #cog_df,

####################
# temp by region
####################

load(file = './data processed/grid_EBS_NBS.RData')

temp_region<-
  aggregate(Temp ~ Year + region,grid.ebs_year,FUN=mean)
names(temp_region)<-c('year','region','temp')

alltemp<-aggregate(temp ~ year, temp_region,FUN=mean)
alltemp$region<-'all'
alltemp<-alltemp[,c("year","region","temp")]

temp_region<-rbind(temp_region,alltemp)

###################
# merge 
###################

cog_df1<-merge(cog_df,temp_region,by=c('year','region'))

#warm<-c(2002:2005,2014:2016)
#cold<-c(2006:2013)
cog_df1$coloryear<-ifelse(cog_df1$year %in% c(2002:2005,2014:2016),'red','blue')

###################
# plot region temp vs depth cog over species and regions 
###################

ggplot()+
  geom_point(data=cog_df1,aes(x=depth,y=temp,color=as.factor(year)))+
  facet_wrap(region~sp,scales='free',nrow=3)+
  labs(x='bathymetric center of gravity',y='mean region SBT')+
  scale_color_manual(values = coloryear, name='year')+
  geom_smooth(data=cog_df1,aes(x=depth,y=temp),method = "lm", se = FALSE,color='black')  # Add a smoothed line with confidence intervals ,formula = y ~ x + I(x^2)
  


ggplot() +
  geom_point(data = cog_df1, aes(x = depth, y = temp, color = as.factor(year)),alpha=0.6) +
  facet_wrap(region ~ sp, scales = 'free', nrow = 4) +
  labs(x = 'bathymetric center of gravity', y = 'mean region SBT') +
  scale_color_manual(values = unique(cog_df1[,c("coloryear","year")])[,'coloryear'], 
                     breaks = unique(cog_df1[,c("coloryear","year")])[,'year'],  # Specify the breaks based on unique year values
                     labels = unique(cog_df1[,c("coloryear","year")])[,'year'],  # Use unique year values as labels
                     name = 'year') +
  theme_bw()+
  #geom_smooth(data = cog_df1, aes(x = depth, y = temp), method = "lm", se = FALSE, color = 'black',alpha=0.9)
  geom_line(data = cog_df1, aes(x = depth, y = temp),stat="smooth",method = "lm", 
            linetype ="dashed",
            color='black')
  
 


library(ggplot2)

# Define plot extent (through trial end error) units km
panel_extent <- data.frame(x = c(-1716559.21, -77636.05), #x = c(-1326559.21, -87636.05),
                           y = c(483099.5, 2194909.7)) #y = c(533099.5, 1894909.7))


#Alaska land shapefile layer
ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")
ak_sppoly<-as(ebs_layers$akland, 'Spatial')

library(sp)

for (s in unique(cog_df1$sp)) {
  
  #s<-unique(cog_df1$sp)[1]
  
  cog_df2<-subset(cog_df1,sp==s)
  
  coordinates(cog_df2) <- ~ lon + lat
  proj4string(cog_df2)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  #proj4string(strata_pol) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs') 
  cog_df2<-spTransform(cog_df2,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))
  cog_df2<-as.data.frame(cog_df2,xy=TRUE)
  cog_df2$year<-as.numeric(cog_df2$year)
  
  print(
  ggplot()+
    geom_point(aes(x=mean(cog_df2$lon),y=mean(cog_df2$lat)),shape=3,size=4)+
    geom_polygon(data=bs_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
    geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),color='black',linewidth=0.2,fill = 'grey80')+
    geom_point(data=cog_df2,aes(x=lon,y=lat,color=as.factor(year)))+
    scale_color_manual(values = sunrise_palette, name='year')+
    geom_path(data=cog_df2,aes(x=lon,y=lat))+
    theme_bw()+
    theme(aspect.ratio = 1,panel.grid.major = element_line(color = rgb(0, 0, 0,20, maxColorValue = 285), linetype = 'dashed', linewidth =  0.5),
          panel.background = element_rect(fill = NA),panel.ontop = TRUE,text = element_text(size=12),
          legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(25, 'points'),
          legend.key.width= unit(25, 'points'),axis.title = element_blank(),legend.position = c(0.12,0.47), #c(0.12,0.47)
          panel.border = element_rect(fill = NA, colour = 'black'),legend.key = element_rect(color="black"),
          axis.text = element_text(color='black'),legend.spacing.y = unit(10, 'points'),
          axis.text.y.right = element_text(hjust= 0.1 ,margin = margin(0,7,0,-25, unit = 'points'),color='black'),
          axis.text.x = element_text(vjust = 6, margin = margin(-7,0,7,0, unit = 'points'),color='black'),
          axis.ticks.length = unit(-5,"points"))+
    scale_x_continuous(expand = c(0,0),position = 'bottom',
                       breaks = c(-180,-175,-170,-165,-160,-155),sec.axis = dup_axis())+
    coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
             xlim = panel_extent$x,
             ylim = panel_extent$y,
             label_axes = "-NE-")+
    labs(title=paste0(s),xlab='',ylab='')
  )
    
  
}


#####################
# with index_gctl
#####################

# Sunrise Palette with 15 colors
sunrise_palette1 <- c("#FFDAC1", "#FFD3B6", "#FFCCAB", "#FFC49F", "#FFBC94", "#FFB48A", "#FFAB7F", "#FFA275", "#FF9A6A", "#FF9260", "#FF8A55", "#FF824B", "#FF7B43", "#FF7338", "#FF6B2E")

# Twilight Palette with 15 colors
twilight_palette1 <- c("#0B0033", "#11224A", "#183B60", "#1D4E72", "#225F82", "#276F91", "#2C7F9F", "#338FAB", "#3B9EB6", "#449CC0", "#4F9AC8", "#5A98CF", "#6596D6", "#7194DC", "#7E92E2")




cog_df1<-data.frame(matrix(NA,nrow = 0,ncol=5))
names(cog_df1)<-c('sp','sim','y','lat','lon')


for (s in sp_shelfslope) {
  
  #s<-sp_shelfslope[1];si<-dimnames(sim_dens1_shelf)[[4]][1];y<-dimnames(sim_dens1_slope)[[3]][1]
  
    cat(paste('#################### ',s," - ",'\n'))
    
  load(paste0('./shelf EBS NBS VAST/',s,'/fit.RData'))
  
  fit_shelf<-fit
  dens_shelf<-fit_shelf$Report$Index_gctl[,,,1]
  
  load(paste0('./slope EBS VAST/',s,'/fit.RData'))
  
  fit_slope<-fit
  dens_slope<-fit_slope$Report$Index_gctl[,,,1]
  
    for (y in colnames(dens_slope)) {
      
      
      all_dens2<-c(dens_shelf[,y],dens_slope[,y])
      
      cog_df<-rbind(cog_df,
                    data.frame('sp'=s,
                               'sim'=si,
                               'year'=y,
                               'lat'=sum(all_dens2*grid$Lat)/sum(all_dens2),
                               'lon'=sum(all_dens2*grid$Lon)/sum(all_dens2)))  
      
    }
  
  cog_df1<-rbind(cog_df,cog_df1)
}
  
list_plot<-list()

for (s in sp_shelfslope) {
  
  #s<-unique(cog_df1$sp)[1]
  
  cog_df2<-subset(cog_df1,sp==s)
  
  coordinates(cog_df2) <- ~ lon + lat
  proj4string(cog_df2)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  #proj4string(strata_pol) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs') 
  cog_df2<-spTransform(cog_df2,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))
  cog_df2<-as.data.frame(cog_df2,xy=TRUE)
  cog_df2$year<-as.numeric(cog_df2$year)
  
  #print(
  pp<-
    ggplot()+
      #geom_path(data=cog_df2,aes(x=lon,y=lat))+
      #geom_point(aes(x=mean(cog_df2$lon),y=mean(cog_df2$lat)),shape=3,size=4)+
      geom_polygon(data=bs_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
      #geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),color='black',linewidth=0.2,fill = 'grey80')+
      geom_point(data=cog_df2,aes(x=lon,y=lat,color=year))+
      scale_x_continuous(expand = c(0,0),position = 'bottom',
                         breaks = c(-180,-175,-170,-165,-160,-155),sec.axis = dup_axis())+
      theme_bw()+
      scale_color_viridis_c(option = 'magma', name='year',direction = -1)+
      coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
               xlim = c(-1422205, -194348.5),
               ylim = c(547256.1, 1813031),
               label_axes = "-NE-")+
      theme(aspect.ratio = 1,panel.grid.major = element_blank(),
            panel.background = element_rect(fill = NA),panel.ontop = TRUE,
            legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(20, 'points'),
            legend.key.width= unit(20, 'points'),axis.title = element_blank(),legend.position = 'none',
            panel.border = element_rect(fill = NA, colour = NA),legend.key = element_rect(color="black"),
            legend.spacing.y = unit(8, 'points'),
            axis.text=element_blank(),axis.ticks = element_blank(),
            plot.margin = margin(0.01,0.01,0.01,0.01), 
            axis.ticks.length = unit(-5,"points"),plot.title = element_text(size=18,hjust = 0.5,vjust=-5, face="bold"))+
      labs(title=paste0(s))#)
           
    list_plot[[s]]<-pp
           
}

legend<-
  ggplot()+
  #geom_path(data=cog_df2,aes(x=lon,y=lat))+
  #geom_point(aes(x=mean(cog_df2$lon),y=mean(cog_df2$lat)),shape=3,size=4)+
  geom_polygon(data=bs_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  #geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),color='black',linewidth=0.2,fill = 'grey80')+
  geom_point(data=cog_df2,aes(x=lon,y=lat,color=year))+
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-180,-175,-170,-165,-160,-155),sec.axis = dup_axis())+
  theme_bw()+
  scale_color_viridis_c(option = 'magma', name='year',direction = -1)+
  coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
           xlim = c(-1422205, -194348.5),
           ylim = c(547256.1, 1813031),
           label_axes = "-NE-")+
  theme(aspect.ratio = 1,panel.grid.major = element_blank(),
        panel.background = element_rect(fill = NA),panel.ontop = TRUE,
        legend.background =  element_rect(fill = "transparent", colour = "transparent"),legend.key.height= unit(20, 'points'),
        legend.key.width= unit(20, 'points'),axis.title = element_blank(),legend.position = 'none',
        panel.border = element_rect(fill = NA, colour = NA),legend.key = element_rect(color="black"),
        legend.spacing.y = unit(8, 'points'),
        axis.text=element_blank(),axis.ticks = element_blank(),
        plot.margin = margin(0.01,0.01,0.01,0.01), 
        axis.ticks.length = unit(-5,"points"),plot.title = element_text(size=18,hjust = 0.5,vjust=-5, face="bold"))+
  labs(title=paste0(s))

legend1 <- cowplot::get_legend( 
  legend + 
    theme(legend.position = "right",legend.key.height =  unit(4, 'cm'),legend.key.width = unit(1, 'cm')) 
) 

pgrid1<-cowplot::plot_grid(plotlist = list_plot, nrow = 2)

#pgridi1<-cowplot::plot_grid(plotlist = plot_list_d[c(1:7)],nrow=1)
pgridi1<-cowplot::plot_grid(pgrid1,legend1, ncol = 2, rel_widths =  c(1, .14))
pgridi1

#save plots
ragg::agg_png(paste0('./figures slope/cog.png'), width = 20, height = 7, units = "in", res = 300)
print(cowplot::plot_grid(pgrid1,legend1, ncol = 2, rel_widths =  c(1, .14)))
dev.off()

######################
# warm vs cold years
######################

load(file = './data processed/grid_EBS_NBS.RData')

temp_region<-
aggregate(Temp ~ Year + region,grid.ebs_year,FUN=mean)

temp_all<-
  aggregate(Temp ~ Year,grid.ebs_year,FUN=mean)


temp_region$scaleTemp<-scale(temp_region$Temp)

temp_region1<-
temp_region %>% group_by(region) %>% mutate(scaleTemp = scale(Temp))

temp_region4<-
  aggregate(scale(Temp) ~ region,grid.ebs_year,FUN=sd)

temp_all$scaleTemp<-scale(temp_all$Temp)

ggplot()+
  geom_rect(aes(xmin=2002,xmax=2016,ymin=-Inf,ymax=Inf),alpha=0.4)+
  geom_hline(yintercept = 0,alpha=0.5)+
  theme_bw()+
  geom_point(data=temp_all,aes(x=Year,y=scaleTemp))

ggplot()+
  geom_rect(aes(xmin=2002,xmax=2016,ymin=temp_region4$V1[1],ymax=-temp_region4$V1[1]),alpha=0.4)+
  geom_point(data=temp_region1,aes(x=Year,y=scaleTemp,color=region))+
  scale_color_manual(values=c('EBSslope'='#B4AF46','EBSshelf'='#B4464B','NBS'='#4682B4'))+
  geom_hline(yintercept = 0,alpha=0.5)+
  theme_bw()+
  geom_hline(yintercept = 0,alpha=0.5)+
  geom_vline(xintercept = c(2003:2005,2015,2016),linetype='dashed',color='red',alpha=0.5)+
  geom_vline(xintercept = c(2008,2009,2012),linetype='dashed',color='blue',alpha=0.5)


temp_region2<-subset(temp_region1,Year %in% (2002:2016) & region == 'EBSshelf')
print(temp_region2,n=30)
temp_region2<-temp_region2[order(temp_region2$scaleTemp),]

# Generate color palette from blue to red
color_palette <- colorRampPalette(c("blue", "red"))

# Add a column of color codes
temp_region2$ColorCode <- color_palette(nrow(temp_region2))
#warm<-c(2002:2005,2014:2016)
#cold<-c(2006:2013)

#####################
# with simdata
#####################

conv_tab<-read.csv('./tables/slope_ebsnbs_convspp.csv')
sp_shelfslope<-conv_tab[which(conv_tab$slope=='There is no evidence that the model is not converged' & conv_tab$EBS_NBS=='There is no evidence that the model is not converged'),'spp']
sp_shelf<-conv_tab[which(conv_tab$EBS_NBS=='There is no evidence that the model is not converged'),'spp']
sp_slope<-conv_tab[which(conv_tab$slope=='There is no evidence that the model is not converged' ),'spp']

load('/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/output/species/ms_sim_dens.RData') #index_hist
sim_dens1_shelf<-sim_dens1
load('/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/output/species/ms_sim_dens_slope.RData') #index_hist
sim_dens1_slope<-sim_dens1
dim(sim_dens1_shelf);dim(sim_dens1_slope)






dimnames(sim_dens1)

dens_df<-data.frame(matrix(NA,nrow = 0,ncol=6))
names(dens_df)<-c('sp','sim','y','dens','depth','region')


for (s in sp_shelfslope) {
  
  #s<-sp_shelfslope[1];si<-dimnames(sim_dens1_shelf)[[4]][1];y<-dimnames(sim_dens1_slope)[[3]][1]
  
  for (si in dimnames(sim_dens1_shelf)[[4]]) {
    
    cat(paste('#################### ',s," - ",si,'\n'))
    
    for (y in dimnames(sim_dens1_slope)[[3]]) {
      
      
      #sim_dens2<-c(sim_dens1_shelf[,s,y,si],sim_dens1_slope[,s,y,si])
      for (r in c('EBSshelf','NBS','EBSslope','all')) {
        
        if (r=='EBSshelf') {
          c<-which(grid.ebs_year1$region=='EBSshelf')
        } else if (r=='NBS') {
          c<-which(grid.ebs_year1$region=='NBS')
        } else if (r=="EBSslope") {
          c<-which(grid.ebs_year1$region=='EBSslope')
        } else if (r=='all') {
          c<-which(!is.na(grid.ebs_year1$region))
        }
        
        igrid<-grid.ebs_year1[c,]
        sim_dens3<-sim_dens2[c]
        
        
      dens_df<-rbind(dens_df,
                    data.frame('sp'=s,
                               'sim'=si,
                               'year'=y,
                               'bio'=sim_dens3,
                               'depth'=grid.ebs_year$DepthGEBCO[c],
                               'region'=grid.ebs_year1$region[c])) 
      }
    }
  }
}

save(dens_df,file='./output/dens_df.RData')
load(file='./output/dens_df.RData') #dens_df




#####################
# with index_gctl
#####################
#distribution abundance over depth by region

setwd('/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/')

#load slope grid
load('./extrapolation grids/bering_sea_slope_grid.rda')
dim(bering_sea_slope_grid)
names(bering_sea_slope_grid)[4]<-'Stratum'
bering_sea_slope_grid$Stratum<-999
#gridslope<-data.frame(bering_sea_slope_grid,region='SLP')

#load EBS+NBS grid
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),
                          data.frame(eastern_bering_sea_grid,region='EBS'),
                          data.frame(bering_sea_slope_grid,region='SLP')))
grid$cell<-1:nrow(grid)
grid$cell<-as.numeric(grid$cell)

library(raster);library(ggplot2);library(VAST)

#read raster GEBCO data
r<-raster('./bathymetry/gebco_2022_n70.0_s50.0_w-180.0_e-155.0.asc')

#extract depth values for each station of grid - using GEBCO data
rr<-extract(r, SpatialPoints(cbind(grid$Lon,grid$Lat)))
grid$DepthGEBCO<--rr

ggplot()+
  geom_histogram(data=grid,aes(y=DepthGEBCO,color=region))+
  facet_wrap(~region+year)


load('./data processed/grid_EBS_NBS.RData')

ggplot()+
  geom_histogram(data=grid.ebs_year,aes(y=DepthGEBCO,color=region))+
  facet_wrap(~region+Year)

grid1<-subset(grid,DepthGEBCO>=1)
nrow(grid1)*length(unique(grid.ebs_year$Year))
nrow(grid.ebs_year)/length(unique(grid.ebs_year$Year))

for (s in sp_shelfslope) {
  
  s<-sp_shelfslope[1];si<-dimnames(sim_dens1_shelf)[[4]][1];y<-dimnames(sim_dens1_slope)[[3]][1]
  
  cat(paste('#################### ',s," - ",'\n'))
  
  load(paste0('./shelf EBS NBS VAST/',s,'/fit.RData'))
  
  fit_shelf<-fit
  dens_shelf<-fit_shelf$Report$Index_gctl[,,,1]
  
  load(paste0('./slope EBS VAST/',s,'/fit.RData'))
  
  fit_slope<-fit
  dens_slope<-fit_slope$Report$Index_gctl[,,,1]
  
  D_d_ebsnbs<-
  cbind(fit_shelf$Report$D_gct,subset(grid,region %in% c('EBS','NBS')))
  D_d_ebsnbs<-drop_units(D_d_ebsnbs)
  names(D_d_ebsnbs)[1:length(1982:2022)]<-c(paste0('y',1982:2022))
  (D_d_ebsnbs)[c(1:length(1982:2022),ncol(D_d_ebsnbs))]
  D_d_ebsnbs1<-reshape2::melt((D_d_ebsnbs)[c(1:length(1982:2022),ncol(D_d_ebsnbs)-2,ncol(D_d_ebsnbs))], id=c('DepthGEBCO','region'))
  
  D_d_slope<-
    cbind(fit_slope$Report$D_gct,subset(grid,region %in% c('SLP')))
  D_d_slope<-drop_units(D_d_slope)
  names(D_d_slope)[1:length(2002:2016)]<-c(paste0('y',2002:2016))
  D_d_slope1<-reshape2::melt((D_d_slope)[c(1:length(2002:2016),ncol(D_d_slope)-2,ncol(D_d_slope))], id=c('DepthGEBCO','region'))
  
  
  D_d1<-rbind(D_d_ebsnbs1,D_d_slope1)
  
  plot_list<-list()
  
  for (reg in unique(D_d1$region)) {

  #reg=='NBS'
  
  if (reg %in% c('NBS','EBS')) {
    D_d2<-subset(D_d1,region==reg & variable %in% paste0('y',2002:2016) & value != 0 & DepthGEBCO <=250 & DepthGEBCO >=1)
    } else if (reg == 'SLP') {
      D_d2<-subset(D_d1,region==reg & variable %in% paste0('y',2002:2016) & value != 0 & DepthGEBCO <=500  & DepthGEBCO >=150)
    }
    
  #D_d2<-subset(D_d1,region==reg & variable %in% paste0('y',2002:2016) & value != 0 & DepthGEBCO )
  names(D_d2)<-c('Depth','Region','Year','Biomass')
  D_d2$Year<-as.numeric(gsub('y','',D_d2$Year))
  
    
  # D_d3<-subset(D_d1,region==reg & variable %in% paste0('y',2002:2016) )
  # D_d3$presence<-ifelse(D_d3$value==0,0,1)
  
  # # Create a scatterplot with a smoothed line
  # ggplot(D_d2, aes(x = Depth, y = Biomass, color = as.factor(Year))) +
  #   geom_point(alpha=0.1) +  # Add points
  #   geom_smooth(method = "lm", se = TRUE,formula = y ~ x + I(x^2)) +  # Add a smoothed line with confidence intervals
  #   labs(x = "Depth", y = "Biomass",title=paste0(s,r)) +  # Label axes
  #   theme_minimal()
  
  #p<-# Create a scatterplot with a smoothed line
  ggplot(D_d2, aes(x = Depth, y = Biomass, color = Year, group=Year)) +
    geom_point(alpha=0.1,color='black') +  # Add points
    geom_smooth(method = "lm", se = FALSE,formula = y ~ x + I(x^2)) +  # Add a smoothed line with confidence intervals
    labs(x = "Depth", y = "Biomass",title=paste0(s,' - ',reg)) +  # Label axes
    theme_minimal()+
    scale_color_viridis_c(option = 'magma', name='year',direction = -1)#+
    #facet_wrap(~Year,scales='free')#+
    #guides(color=guide_legend())
  
  
  ggplot(D_d2, aes(x = Depth, y = Biomass, color = as.factor(Year), group=as.factor(Year))) +
    geom_point(alpha=0.1) +  # Add points
    geom_smooth(method = "lm", se = FALSE,formula = y ~ x + I(x^2)) +  # Add a smoothed line with confidence intervals
    labs(x = "Depth", y = "Biomass",title=paste0(s,' - ',reg)) +  # Label axes
    theme_minimal()+
    scale_color_viridis_d(option = 'magma', name='year',direction = -1)
  
  
  # Convert depth data to factor with depth intervals
  depth_intervals <- cut(D_d2$Depth, breaks = c(0,50,100,150,200,250,300,350,400,450))
  D_d2$depthbin<-depth_intervals

  #p<-
  ggplot(na.omit(D_d2), aes(x = Year, y = Biomass, color=Year,group=Year)) +
    geom_boxplot(position='dodge')+
    facet_wrap(.~depthbin,scales = 'free',nrow=1)+
    scale_color_manual(temp_region2=)
    
  
    
  plot_list[[reg]]<-p

    }
  }

cowplot::plot_grid(plotlist = plot_list,nrow=3)  

}


# Load necessary libraries
library(ggplot2)

# Fit logistic regression model
model <- glm(presence ~ DepthGEBCO, data = D_d3, family = binomial)

# Predict probabilities
D_d3$predicted_prob <- predict(model, type = "response")

# Plotting
ggplot(D_d3, aes(x = depth, y = predicted_prob)) +
  geom_point(aes(color = year)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  labs(x = "Depth", y = "Predicted Probability of Presence", color = "Year") +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal()


#########
# calculate fraction of sp by calculating their abundance index

library(VAST)

#for (s in sp_shelfslope) {
for (s in sp_shelfslope) {
  
  s<-sp_shelfslope[1]
  
  cat(paste('#################### ',s," - ",'\n'))
  

load(paste0('./shelf EBS NBS VAST/',s,'/fit.RData'))
  
fit_shelf<-fit

plot(fit$Report$Index_ctl)

ggplot()+
  #geom_point(aes(y=as.vector(fit_shelf$Report$Index_ctl[,,1]),x=1982:2022),color='black')+
  geom_point(aes(y=as.vector(fit_shelf$Report$Index_ctl[,,2]),x=1982:2022),color='green')+
  geom_point(aes(y=as.vector(fit_shelf$Report$Index_ctl[,,3]),x=1982:2022),color='red')+
  geom_point(aes(y=as.vector(colSums(fit_shelf$Report$Index_gctl[,1,,1])),x=1982:2022),shape=3)+
  theme_bw()+
  ylab(label = 'abundance')+
  xlab(label = '')
  
load(paste0('./slope EBS VAST/',s,'/fit.RData'))

fit_slope<-fit

ggplot()+
  #geom_point(aes(y=as.vector(fit_shelf$Report$Index_ctl[,,1]),x=1982:2022),color='black')+
  geom_point(aes(y=as.vector(fit_shelf$Report$Index_ctl[,,2]),x=1982:2022),color='green')+
  geom_point(aes(y=as.vector(fit_shelf$Report$Index_ctl[,,3]),x=1982:2022),color='red')+
  geom_point(aes(y=as.vector(colSums(fit_slope$Report$Index_gctl[,1,,1])),x=2002:2016),color='blue')+
  theme_bw()+
  geom_vline(xintercept=2016)+
  geom_vline(xintercept=2002)+
  ylab(label = 'abundance')+
  xlab(label = '')

 
abundance_region<-
data.frame(
  'ebs_nbs'=fit_shelf$Report$Index_ctl[,,1],
  'ebs'=fit_shelf$Report$Index_ctl[,,3],
  'nbs'=fit_shelf$Report$Index_ctl[,,2],
  'slope'=c(rep(NA,20),colSums(fit_slope$Report$Index_gctl[,1,,1]),rep(NA,6)))

abundance_region<-drop_units(abundance_region)
abundance_region$ebs_nbs_slope<-abundance_region$ebs_nbs+abundance_region$slope
abundance_region_percent<-abundance_region/abundance_region$ebs_nbs_slope*100
abundance_region_percent1<-reshape2::melt(as.matrix(abundance_region_percent))

ggplot()+
  geom_line(data = abundance_region_percent1,aes(x=Var1,y=value,color=Var2))
}

#########
# corr_df

load('/Users/daniel/Work/Adapting Monitoring to a Changing Seascape/output/full_spearman_1-8.RData')
head(corr_df)
aggregate(rho ~ spp,corr_df,FUN=mean)

library(dplyr)

corr_df %>%
  group_by(spp) %>%
  #summarize(count_high_pvalues = sum(pvalue > 0.05, na.rm = TRUE))
  summarize(percent_above = mean(pvalue > 0.05, na.rm=TRUE) * 100)


