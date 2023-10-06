#setwd
out_dir<-'C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/'
setwd(out_dir)

#load grid of NBS and EBS
load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
x<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
x$cell<-1:nrow(x)

x1<-x[,c('Lon','Lat','cell')]
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



######################################################


# x<-rgdal::readOGR('./shapefiles/regridadjacentcells/bs_grid.shp')
# y<-rgdal::readOGR('./shapefiles/regridadjacentcells/bs_grid_w_corners.shp')
# 
# plot(x)
# plot(y)


#load('./output/baseline_strata.RData')
#baseline_strata$locations
#plot(baseline_strata$locations$longitude,baseline_strata$locations$latitude)

load('./output/species/Gadus macrocephalus/optimization data/optimization_results_scn16.RData')

result_list$sample_allocations
result_list$sol_by_cell

x6<-merge(result_list$sol_by_cell,x5,by.x='ID',by.y='z',all.x=TRUE)
coordinates(x6)<-~x+y
#reproject shapefile
proj4string(x6) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs') 
x7<-spTransform(x6,CRSobj = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
x7<-as.data.frame(x7)

#iterations for samples
n_iter<-100

#to store samples
all_points<-array(NA,dim = list(sum(result_list$sample_allocations),4,n_iter),
                  dimnames = list(c(1:sum(result_list$sample_allocations)),c('Lon','Lat','cell','strata'),1:n_iter))

#loop over iterations
for (iter in 1:n_iter) {
  
  iter<-1

  #print scenario to check progress
  cat(paste(" #############  ", iter, " #############\n"))
  
  #to store points
  dfpoints<-data.frame(matrix(NA,nrow=0,ncol=4))
  colnames(dfpoints)<-c('Lon','Lat','cell','strata')
  
  #random sample for each strata wit distance constrains
  for(istrata in result_list$solution$aggr_strata$STRATO) {
    
    #istrata<-2
    
    #subset cells for strata
    df<-subset(x7,X1==istrata)
    df1<-df
    
    #ratio of available cells and samples to take
    #ratio<-dim(df)[1]/allocations[istrata]
    i<-result_list$sample_allocations[istrata]
      
    dropcell<-c()
    selcell<-c()
    
    for (ii in rep(1,times=i)) {

      cell_i<-sample(df1$ID,ii)
      row<-df1[which(df1$ID==cell_i),c('row','col')][1,1]
      col<-df1[which(df1$ID==cell_i),c('row','col')][1,2]
      adj_cells<-expand.grid(row=c(row,row+(1:200),row-(1:200)),
                            col=c(col,col+(1:200),col-(1:200)))
      adj_cells1<-merge(adj_cells,df,by=c('row', 'col'))[,'ID']
      if (istrata == 2 ) {
        x<-merge(adj_cells,df,by=c('row', 'col'))
      }
      
      dropcell_i<-c(cell_i,adj_cells1)
      selcell<-c(selcell,cell_i)
      dropcell<-c(dropcell,dropcell_i)
      
      df1<-subset(df1, !(ID %in% dropcell))
      

    }
  
    if (length(selcell)!=i) {
      #print scenario to check progress
      cat(paste(" #############  not enough samples  #############\n"))}
    
    points<-subset(df,ID %in% selcell)[,c('x','y','ID','X1')]
    names(points)<-c('Lon','Lat','cell','strata')
  
    #append data
    dfpoints<-rbind(dfpoints,points)
    
  }
  
  #append points
  all_points[,,as.character(iter)]<-unlist(dfpoints)
}


loc_sur<-data.frame(unlist(all_points[,,as.character(1)]))
#change projection of spatial object 
coordinates(loc_sur)<- ~ Lon + Lat
#reproject shapefile
proj4string(loc_sur) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
#                            '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'
loc_sur<-spTransform(loc_sur,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))

#to dataframe
loc_sur<-as.data.frame(loc_sur)

coordinates(x)<- ~ x + y
#reproject shapefile
proj4string(x) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
#                            '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'
x<-spTransform(x,CRSobj = CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'))

#to dataframe
x<-as.data.frame(x)


# #strata by cell  
# strata<-solution$indices
# colnames(strata)<-c('cell','Strata')
# 
# #add a strata value to each cell
# D8<-merge(D6,strata,by='cell',all.x=TRUE)
# D8<-D8[,c("cell","Lat","Lon","Strata")]
# D8$Strata<-as.numeric(D8$Strata)
# D8$Strata<-ifelse(is.na(D8$Strata),999,D8$Strata)
# 
# #df to spatialpoint df
# coordinates(D8) <- ~ Lon + Lat
# crs(D8)<-c(crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
# 
# #reproject coordinates for plotting purposes
# D8_1<-spTransform(D8,'+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
# #D8_2<-data.frame(D8_1)
# 
# #x and y cells
 xycells<-as.integer(sqrt(nrow(x6)))

# create a template raster
r1 <- raster(ext=extent(x6),ncol=xycells, nrow=xycells) #c(15800,15800) 7000

#create raster
r2<-rasterize(x6, r1 ,field='X1')
crs(r2) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
r2[r2==999] <- NA
plot(r2)

gr_sh<-rgdal::readOGR('./shapefiles/regridadjacentcells/bs_grid.shp')

ggplot()+
  geom_raster(data=as.data.frame(r2, xy = TRUE),aes(x=x,y=y,fill=layer))+
  #geom_point(data=D8_2, aes(Lon, Lat, fill=Strata, group=NULL),size=2, stroke=0,shape=21)+
  scale_fill_gradientn(colours=c("#ea5545", "#f46a9b", "#ef9b20", "#edbf33", "#ede15b", "#bdcf32", "#87bc45", "#27aeef", "#b33dc6"),
                       #guide = guide_legend(),na.value = 'white',breaks=sort(unique(D8_1$Strata)),
                       #labels=paste0(sort(unique(D8_1$Strata))," (n=",allocations,')'))+ #,,
                       guide = guide_legend(),na.value = 'white',breaks=sort(na.omit(unique(values(r2)))),
                       labels=paste0(sort(na.omit(unique(values(r2))))," (n=",allocations,')'))+ #,,
  #geom_point(data=st_EBS,aes(x=longitude,y=latitude),shape=4,size=1)+
  #geom_point(data=st_corners1,aes(x=longitude,y=latitude),color='red',shape=20,size=1)+
  geom_polygon(data=ak_sppoly,aes(x=long,y=lat,group=group),fill = 'grey60')+
  geom_polygon(data=gr_sh,aes(x=long,y=lat,group=group),fill=NA,color='grey40',linetype='dashed')+
  scale_x_continuous(expand = c(0,0),position = 'bottom',
                     breaks = c(-175,-170,-165,-160,-155),sec.axis = dup_axis())+
  geom_polygon(data=NBS_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSshelf_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  geom_polygon(data=EBSslope_sh,aes(x=long,y=lat,group=group),fill=NA,col='black')+
  #geom_point(data=data.frame(x),aes(x=x,y=y),color='green',size=1.2,alpha=0.1)+
  geom_point(data=data.frame(loc_sur),aes(x=Lon,y=Lat),color='black',size=1.2)+
  
  
  # scale_shape_manual(values = c('optimization'=21,
  #                               'current design'=4,
  #                               'corner crab'=8),
  #                    breaks=unique(loc_sur$Stations),
  #                    labels=paste0(unique(loc_sur$Stations)," (n=",c(nrow(points1),nrow(st_EBS),nrow(st_corners1)),')'))+
  coord_sf(crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
           xlim = panel_extent$x,
           ylim = c(453099.5,2004909.7),
           label_axes = "-NE-")

