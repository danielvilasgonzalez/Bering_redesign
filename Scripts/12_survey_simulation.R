####################################################################
####################################################################
##    
##    fit single sp VAST model for the EBS shelf and NBS using temp (SBT, BotTemp) with cubic effect 
##    Daniel Vilas (danielvilasgonzalez@gmail.com/dvilasg@uw.edu)
##    evaluate parameters : https://github.com/James-Thorson-NOAA/VAST/blob/main/R/make_parameters.R
##    https://rdrr.io/github/James-Thorson/VAST/man/check_fit.html
##
####################################################################
####################################################################

#clear all objects
rm(list = ls(all.names = TRUE)) 
#free up memrory and report the memory usage
gc() 

#libraries from cran to call or install/load
pack_cran<-c("splines",'dplyr')

#install pacman to use p_load function - call library and if not installed, then install
if (!('pacman' %in% installed.packages())) {
  install.packages("pacman")}

#install coldpool to extract SBT for the EBS
if (!('VAST' %in% installed.packages())) {
  devtools::install_github("james-thorson/VAST@main", INSTALL_opts="--no-staged-install")};library(VAST)

#load/install packages
pacman::p_load(pack_cran,character.only = TRUE)

#setwd
out_dir<-'C:/Users/Daniel.Vilas/Work/Adapting Monitoring to a Changing Seascape/'
setwd(out_dir)

#version VAST (cpp)
version<-'VAST_v13_1_0'

#number of knots
knots<-'500' #1000 

#years
yrs<-1982:2022

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

load('./extrapolation grids/northern_bering_sea_grid.rda')
load('./extrapolation grids/eastern_bering_sea_grid.rda')
grid<-as.data.frame(rbind(data.frame(northern_bering_sea_grid,region='NBS'),data.frame(eastern_bering_sea_grid,region='EBS')))
grid$cell<-1:nrow(grid)
  
  
######################
# GET ANNUAL PREDICTIONS (2023-2027) SBT-
######################

for (sp in spp) {
  
  sp<-"Gadus macrocephalus"
  
  lf<-list.files(path = paste0('./output/species/',sp,'/'),pattern = 'samples_optimization_')
  lf1<-gsub('samples_optimization_','',lf)
  lf2<-gsub('.RData','',lf1)
  load( file = paste0("./output/species/",sp,'/fit_projection.RData')) #pr_list
  
  for (samp in lf2) {
    
    samp_proj<-list()
    load(file=paste0('./output/species/',sp,'/samples_optimization_',samp,'.RData')) #all_points
    for (sbt in names(pr_list)) {
      
      #print scenario to check progress
      cat(paste(" #############     PROJECTING    #############\n",
                " #############   Species", sp, match(sp,spp), 'out of',length(spp),  "  #############\n",
                " #############  Sampling ", samp, " #############\n",
                " #############  SBT ", sbt, " #############\n"))
   
    #sbt<-names(pr_list)[1]  
    #samp<-lf2[1]
    

    bio_proj<-stack(paste0('./output/species/',sp,'/biomass_projection_',sbt,'.grd'))
    names(bio_proj)<-paste0('y',c(2023:2027))
    #all_points[which(all_points[,'strata',]==1),,1:2]
    
    points2<-data.frame(matrix(nrow = 0,ncol = 10))
    colnames(points2)<-c(dimnames(all_points)[[2]],names(bio_proj),'iter')
    
    for (iter in dimnames(all_points)[[3]]) {
      #iter<-'1'
      
      points<-data.frame(all_points[,,iter])

      coordinates(points)<- ~ Lon + Lat
      
      points1<-data.frame(data.frame(all_points[,,iter]),raster::extract(bio_proj,points))
      points1$iter<-iter
      points2<-rbind(points2,points1)
      
      
    }
    
    samp_proj[[paste0('SBT',sbt)]]<-points2
    
    }
    
    save(samp_proj, file = paste0("./output/species/",sp,'/projections_survey_',samp,'.RData'))
    
  }
  
}


#####################
# GET km2 for each cell
#####################

for (sp in spp) {
  
  sp<-"Gadus macrocephalus"
  
  #get samples scenarios
  lf<-list.files(path = paste0('./output/species/',sp,'/'),pattern = 'samples_optimization_')
  lf1<-gsub('samples_optimization_','',lf)
  lf2<-gsub('.RData','',lf1)

  
  for (samp in lf2) {
    
    samp<-lf2[1]
    
    #load survey samples
    load( file = paste0("./output/species/",sp,'/projections_survey_',samp,'.RData')) #samp_proj
    #load optimization results
    load(paste0('./output/species/',sp,'/optimization_results_',samp,'.RData')) #result_list
    
    result_list$solution$
      
    for (sbt in names(samp_proj)) {
      
    x<-samp_proj[[sbt]]
    
    #merge biomass samples of survey with grid
    y<-merge(x,grid,by='cell')
    
    #strata data
    survey_detail <- data.frame("Stratum" = result_list$solution$aggr_strata$STRATO,
                                "nh" = result_list$sample_allocations)
    
    #weight of strata for each
    survey_detail$wh<-survey_detail$nh/sum(survey_detail$nh)
    
    # Calculate stratum area
    strata_areas <- aggregate(cell_areas ~ solution, 
                              FUN = sum,
                              data = with(input, data.frame(solution, 
                                                            cell_areas)))
    strata_areas <- subset(strata_areas, solution %in% survey_detail$Stratum)
    
    coordinates(x)<- ~ Lon + Lat
    proj4string(x) <- CRS('+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs') 
    x<-spTransform(x,CRSobj = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    y<-data.frame(x)
    #merge(x,grid,by='cell')
    
    
    }
  }
}
    
    
grid.ebs_year

mdl <- make_map_info(Region = fit$settings$Region,
                     spatial_list = fit$spatial_list,
                     Extrapolation_List = fit$extrapolation_list)
gridlatlon<-as.data.frame(mdl$PlotDF)
#names(grid)
#names(gridlatlon)
grid1<-merge(gridlatlon,grid,by=c('Lat','Lon'))
grid1<-grid1[order(grid1$x2i),]
dim(grid1)

#get predictions
D_gt<-fit$Report$D_gct[,1,]
#D_gt_proj<-D_gt[,paste0(project_yrs)]
D_gt<-drop_units(D_gt)
D_gt<-data.frame('cell'=c(1:fit$spatial_list$n_g),D_gt)

colnames(D_gt)<-c('cell',c(fit$year_labels))
D_gt1<-reshape2::melt(D_gt,id=c('cell'))

D <- merge(D_gt1, mdl$PlotDF, by.x='cell', by.y='x2i')


index<-fit$Report$Index_ctl

temp_true<-t(apply(X = index,
                   MARGIN = 2,FUN = sum))

sdata<-simulate_data(fit = fit)
sim_data <- array(data = NA, dim = c(nrow(grid1),
                                     length(unique(fit$year_labels)),
                                     1))

sim_data[, , 1] <- matrix(sdata$D_gct*0.001)


temp_density <-sweep(x = sim_data[, , 1] / 0.001,
      MARGIN = 1,
      STATS = grid1$Area_in_survey_km2,
      FUN = "/")


load('./output/species/Gadus macrocephalus/optimization_static_data.RData')

load('./output/species/Gadus macrocephalus/optimization_summary_stats.RData')

input = list(
  "density" = temp_density,
  "cell_areas" = grid1$Area_in_survey_km2,
  "solution" = temp_solution, # "current" = grid_goa$stratum_new_label[cell_idx], "opt" = result_list$sol_by_cell[cell_idx])
  
  result_list <- list(solution = solution,
                      # sum_stats = sum_stats,
                      # cvs = cv_by_boat,
                      # sample_allocations = sample_allocations,
                      # sol_by_cell = temp_ids = solution$indices$X1)

  "allocation" = temp_allocation,
  "true_index_district" = temp_true_index,
  "post_strata" = district_vals[cell_idx]
)


do_STRS <- function(input){
  
  #Some constants
  n_cells <- dim(input$density)[1]
  n_time <-  dim(input$density)[2]
  n_dom <- length(table(input$post_strata))
  
  survey_detail <- 
    data.frame("Stratum" = as.integer(names(table(input$solution))),
               "Nh" = as.integer(table(input$solution)),
               "nh" = input$allocation)
  
  # Remove strata with zero stations
  strata_to_use <- survey_detail$nh > 0
  survey_detail <- survey_detail[strata_to_use, ]
  
  # Assume stratum weights include untrawlabe areas
  survey_detail$Wh <- survey_detail$Nh / sum(survey_detail$Nh)
  survey_detail$wh <- with(survey_detail, nh/Nh)
  
  # Calculate stratum area
  strata_areas <- aggregate(cell_areas ~ solution, 
                            FUN = sum,
                            data = with(input, data.frame(solution, 
                                                          cell_areas)))
  strata_areas <- subset(strata_areas, solution %in% survey_detail$Stratum)
  
  # Create result objects
  mean_density <- cv <- index <- rel_bias <- rel_log_bias <- c()
  index_district <- array(dim = c(n_time, n_dom))
  
  for (iyear in 1:n_time) { ## Loop through years -- start
    
    # Subset density df by year
    sub_df <- input$density[, iyear]
    
    # Take a random sample based on the allocation and stratum
    sample_vec <- c()
    for(istrata in 1:nrow(survey_detail)) {
      sample_vec <- c(sample_vec,
                      sample(x = which(input$solution == survey_detail$Stratum[istrata]),
                             size = survey_detail$nh[istrata]) )
    }
    
    sampled_strata <- rep(x = survey_detail$Stratum, 
                          times = survey_detail$nh)
    
    # Subset sub_df by which cells were chosen
    sample_df <- as.vector(sub_df[sample_vec])
    
    # Calculate STRS mean density
    strata_mean <- tapply(X = sample_df, 
                          INDEX = sampled_strata,
                          FUN = mean)
    STRS_mean <- sum(strata_mean * survey_detail$Wh)
    
    # Calculate STRS variance of mean density
    strata_var <- tapply(X = sample_df, 
                         INDEX = sampled_strata,
                         FUN = var)
    STRS_var <- sum(strata_var * with(survey_detail, Wh^2 * (1 - wh) / nh) )
    
    # Save mean and cv of estimates across species
    cv[iyear] <- sqrt(STRS_var) / STRS_mean
    
    # Calculate index of abundance by district
    index_df <- data.frame(Area_km2 = input$cell_areas,
                           stratum = input$solution,
                           district = input$post_strata,
                           mean_dens = strata_mean[paste(input$solution)])
    
    index_district[iyear, ] <- 
      tapply(X = index_df$mean_dens * index_df$Area_km2,
             INDEX = index_df$district,
             FUN = sum,
             na.rm = TRUE) * 0.001
    
    # Calculate total index
    index[iyear] <- sum(strata_areas * strata_mean) * 0.001
    
  } ## Loop through years -- end
  
  # Calculate Relative bias of index over the entire domain and by districts
  rel_bias <- 100 * (index - rowSums(input$true_index_district)) /
    rowSums(input$true_index_district)
  rel_log_bias <- 
    log10(rowSums(index_district) / rowSums(input$true_index_district))
  
  bias_index_district <- 100 * (index_district - input$true_index_district) /
    input$true_index_district
  log_bias_index_district <- log10(index_district / 
                                     input$true_index_district)
  
  return(list("strs_mean" = STRS_mean,
              "strs_index" = index,
              "cv" = round(cv, 4),
              "rel_bias" = round(rel_bias, 2),
              "rel_log_bias" = round(rel_log_bias, 3),
              "bias_index_district" = round(bias_index_district, 2),
              "log_bias_index_district" = round(log_bias_index_district, 3) )
  )
}



sim_survey <- do_STRS(input = list(
  "density" = temp_density,
  "cell_areas" = grid_goa$Area_km2[cell_idx],
  "solution" = temp_solution,
  "allocation" = temp_allocation,
  "true_index_district" = temp_true_index_district,
  "post_strata" = district_vals[cell_idx]
))



fit$Report$Index_gctl[1,,'1982',1]
fit$Report$Index_gctl[,,'1982',3]
fit$Report$Index_gctl[,,'2020',]

index<-fit$Report$Index_gctl[, , ,1]

#save projection list
load(paste0("./output/species/",sp,'/fit_projection.RData')) #pr_list



sdata$Index_gctl


dim(sdata$Index_gctl)

sim_data <- array(data = NA, dim = c(nrow(grid_goa),
                                     length(unique(data$YEAR)),
                                     1000))

for (isim in 1:1000) {
  Sim1 <- FishStatsUtils::simulate_data(fit = fit_sim, 
                                        type = 1, 
                                        random_seed = isim)
  sim_data[, , isim] <- matrix(data = Sim1$b_i[pred_TF == 1] * 0.001, 
                               nrow = nrow(grid_goa), 
                               ncol = length(unique(data$YEAR)))
  if(isim%%100 == 0) print(paste("Done with", ispp, "Iteration", isim))
}

save(sim_data, file = paste0(result_dir, "/simulated_data.RData"))
