# ---
# title: Wrangling ensemble
# author: Bryan Alpecho
# date: 2026 Feb 12 - 
# output: csv and plot for predictions
# ---

#Aims
#01 to create an ensemble median of surface chlorophyll projections
#02 to plot ensemble at 2015 and 2100
#03 to plot for delta of chlos between 2015 and 2100
#04 to compute for summary statistics of chlos projections

#Load libraries
  library(tidyverse)
  library(tmap)
  library(stars)
  
#setup directory
  # set date
  date <- "16022026"
  # define variable
  var <- "chlos"
  # access ESM output
  base_dir <- "./hotrstuff/" # Just start relative to this location

##01 to plot ensemble at 2015 and 2100
  #select ensemble
  ensemble_median <- list.files(file.path(base_dir, "data", "proc", "ensemble", "median", var), full.names = TRUE)
  
  ens_selected <- read_stars(ensemble_median, quiet = TRUE, proxy = TRUE) %>% setNames("chlos")
  ens_name <- str_extract(basename(ensemble_median), ".*(?=\\.nc)")
  #ens_snapshot_baseline <- ens_selected[,,,1] #2015:1 :: 2100:86
  #ens_snapshot_future <- ens_selected[,,,86] 
  
  #to plot ensemble snapshot function
  plot_esm_snapshot <- function(esm_snapshot){
    
    chla_sqrt <- function(x){ sqrt(x * 1000000) }
    esm_converted <- st_apply(esm_model,1:2, chla_sqrt, rename = TRUE)
    
    #convert star to df
    esm_df <- as.data.frame(esm_converted, xy = TRUE, na.rm = FALSE) #xy coordinates and kept 'NA'
    #simplify esm_df
    esm_df <- esm_df %>% select(c(chla_sqrt,x,y))
    ##

    esm_st <- st_as_stars(esm_df, dims = c("x","y"))

    TG_projection <- tm_shape(esm_st) +
      tm_raster(
        col = "chla_sqrt",      # The name of your data column
        col.scale = tm_scale_continuous(values = c("white", "darkgreen")),
        col.legend = tm_legend("Chlos ensemble ssp585 2015",
                               orientation = "landscape", 
                               frame = FALSE,
                               position = tm_pos_out("center", "bottom"))) +
      tm_layout(legend.outside = TRUE, bg.color = "black") 
    
    tmap_save(TG_projection, filename=paste("Output/map/projections/",ens_name,"_",date,".png",sep=""))
    
  }

  plot_esm_snapshot(ens_snapshot_baseline)
  
##02 to plot ensemble delta function (non-working function)
  plot_esm_delta <- function(esm){
    
    #convert chlos (kg/m^3) to chla_sqrt (sqrt (mg/m^3) )
    chla_sqrt <- function(x){ sqrt(x * 1000000) }
    esm_converted <- st_apply(esm,1:2, chla_sqrt, rename = TRUE)

    ens_snapshot_baseline <- esm_converted[,,,1] #2015:1 :: 2100:86
    ens_snapshot_future <- esm_converted[,,,86]
    
    #convert star to df
    ens_baseline_df <- as.data.frame(ens_snapshot_baseline, xy = TRUE, na.rm = T) %>% rename("chla_sqrt_baseline" = "chla_sqrt")
    #xy coordinates and kept 'NA'
    ens_future_df <- as.data.frame(ens_snapshot_future, xy = TRUE, na.rm = T) %>% rename("chla_sqrt_future" = "chla_sqrt")
    #xy coordinates and kept 'NA'

    #simplify esm_df
    ens_df <- ens_baseline_df %>% 
      left_join(ens_future_df, by = c("x","y")) %>% 
      mutate(chla_sqrt_delta = ((chla_sqrt_future - chla_sqrt_baseline)/chla_sqrt_baseline) * 100)
      
    #save
    saveRDS(ens_df, file=paste("Output/data/projections/chlos_delta_",date,".RData",sep=""))
    print(paste("Output/data/projections/chlos_delta_",date,".RData",sep=""))

    #to plot
    #convert to stars
    st_projections <- st_as_stars(ens_df, dims = c("x","y"))
  
    TG_projection <- tm_shape(st_projections) +
      tm_raster(
        col = "chla_sqrt_delta",      # The name of your data column
        col.scale = tm_scale_continuous(values = c("blue", "white","red")),
        col.legend = tm_legend(paste("Change in surface chlorophyll under climate change (SSP585)",sep=""),
                               orientation = "landscape", 
                               frame = FALSE)) +
      tm_layout(legend.outside = TRUE, bg.color = "black", text.size = 36) 
    
    tmap_save(TG_projection, filename=paste("Output/map/projections/chlos_delta_",date,".png",sep=""))
    
  }
  
  plot_esm_delta(ens_selected)
  
##03 to compute for summary statistics of chlos
  
  summary_stats_chlos <- function(ensemble_list){
    for(i in 1:length(ensemble_list)){
      
      ens_selected <- read_stars(ensemble_list[i], quiet = TRUE, proxy = TRUE) %>% setNames("chlos")
      ens_name <- str_extract(basename(ensemble_list[i]), ".*(?=\\.nc)")
      
      #convert chlos (kg/m^3) to chla_sqrt (sqrt (mg/m^3) )
      chla_sqrt <- function(x){ round(sqrt(x * 1000000), 5) }
      esm_converted <- st_apply(ens_selected,1:2, chla_sqrt, rename = TRUE)
      
      ens_snapshot_baseline <- esm_converted[,,,1] #2015:1 :: 2100:86
      ens_snapshot_future <- esm_converted[,,,86]
      
      #convert star to df
      ens_baseline_df <- as.data.frame(ens_snapshot_baseline, xy = TRUE, na.rm = T) %>% rename("chla_sqrt_baseline" = "chla_sqrt")
      #xy coordinates and kept 'NA'
      ens_future_df <- as.data.frame(ens_snapshot_future, xy = TRUE, na.rm = T) %>% rename("chla_sqrt_future" = "chla_sqrt")
      #xy coordinates and kept 'NA'
      
      options(digits = 10)
      ens_baseline_summary <- ens_baseline_df %>% 
        summarise(mean_delta = mean(chla_sqrt_baseline, na.rm =T), 
                  median_delta = median(chla_sqrt_baseline, na.rm =T), 
                  Q1_delta = quantile(chla_sqrt_baseline, 0.25, na.rm =T), 
                  Q3_delta = quantile(chla_sqrt_baseline, 0.75, na.rm =T), 
                  min_delta = min(chla_sqrt_baseline, na.rm = T),
                  max_delta = max(chla_sqrt_baseline, na.rm = T))
      
      print(paste("Ensemble: ",ens_name,sep=""))
      print("Baseline (yr 2015)")
      print(ens_baseline_summary) 
      
      ens_future_summary <- ens_future_df %>% 
        summarise(mean_delta = mean(chla_sqrt_future, na.rm =T), 
                  median_delta = median(chla_sqrt_future, na.rm =T), 
                  Q1_delta = quantile(chla_sqrt_future, 0.25, na.rm =T), 
                  Q3_delta = quantile(chla_sqrt_future, 0.75, na.rm =T), 
                  min_delta = min(chla_sqrt_future, na.rm = T),
                  max_delta = max(chla_sqrt_future, na.rm = T))
      
      print("Future (yr 2100)")
      print(ens_future_summary) 
      
    }
  }
  
  summary_stats_chlos(ensemble_median)
  
##