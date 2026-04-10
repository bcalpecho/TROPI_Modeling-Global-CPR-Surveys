# ---
# title: Helper functions for handling ensemble data
# author: Bryan Alpecho
# date: 2026 Feb 12 - 
# output: csv and plot for predictions
# ---

#These are additional functions to help in visualizing and summarizing the ensemble for predicting the zooplankton trophic groups (6b_predict_globalCPR)
#01 to plot ensemble at 2015 and 2100
#02 to plot for delta of chlos between 2015 and 2100
#03 to compute for summary statistics of chlos projections

#Load libraries
packages <- c("tidyverse",
              "stars",
              "tmap",
              "marginaleffects")

# Function to download the packages if necessary. Otherwise, these are loaded.
package.check <- lapply(
  packages,
  FUN = function(x)
  {
    if (!require(x, character.only = TRUE))
    {
      install.packages(x, dependencies = TRUE,
                       repos = "http://cran.us.r-project.org")
      library(x, character.only = TRUE)
    }
  }
)
  
#setup directory
  # # set date
  # date <- "25022026"
  # # define variable
  # var <- "chlos"
  # # access ESM output
  # base_dir <- "./hotrstuff/" # Just start relative to this location

##01 to plot ensemble at 2015 and 2100
  # #select ensemble
  # ensemble_median <- list.files(file.path(base_dir, "data", "proc", "ensemble", "median", var), full.names = TRUE)
  # 
  # ens_selected <- read_stars(ensemble_median, quiet = TRUE, proxy = TRUE) %>% setNames("chlos")
  # ens_name <- str_extract(basename(ensemble_median), ".*(?=\\.nc)")
  # #ens_snapshot_baseline <- ens_selected[,,,1] #2015:1 :: 2100:86
  # #ens_snapshot_future <- ens_selected[,,,86] 
  
  #to plot ensemble snapshot function
  plot_esm_snapshot <- function(ens){
    ensemble <- read_stars(ens)
      ens_name <- str_extract(basename(ens), ".*(?=\\_r1i1p1f1_RegriddedAnnual_20150101-21001231.nc)")
    
    #convert chlos (kg/m^3) to chla_sqrt (sqrt (mg/m^3) )
    chla_sqrt <- function(x){ sqrt(x * 1000000) }
    ens_converted <- st_apply(ensemble,1:3, chla_sqrt, rename = TRUE)
    
    ens_snapshot_baseline <- ens_converted[,,,1] #2015:1 :: 2100:86
    ens_snapshot_future <- ens_converted[,,,86]
    
    ens_list <- list(ens_snapshot_baseline, ens_snapshot_future)
    names(ens_list) <- c("baseline_2015","future_2100")
    for(i in 1:length(ens_list)){
      #ens_st <- st_as_stars(ens_list[i], dims = c("x","y"))
      #to plot
      TG_projection <- tm_shape(ens_list[[i]]) +
        tm_raster(
          col = "chla_sqrt",      # The name of your data column
          col.scale = tm_scale_continuous(values = c("white", "darkgreen")),
          col.legend = tm_legend(paste0("Chl os ensemble ",names(ens_list[i])),
                                 orientation = "landscape", 
                                 frame = FALSE,
                                 position = tm_pos_out("center", "bottom"))) +
        tm_layout(legend.outside = TRUE, bg.color = "black") 
      #to save plot
      tmap_save(TG_projection, filename=paste0("output/plots/",ens_name,"_",names(ens_list[i]),"_",date,".png"))
    }
  }

  #plot_esm_snapshot(ens_snapshot_baseline)
  
##02 to plot ensemble delta function 
  plot_esm_delta <- function(ens){
    #read in ensemble as stars object
    ens_selected <- read_stars(ens, quiet = TRUE, proxy = TRUE) %>% setNames("chlos")
    ens_name <- str_extract(basename(ens), ".*(?=\\_r1i1p1f1_RegriddedAnnual_20150101-21001231.nc)")
    print(paste0("SSP scenario in process: ",ens_name))
    
    str_extract(basename(ens), "_ssp\\d{3}_")
    
    #convert chlos (kg/m^3) to chla_sqrt (sqrt (mg/m^3) )
    chla_sqrt <- function(x){ sqrt(x * 1000000) }
    esm_converted <- st_apply(ens_selected,1:3, chla_sqrt, rename = TRUE)
    #separate baseline (2015) and future (2100) snapshot of chlos projections
    ens_snapshot_baseline <- esm_converted[,,,1] #2015:1 :: 2100:86
    ens_snapshot_future <- esm_converted[,,,86]
    
    #convert star to df
    ens_baseline_df <- as.data.frame(ens_snapshot_baseline, xy = TRUE, na.rm = T) %>% rename("chla_sqrt_baseline" = "chla_sqrt")
    #xy coordinates and kept 'NA'
    ens_future_df <- as.data.frame(ens_snapshot_future, xy = TRUE, na.rm = T) %>% rename("chla_sqrt_future" = "chla_sqrt")
    #xy coordinates and kept 'NA'

    #compute for delta
    ens_delta_summary <- ens_baseline_df %>% select("chla_sqrt_baseline", "x", "y") %>% 
      left_join(ens_future_df %>% select("chla_sqrt_future", "x", "y"), by = c("x","y")) %>% 
      mutate("chla_delta" = (((chla_sqrt_future - chla_sqrt_baseline)/chla_sqrt_baseline)*100))
      
    #save
    saveRDS(ens_delta_summary, file=paste0("output/projections/Delta_",ens_name,"_",date,".RData"))
    print(paste0("File (computed delta) saved: output/projections/Delta_",ens_name,"_",date,".RData"))

    #to plot
    #convert to stars
    st_projections <- st_as_stars(ens_delta_summary, dims = c("x","y"))
    
    #print(summary(is.na(st_projections$chla_sqrt_delta)))
    TG_projection <- tm_shape(st_projections) +
      tm_raster(
        col = "chla_delta",      # The name of your data column
        col.scale = tm_scale_continuous(midpoint = 0, values = c("blue", "white","red")),
        col.legend = tm_legend(paste("Change in surface chlorophyll under climate change (SSP585)",sep=""),
                               orientation = "landscape", 
                               frame = FALSE)) +
      tm_layout(legend.outside = TRUE, bg.color = "black", text.size = 36) 
    
    tmap_save(TG_projection, filename=paste0("output/plots/chlos_delta_",date,".png"))
    print(paste0("Plot saved: output/plots/chlos_delta_",date,".png"))
  }
  
  #plot_esm_delta(ens_selected)
  
##03 to compute for summary statistics of chlos
  
  summary_stats_chlos <- function(ensemble_list){
    for(i in 1:length(ensemble_list)){
      
      ens_selected <- read_stars(ensemble_list[i], quiet = TRUE, proxy = TRUE) %>% setNames("chlos")
      ens_name <- str_extract(basename(ensemble_list[i]), ".*(?=\\.nc)")
      
      #convert chlos (kg/m^3) to chla_sqrt (sqrt (mg/m^3) )
      chla_sqrt <- function(x){ sqrt(x * 1000000) }
      esm_converted <- st_apply(ens_selected,1:3, chla_sqrt, rename = TRUE)
      
      ens_snapshot_baseline <- esm_converted[,,,1] #2015:1 :: 2100:86
      ens_snapshot_future <- esm_converted[,,,86]
      
      #convert star to df
      ens_baseline_df <- as.data.frame(ens_snapshot_baseline, xy = TRUE, na.rm = T) %>% rename("chla_sqrt_baseline" = "chla_sqrt")
      #xy coordinates and kept 'NA'
      ens_future_df <- as.data.frame(ens_snapshot_future, xy = TRUE, na.rm = T) %>% rename("chla_sqrt_future" = "chla_sqrt")
      #xy coordinates and kept 'NA'
      
      options(digits = 10)
      ens_baseline_summary <- ens_baseline_df %>% 
        summarise(mean = mean(chla_sqrt_baseline, na.rm =T), 
                  median = median(chla_sqrt_baseline, na.rm =T), 
                  Q1 = quantile(chla_sqrt_baseline, 0.25, na.rm =T), 
                  Q3 = quantile(chla_sqrt_baseline, 0.75, na.rm =T), 
                  min = min(chla_sqrt_baseline, na.rm = T),
                  max = max(chla_sqrt_baseline, na.rm = T))
      
      print(paste("Ensemble: ",ens_name,sep=""))
      print("Baseline (yr 2015)")
      print(ens_baseline_summary) 
      
      ens_future_summary <- ens_future_df %>% 
        summarise(mean = mean(chla_sqrt_future, na.rm =T), 
                  median = median(chla_sqrt_future, na.rm =T), 
                  Q1 = quantile(chla_sqrt_future, 0.25, na.rm =T), 
                  Q3 = quantile(chla_sqrt_future, 0.75, na.rm =T), 
                  min = min(chla_sqrt_future, na.rm = T),
                  max = max(chla_sqrt_future, na.rm = T))
      
      print("Future (yr 2100)")
      print(ens_future_summary) 
      
      ens_delta_summary <- ens_baseline_df %>% select("chla_sqrt_baseline", "x", "y") %>% 
        left_join(ens_future_df %>% select("chla_sqrt_future", "x", "y"), by = c("x","y")) %>% 
        mutate("chla_delta" = (((chla_sqrt_future - chla_sqrt_baseline)/chla_sqrt_baseline))) %>% 
        summarise(mean = mean(chla_delta, na.rm =T), 
                  median = median(chla_delta, na.rm =T), 
                  Q1 = quantile(chla_delta, 0.25, na.rm =T), 
                  Q3 = quantile(chla_delta, 0.75, na.rm =T), 
                  min = min(chla_delta, na.rm = T),
                  max = max(chla_delta, na.rm = T))
        
      print("Change by 2100 relative to 2015 (reported in proportional change)")
      print(ens_delta_summary)
    }
  }
  
  #summary_stats_chlos(ensemble_median)
  
##