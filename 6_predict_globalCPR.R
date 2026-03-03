# ---
# title: Project Global CPR
# author: Bryan Alpecho
# date: 2026 Feb 10  - 
# output: csv for predictions
# ---

#Aims
#01 to predict zooplantkon trophic group
#02 to compute for delta of trophic groups between 2015 and 2100
#03 to provide summary statistics for projected trophic groups' relative abundance

#basic aim: compute for change in zooplankton trophic group by 2100 relative to 2015

##setup
  date <- "24022026"

  #load libraries
    library(marginaleffects)
    library(stars)
    library(tidyverse)
    
  #load ensemble
    var <- "chlos"
    base_dir <- "./hotrstuff/" 
    #ensembles_mean <- list.files(file.path(base_dir, "data", "proc", "ensemble", "mean", var), full.names = TRUE)
    ensembles_median <- list.files(file.path(base_dir, "data", "proc", "ensemble", "median", var), full.names = TRUE)
    
    #process the SSP scenarios one-by-one (One ensemble each SSP scenario)
    
  #load zooplankton GLMs
    load("Output/previousModels/revision/final_zoop_mdls_21022026.RData") #need to revise the directory
    mdl_list <- list(Carni_mdl_zib, Omni_mdl_zib, Filter_mdl_zib)
    names(mdl_list) <- c("Carni","Omni","Filter")

#01 to predict zooplantkon trophic group

  predict_zoop_delta <- function(ensemble, mdls){
    
    #load ensemble
    ens_selected <- read_stars(ensemble, quiet = TRUE, proxy = TRUE) %>% setNames("chlos")
    ens_snapshot_baseline <- ens_selected[,,,1] #2015:1 :: 2100:86
    ens_snapshot_future <- ens_selected[,,,86] 
    #identify ssp scenario
    ssp_scenario <- str_extract(basename(ensemble), "ssp\\d{3}")
    
    #convert chlos (kg/m^3) to chla_sqrt (sqrt (mg/m^3) )
    esm_list <- list(ens_snapshot_baseline, ens_snapshot_future)
    names(esm_list) <- c("baseline","future")
    
    for(i in 1:length(mdls)){
      #identify trophic group
      TG <- names(mdls)[i]
      
      for(j in 1:length(esm_list)){
        esm_level <- names(esm_list)[j]
        chla_sqrt <- function(x){ sqrt(x * 1000000) } #conversion factor from kg m-3 to 
        esm_converted <- st_apply(esm_list[[j]],1:2, chla_sqrt, rename = TRUE)
        
        #convert star to df
        esm_df <- as.data.frame(esm_converted, xy = TRUE, na.rm = FALSE) #xy coordinates and kept 'NA'
        #simplify esm_df
        esm_df <- esm_df %>% select(c(chla_sqrt,x,y))
        ##
        
        esm_pred <- predictions(mdls[[i]],
                                newdata = datagrid(chla_sqrt = esm_df$chla_sqrt),
                                re.form = NA)
        
        esm_pred_merged <- esm_pred %>%
          left_join(esm_df, by = "chla_sqrt")
        
        saveRDS(esm_pred_merged, file=paste("Output/data/projections/",TG,"_",ssp_scenario,"_",esm_level,"_",date,".RData",sep=""))
        print(paste("Output/data/projections/",TG,"_",ssp_scenario,"_",esm_level,"_",date,".RData",sep=""))
      }
    }
  }
  
  #ssp126
  predict_zoop_delta(ensembles_median[1], mdl_list)
  #ssp245
  predict_zoop_delta(ensembles_median[2], mdl_list)
  #ssp370
  predict_zoop_delta(ensembles_median[3], mdl_list)
  #ssp585
  predict_zoop_delta(ensembles_median[4], mdl_list)

#prev versions 14022026

#02 to compute for delta of trophic groups between 2015 and 2100
 
  compute_zoop_delta <- function(mdls){
    ssp_list <- list("ssp126","ssp245","ssp370","ssp585")
    
    for(i in 1:length(ssp_list)){
      
      print(paste("SSP Scenario: ",ssp_list[i],sep=""))
      for(j in 1:length(mdls)){
        
        TG <- names(mdls)[j]
        
        if(file.exists(paste("Output/data/projections/",TG,"_",ssp_list[i],"_baseline_",date,".RData",sep=""))) {
          projection_baseline <- readRDS(file=paste("Output/data/projections/",TG,"_",ssp_list[i],"_baseline_",date,".RData",sep=""))
          
          projection_future <- readRDS(file=paste("Output/data/projections/",TG,"_",ssp_list[i],"_future_",date,".RData",sep=""))
          
          
          projection_delta <- projection_baseline %>% 
            left_join(projection_future %>% 
                        rename("future_estimate" = "estimate") %>% 
                        select(c("future_estimate","x","y")), by = c("x","y")) %>% 
            mutate(delta = ((future_estimate - estimate)/estimate) * 100)
          
          delta_summary <- projection_delta %>% 
            summarise(mean_delta = mean(delta, na.rm =T), 
                      median_delta = median(delta, na.rm =T), 
                      Q1_delta = quantile(delta, 0.25, na.rm =T), 
                      Q3_delta = quantile(delta, 0.75, na.rm =T), 
                      min_delta = min(delta, na.rm = T),
                      max_delta = max(delta, na.rm = T))
          
          print(TG)
          print(delta_summary)
          
        } else {
          print(paste("File does not exist: Output/data/projections/",TG,"_",ssp_list[i],"_baseline_",date,".RData",sep=""))
        }
      }
    }
  }
  
  compute_zoop_delta(mdl_list)

#03 to provide summary statistics for projected trophic groups' relative abundance
  
  summary_stats_TG <- function(mdls){
    
    ssp_list <- list("ssp126","ssp245","ssp370","ssp585")
    for(i in 1:length(mdls)){
      TG <- names(mdls)[i]
      
      for(j in 1:length(ssp_list)){
        
        projection_future <- readRDS(file=paste("Output/data/projections/",TG,"_",ssp_list[j],"_future_",date,".RData",sep=""))
        
        TG_summary <- projection_future %>% 
          summarise(mean = mean(estimate, na.rm =T), 
                    median = median(estimate, na.rm =T), 
                    Q1 = quantile(estimate, 0.25, na.rm =T), 
                    Q3 = quantile(estimate, 0.75, na.rm =T), 
                    min = min(estimate, na.rm = T),
                    max = max(estimate, na.rm = T))
        
        print(paste("Trophic Group: ",TG,"; SSP Scenario: ", ssp_list[j], sep=""))
        print(TG_summary) 
        
      }
    }
    
  }
  
  summary_stats_TG(mdl_list)
  
  baseline_stats_TG <- function(mdls){
    
    ssp_list <- list("ssp126","ssp245","ssp370","ssp585")
    
    print("Baseline (2015) summary statistics")
    for(i in 1:length(mdls)){
      TG <- names(mdls)[i]
      
      for(j in 1:length(ssp_list)){
        
        projection_baseline <- readRDS(file=paste("Output/data/projections/",TG,"_",ssp_list[j],"_baseline_",date,".RData",sep=""))
        
        TG_summary <- projection_baseline %>% 
          summarise(mean = mean(estimate, na.rm =T), 
                    median = median(estimate, na.rm =T), 
                    Q1 = quantile(estimate, 0.25, na.rm =T), 
                    Q3 = quantile(estimate, 0.75, na.rm =T), 
                    min = min(estimate, na.rm = T),
                    max = max(estimate, na.rm = T))
        
        print(paste("Trophic Group: ",TG,"; SSP Scenario: ", ssp_list[j], sep=""))
        print(TG_summary) 
        
      }
    }
    
  }
  
  baseline_stats_TG(mdl_list)
  
##