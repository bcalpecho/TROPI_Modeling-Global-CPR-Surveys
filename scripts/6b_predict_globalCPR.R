# ---
# title: Predict zooplankton trophic groups on ensemble of surface chlorophyll-a (chlos) projections
# author: Bryan Alpecho
# date: 2026 Feb 10  - 
# output: csv for predictions
# ---

#Aims
#01 to visualize and summarize the ensemble of chlos projections
#02 to predict proportions of zooplantkon trophic groups on the ensemble of chlos projections
#03 to compute for delta of proportions of trophic groups between 2015 and 2100
#04 to provide summary statistics for projected trophic groups' relative abundance

#basic aim: compute for change in zooplankton trophic group by 2100 relative to 2015

##setup
  source("helper_ensemble.R")
  date <- "20032026"

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

#01 to visualize and summarize the ensemble of chlos projections
    
    #to plot delta (2100 relative to 2015 chlos)
    #ssp126
    ens_ssp126 <- read_stars(ensembles_median[1])
    plot_esm_delta(ens_ssp126)
      
    #ssp585
    ens_ssp585 <- read_stars(ensembles_median[4])
    plot_esm_delta(ens_ssp585)
    
    #summary statistics
    summary_stats_chlos(ensembles_median)
    
#02 to predict zooplankton trophic group

  predict_zoop_delta <- function(ensemble, mdls){
    
    #1.1 load ensemble
    ens_selected <- read_stars(ensemble, quiet = TRUE, proxy = TRUE) %>% setNames("chlos")
    ens_snapshot_baseline <- ens_selected[,,,1] #2015:1 :: 2100:86
    ens_snapshot_future <- ens_selected[,,,86] 
    #1.2 identify ssp scenario
      ssp_scenario <- str_extract(basename(ensemble), "ssp\\d{3}")
    #1.3 list baseline and future snapshot under ssp scenario
      esm_list <- list(ens_snapshot_baseline, ens_snapshot_future)
      names(esm_list) <- c("baseline","future")
    
    #1.4 perform iteration of predicting outcomes per zooplankton trophic group    
    for(i in 1:length(mdls)){
      #identify trophic group
      TG <- names(mdls)[i]
      #1.5 perform iteration of predicting outcomes for baseline and future snapshot 
      for(j in 1:length(esm_list)){
        esm_level <- names(esm_list)[j]
        #1.6 convert chlos (kg/m^3) to chla_sqrt (sqrt (mg/m^3) )
          chla_sqrt <- function(x){ sqrt(x * 1000000) } #conversion factor from kg m-3 to 
          esm_converted <- st_apply(esm_list[[j]],1:2, chla_sqrt, rename = TRUE)
        
        #convert star to df
        esm_df <- as.data.frame(esm_converted, xy = TRUE, na.rm = FALSE) #xy coordinates and kept 'NA'
        #simplify esm_df
        esm_df <- esm_df %>% select(c(chla_sqrt,x,y))
        ##
        
        #1.7 predict outcomes of models given the Chl-a projections
        esm_pred <- predictions(mdls[[i]], newdata = datagrid(chla_sqrt = esm_df$chla_sqrt), re.form = NA)

        esm_pred_merged <- esm_df %>%
           left_join(esm_pred %>% select("chla_sqrt","estimate"), by = c("chla_sqrt"))
        
        #1.8 Save predictions in an R data output 
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

#03 to compute for delta of trophic groups between 2015 and 2100

  compute_zoop_delta <- function(mdls){
    #stable projections - 24 & 25022026
    date_projections <- "24022026"
    
    ssp_list <- list("ssp126","ssp245","ssp370","ssp585")
    
    for(i in 1:length(ssp_list)){
      
      print(paste("SSP Scenario: ",ssp_list[i],sep=""))
      for(j in 1:length(mdls)){
        
        TG <- names(mdls)[j]
        
        if(file.exists(paste("Output/data/projections/",TG,"_",ssp_list[i],"_baseline_",date_projections,".RData",sep=""))) {
          projection_baseline <- readRDS(file=paste("Output/data/projections/",TG,"_",ssp_list[i],"_baseline_",date_projections,".RData",sep=""))
          
          projection_future <- readRDS(file=paste("Output/data/projections/",TG,"_",ssp_list[i],"_future_",date_projections,".RData",sep=""))
          
          
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
          print(paste("File does not exist: Output/data/projections/",TG,"_",ssp_list[i],"_baseline_",date_projections,".RData",sep=""))
        }
      }
    }
  }
  
  compute_zoop_delta(mdl_list)

#04 to provide summary statistics for projected trophic groups' relative abundance
  
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
  
  summary_stats_TG(mdl_list)
  baseline_stats_TG(mdl_list)
  
##