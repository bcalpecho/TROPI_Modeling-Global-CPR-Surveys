# ---
# title: Extracting OC-CCI data
# author: Bryan Alpecho
# date: 2025-
# output: RData for extracted chl-a
# ---

#Aims
#1. aggregate raster file of 8-day OC-CCI (based on 'wrangle-netcdf-imos' script from 2024 UQ MME Lab R Workshop)
#2. extract 8-day OC-CCI (raster data) at CPR sampling points
#3. fill-up gaps of 8-day by monthly OC-CCI values 

# Load libraries
  library(tidyverse)
  library(stars)
  library(raster) #to be replaced by package(terra)

#1. aggregate raster file of OC-CCI (based on 'wrangle-netcdf-imos' script from 2024 UQ MME Lab R Workshop)
  #set date for version control
  date <- "02032026"
  
  extract_ncdf <- function(survey_list, frequency) {
    
    ## Load data - read in list of netcdf files
      if(frequency == "eightday"){
        #eightday_OCCCI
        ncdf_list <- list.files(path = "C:/Users/power/Desktop/OC-CCI/Global", pattern = "*\\.nc4", full.names = TRUE) #temporary file source due to large file size
        #eightday_OCCCI <- list.files(path = "data_input/OCCCI/eightday/", pattern = "*\\.nc4", full.names = TRUE) #temporary file source due to large file size
        
      }else if(frequency == "monthly"){
        #monthly_OCCCI
        ncdf_list <- list.files(path = "C:/Users/power/Desktop/OC-CCI/Global monthly", pattern = "*\\.nc4", full.names = TRUE) #temporary file source due to large file size
        #monthly_OCCCI <- list.files(path = "data_input/OCCCI/monthly/", pattern = "*\\.nc4", full.names = TRUE) #temporary file source due to large file size
        
      }else{ stop("Wrong frequency indicated", call. = FALSE) }
    
    #determine study area for each survey from the survey_list
    for(i in 1:length(survey_list)){
      #set study area for extraction 
      print(survey_list[i])
      study_area <- read_csv(paste("data_input/CPR_on-process/cpr_",survey_list[i],"_metadata.csv",sep=""), col_names=TRUE) %>% 
        dplyr::select(latitude, longitude) %>% 
        st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
      
      filenames <- basename(ncdf_list)
      
        time.l <- list()
        for(j in 1:length(ncdf_list)){
          dat <- read_ncdf(ncdf_list[j], var = c('chlor_a'), make_time = T, proxy = T)
          freq <- str_extract(filenames[j], "8DAY|MONTHLY")
          
          #extract time dimension from ncdf
          time_values <- data.frame(st_get_dimension_values(dat, "time"))
          time.l <- c(time.l, list(time_values))  
          #print(raster_files[i])
          #dat <- raster(raster_files[i])
        
          # make a grid to chunk for processing 
          chunk <- st_make_grid(st_bbox(dat), cellsize = c(10, 10)) %>% st_as_sf()
          
          # loop through chunks and make coarser res
          tmp <- list()
          for(k in 1:nrow(chunk)){ 
            tmp[[k]] <- dat[st_bbox(chunk[k,])] %>% 
              st_as_stars(proxy = F) %>% 
              st_warp(cellsize = c(0.2, 0.2), crs = st_crs(4326))
          }
           
          # bind processed raster cubes back together
          dat_all <- do.call(st_mosaic, tmp)
            #turns off use of spherical geometry
            sf_use_s2(FALSE)
            #assumes planar coordinates
          
          dat_all2 <- dat_all[st_bbox(study_area)] %>%
            st_set_dimensions("band", values = time_values[,1], names = "time")
          # sequence generation of bands from 1 to 17 by 1 step
          # dat_all2 #shows summary of stars object bydimensions and attribute
         
          #save
          print(paste("file #:",j))
          print(paste("Output/raster/",freq,"/",survey_list[i],"/",sub("*\\.nc4","",filenames[j]),"_",date,".grd",sep=""))
          write_stars(dat_all2, paste("Output/raster/",freq,"/",survey_list[i],"/",sub("*\\.nc4","",filenames[j]),"_",date,".nc",sep=""))
          #write_stars(dat_all2, 'Output/raster/')
        
        }
      print(paste("Timelist saved: Output/raster/",freq,"/",survey_list[i],"/raster_time_list_",date,".RData",sep=""))
      saveRDS(time.l, file=paste("Output/raster/",freq,"/",survey_list[i],"/raster_time_list_",date,".RData",sep=""))
    }  
  }
  
  survey_list <- c("auscpr","npacific","natlantic", "socpr")
  
  #8-day OC-CCI
  extract_ncdf(survey_list, "eightday")
  #monthly OC-CCI
  extract_ncdf(survey_list, "monthly")

  
#2. extract raster data at CPR sampling points
  # set date for version control
  date <- "03032026"
  
  extract_chla <- function(survey_list, frequency){
    
    for(i in 1:length(survey_list)){  
     #read in cpr spatiotemporal coordinates
      if(frequency == "eightday"){
        print("OC-CCI frequency: 8-day")
        rast_list <- list.files(path = paste("Output/raster/8DAY/", survey_list[i], sep=""), pattern = "*\\.grd$", full.names = TRUE)
        #read time list for raster files
        time.l <- read_rds(paste("Output/raster/8DAY/",survey_list[i],"/raster_time_list.RData",sep=""))
        
      }else if(frequency == "monthly"){
        print("OC-CCI frequency: Monthly")
        rast_list <- list.files(path = paste("Output/raster/MONTHLY/", survey_list[i],"/", sep=""), pattern = "*\\.grd$", full.names = TRUE)
        #read time list for raster files
        time.l <- read_rds(paste("Output/raster/MONTHLY/",survey_list[i],"/raster_time_list.RData",sep=""))
        
      }else{ stop("Wrong frequency indicated", call. = FALSE) }
      
      cpr_coord_file <- file.path(paste("data_input/CPR_on-process/cpr_",survey_list[i],"_metadata.csv", sep=""))
      cpr_coords <- read_csv(cpr_coord_file, col_names=T) %>% 
        dplyr::select(sample_id, latitude, longitude, sampleTime_utc) %>% 
        st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
      
 
        
        #create df template for saving extracted chl-a output 
        chla_df <- cpr_coords
        for(j in 1:length(rast_list)){
          
          #insert time coordinates
          chla <- read_stars(rast_list[j]) %>% 
            st_set_dimensions("band", values = time.l[[j]][,1], names = "time")
          
          #progress update
          print(paste("file #:",j," out of ", length(rast_list),sep =""))
          print(basename(rast_list[j]))
          #extract
          chla_ext <- st_extract(chla, cpr_coords, time_column = "sampleTime_utc", interpolate_time = T)
          
          chla_df <- c(chla_df, chla_ext[1]) 
          
          #chla_df <- chla_df %>% left_join(chla_ext[1], by = c("longitude", "latitude"))
        }
      
      chla_df <- cbind(chla_df, cpr_coords %>% dplyr::select("sample_id", "sampleTime_utc"))
      saveRDS(chla_df, file=paste("Output/data/chla/",survey_list[i],"/chla_",frequency,"_",date,".rds",sep=""))
      print(paste("Extracted Chl-a Output: ","Output/data/chla/",survey_list[i],"/chla_",frequency,"_",date,".rds",sep=""))
      
    }
  }
  
  #set which surveys are use to extract chl-a
  survey_list <- c("auscpr", "natlantic", "npacific", "socpr") #can be modified to select the CPR surveys used for chl-a extraction

  #8-day OC-CCI
  extract_chla(survey_list, "eightday") 
  #monthly OC-CCI
  extract_chla(survey_list, "monthly")

#3. Fill-up gaps by monthly OC-CCI values
 
  #to be functional, convert file source from Output to data_input 
  
  fill_up_gaps <- function(survey_list){
    #setup for date and data input
    date_chla <- "03032026"
    
    for(i in 1:length(survey_list)){
      #extract from 8-day (to read output of previous function "extract_chla")
      weekly_chla <- read_rds(paste("Output/data/chla/",survey_list[i],"/chla_eightday_",date_chla,".rds",sep="")) 
      chla_df <- weekly_chla[1:57]
      chla_df$geometry = NULL
      
      chla_df <- chla_df %>%
        pivot_longer(-c("sample_id", "sampleTime_utc"), values_drop_na = T) %>% 
        rename(chla_eightday = value)
      
      #extract from monthly (to read output of previous function "extract_chla")
      monthly_chla <- read_rds(paste("Output/data/chla/",survey_list[i],"/chla_monthly_",date_chla,".rds",sep="")) 
      
      chla_df_m <- monthly_chla[1:9]
      chla_df_m$geometry = NULL
      
      #remove the rows containing only 'NA'
      chla_df_m <- chla_df_m %>%
        pivot_longer(-c("sample_id", "sampleTime_utc"), values_drop_na = T) %>% 
        rename(chla_monthly = value)
      
      #plain chla df
      chla <- weekly_chla %>% 
        dplyr::select(sample_id)
      chla$geometry = NULL
      
      #merge 8-day and monthly
      chla_merged <- chla %>% 
        left_join(chla_df_m, by="sample_id") %>%
        left_join(chla_df, by="sample_id") %>% 
        dplyr::select(-c(name.x, name.y, sampleTime_utc.y)) %>%
        rename(sampleTime_utc = "sampleTime_utc.x") %>% 
        mutate(chla = ifelse(is.na(chla_eightday), chla_monthly, chla_eightday))
      
      print(paste("Survey: ",survey_list[i],sep=""))
      print(summary(is.na(chla_merged$chla_f))) 
      
      file <- paste("Output/data/chla/",survey_list[i],"/chla_merged_",date,".rds",sep="")
      write_rds(chla_merged, file)
    }
  }

  survey_list <- c("auscpr", "natlantic", "npacific", "socpr")
  fill_up_gaps(survey_list)
  
  ##
      