####                  ####
#     FUNCTIONS BANK     #
####                  ####

#Load libraries
  packages <- c("tidyverse",
                "stars",
                "terra",
                "patchwork",
                "ggpointdensity",
                "marginaleffects",
                "glmmTMB",
                "MuMIn",
                "GGally",
                "DHARMa",
                "gstat",
                "sp",
                "sf",
                "viridis",
                "gridGraphics")
  
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

#### 1_generate_traits ####
  #1.1 Integrate Pata & Hunt (2023) assigned trait value to species list
  import_TG.PataHunt <- function(species.list){
    # Read the trait dataset
    zoop.traits <- read_csv("data_input/traits/TG_trait_subset_03-06-2025.csv")
    copepod.traits <- read_csv("data_input/traits/TG_copepods_combined_03-06-2025.csv")
    
    # Select the relevant columns
    zoop.traits <- zoop.traits %>% 
      select(taxonID, scientificName, class, order, family, genus, traitValue) %>% 
      distinct() 
    
    copepod.traits <- copepod.traits %>% 
      select(taxonID, scientificName, class, order, family, genus, traitValue) %>% 
      rename(copepod_trait = traitValue) %>% 
      distinct()
    
    summary(as.factor(species.list$traitValue))
    summary(as.factor(copepod.traits$copepod_trait))
    summary(as.factor(zoop.traits$traitValue))
    
    # Join trait with the species list
    species.list <- species.list %>% 
      left_join(zoop.traits, by = c("taxonID", "scientificName", "class", "order", "family", "genus")) %>% 
      left_join(copepod.traits, by = c("taxonID", "scientificName", "class", "order", "family", "genus")) %>% 
      #rename column of traitValue to "Pata&Hunt.TG"
      mutate(traitValue = ifelse(is.na(traitValue), "not_determined", traitValue))
    
    return(species.list)
  }
  
  #1.2 Assign 'gelatinous filter-feeders' group
  assign_filter_or_not <- function(zoop_list){
    ff.order <- c("Pyrosomida","Salpida","Dolioida","Copelata","Aplousobranchia","Phlebobranchia","Stolidobranchia")
    zoop_list %>% 
      mutate(filter_or_not = case_when(order %in% ff.order ~ "1",
                                       .default = "0")) %>% 
      relocate(filter_or_not, .after = "traitValue")
  }
  
  #1.3 check for duplicates and save finalized trait table
  detect_duplicateAssignedTrait <- function(species_list){
    
    #check which have duplicated 'traitValues'
    if(anyDuplicated(species.list$aphiaID) > 0){
    cat("Need for review!\nDetected duplicate in assigned trait values for: ")
    print(unique(species.list$scientificName[duplicated(species.list$aphiaID)]))
    }else if(anyDuplicated(species.list$aphiaID) == 0){
    cat("Looks good!\nNo detected duplicates")
    }else{ stop("Error in duplicate checking", call. = FALSE)}
    
    ##Save trait table
    write_csv(species.list, paste("output/df/TG_trait-table-",date,".csv",sep=""))
    write_csv(species.list, paste("data_input/traits/TG_trait-table-",date,".csv",sep=""))
    print(paste0("Trait table saved: TG_trait-table-",date,".csv"))
  }
  
  #1.4 compare previous and current trait table 
  compare_traitTables <- function(){
    oldTraitTable <- read_csv(paste0("data_input/traits/TG_trait-table-08082d025.csv"))
    newTraitTable <- read_csv(paste0("output/df/TG_trait-table-",date,".csv"))
    
    comparison_tibble <- oldTraitTable %>% 
      select(c("aphiaID","scientificName","FG_latest")) %>% 
      rename("TG_PreviousTraitTable" = "FG_latest") %>% 
      full_join(newTraitTable %>% select(c("aphiaID","traitValue")), by = c("aphiaID")) %>% 
      rename("TG_newTraitTable" = "traitValue")
    
    view(comparison_tibble)
  }
  
  ## See '1_generate_traits' script for full preparation of trait table.
  
#### 2_extract_chla #### 
  #2.1 aggregate raster file of OC-CCI (based on 'wrangle-netcdf-imos' script from 2024 UQ MME Lab R Workshop)
  aggregate_ncdf <- function(survey_list, frequency, date) {
    
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
      study_area <- read_csv(paste("data_input/CPR/cpr_",survey_list[i],"_metadata.csv",sep=""), col_names=TRUE) %>% 
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
        print(paste("output/raster/",freq,"/",survey_list[i],"/",sub("*\\.nc4","",filenames[j]),"_",date,".grd",sep=""))
        write_stars(dat_all2, paste("output/raster/",freq,"/",survey_list[i],"/",sub("*\\.nc4","",filenames[j]),"_",date,".nc",sep=""))
        #write_stars(dat_all2, 'output/raster/')
        
      }
      print(paste("Timelist saved: output/raster/",freq,"/",survey_list[i],"/raster_time_list_",date,".RData",sep=""))
      saveRDS(time.l, file=paste("output/raster/",freq,"/",survey_list[i],"/raster_time_list_",date,".RData",sep=""))
    }  
  }

  #2.2 extract raster data at CPR sampling points
  extract_chla <- function(survey_list, frequency){
    
    for(i in 1:length(survey_list)){  
      #read in cpr spatiotemporal coordinates
      if(frequency == "eightday"){
        print("OC-CCI frequency: 8-day")
        rast_list <- list.files(path = paste("output/raster/8DAY/", survey_list[i], sep=""), pattern = "*\\.grd$", full.names = TRUE)
        #read time list for raster files
        time.l <- read_rds(paste("output/raster/8DAY/",survey_list[i],"/raster_time_list.RData",sep=""))
        
      }else if(frequency == "monthly"){
        print("OC-CCI frequency: Monthly")
        rast_list <- list.files(path = paste("output/raster/MONTHLY/", survey_list[i],"/", sep=""), pattern = "*\\.grd$", full.names = TRUE)
        #read time list for raster files
        time.l <- read_rds(paste("output/raster/MONTHLY/",survey_list[i],"/raster_time_list.RData",sep=""))
        
      }else{ stop("Wrong frequency indicated", call. = FALSE) }
      
      print(paste0("Survey: ",survey_list[i]))
      cpr_coord_file <- file.path(paste("data_input/CPR/cpr_",survey_list[i],"_metadata.csv", sep=""))
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
        
        chla_df <- c(chla_df, chla_ext[1]) #use left-join instead
        
        #chla_df <- chla_df %>% left_join(chla_ext[1], by = c("longitude", "latitude"))
      }
      
      chla_df <- cbind(chla_df, cpr_coords %>% dplyr::select("sample_id", "sampleTime_utc"))
      saveRDS(chla_df, file=paste("output/chla/",survey_list[i],"/chla_",frequency,"_",date,".rds",sep=""))
      print(paste("Extracted Chl-a Output: ","output/chla/",survey_list[i],"/chla_",frequency,"_",date,".rds",sep=""))
      
    }
  }
  
  #2.3 Fill-up gaps by monthly OC-CCI values
  
  fill_up_gaps <- function(survey_list,date_chla){
    #setup for date and data input

    for(i in 1:length(survey_list)){
      #extract from 8-day (to read output of previous function "extract_chla")
      weekly_chla <- read_rds(paste("output/chla/",survey_list[i],"/chla_eightday_",date_chla,".rds",sep="")) 
      chla_df <- weekly_chla[1:57]
      chla_df$geometry = NULL
      
      chla_df <- chla_df %>%
        pivot_longer(-c("sample_id", "sampleTime_utc"), values_drop_na = T) %>% 
        rename(chla_eightday = value)
      
      #extract from monthly (to read output of previous function "extract_chla")
      monthly_chla <- read_rds(paste("output/chla/",survey_list[i],"/chla_monthly_",date_chla,".rds",sep="")) 
      
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
      #print(summary(is.na(chla_merged$chla_f))) 
      
      print(paste0("Merged Chl-a output: output/chla/",survey_list[i],"/chla_merged_",date,".rds"))
      file <- paste("output/chla/",survey_list[i],"/chla_merged_",date,".rds",sep="")
      write_rds(chla_merged, file)
    }
  }
  
# 3_generate_completeDF :
  #3.1 Compute for proportions of trophic groups
  compute_proportions_perSurvey <- function(abundance_list){
      #01 get trait list (check version of trait table)
      traits <- read_csv("data_input/traits/TG_trait-table-08082025.csv", col_names=T, show_col_types = F) %>% 
        rename(FG_complete = FG_latest)
      print("Trait table version date: 08082025")
      
      #02 Subset the aphiaIDs belonging to each taxonomic group and trophic group
        #zooplankton community
        aphia_noNA <- traits %>% 
          filter(!is.na(FG_complete))
      
        aphia_FF <- traits %>% 
          filter(FG_complete == "filter-feeder")
        
        aphia_omni <- traits %>% 
          filter(FG_complete == "omnivore")
        
        aphia_carni <- traits %>% 
          filter(FG_complete == "carnivore")
        
        #Copepods 
        aphia_copepods <- traits %>% 
          filter(class == "Copepoda") %>% 
          filter(!is.na(aphiaID)) 
          
          aphia_cope_omni <- aphia_copepods %>% 
            filter(FG_complete == "omnivore")
          
          aphia_cope_carni <- aphia_copepods %>% 
            filter(FG_complete == "carnivore")
        
        aphia_copepods_CalAndCyc <- traits %>% #NOTE: Only Calanoids and Cyclopoids are part of further analysis of copepods
          filter(class =="Copepoda") %>% 
          filter(order == "Calanoida" | order == "Cyclopoida")
        
        #Calanoids
        aphia_calanoids <- traits %>% 
          filter(order == "Calanoida")
        
          aphia_calanoid_omni <- aphia_calanoids %>% 
            filter(FG_complete == "omnivore")
          
          aphia_calanoid_carni <- aphia_calanoids %>% 
            filter(FG_complete == "carnivore")
        
        #Cyclopoids
        aphia_cyclopoids <- traits %>% 
          filter(order == "Cyclopoida")
          
          aphia_cyclopoid_omni <- aphia_cyclopoids %>% 
            filter(FG_complete == "omnivore")
          
          aphia_cyclopoid_carni <- aphia_cyclopoids %>% 
            filter(FG_complete == "carnivore")
        
        #Other Copepods
        aphia_otherCopepods <- traits %>%
          filter(class =="Copepoda") %>%
          filter(order != "Calanoida" & order != "Cyclopoida")

          aphia_otherCopepods_omni <- aphia_otherCopepods %>%
            filter(FG_complete == "omnivore")

          aphia_otherCopepods_carni <- aphia_otherCopepods %>%
            filter(FG_complete == "carnivore")
        
      #03 compute for the total abundance for each trophic group using the prior aphiaID subsets
      for(i in 1:length(abundance_list)){
        #get metadata
          filenames <- basename(abundance_list[i])
          survey <- str_extract(filenames, "(?<=_)[^_]+") 
          #date <- str_extract(filenames, "[^_]+$")
        
        #loop progress
          print(paste("CPR survey in process : ",survey,sep=""))
          #print(date[i])
        
        #read in cpr abundance data
          cpr <- read_csv(abundance_list[i], col_names = T, name_repair = "minimal", show_col_types = F)
          cpr_dat <- cpr %>%  select(-c(sample_id))
        
        #compute for overall zooplankton abundance with assigned trophic group or taxa with non-'NA' trophic group
          cpr_dat_noNAaphia <- cpr_dat %>% 
            select(matches(as.character(aphia_noNA$aphiaID)))
          zoopTotal_sum <- rowSums(cpr_dat_noNAaphia, na.rm = T)
       
        # sums within zooplankton community
          # compute for total abundance of filter-feeder
            FF <- cpr_dat_noNAaphia %>% select(.,matches(as.character(aphia_FF$aphiaID)))
            FF_sum <- rowSums(FF, na.rm = T)
          # compute for total abundance of omnivores
            Omni <- cpr_dat_noNAaphia %>%  select(.,matches(as.character(aphia_omni$aphiaID)))
            Omni_sum <- rowSums(Omni, na.rm = T)
          # compute for total abundance of carnivores
            Carni <- cpr_dat_noNAaphia %>% select(.,matches(as.character(aphia_carni$aphiaID)))
            Carni_sum <- rowSums(Carni, na.rm = T)

        # sums within copepods
            # subset copepods using the aphiaID subset for all copepods
            cpr_dat_copepods <- cpr_dat_noNAaphia %>% 
              select(matches(as.character(aphia_copepods$aphiaID))) 
            
            # subset copepods using the aphiaID subset for Calanoid and Cyclopoids
            cpr_dat_copepods_CalAndCyc <- cpr_dat_noNAaphia %>% 
              select(matches(as.character(aphia_copepods_CalAndCyc$aphiaID))) 
            
              # compute for total abundance of omnivorous Calanoid and Cyclopoids
              Cope_omni <- cpr_dat_copepods_CalAndCyc %>%  select(.,matches(as.character(aphia_cope_omni$aphiaID)))
              Cope_omni_sum <- rowSums(Cope_omni, na.rm = T)
              
              # compute for total abundance of carnivorous Calanoid and Cyclopoids
              Cope_carni <- cpr_dat_copepods_CalAndCyc %>% select(.,matches(as.character(aphia_cope_carni$aphiaID)))
              Cope_carni_sum <- rowSums(Cope_carni, na.rm = T)
        
        # sums within copepod orders of Calanoida and Cyclopoida
            #omnivore calanoids
              Cal_omni <- cpr_dat_copepods_CalAndCyc %>%  select(.,matches(as.character(aphia_calanoid_omni$aphiaID)))
              Cal_omni_sum <- rowSums(Cal_omni, na.rm = T)
            #carnivore calanoids
              Cal_carni <- cpr_dat_copepods_CalAndCyc %>%  select(.,matches(as.character(aphia_calanoid_carni$aphiaID)))
              Cal_carni_sum <- rowSums(Cal_carni, na.rm = T)
              
            #omnivore cyclopoids
              Cyc_omni <- cpr_dat_copepods_CalAndCyc %>%  select(.,matches(as.character(aphia_cyclopoid_omni$aphiaID)))
              Cyc_omni_sum <- rowSums(Cyc_omni, na.rm = T)
            #carnivore cyclopoids
              Cyc_carni <- cpr_dat_copepods_CalAndCyc %>%  select(.,matches(as.character(aphia_cyclopoid_carni$aphiaID)))
              Cyc_carni_sum <- rowSums(Cyc_carni, na.rm = T)
          
            #omnivorous "other copepods"
              otherCopepods_omni <- cpr_dat_copepods %>%  select(.,matches(as.character(aphia_otherCopepods_omni$aphiaID)))
              otherCopepods_omni_sum <- rowSums(otherCopepods_omni, na.rm = T)
            #carnivorous "other copepods"
              otherCopepods_carni <- cpr_dat_copepods %>%  select(.,matches(as.character(aphia_otherCopepods_carni$aphiaID)))
              otherCopepods_carni_sum <- rowSums(otherCopepods_carni, na.rm = T) #no records of carnivorous zooplankton
              
        # bind all sum
          trait_sums <- data.frame(rbind(cpr$sample_id,zoopTotal_sum, FF_sum, Omni_sum, Carni_sum, Cope_omni_sum, Cope_carni_sum, Cal_omni_sum, Cal_carni_sum, Cyc_omni_sum, Cyc_carni_sum, otherCopepods_omni_sum, otherCopepods_carni_sum)) %>% 
            t() %>% 
            data.frame() %>% 
            rename(sample_id = V1)
          
          # print(paste("output/data/sums/df_sums_",survey[i],"_",date, sep=""))
          # write_csv(trait_sums, paste("output/data/sums/df_sums_",survey[i],"_",date, sep=""))
          print(paste("Saved file: ",survey,"/df_sums_",date,".csv", sep=""))
          write_csv(trait_sums, paste("output/df/",survey,"/df_sums_",date,".csv", sep=""))
      }
    }
  
  #3.2 Combining variables altogether into a dataframe
  generate_df_perSurvey <- function(date){
    chl_date <- date 
    trait_date <- date
    
    metadata_files <- list.files(path = "data_input/CPR/", pattern = "*\\metadata.csv", full.names = TRUE)
    
    filenames <- basename(metadata_files)
    
    for(i in 1:length(filenames)){
      
      #00 get survey and coordinates
      file_survey <- str_extract(filenames, "(?<=_)[^_]+")
      print(paste("survey: ",file_survey[i],sep=""))
      df <- read_csv(metadata_files[i], show_col_types = F) 
      
      #01 start df with survey name 
      df <- df %>% 
        mutate(survey = file_survey[i]) 
      
      #02 tow_no
      if (file_survey[i] == "auscpr") {
        df <- df %>% mutate(tow_no = paste0("auscpr_",TripCode)) # Adjust to your actual column
      } else if (file_survey[i] == "socpr") {
        df <- df %>% mutate(tow_no = paste0("socpr_",Tow_Number))
      } else if (file_survey[i] == "npacific"){
        df <- df %>% mutate(tow_no = paste0("npacific_",str_extract(sample_id, "^.*?(?=-)")))
      } else if (file_survey[i] == "natlantic"){
        df <- df %>% mutate(tow_no = paste0("natlantic_",str_extract(sample_id, "^.*?(?=-)")))
      } else { stop("Unidentified CPR survey")}
      
      #03 tow_days
      df <- df %>% 
        group_by(tow_no) %>% 
        mutate(tow_days = as.double(difftime(sampleTime_utc, min(sampleTime_utc), units = "days"))) %>% 
        ungroup()
      
      #04 Longhurst
      #load longhurst provinces
      LH <- st_read(dsn = "data_input/LonghurstProvinces/", layer = "Longhurst_world_v4_2010")
      #downloaded from https://marineregions.org/downloads.php
      
      #to reinsert longitude and latitue values
      coords <- read_csv(metadata_files[i], show_col_types = F) %>% 
        dplyr::select(c("sample_id", "latitude", "longitude")) 
      
      #see resolved issue (Github issue #1732 for "sf" package under issues on s2)
      #st_is_valid(LH) #under review
      LH <- st_make_valid(LH)
      
      #set geometry and match crs
      df <- df %>% 
        st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
        st_transform(crs = st_crs(LH)) 
      
      #intersect LH and CPR points
      df <- st_intersection(df, LH) %>% 
        dplyr::select(-c("ProvDescr", "geometry")) %>% 
        left_join(coords, by="sample_id") %>% 
        relocate(c("longitude","latitude"), .after=(sample_id)) %>% 
        rename("longhurst" = "ProvCode")
      
      #05 Chl-a
      ## Extracted chlorophyll matched by time and space
      chl <- read_rds(paste("output/chla/",file_survey[i],"/chla_merged_",chl_date,".rds",sep="")) %>% select("sample_id","chla")
      df <- df %>% 
        left_join(chl, by = "sample_id")
      
      #06 Traits (TG)
      traits <- read_csv(paste("output/df/",file_survey[i],"/df_sums_",trait_date,".csv",sep=""), col_names =T, show_col_types =  F) 
      df <- df %>% 
        left_join(traits, by = "sample_id")
      df$geometry = NULL
      
      #07 export df
      
      print(paste("Filename output: ",file_survey[i],"/df_complete_",date,".csv", sep=""))
      write_rds(df, paste("output/df/",file_survey[i],"/df_complete_",date,".rds",sep=""))
    }
  }
  
  #3.3 Update dataframes (Helper functions)
    # ##3A. Update SURVEY-SPECIFIC dataframe with proportions of trophic groups
    # #NOTE: file list not yet in sync with compute_sums function
    # #set date for version of updated proportions
    # date <- "01072025"
    # 
    # #
    # survey <- c("auscpr","natlantic","npacific","socpr")
    # file.list <- NA
    # for(i in 1:length(survey)){
    #   file.list[i] <- list.files(path = paste("output/data/",survey[i],"/",sep=""), pattern = "*\\complete_04082025.rds", full.names = TRUE)
    # }
    # remove(i)
    
    update_df_traits <- function(previousDF_list){
      for(i in 1:length(previousDF_list)){
        #get metadata
        filenames <- basename(previousDF_list[i])
        survey <- str_extract(previousDF_list[i], "(?<=/data/)[^/]+") 
        date_oldVersion <- str_extract(filenames, "[^_]+$")
        
        #progress
        print(paste("Update DF of ",survey," from version ",date_oldVersion," to ",date,sep=""))
        
        #read in previous df
        df_prev <- paste("output/data/",survey,"/",filenames,sep="")
        df <- read_rds(df_prev)
        
        #remove previous traits
        df <- df %>% select(c(sample_id, longitude, latitude, sampleTime_utc, tow_no, tow_days, longhurst, chla))
        
        #load updated traits
        traits_file <- paste("output/data/",survey,"/df_sums_",date,".csv", sep="")
        traits <- read_csv(traits_file)
        
        #update df
        df_new <- df %>% 
          left_join(traits, by = "sample_id")
        
        #save updated df
        file_output <- paste("output/data/",survey,"/df_complete_",date,".rds",sep="")
        print(file_output)
        write_rds(df_new, file_output)
      } 
    }
    
    
    ##3B. Update SURVEY-SPECIFIC dataframe with chl-a 
    # #set dates for version of chl-a and previous df
    # date <- "16092025"
    # date_previousDF <- "16092025"
    # 
    # #
    # survey <- c("auscpr","natlantic","npacific","socpr")
    # file.list <- NA
    # for(i in 1:length(survey)){
    #   file.list[i] <- list.files(path = paste("output/data/",survey[i],"/",sep=""), pattern = paste("*\\complete_",date_previousDF,".rds",sep=""), full.names = TRUE)
    # }
    # remove(i)
    
    update_df_chla <- function(previousDF_list){
      for(i in 1:length(previousDF_list)){
        #get metadata
        filenames <- basename(previousDF_list[i])
        survey <- str_extract(previousDF_list[i], "(?<=/data/)[^/]+") 
        date_oldVersion <- str_extract(filenames, "[^_]+$")
        
        #progress
        print(paste("Update DF of ",survey," from version ",date_oldVersion," to ",date,sep=""))
        
        #read in previous df
        df_prev <- paste("output/data/",survey,"/",filenames,sep="")
        df <- read_rds(df_prev)
        
        #remove previous chla
        df <- df %>% dplyr::select(-chla)
        
        #load updated chla
        #CHECK the file name
        chla_file <- paste("output/data/",survey,"/chla_monthly_",date,".rds", sep="")
        chla <- read_rds(chla_file) %>% select(c("sample_id","chla"))
        
        #update df
        df_new <- df %>% 
          left_join(chla, by = "sample_id")
        
        #save updated df
        #CHECK the file name
        file_output <- paste("output/data/",survey,"/df_complete_",date,"_monthlychla.rds",sep="")
        print(file_output)
        write_rds(df_new, file_output)
      } 
    }

  
  #3.4 generate GLOBAL dataframe and compute for ratios
    generate_globalCPR_dataframe <- function(date_newVersion){
      survey_list <- c("auscpr","natlantic","npacific","socpr")
      for(i in 1:length(survey_list)){
        print(paste0("Reading dataframe of: ",survey_list[i]))
        assign(paste(survey_list[i]), read_rds(paste("output/df/",survey_list[i],"/df_complete_",date_newVersion,".rds",sep="")) %>%  
                 mutate(survey = survey_list[i]) %>%  
                 select(c("survey","sample_id", "latitude", "longitude", "sampleTime_utc", "tow_no", "tow_days", "longhurst", "chla", "zoopTotal_sum","FF_sum","Omni_sum","Carni_sum","Cope_omni_sum","Cope_carni_sum", "Cal_omni_sum", "Cal_carni_sum","Cyc_omni_sum","Cyc_carni_sum")))
      }
      
      df <- auscpr %>% 
        rows_insert(socpr, by = "sample_id") %>% 
        rows_insert(npacific, by ="sample_id") %>% 
        rows_insert(natlantic, by = "sample_id")
      
      #transform chl-a and compute for trophic group ratios
      df <- df %>% 
        mutate(chla_sqrt = sqrt(chla)) %>%
        mutate(Cope_ROC = Cope_omni_sum/(Cope_omni_sum + Cope_carni_sum)) %>% 
        mutate(Cope_RCO = Cope_carni_sum/(Cope_omni_sum + Cope_carni_sum)) %>% 
        mutate(ROC = Omni_sum/(Omni_sum + Carni_sum)) %>% 
        mutate(RCO = Carni_sum/(Omni_sum + Carni_sum)) %>% 
        mutate(RFF = FF_sum/zoopTotal_sum) %>% 
        mutate(ROC_Cal = Cal_omni_sum/(Cal_omni_sum + Cal_carni_sum)) %>% 
        mutate(RCO_Cal = Cal_carni_sum/(Cal_omni_sum + Cal_carni_sum)) %>% 
        mutate(ROC_Cyc = Cyc_omni_sum/(Cyc_omni_sum + Cyc_carni_sum)) %>% 
        mutate(RCO_Cyc = Cyc_carni_sum/(Cyc_omni_sum + Cyc_carni_sum)) 
      
      df <- df %>% 
        mutate(RFF_SVT_zib = case_when(
          RFF == 1 ~ ((RFF*325075)+0.5)/325076,
          .default = RFF)) %>% 
        mutate(ROC_SVT_zib = case_when(
          ROC == 1 ~ ((ROC*325075)+0.5)/325076,
          .default = ROC)) %>% 
        mutate(RCO_SVT_zib = case_when(
          RCO == 1 ~ ((RCO*325075)+0.5)/325076,
          .default = RCO)) %>% 
        mutate(ROC_Cal_SVT_zib = case_when(
          ROC_Cal == 1 ~ ((ROC_Cal*325075)+0.5)/325076,
          .default = ROC_Cal)) %>% 
        mutate(RCO_Cal_SVT_zib = case_when(
          RCO_Cal == 1 ~ ((RCO_Cal*325075)+0.5)/325076,
          .default = RCO_Cal)) %>% 
        mutate(ROC_Cyc_SVT_zib = case_when(
          ROC_Cyc == 1 ~ ((ROC_Cyc*325075)+0.5)/325076,
          .default = ROC_Cyc)) %>% 
        mutate(RCO_Cyc_SVT_zib = case_when(
          RCO_Cyc == 1 ~ ((RCO_Cyc*325075)+0.5)/325076,
          .default = RCO_Cyc)) 
      
      #finalize df for omnivores and carnivores
      df <- df %>% 
        #relocate(survey, .before= sample_id) %>% 
        relocate(chla_sqrt, .after=chla)
      
      print(paste0("Dataframe for omnivores and carnivores saved: df_complete_",date_newVersion,".rds"))
      write_rds(df, paste0("output/df/global_df_complete_",date_newVersion,".rds"))   
      write_rds(df, paste0("data_input/global_df_complete_",date_newVersion,".rds")) #second copy for later steps 
      
      #finalize df for filter-feeders
      
      LH_modification <- c("ARCT", "BENG", "CNRY", "EAFR", "GFST", "GUIN", "MEDI", "NADR", "NASE", "NASW", "NATR", "NECS", "NPPF", "NWCS", "SARC")
      df_modified <- df %>% filter(!(df$longhurst %in% LH_modification))
      
      print(paste0(c("Modification of dataframe for filter-feeders is the removal of the following Longhurst Provinces: ", LH_modification), collapse = " "))
      print(paste0("Dataframe for filter-feeders saved: global_df_complete_Filter_",date_newVersion,".rds"))
      write_rds(df_modified, paste0("output/df/global_df_complete_Filter_",date_newVersion,".rds"))   
      write_rds(df_modified, paste0("data_input/global_df_complete_Filter_",date_newVersion,".rds")) #second copy for later steps 
      
    }
    
# 4_model_globalCPR :
    
    #4.1 to generate a summary of proportion contributed by fixed and random effects to the total variance of the model
    summary_mdlVariance <- function(mdls){
      tibble_variance_summary <- tibble(mdl = c("Filter_mdl_zib", "Omni_mdl_zib", "Carni_mdl_zib"),
                                        FixedEffect_proportion = c(1,2,3),
                                        RandomEffect_proportion = c(1,2,3),
                                        TowWithinSurvey_proportion = c(1,2,3),
                                        LonghurstProvinces_proportion = c(1,2,3))  
      
      for(i in 1:length(mdls)){
        mdl_summary <- summary(mdls[[i]])
        TG <- names(mdls[i])
        if(TG == "Carni"){
          if(exists("Carni_mdl_zib_R2") == FALSE){
            mdl_Rsquare <- MuMIn::r.squaredGLMM(mdls[[i]])
            print("missing")
          }else{ mdl_Rsquare <- Carni_mdl_zib_R2}
          
        }else if(TG == "Omni"){
          if(exists("Omni_mdl_zib_R2") == FALSE){
            mdl_Rsquare <- MuMIn::r.squaredGLMM(mdls[[i]])
            print("missing")
          }else{ mdl_Rsquare <- Omni_mdl_zib_R2}
          
        }else if(TG == "Filter"){
          if(exists("Filter_mdl_zib_R2") == FALSE){
            mdl_Rsquare <- MuMIn::r.squaredGLMM(mdls[[i]])
            print("missing")
          }else{ mdl_Rsquare <- Filter_mdl_zib_R2}
          
        }else(stop ("Error in determining R^2: misidentified trophic group"))
        
        #(proportion of variation) overall contribution of fixed effects to mdl variance
        var_FE <- mdl_Rsquare[1]
        #(proportion of variation) overall contribution of REs to mdl variance
        var_RE <- mdl_Rsquare[2] - mdl_Rsquare[1]
        #Conditional R2
        var_total <- mdl_Rsquare[2]
        
        #get variance
        #Tow within Survey
        var_survey <- mdl_summary$varcor$cond$`survey:tow_no`[1] #gets the variance 
        #Longhurst Provinces
        var_lh <- mdl_summary$varcor$cond$`longhurst`[1]      
        
        #contribution of each RE in relative terms
        #Tow within Survey
        rel_var_survey <- (var_survey)/(var_survey + var_lh) * var_RE
        #Longhurst Provinces
        rel_var_LH <- (var_lh)/(var_survey + var_lh) * var_RE
        
        tibble_variance_summary$FixedEffect_proportion[i] <- var_FE
        tibble_variance_summary$RandomEffect_proportion[i] <- var_RE
        tibble_variance_summary$TowWithinSurvey_proportion[i] <- rel_var_survey
        tibble_variance_summary$FixedEffect_proportion[i] <- rel_var_LH
      }        
      view(tibble_variance_summary)
    } 
    
    
# 5_assess_models : 
  
  # Create a publication-ready theme (Adapted from 2025 UQ MME Lab Winter R Workshop)
  pub_theme <- theme_classic(base_size = 10, base_family = "sans") + # Family including Arial, Helvetica, Futura, Verdana, and Calibri
    theme(
      # Axis styling
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.5),
      
      # Legend styling
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 8),
      legend.position = "inside",
      #legend.position.inside = c(0.85, 0.25),
      legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
      legend.margin = margin(2, 2, 2, 2),
      
      # Panel styling
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      
      # Remove grid lines for cleaner look
      panel.grid = element_blank()
    )

  #5.1 quantile-quantile plot to assess normality of residuals
  plot_QQ <- function(mdl_list) {
    
    for(i in 1:length(mdl_list)){
    
    TG <- names(mdl_list[i])  
    print(paste0("QQ plot in process for: ", TG))
    simulationOutput <- simulateResiduals(fittedModel = mdl_list[[i]], plot = FALSE)
    
    qq_plot <- wrap_elements(panel = ~plotQQunif(
      simulationOutput, 
      testUniformity = FALSE, 
      testOutliers = FALSE, 
      testDispersion = FALSE
    ))
    
    if(TG == "Omni"){
      Omni_qq <- qq_plot
    }else if(TG == "Carni"){
      Carni_qq <- qq_plot
    }else if(TG == "Filter"){
      Filter_qq <- qq_plot
    }else{ stop("Error in generating QQ plot") }
    # To save individual plots
    # ggsave(paste0("output/plots/",TG,"_QQplot_",date,".png"), plot = qq_plot,
    #        width = 8, height = 10, dpi = 300)
    # print(paste0("QQ plot saved: ",TG,"_QQplot_",date,".png"))
    
    }
    # To save merged QQ plots 
    print(paste0("Plotting and combining the plots"))
    QQ_patch <- Omni_qq + Carni_qq + Filter_qq +
      plot_layout(design = "ABC") +
      plot_annotation(tag_levels = "A", 
                      theme = theme(plot.title = element_text(size = 16, face = "bold"),
                                    plot.margin = margin(2,2,2,2)))
    
    ggsave(paste("output/plots/QQplot_",date,".png",sep=""), plot = QQ_patch,
           width = 10, height = 4, dpi = 600)
    print(paste0("Plot saved: QQplot_",date,".png"))
    
  }
  
  #5.2 mean variance plot to assess homogeneity of variance
  plot_meanVariance <- function(mdl_list){
    date_newVersion <- date
    df <- read_rds(paste0("data_input/global_df_complete_",date_newVersion,".rds"))
    df_filter <- read_rds(paste0("data_input/global_df_complete_Filter_",date_newVersion,".rds"))
    
    for(i in 1:length(mdl_list)){
      TG <- names(mdl_list[i])  
      print(paste0("Model in process: ",TG))
      
      if(TG == "Omni"){
        df_Omni_noNA <- df %>% filter(!is.na(chla_sqrt)) %>% filter(!is.na(ROC_SVT_zib)) %>% filter(!is.na(tow_days))
        Omni_meanVariance <- ggplot(df_Omni_noNA) + aes(fitted(Omni_mdl_zib), resid(Omni_mdl_zib, type = "pearson")) +
          geom_hex(bins = 80) + geom_smooth(method = "loess") + pub_theme + theme(legend.position = c(.9, .8)) +
          scale_fill_viridis(begin = 0, end = .9, option = "C", limit = range(c(1,1200))) +
          theme(plot.title = element_text(hjust = 0.5), ) +
          labs(fill = "Frequency") + ylab("Residual") + xlab("Fitted Value")
        
      }else if(TG == "Carni"){
        df_Carni_noNA <- df %>% filter(!is.na(chla_sqrt)) %>% filter(!is.na(RCO_SVT_zib)) %>% filter(!is.na(tow_days))
        Carni_meanVariance <- ggplot(df_Carni_noNA) + aes(fitted(Carni_mdl_zib), resid(Carni_mdl_zib, type = "pearson")) +
          geom_hex(bins = 80) + geom_smooth(method = "loess") + pub_theme + theme(legend.position = c(.9, .8)) +
          scale_fill_viridis(begin = 0, end = .9, option = "C", limit = range(c(1,1200))) +
          theme(plot.title = element_text(hjust = 0.5), ) +
          labs(fill = "Frequency") + ylab("Residual") + xlab("Fitted Value")
        
      }else if(TG == "Filter"){
        df_FF_noNA <- df_filter %>% filter(!is.na(chla_sqrt)) %>% filter(!is.na(RFF_SVT_zib)) %>% filter(!is.na(tow_days))
        Filter_meanVariance <- ggplot(df_FF_noNA) + aes(fitted(Filter_mdl_zib), resid(Filter_mdl_zib, type = "pearson")) +
          geom_hex(bins = 80) + geom_smooth(method = "loess") + pub_theme + theme(legend.position = c(.9, .8)) +
          scale_fill_viridis(begin = 0, end = .9, option = "C", limit = range(c(1,1200))) +
          theme(plot.title = element_text(hjust = 0.5), ) +
          labs(fill = "Frequency") + ylab("Residual") + xlab("Fitted Value")
        
      }else{stop("Error in designated trophic group: possible misidentified model")}
    }
    
    ## Residual plot
    # #to save individual plots
    # ggsave(paste("output/illustrations/Filter_meanVariance_",date,".png",sep=""), plot = Filter_meanVariance,
    #        width = 9, height = 6, dpi = 300)
    #
    # ggsave(paste("output/illustrations/Omni_meanVariance_",date,".png",sep=""), plot = Omni_meanVariance,
    #        width = 9, height = 6, dpi = 300)
    #
    # ggsave(paste("output/illustrations/Carni_meanVariance_",date,".png",sep=""), plot = Carni_meanVariance,
    #        width = 9, height = 6, dpi = 300)
    # #
    
    # design <- "
    #             A
    #             B
    #             C
    #           "
    # 
    # (meanVariance_patch <- Filter_meanVariance + Omni_meanVariance + Carni_meanVariance +
    #    plot_layout(design = design) +
    #    plot_annotation(
    #      tag_levels = "A",
    #      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    #    )
    # )
    
    print(paste0("Plotting and combining the plots"))
    meanVariance_patch <- Omni_meanVariance / Carni_meanVariance / Filter_meanVariance +
      plot_annotation(tag_levels = "A", 
                      theme = theme(plot.title = element_text(size = 16, face = "bold"),
                                    plot.margin = margin(10,10,10,10)))
    
    ggsave(paste("output/plots/meanVariance_",date,".png",sep=""), plot = meanVariance_patch,
           width = 10, height = 3, dpi = 300)
    print(paste0("Plot saved: meanVariance_",date,".png"))
  }
  
  #5.3 to plot the intercepts of Longhurst Provinces
  plot_Longhurst <- function(model_list){
    for(i in 1:length(model_list)){
      REs <- ranef(model_list[[i]], condVar = TRUE) #Update the model input
      TG <- names(model_list[i])
      print(TG)
      
      REs_Longhurst <- REs$cond$longhurst
      
      ### Longhurst Province ###
      qq <- attr(REs_Longhurst, "condVar")
      
      # Extract intercepts
      rand.interc <- REs$cond$longhurst
      
      # Make a dataframe for plotting
      df_plot <- data.frame(Intercepts = REs$cond$longhurst[,1],
                            sd.interc = 2*sqrt(qq[,,1:length(qq)]),
                            lev.names = factor(rownames(rand.interc))) %>% 
        arrange(Intercepts) %>% 
        within({  # Reorder levels
          lev.names <- factor(as.character(lev.names),
                              as.character(lev.names))
        })
      
      re_lh_plot <- ggplot(df_plot, aes(lev.names, Intercepts)) + pub_theme + 
        geom_hline(yintercept=0, linetype = "dashed", color = "black") +
        geom_errorbar(aes(ymin=Intercepts-sd.interc, 
                          ymax=Intercepts+sd.interc),
                      width = 0,color="black") +
        geom_point(color = "black", size = 2) +
        guides(size="none",shape="none") + 
        theme(axis.text.x=element_text(size=10), 
              axis.title.x=element_text(size=13),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.x = element_blank()) +
        coord_flip() + 
        labs(y = "Intercept", x = "Longhurst Provinces")
      
      #annotate("text", y = -2.3, x = 28, label = "(D)", size = 5)
      
      ggsave(paste0("output/plots/",TG,"_Longhurst-",date,".png"), plot = re_lh_plot,
             width = 5, height = 6, dpi = 300)
      print(paste0("File saved: ",TG,"_Longhurst-",date,".png"))
    }
  }
  
  #5.4 to plot point density plots of Tow within Survey slope and intercept
  ### Tow slope and intercept ###
  
  plot_TowSlopeAndIntercept <- function(model_list){
    for(i in 1:length(model_list)){
      REs <- ranef(model_list[[i]], condVar = TRUE) #Update the model input
      TG <- names(model_list[i])
      print(TG)
      
      # Create a dataframe for plotting
      REs_towNo <- REs$cond$`survey:tow_no`
      towDf <- as.data.frame(REs_towNo) %>% setNames(c("Intercept", "Slope"))
      
      # Specify frequency breaks
      my_breaks <- c(1,5,30,175,1000)
      
      #
      library(ggpointdensity)
      re_tow_plot_point_density <- ggplot(towDf, aes(Intercept, Slope)) +
        pub_theme +
        geom_pointdensity() +
        scale_colour_viridis_c() +
        labs(fill = "Density", x = "Tow within Survey Intercept", y = "Tow within Survey Slope") +
        theme(
          axis.text = element_text(size = 14),    # Adjust the size for axis numbers/text
          legend.position = c(0.9,0.9)
          ) 
      
      #update file name
      ggsave(paste0("output/plots/",TG,"_Tow-Survey-",date,".png"), plot = re_tow_plot_point_density,
             width = 7, height = 8, dpi = 300)
      print(paste0("File saved: ",TG,"_Tow-Survey-",date,".png"))
    }  
  }
  
  #5.5 to plot residuals of models with and without the random effects
  plot_residuals <- function(model_list){
    date_newVersion <- date
    df <- read_rds(paste0("data_input/global_df_complete_",date_newVersion,".rds"))
    df_filter <- read_rds(paste0("data_input/global_df_complete_Filter_",date_newVersion,".rds"))
    
    for(i in 1:length(model_list)){
      vc <- VarCorr(model_list[[i]])
      TG <- names(model_list[i])
      print(TG) 
      
      if(TG == "Filter"){
        df_noNA <- df_filter %>% filter(!is.na(chla_sqrt)) %>% filter(!is.na(RFF_SVT_zib)) %>% filter(!is.na(tow_days))
      }else if(TG == "Carni"){ 
        df_noNA <- df %>% filter(!is.na(chla_sqrt)) %>% filter(!is.na(RCO_SVT_zib)) %>% filter(!is.na(tow_days))
      }else if(TG == "Omni"){ 
        df_noNA <- df %>% filter(!is.na(chla_sqrt)) %>% filter(!is.na(ROC_SVT_zib)) %>% filter(!is.na(tow_days))
      }else{stop("Error in reading dataframe")}
      
      #print(summary(df_noNA))
      
      # Data format for variograms
      autocorData <- data.frame(Longitude = df_noNA$longitude, 
                                Latitude = df_noNA$latitude, 
                                resids = resid(model_list[[i]])) %>% # Extract residuals
        within({
          signres <- sign(resids)
        })
      
      ## Plot the residuals in space
      world <- map_data("world")
      
      # Change longitude so it matches up with the world map
      autocorData$Longitude[autocorData$Longitude < (-170)] <- autocorData$Longitude[autocorData$Longitude < (-170)] + 360
      
      # Bubble plot WITH random effects
      residuals_withRandom <- ggplot(data = autocorData, aes(x = Longitude, y = Latitude)) + 
        geom_map(data=world, map=world, aes(x = long, y = lat, map_id=region),
                 color="white", fill="gray94", size=0.08) + 
        geom_point(aes(size = abs(resids)/4, color = sign(resids)), shape = 1,
                   alpha = 0.4) + 
        scale_size_continuous(range=c(.1,4)) + 
        scale_colour_gradient(low = "springgreen3", high = "magenta3") + 
        ylab(NULL) + xlab(NULL) + 
        #annotate("text", x = -190, y = 90, label = "(b)", size = 9) +
        guides(colour = "none", size = guide_legend(title = "Magnitude"))
      
      #to save individual plot
      # ggsave(paste("output/plots/",TG,"_residualMap_withRandom_",date,".png",sep=""), plot = residuals_withRandom,
      #        width = 9, height = 6, dpi = 300)
      # print(paste("File saved (residual map with random effects): output/plots/",TG,"_residualMap_withRandom_",date,".png",sep=""))
      
      # Model without random effects 
      if(TG == "Filter"){ 
             Mdl_noRandom <- glmmTMB(RFF_SVT_zib ~ chla_sqrt + survey, 
                                     ziformula = ~1,
                                     data = df_noNA, family = beta_family(link = "logit"))
      }else if(TG == "Carni"){ 
                    Mdl_noRandom <- glmmTMB(RCO_SVT_zib ~ chla_sqrt + survey, 
                                            ziformula = ~1,
                                            data = df_noNA, family = beta_family(link = "logit"))
      }else if(TG == "Omni"){
                    Mdl_noRandom <- glmmTMB(ROC_SVT_zib ~ chla_sqrt + survey, 
                                            ziformula = ~1,
                                            data = df_noNA, family = beta_family(link = "logit"))
      }else{ stop("Error in generating model without random effects")}
      
      autocorData$resids <- resid(Mdl_noRandom)
      autocorData$signres<- sign(autocorData$resids)
      
      # Bubble plot WITHOUT random effects
      residuals_noRandom <- ggplot(data = autocorData, aes(x = Longitude, y = Latitude)) + 
        geom_map(data=world, map=world, aes(x = long, y = lat, map_id=region),
                 color="white", fill="gray94", size=0.08) + 
        geom_point(aes(size = abs(resids)/4, color = sign(resids)), shape = 1,
                   alpha = 0.4)  + 
        scale_size_continuous(range=c(.1,4)) + 
        scale_colour_gradient(low = "springgreen3", high = "magenta3") + 
        ylab(NULL) + xlab(NULL) + 
        #annotate("text", x = -190, y = 90, label = "(a)", size = 9) +
        guides(colour = "none", size = guide_legend(title = "Magnitude"))
      
      # #save plot of residuals without random effects
      # ggsave(paste("output/plots/",TG,"_residualMap_noRandom_",date,".png",sep=""), plot = residuals_noRandom,
      #        width = 9, height = 6, dpi = 300)
      # print(paste("File saved (residual map without random effects): output/plots/",TG,"_residualMap_noRandom_",date,".png",sep=""))
      
      #combine
      design <- "
                  A
                  B
                "
      
      #update file names
      (residualMap <- residuals_withRandom + residuals_noRandom + 
         plot_layout(design = design) +
         plot_annotation(
           tag_levels = "A",
           theme = theme(plot.title = element_text(size = 16, face = "bold"))
         )
      )   
      ggsave(paste("output/plots/residualMap_",TG,"_",date,".png",sep=""), plot = residualMap,
             width = 8, height = 10, dpi = 300)
      print(paste("File saved: residualMap_",TG,"_",date,".png",sep=""))
      print("Map of residuals for model (A) with and (B) without random effects")
      
    }
  }

# 6b_predict_globalCPR : 

  #6.1 to predict zooplankton trophic group
  
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
    
    print(paste0("Generating predictions under ",ssp_scenario," scenario"))
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
        saveRDS(esm_pred_merged, file=paste("Output/projections/",TG,"_",ssp_scenario,"_",esm_level,"_",date,".RData",sep=""))
        print(paste("Output saved: projections/",TG,"_",ssp_scenario,"_",esm_level,"_",date,".RData",sep=""))
      }
    }
  }
  
  
  #6.2 to compute for delta of trophic groups between 2015 and 2100
  
  compute_zoop_delta <- function(mdls){
    #stable projections - 24 & 25022026
    date_projections <- date
    
    ssp_list <- list("ssp126","ssp245","ssp370","ssp585")
    
    for(i in 1:length(ssp_list)){
      
      print(paste("SSP Scenario: ",ssp_list[i],sep=""))
      for(j in 1:length(mdls)){
        
        TG <- names(mdls)[j]
        
        if(file.exists(paste("output/projections/",TG,"_",ssp_list[i],"_baseline_",date_projections,".RData",sep=""))) {
          projection_baseline <- readRDS(file=paste("output/projections/",TG,"_",ssp_list[i],"_baseline_",date_projections,".RData",sep=""))
          
          projection_future <- readRDS(file=paste("output/projections/",TG,"_",ssp_list[i],"_future_",date_projections,".RData",sep=""))
          
          
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
          print(paste("File does not exist: output/projections/",TG,"_",ssp_list[i],"_baseline_",date_projections,".RData",sep=""))
        }
      }
    }
  }
  
  #6.3 to provide summary statistics for projected trophic groups' relative abundance
  
  summary_stats_TG <- function(mdls){
    
    ssp_list <- list("ssp126","ssp245","ssp370","ssp585")
    for(i in 1:length(mdls)){
      TG <- names(mdls)[i]
      
      for(j in 1:length(ssp_list)){
        
        projection_future <- readRDS(file=paste("output/projections/",TG,"_",ssp_list[j],"_future_",date,".RData",sep=""))
        
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
        
        projection_baseline <- readRDS(file=paste("output/projections/",TG,"_",ssp_list[j],"_baseline_",date,".RData",sep=""))
        
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

# 7_plot_modelsummary : to plot visual summary of glmmTMB model

    plot_model_summary_omnivores <- function(date){
      
      #01 read in df and model
      # df <- read_csv(paste0("data_input/global_df_complete_",date,".rds"))
      df_noNA <- df %>% filter(!is.na(chla_sqrt)) %>% filter(!is.na(ROC_SVT_zib)) %>% filter(!is.na(tow_days))

      #FIGURE 3
      print("Visual summary of glm for omnivorous zooplankton in preparation")
      #02 plot model predictions by (A) chl-a and (B) CPR Survey
      #Omnivores proportion vs. Chl-a
      pop_preds_omni_chla <- predictions(Omni_mdl_zib, 
                                         newdata = datagrid(chla_sqrt = seq(floor(min(df_noNA$chla_sqrt)),ceiling(max(df_noNA$chla_sqrt)),0.5)), 
                                         re.form = NA) # Zeros out random effects
      
      omni_plot_chla <- ggplot(data = pop_preds_omni_chla) + pub_theme +
        geom_point(data = df_noNA, 
                   aes(x = chla_sqrt, y = ROC_SVT_zib), alpha = 0.1) +
        geom_ribbon(aes(x = chla_sqrt, y = estimate, 
                        ymin = conf.low, ymax = conf.high), alpha = 0.3, colour = "blue", linetype=2) + 
        geom_line(aes(x = chla_sqrt, y = estimate),colour="blue") + 
        labs(y = "Proportion of omnivores", x = expression(bold(sqrt("Chl-a")))) 
      #annotate("text", y = 0.8, x = 3.8, label = expression("Estimate= 0.055; R"^{2}*"= 0.26; p-value=0.006 "), size = 5) 
      
      #Omnivores proportion vs. Surveys
      pop_preds_omni_surveys <- predictions(Omni_mdl_zib, 
                                            newdata = datagrid(survey = unique(df_noNA$survey), 
                                                               chla_sqrt = seq(floor(min(df_noNA$chla_sqrt)),ceiling(max(df_noNA$chla_sqrt)),0.5)), 
                                            re.form = NA) # Zeros out random effects
      
      omni_plot_surveys <- ggplot(data = pop_preds_omni_surveys) + pub_theme +
        geom_point(data = df_noNA, 
                   aes(x = chla_sqrt, y = ROC_SVT_zib), alpha = 0.1) +
        geom_ribbon(aes(x = chla_sqrt, y = estimate, 
                        ymin = conf.low, ymax = conf.high, 
                        colour = survey, fill = survey), alpha = 0.3, colour = NA) + 
        geom_line(aes(x = chla_sqrt, y = estimate, colour = survey), show.legend = F) + 
        labs(fill = "Survey", y = "Proportion of omnivores", x = expression(bold(sqrt("Chl-a")))) +
        theme(legend.position.inside = c(0.65, 0.25)) +
        scale_fill_discrete(labels=c('Australian CPR', 'Atlantic CPR','North Pacific CPR','SCAR Southern Ocean CPR'))
      
      #
      print("Plotting intercepts and slope of Longhurst Provinces and Tow within Survey")
      REs <- ranef(Omni_mdl_zib, condVar = TRUE) 
      TG <- "Omni"
      
      re_lh <- REs$cond$longhurst
      
      ### Longhurst Province ###
      qq <- attr(re_lh, "condVar")
      
      # Extract intercepts
      rand.interc <- REs$cond$longhurst
      
      # Make a dataframe for plotting
      df_plot <- data.frame(Intercepts = REs$cond$longhurst[,1],
                            sd.interc = 2*sqrt(qq[,,1:length(qq)]),
                            lev.names = factor(rownames(rand.interc))) %>% 
        arrange(Intercepts) %>% 
        within({  # Reorder levels
          lev.names <- factor(as.character(lev.names),
                              as.character(lev.names))
        })
      
      re_lh_plot <- ggplot(df_plot, aes(lev.names, Intercepts)) + pub_theme + 
        geom_hline(yintercept=0, linetype = "dashed", color = "black") +
        geom_errorbar(aes(ymin=Intercepts-sd.interc, 
                          ymax=Intercepts+sd.interc),
                      width = 0,color="black") +
        geom_point(color = "black", size = 2) +
        guides(size=FALSE,shape=FALSE) + 
        theme(axis.text.x=element_text(size=10), 
              axis.title.x=element_text(size=13),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.x = element_blank()) +
        coord_flip() + 
        labs(y = "Intercept", x = "Longhurst Provinces")
      
      #annotate("text", y = -2.3, x = 28, label = "(D)", size = 5)
      
      ### Tow slope and intercept ###
      
      # Create a dataframe for plotting
      re_tow <- REs$cond$`survey:tow_no`
      
      towDf <- as.data.frame(re_tow) %>% setNames(c("Intercept", "Slope"))
      
      # Specify frequency breaks
      my_breaks <- c(1,5,30,175,1000)
      
      #
      re_tow_plot_point_density <- ggplot(towDf, aes(Intercept, Slope)) +
        pub_theme +
        geom_pointdensity() +
        scale_colour_viridis_c() +
        labs(colour = "Density", x = "Tow within survey intercept", y = "Tow within survey slope") +
        theme(
          axis.text = element_text(size = 14), legend.position.inside = c(0.85, 0.85)    # Adjust the size for axis numbers/text
        ) 
      
      #04 Combine plots per trophic group. 
      design <- "
                    AABB
                    AABB
                    CCDD
                    CCDD
                    CCDD
                    CCDD
                  "
      
      #Figure 3. Omnivores
      print("Patching plots together")
      (final_patch <- omni_plot_chla + omni_plot_surveys + re_lh_plot + re_tow_plot_point_density +
          plot_layout(design = design) +
          plot_annotation(
            tag_levels = "A",
            theme = theme(plot.title = element_text(size = 12, face = "bold"))
          )
      )    
      print(paste("Saved plot: output/plots/Omni_",date,".png",sep=""))
      ggsave(paste("output/plots/Omni_",date,".png",sep=""), plot = final_patch,
             width = 8, height = 10, dpi = 300)
      
    }
    
    
    plot_model_summary_carnivores <- function(date){
      # Create a publication-ready theme (Adapted from 2025 UQ MME Lab Winter R Workshop)
      # #01 read in model
      # load("output/previousModels/revision/fMdl_carni.RData") 
      # df <- read_csv(paste0("data_input/global_df_complete_",date,".rds"))
      df_noNA <- df %>% filter(!is.na(chla_sqrt)) %>% filter(!is.na(RCO_SVT_zib)) %>% filter(!is.na(tow_days))

      
      #FIGURE 4
      print("Visual summary of glm for carnivorous zooplankton in preparation")
      #B. Carnivores
      #Carnivores proportion vs. Chl-a
      
      pop_preds_carni_chla <- predictions(Carni_mdl_zib, 
                                          newdata = datagrid(chla_sqrt = seq(floor(min(df_noNA$chla_sqrt)),ceiling(max(df_noNA$chla_sqrt)),0.05)), 
                                          re.form = NA) # Zeros out random effects
      
      carni_plot_chla <- ggplot(data = pop_preds_carni_chla) + pub_theme + 
        geom_point(data = df_noNA, 
                   aes(x = chla_sqrt, y = RCO_SVT_zib), alpha = 0.1) +
        geom_ribbon(aes(x = chla_sqrt, y = estimate, 
                        ymin = conf.low, ymax = conf.high), alpha = 0.3, colour = "blue", linetype = 2) + 
        geom_line(aes(x = chla_sqrt, y = estimate),colour="blue") + 
        labs(y = "Proportion of carnivores", x = expression(bold(sqrt("Chl-a")))) 
      #annotate("text", y = 0.6, x = 3.8, label = expression("Estimate= -0.22; R"^{2}*"= 0.20; p-value<1e-21 "), size = 5) 
      
      #Carnivores proportion vs. Surveys
      pop_preds_carni_surveys <- predictions(Carni_mdl_zib, 
                                             newdata = datagrid(survey = unique(df_noNA$survey), 
                                                                chla_sqrt = seq(floor(min(df_noNA$chla_sqrt)),ceiling(max(df_noNA$chla_sqrt)),0.05)), 
                                             re.form = NA) # Zeros out random effects
      
      carni_plot_surveys <- ggplot(data = pop_preds_carni_surveys) + pub_theme +
        geom_point(data = df_noNA, 
                   aes(x = chla_sqrt, y = RCO_SVT_zib), alpha = 0.1) +
        geom_ribbon(aes(x = chla_sqrt, y = estimate, 
                        ymin = conf.low, ymax = conf.high, 
                        colour = survey, fill = survey), alpha = 0.3, colour = NA, linetype = 2) + 
        geom_line(aes(x = chla_sqrt, y = estimate, colour = survey), show.legend = F) + 
        labs(fill = "Survey", y = "Proportion of carnivores", x = expression(bold(sqrt("Chl-a")))) +
        theme(legend.position.inside = c(0.65, 0.85)) +
        scale_fill_discrete(labels=c('Australian CPR', 'Atlantic CPR','North Pacific CPR','SCAR Southern Ocean CPR'))
      
      
      REs <- ranef(Carni_mdl_zib, condVar = TRUE) 
      TG <- "Carni"
      
      re_lh <- REs$cond$longhurst
      
      ### Longhurst Province ###
      qq <- attr(re_lh, "condVar")
      
      # Extract intercepts
      rand.interc <- REs$cond$longhurst
      
      # Make a dataframe for plotting
      df_plot <- data.frame(Intercepts = REs$cond$longhurst[,1],
                            sd.interc = 2*sqrt(qq[,,1:length(qq)]),
                            lev.names = factor(rownames(rand.interc))) %>% 
        arrange(Intercepts) %>% 
        within({  # Reorder levels
          lev.names <- factor(as.character(lev.names),
                              as.character(lev.names))
        })
      
      re_lh_plot <- ggplot(df_plot, aes(lev.names, Intercepts)) + pub_theme + 
        geom_hline(yintercept=0, linetype = "dashed", color = "black") +
        geom_errorbar(aes(ymin=Intercepts-sd.interc, 
                          ymax=Intercepts+sd.interc),
                      width = 0,color="black") +
        geom_point(color = "black", size = 2) +
        guides(size=FALSE,shape=FALSE) + 
        theme(axis.text.x=element_text(size=10), 
              axis.title.x=element_text(size=13),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.x = element_blank()) +
        coord_flip() + 
        labs(y = "Intercept", x = "Longhurst Provinces")
      
      #annotate("text", y = -2.3, x = 28, label = "(D)", size = 5)
      
      ### Tow slope and intercept ###
      print("Plotting intercepts and slope of Longhurst Provinces and Tow within Survey")
      # Create a dataframe for plotting
      re_tow <- REs$cond$`survey:tow_no`
      
      towDf <- as.data.frame(re_tow) %>% setNames(c("Intercept", "Slope"))
      
      # Specify frequency breaks
      my_breaks <- c(1,5,30,175,1000)
      
      #
      re_tow_plot_point_density <- ggplot(towDf, aes(Intercept, Slope)) +
        pub_theme +
        geom_pointdensity() +
        scale_colour_viridis_c() +
        labs(colour = "Density", x = "Tow within survey intercept", y = "Tow within survey slope") +
        theme(
          axis.text = element_text(size = 14), legend.position.inside = c(0.85, 0.85)    # Adjust the size for axis numbers/text
        ) 
      
      #04 Combine plots per trophic group. 
      design <- "
                    AABB
                    AABB
                    CCDD
                    CCDD
                    CCDD
                    CCDD
                  "
      
      #Figure 4. Carnivores
      print("Patching plots together")
      (final_patch <- carni_plot_chla + carni_plot_surveys + re_lh_plot + re_tow_plot_point_density +
          plot_layout(design = design) +
          plot_annotation(
            tag_levels = "A",
            theme = theme(plot.title = element_text(size = 12, face = "bold"))
          )
      )    
      
      print(paste("Plot saved: output/plots/Carni_",date,".png",sep=""))
      ggsave(paste("output/plots/Carni_",date,".png",sep=""), plot = final_patch,
             width = 8, height = 10, dpi = 300)
      
    }
    
    
    plot_model_summary_filterfeeders <- function(date){
      # #01 read in model
      # load("output/previousModels/revision/fMdl_filterfeeder.RData") 
      # df_filter <- read_csv(paste0("data_input/global_df_complete_Filter_",date,".rds"))
      df_noNA <- df_filter %>% filter(!is.na(chla_sqrt)) %>% filter(!is.na(RFF_SVT_zib)) %>% filter(!is.na(tow_days))

      #FIGURE 5
      print("Visual summary of glm for gelatinous filter-feeders in preparation")
      #C. Filter-feeders
      #Filter-feeders proportion vs. Chl-a
      pop_preds_FF_chla <- predictions(Filter_mdl_zib, 
                                       newdata = datagrid(chla_sqrt = seq(floor(min(df_noNA$chla_sqrt)),ceiling(max(df_noNA$chla_sqrt)),0.05)), 
                                       re.form = NA) # Zeros out random effects
      
      ff_plot_chla <- ggplot(data = pop_preds_FF_chla) + pub_theme + 
        geom_point(data = df_noNA, 
                   aes(x = chla_sqrt, y = RFF_SVT_zib), alpha = 0.1) +
        geom_ribbon(aes(x = chla_sqrt, y = estimate, 
                        ymin = conf.low, ymax = conf.high), alpha = 0.3, colour = "blue", linetype=2) + 
        geom_line(aes(x = chla_sqrt, y = estimate),colour="blue") + 
        labs(y = "Proportion of filter-feeders", x = expression(bold(sqrt("Chl-a")))) 
      #annotate("text", y = 0.95, x = 3.5, label = expression("Estimate= -0.49; R"^{2}*"= 0.41; p-value<1e-15 "), size = 5)  
      
      #Filter-feeders proportion vs. Surveys
      pop_preds_FF_surveys <- predictions(Filter_mdl_zib, 
                                          newdata = datagrid(survey = unique(df_noNA$survey), 
                                                             chla_sqrt = seq(floor(min(df_noNA$chla_sqrt)),ceiling(max(df_noNA$chla_sqrt)),0.05)), 
                                          re.form = NA) # Zeros out random effects
      
      ff_plot_surveys <- ggplot(data = pop_preds_FF_surveys) + pub_theme +
        geom_point(data = df_noNA, 
                   aes(x = chla_sqrt, y = RFF_SVT_zib), alpha = 0.1) +
        geom_ribbon(aes(x = chla_sqrt, y = estimate, 
                        ymin = conf.low, ymax = conf.high, 
                        colour = survey, fill = survey), alpha = 0.3, colour = NA) + 
        geom_line(aes(x = chla_sqrt, y = estimate, colour = survey), show.legend = F) + 
        labs(fill = "Survey", y = "Proportion of filter-feeders", x = expression(bold(sqrt("Chl-a")))) +
        theme(legend.position.inside = c(0.65, 0.85)) +
        scale_fill_discrete(labels=c('Australian CPR', 'Atlantic CPR','North Pacific CPR','SCAR Southern Ocean CPR'))
      
      #03 plot intercepts and slope of (A) Longhurst Provinces and (B) Tow No. and Days
      print("Plotting intercepts and slope of Longhurst Provinces and Tow within Survey")
      ###Longhurst Province
      # Extract random effects
      REs <- ranef(Filter_mdl_zib, condVar = TRUE) #Update the model input
      TG <- "Filter"
      
      
      
      re_lh <- REs$cond$longhurst
      
      ### Longhurst Province ###
      qq <- attr(re_lh, "condVar")
      
      # Extract intercepts
      rand.interc <- REs$cond$longhurst
      
      # Make a dataframe for plotting
      df_plot <- data.frame(Intercepts = REs$cond$longhurst[,1],
                            sd.interc = 2*sqrt(qq[,,1:length(qq)]),
                            lev.names = factor(rownames(rand.interc))) %>% 
        arrange(Intercepts) %>% 
        within({  # Reorder levels
          lev.names <- factor(as.character(lev.names),
                              as.character(lev.names))
        })
      
      re_lh_plot <- ggplot(df_plot, aes(lev.names, Intercepts)) + pub_theme + 
        geom_hline(yintercept=0, linetype = "dashed", color = "black") +
        geom_errorbar(aes(ymin=Intercepts-sd.interc, 
                          ymax=Intercepts+sd.interc),
                      width = 0,color="black") +
        geom_point(color = "black", size = 2) +
        guides(size=FALSE,shape=FALSE) + 
        theme(axis.text.x=element_text(size=10), 
              axis.title.x=element_text(size=13),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.x = element_blank()) +
        coord_flip() + 
        labs(y = "Intercept", x = "Longhurst Provinces")
      
      #annotate("text", y = -2.3, x = 28, label = "(D)", size = 5)
      
      ### Tow slope and intercept ###
      
      # Create a dataframe for plotting
      re_tow <- REs$cond$`survey:tow_no`
      
      towDf <- as.data.frame(re_tow) %>% setNames(c("Intercept", "Slope"))
      
      # Specify frequency breaks
      my_breaks <- c(1,5,30,175,1000)
      
      #
      re_tow_plot_point_density <- ggplot(towDf, aes(Intercept, Slope)) +
        pub_theme +
        geom_pointdensity() +
        scale_colour_viridis_c() +
        labs(colour = "Density", x = "Tow within survey intercept", y = "Tow within survey slope") +
        theme(
          axis.text = element_text(size = 14), legend.position.inside = c(0.85, 0.85)    # Adjust the size for axis numbers/text
        ) 
      
      #04 Combine plots per trophic group. 
      design <- "
                    AABB
                    AABB
                    CCDD
                    CCDD
                    CCDD
                    CCDD
                  "
      
      #Figure 5. Filter-feeders
      print("Patching plots together")
      
      (final_patch <- ff_plot_chla + ff_plot_surveys + re_lh_plot + re_tow_plot_point_density +
          plot_layout(design = design) +
          plot_annotation(
            tag_levels = "A",
            theme = theme(plot.title = element_text(size = 12, face = "bold"))
          )
      )    
      
      print(paste("Plot saved: output/plots/Filter_",date,".png",sep=""))
      ggsave(paste("output/plots/Filter_",date,".png",sep=""), plot = final_patch,
             width = 8, height = 10, dpi = 300)
      
    }  