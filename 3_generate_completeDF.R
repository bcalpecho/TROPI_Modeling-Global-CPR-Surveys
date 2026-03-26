# ---
# title: generate_CompleteDF
# author: Bryan Alpecho
# date: 2025-
# output: csv file
# ---

#Aims
##1. Compute for proportions of trophic groups                                     
##2. Generate values of random effects and combine variables altogether into a dataframe
##3. Update df with proportions of (3A) trophic groups and (3B) chla 
##4. generate complete DF (global) inc. compute for ratios

# Load libraries
library(tidyverse)
library(sf)
library(stars)

#set date for version control
date <- "21032026"

##1. Compute for proportions of trophic groups
    #00 get abundance list 
    file.list <- list.files(path = "data_input/CPR_on-process/", pattern = "*\\complete.csv", full.names = TRUE)
    
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
          
          # print(paste("Output/data/sums/df_sums_",survey[i],"_",date, sep=""))
          # write_csv(trait_sums, paste("Output/data/sums/df_sums_",survey[i],"_",date, sep=""))
          print(paste("Saved file: ",survey,"/df_sums_",date,".csv", sep=""))
          write_csv(trait_sums, paste("Output/data/",survey,"/df_sums_",date,".csv", sep=""))
      }
    }
    compute_proportions_perSurvey(file.list)
    

##2. Combining variables altogether into a dataframe
    #set dates
 
    generate_df_perSurvey <- function(){
      chl_date <- "15092025"
      trait_date <- "21032026"
      
      metadata_files <- list.files(path = "data_input/CPR_on-process/", pattern = "*\\metadata.csv", full.names = TRUE)
      
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
             chl <- read_rds(paste("Output/data/",file_survey[i],"/chla_merged_",chl_date,".rds",sep="")) %>% select("sample_id","chla")
             df <- df %>% 
               left_join(chl, by = "sample_id")
             
          #06 Traits (TG)
             traits <- read_csv(paste("Output/data/",file_survey[i],"/df_sums_",trait_date,".csv",sep=""), col_names =T, show_col_types =  F) 
             df <- df %>% 
               left_join(traits, by = "sample_id")
             df$geometry = NULL
             
          #07 export df
             
             print(paste("Filename output: ",file_survey[i],"/df_complete_",date,".csv", sep=""))
             write_rds(df, paste("Output/data/",file_survey[i],"/df_complete_",date,".rds",sep=""))
        }
          #return(df)
    }
    generate_df_perSurvey()
   
    
##3A. Update SURVEY-SPECIFIC dataframe with proportions of trophic groups
    #NOTE: file list not yet in sync with compute_sums function
    #set date for version of updated proportions
    date <- "01072025"
    
    #
    update_df_traits <- function(date_){
        survey <- c("auscpr","natlantic","npacific","socpr")
        previousDF_list <- NA
        for(i in 1:length(survey)){
          previousDF_list[i] <- list.files(path = paste("Output/data/",survey[i],"/",sep=""), pattern = "*\\complete_04082025.rds", full.names = TRUE)
        }
        remove(i)
      
      for(i in 1:length(previousDF_list)){
        #get metadata
        filenames <- basename(previousDF_list[i])
        survey <- str_extract(previousDF_list[i], "(?<=/data/)[^/]+") 
        date_oldVersion <- str_extract(filenames, "[^_]+$")
        
        #progress
        print(paste("Update DF of ",survey," from version ",date_oldVersion," to ",date,sep=""))
        
        #read in previous df
        df_prev <- paste("Output/data/",survey,"/",filenames,sep="")
        df <- read_rds(df_prev)
        
        #remove previous traits
        df <- df %>% select(c(sample_id, longitude, latitude, sampleTime_utc, tow_no, tow_days, longhurst, chla))
        
        #load updated traits
        traits_file <- paste("Output/data/",survey,"/df_sums_",date,".csv", sep="")
        traits <- read_csv(traits_file)
        
        #update df
        df_new <- df %>% 
          left_join(traits, by = "sample_id")
        
        #save updated df
        file_output <- paste("Output/data/",survey,"/df_complete_",date,".rds",sep="")
        print(file_output)
        write_rds(df_new, file_output)
      } 
    }

  update_df_traits(file.list)

  
##3B. Update SURVEY-SPECIFIC dataframe with chl-a 
  #set dates for version of chl-a and previous df
  date <- "16092025"
  date_previousDF <- "16092025"
  
  #
  survey <- c("auscpr","natlantic","npacific","socpr")
  file.list <- NA
  for(i in 1:length(survey)){
    file.list[i] <- list.files(path = paste("Output/data/",survey[i],"/",sep=""), pattern = paste("*\\complete_",date_previousDF,".rds",sep=""), full.names = TRUE)
  }
  remove(i)
  
  update_df_chla <- function(previousDF_list){
    for(i in 1:length(previousDF_list)){
      #get metadata
      filenames <- basename(previousDF_list[i])
      survey <- str_extract(previousDF_list[i], "(?<=/data/)[^/]+") 
      date_oldVersion <- str_extract(filenames, "[^_]+$")
      
      #progress
      print(paste("Update DF of ",survey," from version ",date_oldVersion," to ",date,sep=""))
      
      #read in previous df
      df_prev <- paste("Output/data/",survey,"/",filenames,sep="")
      df <- read_rds(df_prev)
      
      #remove previous chla
      df <- df %>% dplyr::select(-chla)
      
      #load updated chla
        #CHECK the file name
        chla_file <- paste("Output/data/",survey,"/chla_monthly_",date,".rds", sep="")
        chla <- read_rds(chla_file) %>% select(c("sample_id","chla"))
        
      #update df
      df_new <- df %>% 
        left_join(chla, by = "sample_id")
      
      #save updated df
        #CHECK the file name
        file_output <- paste("Output/data/",survey,"/df_complete_",date,"_monthlychla.rds",sep="")
        print(file_output)
        write_rds(df_new, file_output)
    } 
  }
  update_df_chla(file.list)
  

##4. generate GLOBAL dataframe and compute for ratios
    #set date for version of updated SURVEY-SPECIFIC dataframes
    date_newVersion <- "21032026"
    
    generate_globalCPR_dataframe <- function(date_newVersion){
      survey_list <- c("auscpr","natlantic","npacific","socpr")
      for(i in 1:length(survey_list)){
        print(paste0("Reading dataframe of: ",survey_list[i]))
        assign(paste(survey_list[i]), read_rds(paste("Output/data/",survey_list[i],"/df_complete_",date_newVersion,".rds",sep="")) %>%  
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
        # mutate(RFF_SVT = case_when(
        #   RFF == 0 | RFF == 1 ~ ((RFF*325075)+0.5)/325076,
        #   .default = RFF)) %>%
        # mutate(ROC_SVT = case_when(
        #   ROC == 0 | ROC == 1 ~ ((ROC*325075)+0.5)/325076,
        #   .default = ROC)) %>%
        # mutate(RCO_SVT = case_when(
        #   RCO == 0 | RCO == 1 ~ ((RCO*325075)+0.5)/325076,
        #   .default = RCO)) %>%
        # mutate(ROC_Cal_SVT = case_when(
        #   ROC_Cal == 0 | ROC_Cal == 1 ~ ((ROC_Cal*325075)+0.5)/325076,
        #   .default = ROC_Cal)) %>%
        # mutate(RCO_Cal_SVT = case_when(
        #   RCO_Cal == 0 | RCO_Cal == 1 ~ ((RCO_Cal*325075)+0.5)/325076,
        #   .default = RCO_Cal)) %>%
        # mutate(ROC_Cyc_SVT = case_when(
        #   ROC_Cyc == 0 | ROC_Cyc == 1 ~ ((ROC_Cyc*325075)+0.5)/325076,
        #   .default = ROC_Cyc)) %>%
        # mutate(RCO_Cyc_SVT = case_when(
        #   RCO_Cyc == 0 | RCO_Cyc == 1 ~ ((RCO_Cyc*325075)+0.5)/325076,
        #   .default = RCO_Cyc)) %>%
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
      
      #finalize df without Group
      df <- df %>% 
        #relocate(survey, .before= sample_id) %>% 
        relocate(chla_sqrt, .after=chla)
      
      print(paste0("File saved: df_complete_",date_newVersion,".rds"))
      write_rds(df, paste0("Output/data/global/df_complete_",date_newVersion,".rds"))   
      write_rds(df, paste0("data_input/global/df_complete_",date_newVersion,".rds")) 
      
    }

    generate_globalCPR_dataframe(survey_list)
    
    #Combining the survey-specific CPR datasets (selected parameters; no "other copepods" for 04.08.2025)
    #CHECK file name of updated survey-specific df
  
#latest revision:
#30.06.2025  sync functions; generate_df not yet working
#04.08.2025 generate_df function working
#16.09.2025 update_df_chla working
#21.03.2026 generate_globalCPR_dataframe working    