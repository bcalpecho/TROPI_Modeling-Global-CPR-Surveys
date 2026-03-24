# ---
# title: Wrangling CPR data
# author: Bryan Alpecho
# date: 2025-
# output: rds or csv files
# ---

# Aims
#These are functions to help in pre-processing CPR data for the data analyses
# 01. to map global cpr
# 02. to update global cpr taxon list
# 03. to update the column names of cpr tables from taxon name into aphia id
# 04. to separate mba data into north pacific and north atlantic data
# 05. to compute abundances per taxon

# setup
library(tidyverse)
library(stars)
library(tmap)
library(rnaturalearth)

#Current iterations of CPR data (DOI or latest date)
#auscpr - coverage 05/10/2024 
#mba (npacific + natlantic) - 10.17031/685d5b713a41b
#socpr - 10.26179/r7sm-9y85

# set date for version control of output
date <- "09112025" #format DDMMYYYY

#01 to map global CPR

cpr_meta <- list.files(path = "data_input/CPR_on-process/", pattern = "*\\metadata.csv", full.names = TRUE)

map_globalcpr <- function(file.list){
  #extract metadata  
  filenames <- basename(file.list)
  cpr <- data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("survey", "sample_id", "latitude", "longitude", "sampleTime_utc"))))
  cpr[,1:2] <- lapply(cpr[,1:2], as.character) # survey and sample id
  cpr[,3:4] <- lapply(cpr[,3:4], as.double) # latitude and longitude
  cpr[,5] <- as.POSIXct(cpr[,5]) #sample time UTC
  
  
  #loop through each survey
    for(i in 1:length(filenames)){
    file_survey <- str_extract(filenames[i], "(?<=_)[^_]+")
    print(paste("Survey: ",file_survey, sep=""))
    
    cpr_metadata <- read_csv(file.list[i]) %>% 
      mutate(survey = file_survey) %>% 
      select("survey","sample_id","latitude","longitude","sampleTime_utc")
    
    print(paste("File No.", i, sep=""))
    #integrate into a global cpr
    cpr <- cpr %>% 
      rows_insert(cpr_metadata, by = "survey") 
    }
    
    #identify coordinates 
    cpr <- cpr %>% 
      st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
      rename("Survey" = "survey")
    
    #set projection 
    world_projection <- '+proj=eqearth +lon_0=0 +datum=WGS84 +units=m +no_defs'
    #get world coastline
    world <- ne_coastline(scale = "medium")
    
    #to plot
    globalcpr_map <- tm_shape(world) +
      tm_polygons(col = "white") +
      tm_shape(cpr, crs = world_projection, raster.warp = TRUE) +
      tm_dots(col = "Survey", fill_alpha = 0.5, size = 0.3,
              labels = c("Australian CPR", "Atlantic CPR", "North Pacific CPR", "SCAR Southern Ocean CPR")) + 
      tm_graticules(alpha = 0.5, 
                    x = c(-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180), 
                    y = c(-90, -60, -30, 0, 30, 60, 90), 
                    labels.size = 1) +
      tm_legend(position = c("left", "center"),
                text.size = 1.2, title.size = 1.4)
    
    #to export plot
    tmap_save(globalcpr_map, filename=paste("Output/map/Global/Global-CPR-map_",date,".png",sep=""),
              width = 400,
              height = 200,
              units = "mm",
              dpi = 400)
}
map_globalcpr(cpr_meta)

#01.B to map survey-specific CPR
#Derivative of global cpr into survey-specific CPR

    #import cpr data
    mba <- read_csv("data_input/CPR_on-process/")
    mba_npacific <- read_csv("data_input/CPR_on-process/cpr_npacific_metadata.csv") %>% 
      mutate(Survey = "North Pacific CPR") %>% 
      select("Survey","Sample_ID","latitude","longitude") 
    
    mba_natlantic <- read_csv("data_input/CPR_on-process/cpr_natlantic_metadata.csv") %>% 
      mutate(Survey = "Atlantic CPR") %>% #update North Atlantic to Atlantic 05.11.2025
      select("Survey","sample_id","latitude","longitude") 
    
    auscpr <- read_csv("data_input/CPR_on-process/cpr_auscpr_metadata.csv") %>% 
      mutate(Survey = "Australian CPR") %>% 
      select("Survey","Sample_ID","Latitude","Longitude") %>% 
      rename(c(sample_id = Sample_ID,latitude = Latitude, longitude = Longitude))
    
    socpr <- read_csv("data_input/CPR_on-process/cpr_socpr_metadata.csv") %>% 
      mutate(Survey = "SCAR Southern Ocean CPR") %>% 
      select("Survey","Sample_ID","Latitude","Longitude") %>% 
      rename(c(sample_id = Sample_ID,latitude = Latitude, longitude = Longitude))
    
    #plot to double check
    #Global
    mba_cpr <-  mba_npacific %>%
      rows_insert(mba_natlantic, by="Survey") %>% 
      rows_insert(auscpr, by="Survey") %>% 
      rows_insert(socpr, by="Survey") %>% 
      #select(Latitude, Longitude, SampleTime_UTC) %>% 
      st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
    
    plot(mba_cpr["geometry"], axes = TRUE)
    world_projection <- '+proj=eqearth +lon_0=0 +datum=WGS84 +units=m +no_defs'
    world <- ne_coastline(scale = "medium")
    
    globalcpr_map <- tm_shape(world) + 
      tm_polygons(col = "white") +
      tm_shape(mba_cpr, projection = world_projection, raster.warp = TRUE) +
      tm_dots(col = "Survey", color = "Survey",
              title = "CPR Survey", alpha = 0.5, icon.scale = 3, labels.size = 20) + 
      tm_graticules(alpha = 0.5, 
                    x = c(-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180), 
                    y = c(-90, -60, -30, 0, 30, 60, 90), 
                    labels.size = 1) 
      #tm_compass(type = "8star", position = c("left", "center")) 
    
    tmap_save(globalcpr_map, filename=paste("Output/map/Global/Global-CPR-map_",date,".png",sep=""),
              width = 400,
              height = 200,
              units = "mm",
              dpi = 216)


    
#02. to update species list
  #set date for version control
  old_version <- "02062025"
  date <- "26072025"
    #old list
    taxon_list_old <- read_csv(paste("data_input/CPR_on-process/cpr_all_aphia-taxon_",date,".csv", sep=""))
    #new list
    taxon_list_new <- read_csv(paste("data_input/CPR_on-process/cpr_all_aphia-taxon_",date,".csv", sep="")) #update file name
    #combine
    taxon_list <- taxon_list_new %>% left_join(taxon_list_old, by = ("aphiaID"))
    #add new to old list
    taxon_list <- taxon_list %>% mutate(taxon = case_when())
    #export
    write_csv(taxon_list, paste("data_input/CPR_on-process/cpr_all_aphia-taxon_",date,".csv", sep=""))
    
#03. to update the column names from taxon name into aphia id
    #get taxon list with aphiaid
    #taxon_list <- read_csv("data_input/CPR_on-process/cpr_all_aphia-taxon_26072025.csv") 
    #taxon_list$aphiaID <- as.factor(species.list$aphiaID)
    mba_taxonlist <- read_csv("data_input/CPR_on-process/cpr_mba_taxonlist.csv")
    
    #update column names to aphia id 
    #for mba data
    mba <- read_csv("data_input/CPR_on-process/cpr_mba_complete.csv")
    mba_meta <- read_csv("data_input/CPR_on-process/cpr_mba_metadata.csv")
    ## remove samples prior to september 01 1997 (start date of OC-CCI coverage)
    mba <- mba %>% filter(date(midpoint_date_gmt) >= "1997-09-01") 
    
    #update column names to aphia id
    names(mba)[match(mba_taxonlist$taxa_name, names(mba))] <- mba_meta$aphia_id
    
    #combine columns of the same name
    mba_long <- mba[,c(2,6:559)] %>%
      pivot_longer(cols = -sample_id, names_to = "AphiaID", values_to = "value")
    
    DFsum <- mba_long %>% group_by(sample_id, AphiaID) %>% summarize(Total=sum(value))
    
    DFfinal <- DFsum %>% pivot_wider(names_from = "AphiaID",values_from="Total")
    
    cpr_dat <- DFfinal %>% ungroup()
    
    #combine columns of same aphiaIDs
    # DFlong <- mba %>% pivot_longer(cols = -sample_id,names_to="aphiaID",values_to="value")
    # DFsum <- DFlong %>% group_by(Sample_ID,aphiaID) %>% summarize(Total=sum(value))
    # DFfinal <- DFsum %>% pivot_wider(names_from = "aphiaID",values_from="Total")
    # cpr_dat <- DFfinal %>% ungroup() %>% select(-c("Sample_ID"))
    
    write_csv(cpr_dat, "data_input/CPR_on-process/cpr_mba_complete.csv")
    
    #return the metadata
    mba_final <- DFfinal %>% 
      left_join(mba_meta, by="sample_id")
    
    
#04. to separate mba data into north pacific and north atlantic data
    #for mba data
    mba <- read_csv("data_input/CPR_on-process/mba/cpr_mba_complete.csv")
    mba_meta <- read_csv("data_input/CPR_on-process/mba/cpr_mba_metadata.csv") %>% 
      mutate(sampleTime_utc = as.POSIXct(midpoint_date_gmt, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
    mba_final <- mba %>% 
      left_join(mba_meta, by = "sample_id")
    
    #to remove taxa for avoidance of double counting (Richardson et al., 2006)
    
    #identify cpr belonging to North Atlantic and North Pacific CPR
    mba_npacific <- mba_final %>% 
      filter(longitude >= 100 | longitude <= -100) %>%
      mutate(survey = "North Pacific CPR") 
    
    mba_natlantic <- mba_final %>% 
      filter(longitude < 100 & longitude > -100) %>% 
      mutate(survey = "North Atlantic CPR") 
      #North Atlantic CPR: remove freshwater (Lake Tanganyika)
      mba_natlantic <- mba_natlantic %>% 
        filter(!sample_id %in% c("2ALT-1","2ALT-5","3ALT-2","3ALT-5"))
      
    #get metadata
    mba_npacific_metadata <- mba_npacific %>% 
      select(c("survey","latitude","longitude","sample_id","sampleTime_utc")) 
    
    mba_natlantic_metadata <- mba_natlantic %>% 
      select(c("survey","latitude","longitude","sample_id","sampleTime_utc"))
    
    #export 
    write_csv(mba_natlantic %>% select(-c(survey, midpoint_date_gmt, latitude, longitude, chlorophyll_index, sampleTime_utc)), "data_input/CPR_on-process/cpr_natlantic_complete.csv")
    write_csv(mba_npacific %>% select(-c(survey, midpoint_date_gmt, latitude, longitude, chlorophyll_index, sampleTime_utc)), "data_input/CPR_on-process/cpr_npacific_complete.csv")
    
    write_csv(mba_natlantic_metadata, "data_input/CPR_on-process/cpr_natlantic_metadata.csv")
    write_csv(mba_npacific_metadata, "data_input/CPR_on-process/cpr_npacific_metadata.csv")
    
#05. to compute abundances per taxon  

    #zoop_list is any list with aphiaID
    #00 get abundance list 
    file.list <- list.files(path = "data_input/CPR_on-process/", pattern = "*\\complete.csv", full.names = TRUE)
    #set date for version control
    date <- "20082025"
    
    #apply function compute_abundance to file.list
    compute_abundance <- function(zoop_list){
      cpr_aphia_list <- read_csv("data_input/traits/TG_trait-table_20082025.csv") %>% select(c(aphiaID, scientificName)) %>%  rename(taxon = scientificName)
      cpr_aphia_list$aphiaID <- as.factor(cpr_aphia_list$aphiaID)
      survey <- c("auscpr","natlantic","npacific","socpr")
      
      for(i in 1:length(survey)){
        #loop progress
        print(survey[i])
        
        #read in cpr abundance data
        cpr <- read_csv(zoop_list[i], col_names = T, name_repair = "minimal")
        #select for abundances (columns with aphiaID)
        cpr_dat <- cpr %>%  select(-c(sample_id)) #double check the metadata in each csv file
        
        # #OPTIONAL: remove CPR samples pre-dating 1997 September
        # cpr_dat <- cpr_dat %>% 
        #   filter(Year>=1998)
        # natlantic_dat_1997 <- cpr %>% 
        #   filter(Year == 1997 & Month >= 9)
        
        cpr_sums <- colSums(cpr_dat, na.rm=T) %>% 
          data.frame() %>% rownames_to_column(var = "aphiaID") %>% rename("totalAbundance"=".") %>% 
          left_join(cpr_aphia_list, by="aphiaID") %>% 
          group_by(aphiaID, taxon) %>%
          summarise(totalAbundance = sum(totalAbundance)) %>%
          ungroup() %>%
          mutate(relativeAbundance = totalAbundance/sum(totalAbundance))%>%
          select(-c("taxon"))
        
        write_rds(cpr_sums, paste("Output/data/abundance/",survey[i],"_",date,".rds",sep=""))
      }
      
      traits_df <- read_csv("data_input/traits/TG_trait-table_20082025.csv")
      traits_df$aphiaID <- as.factor(traits_df$aphiaID)
      
      #combine each survey-specific abundance to the trait table
      for(i in 1:length(survey)){
        survey_abundance <- read_rds(paste("Output/data/abundance/",survey[i],"_",date,".rds",sep=""))
        survey_abundance <- survey_abundance %>% 
          rename_with(~c("aphiaID", paste(survey[i],"_totalAbundance",sep=""),paste(survey[i],"_relativeAbundance",sep="")))
        traits_df <- traits_df %>% left_join(survey_abundance, by="aphiaID")
      }
      
      #compute for the global abundances in the trait table
      traits_df <- traits_df %>%
        rowwise() %>%
        mutate(global_totalAbundance = sum(natlantic_totalAbundance, auscpr_totalAbundance, npacific_totalAbundance, socpr_totalAbundance, na.rm=T))  %>%
        ungroup()
      traits_df <- traits_df %>%
        mutate(global_relativeAbundance = global_totalAbundance/sum(global_totalAbundance)) %>%
        relocate(c("global_totalAbundance", "global_relativeAbundance"),.before = "auscpr_totalAbundance")
      
      write_csv(traits_df, paste("Output/traits/TG_trait-table_withAbundances_",date,".csv",sep=""))
    }  
    
    compute_abundance(file.list)  
    
#Latest revision 06.08.2025