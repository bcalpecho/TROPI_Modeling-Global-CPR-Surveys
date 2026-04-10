# ---
# title: Helper functions for handling CPR data
# author: Bryan Alpecho
# date: 2025-
# output: rds or csv files
# ---

#These are functions to help in pre-processing CPR data for generating the final dataframe ('3_generate_completeDF' script) 
# 01. to pre-process raw CPR csv files (Identify AphiaID) 
# 02. to compute proportions of trophic groups
# 02. to map global CPR
# 03. to update global cpr taxon list
# 04. to compute abundances per taxon

# load libraries
#Load libraries
packages <- c("tidyverse",
              "stars",
              "sf",
              "tmap",
              "rnaturalearth",
              "fs")

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

#Current iterations of CPR data (latest as of 27 March 2026)
#auscpr - downloaded 27 March 2026 from the Australian Ocean Data Network
#mba (npacific + natlantic) - latest as of 2025 from the Archive for Marine Species and Habitats Data
#socpr - latest as of 03 Sept 2025 from the Australian Antarctic Data Centre

# # set date for version control of output
# date <- "27032026" #format DDMMYYYY

#01 pre-processing CPR raw files 
  #pre-process Australian CPR survey
    #load auscpr raw csv file (based on latest available df as of 27 March 2026)
    
    preprocess_auscpr <- function(auscpr_rawfile){
      auscpr_rawfile <- auscpr_rawfile %>% 
        mutate(survey = "Australian CPR") 
      #list metadata content
      auscpr_meta <- c("survey","TripCode","Sample_ID","Region","Latitude","Longitude","SampleTime_UTC","SampleTime_Local","Year_Local","Month_Local","Day_Local","Time_Local24hr","SampleVolume_m3")
      #subset metadata
      df_metadata <- auscpr_rawfile %>% 
        select(all_of(auscpr_meta)) %>% 
        rename(c("sample_id" = "Sample_ID", "latitude" = "Latitude", "longitude" = "Longitude", "sampleTime_UTC" = "SampleTime_UTC"))
      
      #save metadata file
      metadata_filename <- "cpr_auscpr_metadata.csv"
      write_csv(df_metadata, paste0("data_input/CPR/",metadata_filename))
      print(paste0("Metadata dataframe saved: ", metadata_filename))
      
      df_zoop <- auscpr_rawfile %>% 
        select(c("Sample_ID",!all_of(auscpr_meta))) %>% 
        select(-c("SatSST_degC","PCI","SatChlaSurf_mgm3","BiomassIndex_mgm3")) %>% 
        rename("sample_id" = "Sample_ID")
      
      #load taxon list
      auscpr_taxonlist <- read_csv("data_input/CPR/CPR_raw_data/auscpr_taxonlist.csv", show_col_types = F)
      count_recognizedAphiaID <- sum(names(df_zoop) %in% auscpr_taxonlist$taxon_auscpr, na.rm = T)
      print(paste0("Taxon with AphiaID: ", count_recognizedAphiaID, " out of ", (length(names(df_zoop)) - 1)))
      
      #remove unidentified / general taxa
      unidentified_taxa <- c("Egg","Egg mass","Fish egg","Fish j","Fish larvae","Fish scales","Nauplii zooplankton","Trochophore larvae","Unid invert larvae") #auscpr-specific list
      df_noUnid <- df_zoop %>% select(!any_of(unidentified_taxa)) 
      
      #rename columns from taxonomic name to AphiaID 
        #prior step: remove unidentified / general taxa from taxon list and keep only taxon presently identified in auscpr abundance raw file
        auscpr_taxonlist <- auscpr_taxonlist %>% 
          filter(., !(taxon_auscpr %in% unidentified_taxa)) %>% 
          filter(., (taxon_auscpr %in% names(df_noUnid)))
        #to select columns with matched AphiaID
        df_withmatch <- df_noUnid %>% 
          select(-c("sample_id")) %>% 
          select(match(names(df_noUnid), auscpr_taxonlist$taxon_auscpr, nomatch = 0))
        #to rename based on the match
        names(df_withmatch)[match(names(df_withmatch),auscpr_taxonlist$taxon_auscpr)] <- auscpr_taxonlist$aphiaID
      
      #combine columns of the same name
      auscpr_long <- df_withmatch %>%
        add_column("sample_id" = df_noUnid$sample_id, .name_repair = "minimal") %>%
        pivot_longer(cols = !sample_id, names_to = "AphiaID", values_to = "value")
      
      #sum up abundance per taxon
      DF_sum <- auscpr_long %>% group_by(sample_id, AphiaID) %>% summarize(Total=sum(value))
      
      #revert to format of having AphiaID as columns
      DF_reverted <- DF_sum %>% pivot_wider(names_from = "AphiaID",values_from="Total")
      abundance_dataframe <- DF_reverted %>% ungroup()
      
      #save abundance dataframe
      abundance_dataframe_filename <- "cpr_auscpr_abundance.csv"
      write_csv(abundance_dataframe, paste0("data_input/CPR/",abundance_dataframe_filename))
      print(paste0("Abundance dataframe saved: data_input/CPR/", abundance_dataframe_filename))
      
    }
    
  #pre-process Southern Ocean CPR survey
    preprocess_socpr <- function(socpr_rawfile){
      #assign sample_id
      socpr_rawfile <- socpr_rawfile %>% 
        mutate(survey = "Southern Ocean CPR", sample_id = paste("SOCPR", cur_group_rows(), sep = "")) %>% 
        relocate(sample_id, .before=Tow_Number) #primary key = Sample_ID 
      
      
      #list metadata content
      socpr_meta <- c("survey","sample_id","Tow_Number","Ship_Code","Time","Date","Month","Year","Season","Latitude","Longitude","Segment_No.","Segment_Length")
      #subset metadata
      df_metadata <- socpr_rawfile %>% 
        select(all_of(socpr_meta)) %>% 
        mutate(sampleTime_UTC = as.POSIXct(paste(Date, Time, sep="T"), format = "%d-%b-%YT%H:%M:%S", tz = "UTC")) %>% 
        rename(c("latitude" = "Latitude", "longitude" = "Longitude"))
      
      #save metadata file
      metadata_filename <- "cpr_socpr_metadata.csv"
      write_csv(df_metadata, paste0("data_input/CPR/",metadata_filename))
      print(paste0("Metadata dataframe saved: ", metadata_filename))
      
      df_zoop <- socpr_rawfile %>% 
        select(c("sample_id",!all_of(socpr_meta))) %>% 
        select(-c("Total abundance","Phytoplankton_Colour_Index","Fluorescence","Salinity","Water_Temperature","Photosynthetically_Active_Radiation")) 
      
      #load taxon list
      socpr_taxonlist <- read_csv("data_input/CPR/CPR_raw_data/socpr_taxonlist.csv", show_col_types = F)
      
      #remove unidentified / general taxa
      unidentified_taxa <- c("Egg indet","Egg mass","Nauplius indet") #socpr-specific list
      df_noUnid <- df_zoop %>% select(!any_of(unidentified_taxa))
      
      #report number of taxon with recognized aphiaID
      count_recognizedAphiaID <- sum(names(df_noUnid) %in% socpr_taxonlist$taxon, na.rm = T)
      print(paste0("Taxon with AphiaID: ", count_recognizedAphiaID, " out of ", (length(names(df_noUnid)) - 1)))
      
      #rename columns from taxonomic name to AphiaID 
        #prior step: remove unidentified / general taxa from taxon list 
        socpr_taxonlist <- socpr_taxonlist %>% 
          filter(., !(taxon %in% unidentified_taxa)) %>% 
          filter(., (taxon %in% names(df_noUnid)))
        ##to select columns with matched AphiaID
        df_withmatch <- df_noUnid %>% 
          select(match(socpr_taxonlist$taxon,names(df_noUnid)))
        #to rename based on the match    
        names(df_withmatch)[match(names(df_withmatch),socpr_taxonlist$taxon)] <- socpr_taxonlist$aphiaID
      
      #combine columns of the same name
      socpr_long <- df_withmatch %>%  
        add_column("sample_id" = df_noUnid$sample_id, .name_repair = "minimal") %>% #return sample_id column
        pivot_longer(cols = !sample_id, names_to = "AphiaID", values_to = "value")
      
      #sum up abundance per taxon
      DF_sum <- socpr_long %>% group_by(sample_id, AphiaID) %>% summarize(Total=sum(value))
      
      #revert to format of having AphiaID as columns
      DF_reverted <- DF_sum %>% pivot_wider(names_from = "AphiaID",values_from="Total")
      abundance_dataframe <- DF_reverted %>% ungroup()
      
      #save abundance dataframe
      abundance_dataframe_filename <- "cpr_socpr_abundance.csv"
      write_csv(abundance_dataframe, paste0("data_input/CPR/",abundance_dataframe_filename))
      print(paste0("Abundance dataframe saved: data_input/CPR/", abundance_dataframe_filename))
    }
    
  #pre-process MBA data for North Atlantic and North Pacific CPR Surveys
    preprocess_mba_cpr <- function(mbacpr_rawfile){
      
      mbacpr_rawfile <- mbacpr_rawfile %>% 
        mutate(sampleTime_UTC = as.POSIXct(midpoint_date_gmt, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),.after = "sample_id") %>% 
        filter(sampleTime_UTC >= "1997-09-01")
      #to remove taxa for avoidance of double counting (Richardson et al., 2006)
      
      #to separate mba data into north pacific and north atlantic data
      #identify cpr belonging to North Atlantic and North Pacific CPR
      mba_npacific <- mbacpr_rawfile %>% 
        filter(longitude >= 100 | longitude <= -100) %>%
        mutate(survey = "North Pacific CPR") 
          
      mba_natlantic <- mbacpr_rawfile %>% 
        filter(longitude < 100 & longitude > -100) %>% 
        mutate(survey = "Atlantic CPR") 
      
        #North Atlantic CPR: remove freshwater (Lake Tanganyika)
        mba_natlantic <- mba_natlantic %>% 
          filter(!sample_id %in% c("2ALT-1","2ALT-5","3ALT-2","3ALT-5"))
        
      mba_meta <- c("survey","sample_id","midpoint_date_gmt","sampleTime_UTC","latitude","longitude","chlorophyll_index")
      #get metadata
      mba_npacific_metadata <- mba_npacific %>% 
        select(all_of(mba_meta)) 
      
        #save metadata file
        npacific_metadata_filename <- "cpr_npacific_metadata.csv"
        write_csv(mba_npacific_metadata, paste0("data_input/CPR/",npacific_metadata_filename))
        print(paste0("Metadata dataframe saved: ", npacific_metadata_filename))
      
      mba_natlantic_metadata <- mba_natlantic %>% 
        select(all_of(mba_meta)) 
        
        #save metadata file
        natlantic_metadata_filename <- "cpr_natlantic_metadata.csv"
        write_csv(mba_natlantic_metadata, paste0("data_input/CPR/",natlantic_metadata_filename))
        print(paste0("Metadata dataframe saved: ", natlantic_metadata_filename))
      
      survey <- c("npacific","natlantic")
                  
      for(i in 1:2){
        if(survey[i] == "natlantic"){
          df_zoop <- mba_natlantic %>% select(c("sample_id",!all_of(mba_meta)))
        }else if(survey[i] == "npacific"){
          df_zoop <- mba_npacific %>% select(c("sample_id",!all_of(mba_meta)))
        }else{ stop("Error in reading dataframe of zooplankton abundance")}
        
        print(paste0("CPR survey in-process: ",survey[i]))
        #load taxon list
        mba_taxonlist <- read_csv("data_input/CPR/CPR_raw_data/mba_cpr_taxonlist.csv", show_col_types = F)
        
        #remove unidentified / general taxa 
        unidentified_taxa <- c("Egg indet","Egg mass","Nauplius indet") 
        df_noUnid <- df_zoop %>% select(!any_of(unidentified_taxa))
        
        #report number of taxon with recognized aphiaID
        count_recognizedAphiaID <- sum(names(df_noUnid) %in% mba_taxonlist$taxa_name, na.rm = T)
        print(paste0("Taxon with AphiaID: ", count_recognizedAphiaID, " out of ", (length(names(df_noUnid)) - 1))) #account for column "sample_id"
        
        #rename columns from taxonomic name to AphiaID 
        #prior step: remove unidentified / general taxa from taxon list 
        mba_taxonlist <- mba_taxonlist %>% 
          filter(., !(taxa_name %in% unidentified_taxa)) %>% 
          filter(., (taxa_name %in% names(df_noUnid)))
        ##to select columns with matched AphiaID
        df_withmatch <- df_noUnid %>% 
          select(match(mba_taxonlist$taxa_name,names(df_noUnid)))
        #to rename based on the match    
        names(df_withmatch)[match(names(df_withmatch),mba_taxonlist$taxa_name)] <- mba_taxonlist$aphia_id
        
        # #return "sample_id" column
        # abundance_dataframe <- df_withmatch %>% 
        #   add_column("sample_id" = df_noUnid$sample_id, .before = 1, .name_repair = "minimal") 
          
        #combine columns of the same name
        mbacpr_long <- df_withmatch %>%  
          add_column("sample_id" = df_noUnid$sample_id, .name_repair = "minimal") %>% #return sample_id column
          pivot_longer(cols = !sample_id, names_to = "AphiaID", values_to = "value")
        
        #sum up abundance per taxon
        DF_sum <- mbacpr_long %>% group_by(sample_id, AphiaID) %>% summarize(Total=sum(value))
        
        #revert to format of having AphiaID as columns
        DF_reverted <- DF_sum %>% pivot_wider(names_from = "AphiaID",values_from="Total")
        abundance_dataframe <- DF_reverted %>% ungroup()
        
        #save abundance dataframe
        abundance_dataframe_filename <- paste0("cpr_",survey[i],"_abundance.csv")
        write_csv(abundance_dataframe, paste0("data_input/CPR/",abundance_dataframe_filename))
        print(paste0("Abundance dataframe saved: data_input/CPR/", abundance_dataframe_filename))
        
      }
    } 

#02 to compute proportions of trophic groups
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
    
#02 to map global CPR

#cpr_meta <- list.files(path = "data_input/CPR_on-process/", pattern = "*\\metadata.csv", full.names = TRUE)

map_globalcpr <- function(file.list){
  #extract metadata  
  filenames <- basename(file.list)
  cpr <- data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("survey", "sample_id", "latitude", "longitude", "sampleTime_UTC"))))
  cpr[,1:2] <- lapply(cpr[,1:2], as.character) # survey and sample id
  cpr[,3:4] <- lapply(cpr[,3:4], as.double) # latitude and longitude
  cpr[,5] <- as.POSIXct(cpr[,5]) #sample time UTC
  
  
  #loop through each survey
    for(i in 1:length(filenames)){
    file_survey <- str_extract(filenames[i], "(?<=_)[^_]+")
    print(paste("Survey: ",file_survey, sep=""))
    
    cpr_metadata <- read_csv(file.list[i]) %>% 
      mutate(survey = file_survey) %>% 
      select("survey","sample_id","latitude","longitude","sampleTime_UTC")
    
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
    tmap_save(globalcpr_map, filename=paste0("output/plots/Global-CPR-map_",date,".png"),
              width = 400,
              height = 200,
              units = "mm",
              dpi = 400)
}

# 03. to update species list
#set date for version control
  generate_globalCPR_taxonlist <- function(){
    auscpr_taxonlist <- read_csv("data_input/CPR/CPR_raw_data/auscpr_taxonlist.csv")
    mba_taxonlist <- read_csv("data_input/CPR/CPR_raw_data/socpr_taxonlist.csv")
    socpr_taxonlist <- read_csv("data_input/CPR/CPR_raw_data/mba_cpr_taxonlist.csv")
    
    #unique_taxonlist <- merged_taxonlist %>%
      
    pattern <- "^\\w+(?:[ .]+(?!(?:[mjf]|zoea|agg|sol|Grp[135]|larvae|megalopa|CIV|CV)\\b)(?:spp\\.?|\\w+))?"
    
    auscpr_taxonlist$taxon_auscpr <- str_extract(auscpr_taxonlist$taxon_auscpr, pattern)
    auscpr_taxonlist_merged <- auscpr_taxonlist %>% 
      group_by(aphiaID) %>% 
      summarise(taxa_name = first(taxon_auscpr))
    
    merged_taxonlist <- auscpr_taxonlist_merged %>% 
      full_join(socpr_taxonlist %>% select("aphia_id","taxa_name"), by = c("aphiaID" = "aphia_id","taxa_name")) %>% 
      full_join(mba_taxonlist, by = c("aphiaID", "taxa_name" = "taxon"))    
      
    merged_taxonlist <- merged_taxonlist %>% 
      group_by(aphiaID) %>% 
      summarise(taxa_name = first(taxa_name))
    
    #output csv file of merged taxa list
    write_csv(merged_taxonlist, "output/CPR/cpr_merged_taxonlist.csv")
    write_csv(merged_taxonlist, "data_input/CPR/cpr_merged_taxonlist.csv")
  }

#04. to compute abundances per taxon  

    # #zoop_list is any list with aphiaID
    # #00 get abundance list 
    # file.list <- list.files(path = "data_input/CPR_on-process/", pattern = "*\\complete.csv", full.names = TRUE)
    # #set date for version control
    # date <- "20082025"
    
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
    
#Latest revision 06.08.2025