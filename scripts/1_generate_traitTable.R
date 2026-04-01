# ---
# title: generate_traitTable
# author: Bryan Alpecho
# date: 2025-
# output: csv file
# ---

##Aims
#1. Extract from Pata & Hunt (2023) trait table by species and trait list
#2. Fill-up gaps by deriving trait values at genus level
#3. Generate trait table
#4. Compute for taxon abundances

###Code for aims 1 and 2 adapted from
#'7-Database_subset_template'
# ---
# title: "Template for subsetting the database"
# author: "Patrick Pata"
# date: "06/30/2023"
# output: html_document
# Repository: https://github.com/Pelagic-Ecosystems/Zooplankton_trait_database
# ---
# Article Citation: Pata, P. R., & Hunt, B. P. V. (2023). Harmonizing marine zooplankton trait data toward a mechanistic understanding of ecosystem functioning. Limnology and Oceanography, 70(S1), S8–S27. https://doi.org/10.1002/lno.1247
###


#Input 
##Set date (version control of output)
date <- "31032026" 

#### Load libraries and data ####
    packages <- c("tidyverse",
                  "openxlsx")
    
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
    
    source("Pata & Hunt (2023) toolkit.R")
    
    s.format <- read_csv("data_input/PataHunt_traits/trait_dataset_standard_format_20230628.csv")[-1,]
    
    trait.directory <- read_csv("data_input/PataHunt_traits/trait_directory_20230628.csv") %>% 
      distinct(traitID, .keep_all = TRUE)
    
    # taxonomy table
    taxonomy <- read_csv("data_input/PataHunt_traits/taxonomy_table_20230628.csv") 
    
    # stage table
    lifestagelist <- read_csv("data_input/PataHunt_traits/lifestage_directory_20230628.csv") %>% 
      select(-c(majorgroup, notes))
    
    # # Load Level 1 dataset
    # traits.lvl1 <- read_csv("data_input/PataHunt_traits/Trait_dataset_level1/trait_dataset_level1-2023-06-28.csv")
    
    # Load Level 2 dataset
    traits.lvl2 <- read_csv("data_input/PataHunt_traits/Trait_dataset_level2/trait_dataset_level2-2023-09-14.csv")
    
    # Change assocTemperature column class to allow row binds
    traits.lvl2$assocTemperature <- as.character(traits.lvl2$assocTemperature)
    traits.lvl2$verbatimTemperature <- as.character(traits.lvl2$verbatimTemperature)

#1. Extract from Pata & Hunt (2023) trait table by species and trait list
  #### List selected species and traits ####
    # List of traits to extract
      trait.list <- c("trophicGroup")
    
    # List of species to extract
    # Need to know the species name that matches the species names in this database. Alternatively, can do this matching using AphiaIDs and the acceptedNameUsageID field.
    
    #unique aphia ID (836 total unique aphia)
      species.list <- read_csv("data_input/CPR_on-process/cpr_all_aphia-taxon.csv")
      
    # Match species list with taxonomy file to get taxonID AND keep aphiaID without taxonID
      species.list <- species.list  %>%
        left_join(taxonomy, by = c("aphiaID")) %>% 
        distinct(aphiaID, .keep_all = TRUE)
      
      copepod.list <- species.list %>% 
        filter(class == "Copepoda")
          
  # Subset the species and traits 
  # This will return the trait values based on literature or derived using related traits.
     
    # Create a subset of records based on the species and traits list.
    trait.subset <- traits.lvl2 %>% 
      filter(traitName %in% trait.list) %>% 
      filter(taxonID %in% species.list$taxonID) %>% 
      arrange(scientificName)
    
    # Evaluate how many species have trait information
    trait.subset.perc <- trait.subset %>% 
      distinct(traitName, taxonID) %>% 
      group_by(traitName) %>% 
      summarise(Nrecords = n(), .groups = "drop") %>% 
      mutate(Perc.sp = Nrecords / nrow(species.list) * 100)
    
    #TG.trait.subset
    TG.trait.subset <- trait.subset %>% 
      filter(traitName == "trophicGroup") %>% 
      select(taxonID, stageID, scientificName, class, order, family, genus, traitValue) %>% 
      distinct() %>% 
      mutate(traitValue = as.character(traitValue))
    
    TG.trait.subset.copepods <- trait.subset %>% 
      filter(traitName == "trophicGroup") %>% 
      filter(class == "Copepoda")
    
    #export trait subset for trophicGroup from Pata & Hunt (2023) trait table
    write_csv(TG.trait.subset, paste("Output/traits/TG_trait_subset_",date,".csv", sep=""))
    write_csv(TG.trait.subset.copepods, paste("Output/traits/TG_trait_subset_copepods_",date,".csv",sep=""))
    
#2. Fill-up gaps by deriving trait values at genus level 
  ### Filling up the gaps ###
  #### Using Method 2: Generalize at broader taxonomic level ####
  cat.missing <- species.list %>% 
    filter(taxonID %notin% filter(trait.subset, traitName == "trophicGroup")$taxonID) %>% 
    mutate(traitName = "trophicGroup")
  
  # Generalize trophic group at genus level. 
  genus.bin.tg <- traits.lvl2 %>% 
    # Need to know the prefix of the binary trait.
    filter(grepl("TG.",traitName)) %>% 
    mutate(traitValue = as.numeric(traitValue)) %>% 
    # TODO: Choose the taxonomic level
    group_by(genus, traitName) %>% 
    # Collect all instances of a categorical level per group
    summarise(traitValue = sum(traitValue),  
              .groups = "drop") %>% 
    group_by(genus) %>% 
    #sum up trait
    mutate(ntaxa = sum(traitValue)) %>% 
    ungroup() %>% 
    #convert to proportions
    mutate(trait_prop = traitValue / ntaxa)
    #write_csv(genus.bin.tg, "data_output/genus_bin_tg_30052025.csv")
  
    # Generalize trophic group at genus level for copepods only
    genus.bin.tg.copepods <- traits.lvl2 %>% 
      # Need to know the prefix of the binary trait.
      filter(grepl("TG.",traitName)) %>% 
      filter(class == "Copepoda") %>% 
      mutate(traitValue = as.numeric(traitValue)) %>% 
      # TODO: Choose the taxonomic level
      group_by(genus, traitName) %>% 
      # Collect all instances of a categorical level per group
      summarise(traitValue = sum(traitValue),  
                .groups = "drop") %>% 
      group_by(genus) %>% 
      mutate(ntaxa = sum(traitValue)) %>% 
      ungroup() %>% 
      
      # convert to proportions
      mutate(trait_prop = traitValue / ntaxa)
      # write_csv(genus.bin.tg, "data_output/genus_bin_tg_30052025.csv")
  
  #categorize trophicGroup to omnivore and carnivore
  genus.bin.tg.copepods <- genus.bin.tg.copepods %>% 
    mutate(FG_derived = case_when(
      traitName == "TG.herbivore" | traitName == "TG.omnivore" | traitName == "TG.detritivore" ~ "omnivore",
      traitName == "TG.carnivore" ~ "carnivore"
    )) %>% 
    select(c(genus, FG_derived)) %>% 
    filter(!is.na(genus)) %>% 
    distinct()
  
  #determine genus classified as both omnivore and carnivore
  nondistinct <- c(genus.bin.tg.copepods %>%
      group_by(genus) %>% 
      filter(n() != 1) %>% 
      distinct(genus))
  
  #classify omnivore-carnivore to carnivore
  genus.bin.tg.copepods <- genus.bin.tg.copepods %>% 
    mutate(distinct_or_not = genus %in% nondistinct$genus) %>% 
    mutate(FG_distinct = case_when(
      genus %in% nondistinct$genus ~ "carnivore",
      .default = FG_derived
    )) %>% 
    select(c(genus,FG_distinct)) %>% 
    distinct()
  
  #filter for copepods with missing trait TG
  cat.missing.copepods <- cat.missing %>%  filter(class=="Copepoda")
  #create df for missing copepods with filled-up traits 
  cat.missing.copepods <- cat.missing.copepods %>% 
    left_join(genus.bin.tg.copepods, by = "genus") %>% 
    distinct(genus, FG_distinct) %>% 
    filter(!is.na(FG_distinct))
  
  write_csv(cat.missing.copepods, paste("Output/traits/TG_copepods_derivedTrait_",date,".csv",sep=""))
  summary(as.factor(cat.missing.copepods$FG_distinct)) #19 carnivore, 32 omnivore 

  #for Copepods
    #integrate the values of generalized traitValues to the trait subset
    trait.subset.missing.TG <- traits.lvl2 %>% 
      filter(taxonID %in% cat.missing$taxonID) %>% 
      filter(genus %in% cat.missing.copepods$genus) %>%  
      distinct(scientificName, .keep_all=T) %>% 
      mutate(traitName = "trophicGroup") %>% 
      mutate(traitValue = NULL) %>% 
      mutate(traitUnit = NA) %>% 
      mutate(valueType = "categorical") %>% 
      mutate(basisOfRecordDescription = "Trophic group derived at genus level (30.05.2025)") %>% 
      left_join(cat.missing.copepods[c("genus","FG_distinct")], by = "genus") %>% 
      rename(traitValue = FG_distinct) %>% 
      relocate(traitValue, .after = traitName)
  
  #combine original trait table (TG-specific) with traits generalized at genus level
  TG.copepods <-   TG.trait.subset.copepods %>% 
    bind_rows(trait.subset.missing.TG) 
  
  # Evaluate progress of data gap estimation
  TG.copepods %>% 
    distinct(traitName, taxonID) %>% 
    group_by(traitName) %>% 
    summarise(Nrecords = n(), .groups = "drop") %>% 
    mutate(Perc.sp = Nrecords / nrow(copepod.list) * 100) #86.2 % (344 out of 399) out of all copepods
  
  #export trait table for copepods
  write_csv(TG.copepods, paste("Output/traits/TG_copepods_combined_",date,".csv",sep=""))
  
  
#3. Generate trait table
  #Get species list
    #unique aphia ID (836 total unique aphia)
    species.list <- read_csv("data_input/CPR_on-process/cpr_all_aphia-taxon.csv")
    # Match species list with taxonomy file to get taxonID AND keep aphiaID without taxonID
    species.list <- species.list  %>%
      left_join(taxonomy, by = c("aphiaID")) %>% 
      distinct(aphiaID, .keep_all = TRUE)
  
  #Integrate Pata & Hunt (2023) to species list

    import_TG.PataHunt <- function(species.list){
      # Read the trait dataset
      zoop.traits <- read_csv("Output/traits/TG_trait_subset_03-06-2025.csv")
      copepod.traits <- read_csv("Output/traits/TG_copepods_combined_03-06-2025.csv")
      
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
        mutate(traitValue = ifelse(is.na(traitValue), "not determined", traitValue))
    
      return(species.list)
    }

  #filter-or-not
    assign_filter_or_not <- function(zoop_list){
      ff.order <- c("Pyrosomida","Salpida","Dolioida","Copelata","Aplousobranchia","Phlebobranchia","Stolidobranchia")
      zoop_list %>% 
        mutate(filter_or_not = case_when(order %in% ff.order ~ "1",
                                         .default = "0")) %>% 
        relocate(filter_or_not, .after = "traitValue")
    }
    
  #Integrate Campbell (2023) trait table 
    
    
  #Other trait list for comparison
    

  species.list <- species.list %>% 
    import_TG.PataHunt() %>% 
    assign_filter_or_not()

  #check which have duplicated traitValues
    check_duplicate <- function(species_list){
      
      #check which have duplicated 'traitValues'
      if(anyDuplicated(species.list$aphiaID) > 0){
        cat("Need for review!\nDetected duplicate in assigned trait values for: ")
        print(unique(species.list$scientificName[duplicated(species.list$aphiaID)]))
      }else if(anyDuplicated(species.list$aphiaID) == 0){
        cat("Looks good!\nNo detected duplicates")
      }else{ stop("Error in duplicate checking", call. = FALSE)}
      
      ##Save trait table
      print(paste0("Trait table saved: Output/traits/TG_trait-table-",date,".csv"))
      write_csv(species.list, paste("Output/traits/TG_trait-table-",date,".csv",sep=""))
    }

    
  #20.08.2025
    
    #NOT WORKING

    update_tg <- function(tg, taxon, taxonomicRank, trait_table){
        trait_table_rev <- trait_table %>% 
          rename(TG_old = FG_complete)  %>% #need to rename FG_complete globally
          mutate(FG_complete = case_when(
            {{taxonomicRank}} == {{taxon}} ~ {{tg}},
            .default = TG_old
            )) %>% 
          select(-c("TG_old")) %>% 
          relocate(FG_complete, .before = "References")
      
        assign(traits_rev, trait_table_rev)
        return(traits_rev)
    }
    update_tg("omnivore", "Euphausiacea", order, traits) 
      
    #working for Euphasiaceae 
    trait_table_rev <- traits %>% 
      rename(TG_old = FG_complete)  %>% #need to rename FG_complete globally
      mutate(FG_complete = case_when(
        order == "Euphausiacea" ~ "omnivore",
        .default = TG_old
      )) %>% 
      select(-c("TG_old")) %>% 
      relocate(FG_complete, .before = "References")
    
    
    write_csv(trait_table_rev, paste("Output/data/traits/trait-table-revised_",date,".csv",sep=""))
    
  #Remaining issue: duplicated traitValues are not yet removed resulting to many-to-many relationships 
    #Current sol'n: manually review the trait table
  
  #Latest revision for this code: 
  #26.06.2025 add abundance calculation to aim #4
  #28.07.2025 Moved aim #4 to 0_wrangling_CPR
  #future: integrate other trait lists