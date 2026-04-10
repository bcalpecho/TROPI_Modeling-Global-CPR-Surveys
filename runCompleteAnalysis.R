####                          ####
# DATA PREPARATION AND ANALYSIS #
####                          ####

# 0 Setup ####
  # Load functions and set directories
  
  source("functions_bank.R")
  input_dir <- "data_input/" # where files are read from
  output_dir <- "output/" # where files are written to
  
  #### set date and list of models for version control and coverage of CPR surveys
  date <- "08042026"
  survey_list <- c("auscpr", "socpr","npacific","natlantic")

########## Pre-process raw CPR data #########
  # Read in raw CPR csv files (includes translation from taxa_name to AphiaID)
  # Outputs: metadata and abundance dataframe for each CPR survey
  
  #auscpr
    auscpr_rawfile <- read_csv("data_input/CPR/CPR_raw_data/IMOS_-_Zooplankton_Abundance_and_Biomass_Index_(CPR)-raw_data.csv")
    preprocess_auscpr(auscpr_rawfile)

  #socpr
    socpr_rawfile <- read_csv("data_input/CPR/CPR_raw_data/AADC-00099_29August2025.csv")
    preprocess_socpr(socpr_rawfile)

  #mba natlantic and npacific cpr
    mba_rawfile <- read_csv("data_input/CPR/CPR_raw_data/NPacifc_Atlantic_CPR_1958_2021.csv")
      #above dataset is not uploaded in the repository, but earlier versions are available online.
    preprocess_mba_cpr(mba_rawfile)
  
  #to map the global cpr (Output: Figure 1)
    cpr_meta <- list.files(path = "data_input/CPR/", pattern = "*\\metadata.csv", full.names = TRUE)
    map_globalcpr(cpr_meta)
      
########## 1 Generate trait table ##########

  #read in trait table
    fTraitTable <- read_csv("data_input/traits/fTraitTable.csv")
  
  # #Steps in generating trait table for a given CPR taxa list 
  # #Get taxa list
  # # list of unique AphiaID in four CPR surveys (836 total unique AphiaID)
  #   taxa.list <- read_csv("data_input/CPR/cpr_merged_taxonlist.csv") #generated through 'generate_globalCPR_taxonlist()' function in 'helper_CPR' script
  # 
  # # taxonomy table (Pata & Hunt, 2023)
  #   # Citation: Pata, P. R., & Hunt, B. P. V. (2023). Harmonizing marine zooplankton trait data toward a mechanistic understanding of ecosystem functioning. Limnology and Oceanography, 70(S1), S8–S27. https://doi.org/10.1002/lno.12478
  #   taxonomy <- read_csv("data_input/PataHunt_taxonomy_table_20230628.csv")
  # 
  # # match species list with taxonomy file to get taxonID AND keep AphiaID without taxonID
  #   taxa.list <- taxa.list  %>%
  #     left_join(taxonomy, by = c("aphiaID")) %>%
  #     distinct(aphiaID, .keep_all = TRUE) 
  # 
  #   taxa.list <- taxa.list %>%
  #     import_TG.PataHunt() %>%  #Function 1.1 Integrate Pata & Hunt (2023) assigned trait value to species list
  #     assign_filter_or_not() #Function 1.2 Assign 'gelatinous filter-feeders' group
  # 
  #   adapt_reviewedTraitTable(taxa.list) #Function 1.3 check for duplicates and save trait table
  # #To finalize the trait table, a manual review of assigned traits is still needed.
    
########## 2 Extract Chl-a ##########
  
  #Steps in extracting Chl-a data given spatial and temporal coordinates of CPR sampling points
  #2.1 aggregate raster of OC-CCI
    #8-day OC-CCI
    aggregate_ncdf(survey_list, "eightday", date)
    #monthly OC-CCI
    aggregate_ncdf(survey_list, "monthly", date)

  #2.2 extract raster data at CPR sampling points
    #8-day OC-CCI
    extract_chla(survey_list, "eightday")
    #monthly OC-CCI
    extract_chla(survey_list, "monthly")

  #2.3 Fill-up gaps of eight-day by monthly OC-CCI values if available
    fill_up_gaps(survey_list, date)

########## 3 Generate Global data frame ##########
  
  #3.1 get abundance list of CPR surveys (list of csv files of zooplankton abundances from each CPR survey)
    file.list <- list.files(path = "data_input/CPR/", pattern = "*\\_abundance.csv", full.names = TRUE)
  
  #3.2 Compute for proportions of zooplankton trophic groups
    compute_proportions_perSurvey(file.list) 
  
  #3.3 Combining variables altogether into a dataframe for each CPR survey
    generate_df_perSurvey(date)
  
  #3.4 generate GLOBAL dataframe and compute for ratios
    generate_globalCPR_dataframe(date)
    
########## 4 Model the Global CPR ##########

  # 4.1 Import data
    #for omnivores and carnivores
    df <- read_rds(paste0("data_input/global_df_complete_",date,".rds"))
    #for filter-feeders
    #df_filter <- read_rds(paste0("data_input/global_df_complete_Filter_",date,".rds"))
  
  # 4.2 Selected models (see '4_model_globalCPR' script for competing models)
    Omni_mdl_zib <- glmmTMB(ROC_SVT_zib ~ chla_sqrt + survey + (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                            ziformula = ~1,
                            data = df, family = beta_family(link = "logit"))
    
    Carni_mdl_zib <- glmmTMB(RCO_SVT_zib ~ chla_sqrt + survey + (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                             ziformula = ~1,
                             data = df, family = beta_family(link = "logit"))
    
    Filter_mdl_zib <- glmmTMB(RFF_SVT_zib ~ chla_sqrt + survey + (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                              ziformula = ~1,
                              data = df, family = beta_family(link = "logit"))
  
    #list the models
    mdl_list <- list(Carni_mdl_zib, Omni_mdl_zib, Filter_mdl_zib)
    names(mdl_list) <- c("Carni","Omni","Filter")
    
  # Model summary
    summary(Omni_mdl_zib)
    summary(Carni_mdl_zib)
    summary(Filter_mdl_zib)
    
  # 4.3 to summarize the coefficients of the models
     #determine R^2 (coefficient of determination)
      Omni_mdl_zib_R2 <- MuMIn::r.squaredGLMM(Omni_mdl_zib)
      Carni_mdl_zib_R2 <- MuMIn::r.squaredGLMM(Carni_mdl_zib)
      Filter_mdl_zib_R2 <- MuMIn::r.squaredGLMM(Filter_mdl_zib)
    
     #extract coefficients
      Omni_mdl_zib_coeff <- summary(Omni_mdl_zib)$coefficients$cond %>% as.data.frame()
      Carni_mdl_zib_coeff<- summary(Carni_mdl_zib)$coefficients$cond %>% as.data.frame()
      Filter_mdl_zib_coeff <- summary(Filter_mdl_zib)$coefficients$cond %>% as.data.frame()
    
    #to summarize the estimated proportions of zooplankton trophic groups with Chl-a as a predictor
    summary_predictionsChla()
    
    #to summarize the contributions of fixed and random effects to variance in the model
    summary_mdlVariance(mdl_list)
    
    #to summarize the change in proportion of zooplankton trophic groups per unit of Chl-a decline
    summary_Delta_Per_Unit_Chla_Decline(mdl_list)
    
    #to compute for difference in estimates across CPR surveys
    summary_predictionsPerSurvey(mdl_list)
    
   # 4.5 to save/export selected models
    {for(i in length(mdl_list)){
      write_rds(mdl_list[[i]], paste0("output/mdls/",names(mdl_list[i]),"_mdl_zib.rds"))
    }
    remove(i)}

########## 5 Assess the Model ##########
  #to load back the models
    Carni_mdl_zib <- read_rds("output/mdls/Carni_mdl_zib.rds")
    Omni_mdl_zib <- read_rds("output/mdls/Omni_mdl_zib.rds")
    Filter_mdl_zib <- read_rds("output/mdls/Filter_mdl_zib.rds")
    
  #5.1 quantile-quantile plot to assess normality of residuals
    plot_QQ(mdl_list)
    
  #5.2 mean variance plot to assess homogeneity of variance
    plot_meanVariance(mdl_list)
    
  #5.3 to plot the intercepts of Longhurst Provinces
    plot_Longhurst(mdl_list)
    
  #5.4 to plot point density plots of Tow within Survey slope and intercept
    plot_TowSlopeAndIntercept(mdl_list)
    
  #5.5 to plot residuals of models with and without the random effects
    plot_residuals(mdl_list)
  
########## 6b Predict Global CPR ##########
  
  #6. Compute change in proportion of zooplankton trophic groups per unit of Chl-a decline
    
    
  #6b.1 to predict zooplankton trophic group 
    #load ensemble
    ensembles_median <- list.files(file.path("data_input","ensemble_chlos","median","chlos"), full.names = TRUE)
    
    #process the SSP scenarios one-by-one (One ensemble each SSP scenario)
    
  #6b.2 list the selected models (double check the order of mdl_list and names)
    mdl_list <- list(Carni_mdl_zib, Omni_mdl_zib, Filter_mdl_zib)
    names(mdl_list) <- c("Carni","Omni","Filter")
    
    #ssp126
    predict_zoop_delta(ensembles_median[1], mdl_list)
    #ssp245
    predict_zoop_delta(ensembles_median[2], mdl_list)
    #ssp370
    predict_zoop_delta(ensembles_median[3], mdl_list)
    #ssp585
    predict_zoop_delta(ensembles_median[4], mdl_list)
  
   #6b.3 to compute for delta of trophic groups between 2015 and 2100
    compute_zoop_delta(mdl_list)
  
  #to see summary of model predictions per model
    summary_stats_TG(mdl_list)
    baseline_stats_TG(mdl_list)
  
########## 7 Plot visual summary of model ##########
  
  #7.1 summary for model of omnivores
    plot_model_summary_omnivores(date)
    
  #7.2 summary for model of carnivores
    plot_model_summary_carnivores(date)
  
  #7.3 summary for model of gelatinous filter-feeders
    plot_model_summary_filterfeeders(date)
  
## End ##