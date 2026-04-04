####                          ####
# DATA PREPARATION AND ANALYSIS #
####                          ####

# 0 Setup ####
  # Load functions and set directories
  
  source("functions_bank.R")
  input_dir <- "data_input/" # where files are read from
  output_dir <- "output/" # where files are written to
  
  #### set date and list of models for version control and coverage of CPR surveys
  date <- "29032026"
  tsurvey_list <- c("auscpr", "socpr","npacific","natlantic")
  
# 1 Generate trait table ####

  #Get species list
  # list of unique AphiaID in four CPR surveys (836 total unique AphiaID)
    species.list <- read_csv("data_input/CPR/cpr_merged_taxonlist.csv") #generated through 'generate_globalCPR_taxonlist()' function in 'helper_CPR' script

  # taxonomy table (Pata & Hunt, 2023)
    # Citation: Pata, P. R., & Hunt, B. P. V. (2023). Harmonizing marine zooplankton trait data toward a mechanistic understanding of ecosystem functioning. Limnology and Oceanography, 70(S1), S8–S27. https://doi.org/10.1002/lno.12478
    taxonomy <- read_csv("data_input/PataHunt_taxonomy_table_20230628.csv")

  # match species list with taxonomy file to get taxonID AND keep AphiaID without taxonID
    species.list <- species.list  %>%
      left_join(taxonomy, by = c("aphiaID")) %>%
      distinct(aphiaID, .keep_all = TRUE)
  
    species.list <- species.list %>%
      import_TG.PataHunt() %>%  #Function 1.1 Integrate Pata & Hunt (2023) assigned trait value to species list
      assign_filter_or_not() #Function 1.2 Assign 'gelatinous filter-feeders' group
  
    detect_duplicateAssignedTrait(species.list) #Function 1.3 check for duplicates and save trait table
    #To finalize trait table, perform a manual review of assigned traits to verify trait values.
    
# 2 Extract Chl-a ####
  
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

# 3 Generate Global data frame ####
  
  #3.1 get abundance list of CPR surveys (list of csv files of zooplankton abundances from each CPR survey)
    file.list <- list.files(path = "data_input/CPR/", pattern = "*\\complete.csv", full.names = TRUE)
  
  #3.2 Compute for proportions of zooplankton trophic groups
    compute_proportions_perSurvey(file.list)
  
  #3.3 Combining variables altogether into a dataframe for each CPR survey
    generate_df_perSurvey(date)
  
  #3.4 generate GLOBAL dataframe and compute for ratios
    generate_globalCPR_dataframe(date)
    
# 4 Model the Global CPR ####

  # 4.1 Import data
    #for omnivores and carnivores
    df <- read_rds(paste0("data_input/global_df_complete_",date,".rds"))
    #for filter-feeders
    df_filter <- read_rds(paste0("data_input/global_df_complete_Filter_",date,".rds"))
  
  # 4.2 Selected models (see '4_model_globalCPR' script for competing models)
    Omni_mdl_zib <- glmmTMB(ROC_SVT_zib ~ chla_sqrt + survey + (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                            ziformula = ~1,
                            data = df, family = beta_family(link = "logit"))
    
    Carni_mdl_zib <- glmmTMB(RCO_SVT_zib ~ chla_sqrt + survey + (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                             ziformula = ~1,
                             data = df, family = beta_family(link = "logit"))
    
    Filter_mdl_zib <- glmmTMB(RFF_SVT_zib ~ chla_sqrt + survey + (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                              ziformula = ~1,
                              data = df_filter, family = beta_family(link = "logit"))
  
  # 4.3 to summarize the coefficients of the models
    #determine R^2 (coefficient of determination)
    Omni_mdl_zib_R2 <- MuMIn::r.squaredGLMM(Omni_mdl_zib)
    Carni_mdl_zib_R2 <- MuMIn::r.squaredGLMM(Carni_mdl_zib)
    Filter_mdl_zib_R2 <- MuMIn::r.squaredGLMM(Filter_mdl_zib)
    
    #extract coefficients
    Omni_mdl_zib_coeff <- summary(Omni_mdl_zib)$coefficients$cond %>% as.data.frame()
    Carni_mdl_zib_coeff<- summary(Carni_mdl_zib)$coefficients$cond %>% as.data.frame()
    Filter_mdl_zib_coeff <- summary(Filter_mdl_zib)$coefficients$cond %>% as.data.frame()
    
    #summarize
    tibble_models_zib <- tibble(
      mdl = c("Omni_mdl_zib", "Carni_mdl_zib", "Filter_mdl_zib"),
      Estimate = c(Omni_mdl_zib_coeff$Estimate[2], Carni_mdl_zib_coeff$Estimate[2], Filter_mdl_zib_coeff$Estimate[2]),
      StdError = c(Filter_mdl_zib_coeff$`Std. Error`[2],  Omni_mdl_zib_coeff$`Std. Error`[2], Carni_mdl_zib_coeff$`Std. Error`[2]),
      z_value = c(Filter_mdl_zib_coeff$`z value`[2], Omni_mdl_zib_coeff$`z value`[2], Carni_mdl_zib_coeff$`z value`[2]),
      p_value = c(Filter_mdl_zib_coeff$`Pr(>|z|)`[2], Omni_mdl_zib_coeff$`Pr(>|z|)`[2], Carni_mdl_zib_coeff$`Pr(>|z|)`[2]),
      PseudoR2_marginal = c(Filter_mdl_zib_R2[1], Omni_mdl_zib_R2[1], Carni_mdl_zib_R2[1]),
      PseudoR2_conditional = c(Filter_mdl_zib_R2[2], Omni_mdl_zib_R2[2], Carni_mdl_zib_R2[2]))
    
  # 4.4 to generate a summary of proportion contributed by fixed and random effects to the total variance of the model
    
    #prior step
    mdl_list <- list(Carni_mdl_zib, Omni_mdl_zib, Filter_mdl_zib)
    names(mdl_list) <- c("Carni","Omni","Filter")
    
    #to create a table summary of contributions of fixed and random effects to variance in the model
    summary_mdlVariance(mdl_list)
    
  # 4.5 to save/export selected models
    write_rds(Filter_mdl_zib, "output/mdls/Filter_mdl_zib.rds")
    write_rds(Carni_mdl_zib, "output/mdls/Carni_mdl_zib.rds")
    write_rds(Omni_mdl_zib, "output/mdls/Omni_mdl_zib.rds")
  
# 5 Assess the Model ####
  #to load back the models
    
    Carni_mdl_zib <- read_rds("output/mdls/Carni_mdl_zib.rds")
    Omni_mdl_zib <- read_rds("output/mdls/Omni_mdl_zib.rds")
    Filter_mdl_zib <- read_rds("output/mdls/Filter_mdl_zib.rds")
    
  #5.0 list the selected models (double check the order of mdl_list and names)
    mdl_list <- list(Carni_mdl_zib, Omni_mdl_zib, Filter_mdl_zib)
    names(mdl_list) <- c("Carni","Omni","Filter")
    
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
  
# 6b Predict Global CPR ####
  
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
  
# 7 Plot visual summary of model ####
  
  #7.1 summary for model of omnivores
    plot_model_summary_omnivores(date)
    
  #7.2 summary for model of carnivores
    plot_model_summary_carnivores(date)
  
  #7.3 summary for model of gelatinous filter-feeders
    plot_model_summary_filterfeeders(date)
  
## End ##