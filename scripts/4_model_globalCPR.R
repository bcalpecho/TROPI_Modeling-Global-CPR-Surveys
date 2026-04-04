# ---
# title: Modeling (Global CPR)
# author: Bryan Alpecho
# date: 2026-
# output: glmmTMB models, tibble for summary of coefficients
# ---

#### Load libraries and data ####
  packages <- c("tidyverse",
                "glmmTMB",
                "marginaleffects",
                "MuMIn",
                "GGally",
                "DHARMa",
                "reprex",
                "viridis",
                "sp",
                "gstat")

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
  remove(packages)

#set date for version control
  date <- "21032026" #07112025 prior to 26022026
  df_date <- "21032026" #data frame version (verify with latest df in "data input" folder)

# Import data
  df <- read_rds(paste("data_input/global/df_complete_",df_date,".rds",sep=""))
  
  #set ceiling for chl-a values
  df <- df %>%
    mutate(chla_withCeilingAt3 = case_when(chla > 3 ~ 3,
                                            .default = chla)) 
  #to remove NA values,
  df_noNA <- df %>% filter(!is.na(chla)) 
  
  
## Exploratory data analysis
  # Look at spearman correlations and distributions
    ggpairs_plot <- GGally::ggpairs(df, columns = c("ROC","RCO","RFF","chla","chla_sqrt"), progress = FALSE, 
                                    upper = list(continuous = GGally::wrap("cor", method = "spearman", stars = FALSE)),
                                    columnLabels = c("Relative abundance of\nomnivores","Relative abundance of\ncarnivores","Relative abundance of\ngelatinous filter-feeders","Chlorophyll-a","Sqrt-transformed\nchlorophyll-a"))
    ggsave(paste("Output/data/distributions/ggally-",date,".png",sep=""), plot = ggpairs_plot, 
           width = 8, height = 8, dpi = 300)
    
# Run competing models 
  #zero-inflated beta regression model 
    #Predictors: Sqrt transformed Chl-a and Survey 
    #Note: SVT transformation for 1s only; non-transformed [0,1) 
    Filter_mdl_zib <- glmmTMB(RFF_SVT_zib ~ chla_sqrt + survey + (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                          ziformula = ~1,
                          data = df_modified, family = beta_family(link = "logit"))
    
    Omni_mdl_zib <- glmmTMB(ROC_SVT_zib ~ chla_sqrt + survey + (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                        ziformula = ~1,
                        data = df, family = beta_family(link = "logit"))
    
    Carni_mdl_zib <- glmmTMB(RCO_SVT_zib ~ chla_sqrt + survey + (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                         ziformula = ~1,
                         data = df, family = beta_family(link = "logit"))
    
    ###
    #OTHER MODELS
        #Models with chlorophyll-a only
        #Filter chl-a only (non-transformed, sqrt transformed, with Ceiling at 3)
        Filter_mdl_zib_chla_sqrt <- glmmTMB(RFF_SVT_zib ~ chla_sqrt + (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                                  ziformula = ~1,
                                  data = df_modified, family = beta_family(link = "logit"))
 
        Filter_mdl_zib_chla <- glmmTMB(RFF_SVT_zib ~ chla + (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                                       ziformula = ~1,
                                       data = df_modified, family = beta_family(link = "logit"))
        
        Filter_mdl_zib_chla_withCeilingAt3 <- glmmTMB(RFF_SVT_zib ~ chla_withCeilingAt3 + (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                                       ziformula = ~1,
                                       data = df_modified, family = beta_family(link = "logit"))
        
        #Omni chl-a only (non-transformed, sqrt transformed, with Ceiling at 3)
        Omni_mdl_zib_chla <- glmmTMB(ROC_SVT_zib ~ chla + (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                                ziformula = ~1,
                                data = df, family = beta_family(link = "logit"))
        
        Omni_mdl_zib_chla_sqrt <- glmmTMB(ROC_SVT_zib ~ chla_sqrt + (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                                     ziformula = ~1,
                                     data = df, family = beta_family(link = "logit"))
        
        Omni_mdl_zib_chla_withCeilingAt3 <- glmmTMB(ROC_SVT_zib ~ chla_withCeilingAt3 + (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                                     ziformula = ~1,
                                     data = df, family = beta_family(link = "logit"))
        
        #Carni chl-a only (non-transformed, sqrt transformed, with Ceiling at 3)
        Carni_mdl_zib_chla <- glmmTMB(RCO_SVT_zib ~ chla + (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                                 ziformula = ~1,
                                 data = df, family = beta_family(link = "logit"))
        
        Carni_mdl_zib_chla_sqrt <- glmmTMB(RCO_SVT_zib ~ chla_sqrt + (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                                      ziformula = ~1,
                                      data = df, family = beta_family(link = "logit"))
        
        Carni_mdl_zib_chla_withCeilingAt3 <- glmmTMB(RCO_SVT_zib ~ chla_withCeilingAt3 + (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                                      ziformula = ~1,
                                      data = df, family = beta_family(link = "logit"))
        
        #Models with latitude
        Filter_mdl_zib_latitude <- glmmTMB(RFF_SVT_zib ~ chla_sqrt + latitude_Abs + (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                                       ziformula = ~1,
                                       data = df, family = beta_family(link = "logit"))
        Omni_mdl_zib_latitude <- glmmTMB(ROC_SVT_zib ~ chla_sqrt + latitude_Abs + (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                                           ziformula = ~1,
                                           data = df, family = beta_family(link = "logit"))
        Carni_mdl_zib_latitude <- glmmTMB(RCO_SVT_zib ~ chla_sqrt + latitude_Abs + (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                                           ziformula = ~1,
                                           data = df, family = beta_family(link = "logit"))
        #Models with survey and latitude
        Filter_mdl_zib_surveyLat <- glmmTMB(RFF_SVT_zib ~ chla_sqrt + survey + latitude_Abs +  (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                                           ziformula = ~1,
                                           data = df, family = beta_family(link = "logit"))
        
        Carni_mdl_zib_surveyLat <- glmmTMB(RCO_SVT_zib ~ chla_sqrt + survey + latitude_Abs +  (1 + tow_days | survey: tow_no) + (1 | longhurst),  
                                            ziformula = ~1,
                                            data = df, family = beta_family(link = "logit"))
        
        Omni_mdl_zib_surveyLat <- glmmTMB(ROC_SVT_zib ~ chla_sqrt + survey + latitude_Abs +  (1 + tow_days | survey: tow_no) + (1 | longhurst), 
                                            ziformula = ~1,
                                            data = df, family = beta_family(link = "logit"))
        

#Visualize the predictions    
    ##prior step
        #pub_theme from UQ MME Lab R Workshop Winter 2025
        # Create a publication-ready theme
        pub_theme <- theme_classic(base_size = 10, base_family = "sans") + # Family including Arial, Helvetica, Futura, Verdana, and Calibri
          theme(
            # Axis styling
            axis.title = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12, colour = "black"),
            axis.line = element_line(colour = "black", linewidth = 0.5),
            axis.ticks = element_line(colour = "black", linewidth = 0.5),
            
            # Legend styling
            legend.title = element_text(size = 11, face = "bold"),
            legend.text = element_text(size = 10),
            legend.position = "inside",
            legend.position.inside = c(0.9, 0.85),
            legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
            legend.margin = margin(4, 4, 4, 4),
            
            # Panel styling
            panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white"),
            
            # Remove grid lines for cleaner look
            panel.grid = element_blank()
          )    
    ##
        
    #Omnivores proportion vs. Chl-a
        pop_preds_omni_chla <- predictions(Omni_mdl_zib, 
                                     newdata = datagrid(chla_sqrt = seq(floor(min(df_noNA$chla_sqrt)),ceiling(max(df_noNA$chla_sqrt)),0.05)), 
                                     re.form = NA) # Zeros out random effects
          
        omni_plot_chla <- ggplot(data = pop_preds_omni_chla) + pub_theme +
            geom_point(data = df_noNA, 
                       aes(x = chla_sqrt, y = ROC_SVT_zib), alpha = 0.1) +
            geom_ribbon(aes(x = chla_sqrt, y = estimate, 
                            ymin = conf.low, ymax = conf.high), alpha = 0.3, colour = "blue", linetype=2) + 
            geom_line(aes(x = chla_sqrt, y = estimate),colour="blue") + 
          labs(y = "Omnivores proportion", x = expression(bold(sqrt("chl-a")))) 
          #annotate("text", y = 0.8, x = 3.8, label = expression("Estimate= 0.055; R"^{2}*"= 0.26; p-value=0.006 "), size = 5) 
          
        ggsave(paste("Output/predictions/global/Omnivore-prop-zib_ChlaPlot_",date,".png",sep=""), plot = omni_plot_chla, 
                 width = 8, height = 5, dpi = 300)
    
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
          geom_line(aes(x = chla_sqrt, y = estimate, colour = survey)) + 
          labs(y = "Omnivores proportion", x = expression(bold(sqrt("chl-a")))) 
        
        ggsave(paste("Output/predictions/global/Omnivore-prop-zib_SurveyPlot_",date,".png",sep=""), plot = omni_plot_surveys, 
               width = 8, height = 5, dpi = 300)  
          
    #Carnivores proportion vs. Chl-a
          pop_preds_carni_chla <- predictions(Carni_mdl_zib, 
                                       newdata = datagrid(chla_sqrt = seq(floor(min(df_noNA$chla_sqrt)),ceiling(max(df_noNA$chla_sqrt)),0.05)), 
                                       re.form = NA) # Zeros out random effects
          
          carni_chla_plot <- ggplot(data = pop_preds_carni_chla) + pub_theme + 
            geom_point(data = df_noNA, 
                       aes(x = chla_sqrt, y = RCO_SVT_zib), alpha = 0.1) +
            geom_ribbon(aes(x = chla_sqrt, y = estimate, 
                            ymin = conf.low, ymax = conf.high), alpha = 0.3, colour = "blue", linetype = 2) + 
            geom_line(aes(x = chla_sqrt, y = estimate),colour="blue") + 
            labs(y = "Carnivores proportion", x = expression(bold(sqrt("chl-a")))) 
            #annotate("text", y = 0.6, x = 3.8, label = expression("Estimate= -0.22; R"^{2}*"= 0.20; p-value<1e-21 "), size = 5) 
          
          ggsave(paste("Output/predictions/global/Carnivore-prop-zib_ChlaPlot_",date,".png",sep=""), plot = carni_chla_plot, 
                 width = 8, height = 5, dpi = 300)
            
    #Carnivores proportion vs. Surveys
          pop_preds_carni_surveys <- predictions(Carni_mdl_zib, 
                                                 newdata = datagrid(survey = unique(df_noNA$survey), 
                                                                    chla_sqrt = seq(floor(min(df_noNA$chla_sqrt)),ceiling(max(df_noNA$chla_sqrt)),0.05)), 
                                                 re.form = NA) # Zeros out random effects
          
          carni_survey_plot <- ggplot(data = pop_preds_carni_surveys) + pub_theme +
            geom_point(data = df_noNA, 
                       aes(x = chla_sqrt, y = RCO_SVT_zib), alpha = 0.1) +
            geom_ribbon(aes(x = chla_sqrt, y = estimate, 
                            ymin = conf.low, ymax = conf.high, 
                            colour = survey, fill = survey), alpha = 0.3, colour = NA, linetype = 2) + 
            geom_line(aes(x = chla_sqrt, y = estimate, colour = survey)) + 
            labs(y = "Carnivores proportion", x = expression(bold(sqrt("chl-a"))))
          
          ggsave(paste("Output/predictions/global/Carnivore-prop-zib_SurveyPlot_",date,".png",sep=""), plot = carni_survey_plot, 
                 width = 8, height = 5, dpi = 300)
    
    #Filter-feeders proportion vs. Chl-a
    pop_preds_FF_chla <- predictions(Filter_mdl_zib, 
                                     newdata = datagrid(chla_sqrt = seq(floor(min(df_noNA$chla_sqrt)),ceiling(max(df_noNA$chla_sqrt)),0.05)), 
                                     re.form = NA) # Zeros out random effects
    
    ff_chla_plot <- ggplot(data = pop_preds_FF_chla) + pub_theme + 
      geom_point(data = df, 
                 aes(x = chla_sqrt, y = RFF_SVT_zib), alpha = 0.1) +
      geom_ribbon(aes(x = chla_sqrt, y = estimate, 
                      ymin = conf.low, ymax = conf.high), alpha = 0.3, colour = "blue", linetype=2) + 
      geom_line(aes(x = chla_sqrt, y = estimate),colour="blue") + 
      labs(y = "Filter-feeders proportion", x = expression(bold(sqrt("chl-a")))) 
    #annotate("text", y = 0.95, x = 3.5, label = expression("Estimate= -0.49; R"^{2}*"= 0.41; p-value<1e-15 "), size = 5)  
    
    ggsave(paste("Output/predictions/global/Filter-feeder-prop-zib_ChlaPlot_",date,".png",sep=""), plot = ff_chla_plot, 
           width = 8, height = 5, dpi = 300)
    
    #Filter-feeders proportion vs. Surveys
    pop_preds_FF_surveys <- predictions(Filter_mdl_zib, 
                                        newdata = datagrid(survey = unique(df_noNA$survey), 
                                                           chla_sqrt = seq(floor(min(df_noNA$chla_sqrt)),ceiling(max(df_noNA$chla_sqrt)),0.05)), 
                                        re.form = NA) # Zeros out random effects
    
    ff_survey_plot <- ggplot(data = pop_preds_FF_surveys) + pub_theme +
      geom_point(data = df_noNA, 
                 aes(x = chla_sqrt, y = RFF_SVT_zib), alpha = 0.1) +
      geom_ribbon(aes(x = chla_sqrt, y = estimate, 
                      ymin = conf.low, ymax = conf.high, 
                      colour = survey, fill = survey), alpha = 0.3, colour = NA) + 
      geom_line(aes(x = chla_sqrt, y = estimate, colour = survey)) + 
      labs(y = "Filter-feeders proportion", x = expression(bold(sqrt("chl-a")))) 
    
    ggsave(paste("Output/predictions/global/Filter-feeder-prop-zib_SurveyPlot_",date,".png",sep=""), plot = ff_survey_plot, 
           width = 8, height = 5, dpi = 300)            

#Assess model fit
    #Filter_mdls
      Filter_mdl_zib_chla_R2 <- MuMIn::r.squaredGLMM(Filter_mdl_zib_chla)
      Filter_mdl_zib_chla_sqrt_R2 <- MuMIn::r.squaredGLMM(Filter_mdl_zib_chla_sqrt)
      Filter_mdl_zib_chla_withCeilingAt3_R2 <- MuMIn::r.squaredGLMM(Filter_mdl_zib_chla_withCeilingAt3)
      Filter_mdl_zib_R2 <- MuMIn::r.squaredGLMM(Filter_mdl_zib)
      Filter_mdl_zib_surveyLat_R2 <- MuMIn::r.squaredGLMM(Filter_mdl_zib_surveyLat)
      Filter_mdl_zib_Lat_R2 <- MuMIn::r.squaredGLMM(Filter_mdl_zib_latitude)
      Filter_mdl_zib_coeff <- summary(Filter_mdl_zib)$coefficients$cond %>% as.data.frame()
      
    #Omni_mdls
      Omni_mdl_zib_chla_R2 <- MuMIn::r.squaredGLMM(Omni_mdl_zib_chla)
      Omni_mdl_zib_chla_sqrt_R2 <- MuMIn::r.squaredGLMM(Omni_mdl_zib_chla_sqrt)
      Omni_mdl_zib_chla_withCeilingAt3_R2 <- MuMIn::r.squaredGLMM(Omni_mdl_zib_chla_withCeilingAt3)
      Omni_mdl_zib_R2 <- MuMIn::r.squaredGLMM(Omni_mdl_zib)
      Omni_mdl_zib_surveyLat_R2 <- MuMIn::r.squaredGLMM(Omni_mdl_zib_surveyLat)
      Omni_mdl_zib_Lat_R2 <- MuMIn::r.squaredGLMM(Omni_mdl_zib_latitude)
      Omni_mdl_zib_coeff <- summary(Omni_mdl_zib)$coefficients$cond %>% as.data.frame()
    
    #Carni mdls
      Carni_mdl_zib_chla_R2 <- MuMIn::r.squaredGLMM(Carni_mdl_zib_chla)
      Carni_mdl_zib_chla_sqrt_R2 <- MuMIn::r.squaredGLMM(Carni_mdl_zib_chla_sqrt)
      Carni_mdl_zib_chla_withCeilingAt3_R2 <- MuMIn::r.squaredGLMM(Carni_mdl_zib_chla_withCeilingAt3)
      Carni_mdl_zib_R2 <- MuMIn::r.squaredGLMM(Carni_mdl_zib)
      Carni_mdl_zib_surveyLat_R2 <- MuMIn::r.squaredGLMM(Carni_mdl_zib_surveyLat)
      Carni_mdl_zib_Lat_R2 <- MuMIn::r.squaredGLMM(Carni_mdl_zib_latitude)
      Carni_mdl_zib_coeff <- summary(Carni_mdl_zib)$coefficients$cond %>% as.data.frame()
                            
## SUMMARY OF FINAL MODELS
    #R2
      Omni_mdl_zib_R2 <- MuMIn::r.squaredGLMM(Omni_mdl_zib)
      Filter_mdl_zib_R2 <- MuMIn::r.squaredGLMM(Filter_mdl_zib)
      Carni_mdl_zib_R2 <- MuMIn::r.squaredGLMM(Carni_mdl_zib)
      
    #Coefficients
      Filter_mdl_zib_coeff <- summary(Filter_mdl_zib)$coefficients$cond %>% as.data.frame()
      Carni_mdl_zib_coeff<- summary(Carni_mdl_zib)$coefficients$cond %>% as.data.frame()
      Omni_mdl_zib_coeff <- summary(Omni_mdl_zib)$coefficients$cond %>% as.data.frame()
      
    tibble_models_zib <- tibble(
      mdl = c("Filter_mdl_zib", "Omni_mdl_zib", "Carni_mdl_zib"),
      Estimate = c(Filter_mdl_zib_coeff$Estimate[2],  Omni_mdl_zib_coeff$Estimate[2], Carni_mdl_zib_coeff$Estimate[2]),
      StdError = c(Filter_mdl_zib_coeff$`Std. Error`[2],  Omni_mdl_zib_coeff$`Std. Error`[2], Carni_mdl_zib_coeff$`Std. Error`[2]),
      z_value = c(Filter_mdl_zib_coeff$`z value`[2], Omni_mdl_zib_coeff$`z value`[2], Carni_mdl_zib_coeff$`z value`[2]),
      p_value = c(Filter_mdl_zib_coeff$`Pr(>|z|)`[2], Omni_mdl_zib_coeff$`Pr(>|z|)`[2], Carni_mdl_zib_coeff$`Pr(>|z|)`[2]),
      PseudoR2_marginal = c(Filter_mdl_zib_R2[1], Omni_mdl_zib_R2[1], Carni_mdl_zib_R2[1]),
      PseudoR2_conditional = c(Filter_mdl_zib_R2[2], Omni_mdl_zib_R2[2], Carni_mdl_zib_R2[2]))
    
   #to compute for variance % contribution by random effects
    
    mdl_list <- list(Filter_mdl_zib, Omni_mdl_zib, Carni_mdl_zib)
    
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

    summary_mdlVariance(mdl_list)    
    
    #Omni
      Omni_summary <- summary(Omni_mdl_zib)
      #Tow within Survey
      var_survey <- Omni_summary$varcor$cond$`survey:tow_no`[1] + Omni_summary$varcor$cond$`survey:tow_no`[4] #gets the variance 
      #Longhurst Provinces
      var_lh <- Omni_summary$varcor$cond$`longhurst`[1]    
      
      #contribution of each RE in relative terms
      #Tow within Survey
      rel_var_survey <- (var_survey)/(var_survey + var_lh) * (Omni_mdl_zib_R2[2] - Omni_mdl_zib_R2[1])
      #Longhurst Provinces
      rel_var_LH <- (var_lh)/(var_survey + var_lh) * (Omni_mdl_zib_R2[2] - Omni_mdl_zib_R2[1])
    
    
    #Carni
      Carni_summary <- summary(Carni_mdl_zib)
      #Tow within Survey
      var_survey <- Carni_summary$varcor$cond$`survey:tow_no`[1] + Carni_summary$varcor$cond$`survey:tow_no`[4] #gets the variance 
      #Longhurst Provinces
      var_lh <- Carni_summary$varcor$cond$`longhurst`[1]    
      
      #contribution of each RE in relative terms
      #Tow within Survey
      rel_var_survey <- (var_survey)/(var_survey + var_lh) * (Carni_mdl_zib_R2[2] - Carni_mdl_zib_R2[1])
      #Longhurst Provinces
      rel_var_LH <- (var_lh)/(var_survey + var_lh) * (Carni_mdl_zib_R2[2] - Carni_mdl_zib_R2[1])
    
    #Filter
      Filter_summary <- summary(Filter_mdl_zib)
      #Tow within Survey
      var_survey <- Filter_summary$varcor$cond$`survey:tow_no`[1] + Filter_summary$varcor$cond$`survey:tow_no`[4] #gets the variance 
      #Longhurst Provinces
      var_lh <- Filter_summary$varcor$cond$`longhurst`[1]    
      
      #contribution of each RE in relative terms
      #Tow within Survey
      rel_var_survey <- (var_survey)/(var_survey + var_lh) * (Filter_mdl_zib_R2[2] - Filter_mdl_zib_R2[1])
      #Longhurst Provinces
      rel_var_LH <- (var_lh)/(var_survey + var_lh) * (Filter_mdl_zib_R2[2] - Filter_mdl_zib_R2[1])
      