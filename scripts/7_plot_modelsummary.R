# ---
# title: Plot visual summary of final models
# author: Bryan Alpecho
# date: 2025-
# output: Figure 3, 4, and 5
# ---

#setup
library(tidyverse)
library(patchwork)
library(ggpointdensity)
library(marginaleffects)
library(glmmTMB)

date <- "23032026"

plot_model_summary_omnivores <- function(date){
  
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
  
  #01 read in model
  load("Output/previousModels/revision/fMdl_omni.RData") 
  
  #FIGURE 3
  print("Figure 3 (omnivores) in preparation: predictions")
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
  print(paste("Saved plot: Output/illustrations/Omni_",date,".png",sep=""))
  ggsave(paste("Output/illustrations/Omni_",date,".png",sep=""), plot = final_patch,
         width = 8, height = 10, dpi = 300)
  
}
 
  
plot_model_summary_carnivores <- function(date){
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
  
  #01 read in model
  load("Output/previousModels/revision/fMdl_carni.RData") 
  
  #FIGURE 4
  print("Figure 4 (carnivores) in preparation: predicting outcomes")
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
  
  print(paste("Plot saved: Output/illustrations/Carni_",date,".png",sep=""))
  ggsave(paste("Output/illustrations/Carni_",date,".png",sep=""), plot = final_patch,
         width = 8, height = 10, dpi = 300)
  
}
  

plot_model_summary_filterfeeders <- function(date){
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
  
  #01 read in model
  load("Output/previousModels/revision/fMdl_filterfeeder.RData") 
  
  #FIGURE 5
  print("Figure 5 (gelatinous filter-feeders) in preparation: predicting outcomes")
  #C. Filter-feeders
  #Filter-feeders proportion vs. Chl-a
  pop_preds_FF_chla <- predictions(Filter_mdl_zib, 
                                   newdata = datagrid(chla_sqrt = seq(floor(min(df_noNA$chla_sqrt)),ceiling(max(df_noNA$chla_sqrt)),0.05)), 
                                   re.form = NA) # Zeros out random effects
  
  ff_plot_chla <- ggplot(data = pop_preds_FF_chla) + pub_theme + 
    geom_point(data = df, 
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
  
  print(paste("Plot saved: Output/illustrations/Filter_",date,".png",sep=""))
  ggsave(paste("Output/illustrations/Filter_",date,".png",sep=""), plot = final_patch,
         width = 8, height = 10, dpi = 300)
  
}  
  
plot_model_summary_omnivores(date)

plot_model_summary_carnivores(date)

plot_model_summary_filterfeeders(date)

##