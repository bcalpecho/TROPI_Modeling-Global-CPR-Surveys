# ---
# title: Assessing models
# author: Bryan Alpecho
# date: 2026-26-03
# output: QQ plot, residual plots, density plot of random effect intercepts and slope
# ---

library(DHARMa)
library(tidyverse)
library(viridis)
library(patchwork)
library(glmmTMB)
library(sp)
library(gstat)
library(gridGraphics)

#aims
# 01 to check for normality of model residuals using QQ plot
# 02 mean variance plot to assess homogeneity of variance
# 03 to plot the intercepts of Longhurst Provinces
# 04 to plot point density plots of Tow within Survey slope and intercept
# 05 to plot the residuals with and without the random effects

#setup
date <- "26032026"

# Create a publication-ready theme (Adapted from 2025 UQ MME Lab Winter R Workshop)
    pub_theme <- theme_classic(base_size = 11, base_family = "sans") + # Family including Arial, Helvetica, Futura, Verdana, and Calibri
      theme(
        # Axis styling
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 11, colour = "black"),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        
        # Legend styling
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
        legend.margin = margin(4, 4, 4, 4),
        
        # Panel styling
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        
        # Remove grid lines for cleaner look
        panel.grid = element_blank()
      )

  #read in models
    #load("Output/previousModels/revision/merged_chla_mdls.RData") 
    load("Output/previousModels/revision/final_zoop_mdls_21022026.RData")
    mdl_list <- list(Carni_mdl_zib, Omni_mdl_zib, Filter_mdl_zib)
    names(mdl_list) <- c("Carni","Omni","Filter")
  
  #01 check for normality of model residuals using QQ plot
    dharma <- function(model) {
      simulationOutput <- simulateResiduals(fittedModel = model, plot = FALSE)
      
      # Wrap the base R plot call in a formula using wrap_elements
      # This captures the plot so patchwork can treat it like an object
      wrap_elements(panel = ~plotQQunif(
        simulationOutput, 
        testUniformity = FALSE, 
        testOutliers = FALSE, 
        testDispersion = FALSE
      ))
    }
    
    # Note: Ensure your variable names match (FF_qq_plot vs Filter_qq_plot)
    FF_qq_plot    <- dharma(Filter_mdl_zib)
    Omni_qq_plot  <- dharma(Omni_mdl_zib)
    Carni_qq_plot <- dharma(Carni_mdl_zib)
    
    # Combine using patchwork
    design <- "ABC "
    
    qq_patch <- (Omni_qq_plot + Carni_qq_plot + FF_qq_plot) + 
      plot_layout(design = design) + 
      plot_annotation(
        tag_levels = "A",
        theme = theme(plot.title = element_text(size = 16, face = "bold"))
      )
    
    # View the result
    qq_patch
      
      #Save
    ggsave(paste("Output/plots/qq_plot/QQ2_",date,".png",sep=""), plot = qq_patch,
             width = 15, height = 6, dpi = 300)

#02 mean variance plot to assess homogeneity of variance
    
    plot_meanVariance <- function(mdl_list){
      
      load("data_input/global/global_df.RData")
      
      for(i in 1:length(mdl_list)){
        TG <- names(mdl_list[i])  
        print(paste0("Model in process: ",TG))
        
        if(TG == "Omni"){
          print("Success 1")
          df_Omni_noNA <- df %>% filter(!is.na(chla_sqrt)) %>% filter(!is.na(ROC_SVT_zib)) %>% filter(!is.na(tow_days))
          Omni_meanVariance <- ggplot(df_Omni_noNA) + aes(fitted(Omni_mdl_zib), resid(Omni_mdl_zib, type = "pearson")) +
            geom_hex(bins = 80) + geom_smooth(method = "loess") + pub_theme + theme(legend.position = c(.9, .8)) +
            scale_fill_viridis(begin = 0, end = .9, option = "C", limit = range(c(1,1200))) +
            theme(plot.title = element_text(hjust = 0.5), ) +
            labs(fill = "Frequency") + ylab("Residual") + xlab("Fitted Value")
          
        }else if(TG == "Carni"){
          print("Success 2")
          df_Carni_noNA <- df %>% filter(!is.na(chla_sqrt)) %>% filter(!is.na(RCO_SVT_zib)) %>% filter(!is.na(tow_days))
          Carni_meanVariance <- ggplot(df_Carni_noNA) + aes(fitted(Carni_mdl_zib), resid(Carni_mdl_zib, type = "pearson")) +
            geom_hex(bins = 80) + geom_smooth(method = "loess") + pub_theme + theme(legend.position = c(.9, .8)) +
            scale_fill_viridis(begin = 0, end = .9, option = "C", limit = range(c(1,1200))) +
            theme(plot.title = element_text(hjust = 0.5), ) +
            labs(fill = "Frequency") + ylab("Residual") + xlab("Fitted Value")
          
        }else if(TG == "Filter"){
          df_FF_noNA <- df_modified %>% filter(!is.na(chla_sqrt)) %>% filter(!is.na(RFF_SVT_zib)) %>% filter(!is.na(tow_days))
          Filter_meanVariance <- ggplot(df_FF_noNA) + aes(fitted(Filter_mdl_zib), resid(Filter_mdl_zib, type = "pearson")) +
            geom_hex(bins = 80) + geom_smooth(method = "loess") + pub_theme + theme(legend.position = c(.9, .8)) +
            scale_fill_viridis(begin = 0, end = .9, option = "C", limit = range(c(1,1200))) +
            theme(plot.title = element_text(hjust = 0.5), ) +
            labs(fill = "Frequency") + ylab("Residual") + xlab("Fitted Value")
          
          print("Success 3")
        }else{stop("Error in generating residual map")}
      }

      ## Residual plot
      # #to save individual plots
      # ggsave(paste("Output/illustrations/Filter_meanVariance_",date,".png",sep=""), plot = Filter_meanVariance,
      #        width = 9, height = 6, dpi = 300)
      #
      # ggsave(paste("Output/illustrations/Omni_meanVariance_",date,".png",sep=""), plot = Omni_meanVariance,
      #        width = 9, height = 6, dpi = 300)
      #
      # ggsave(paste("Output/illustrations/Carni_meanVariance_",date,".png",sep=""), plot = Carni_meanVariance,
      #        width = 9, height = 6, dpi = 300)
      # #

      design <- "
                A
                B
                C
              "

      (meanVariance_patch <- Filter_meanVariance + Omni_meanVariance + Carni_meanVariance +
         plot_layout(design = design) +
         plot_annotation(
           tag_levels = "A",
           theme = theme(plot.title = element_text(size = 16, face = "bold"))
         )
      )
      ggsave(paste("Output/illustrations/meanVariance_",date,".png",sep=""), plot = meanVariance_patch,
             width = 8, height = 10, dpi = 300)
      print(paste0("Combined plot saved: meanVariance_",date,".png"))
  }
    plot_meanVariance(mdl_list)
    
#03 to plot the intercepts of Longhurst Provinces
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
        guides(size=FALSE,shape=FALSE) + 
        theme(axis.text.x=element_text(size=10), 
              axis.title.x=element_text(size=13),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.x = element_blank()) +
        coord_flip() + 
        labs(y = "Intercept", x = "Longhurst Provinces")
      
      #annotate("text", y = -2.3, x = 28, label = "(D)", size = 5)
      
      ggsave(paste("Output/predictions/global/randomEffects/",TG,"_Longhurst-",date,".png",sep=""), plot = re_lh_plot,
             width = 5, height = 6, dpi = 300)
      print(paste0("Output: ",TG,"_Longhurst-",date,".png"))
    }
  }
  
  plot_Longhurst(mdl_list)

#04 to plot point density plots of Tow within Survey slope and intercept

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
          axis.text = element_text(size = 14)    # Adjust the size for axis numbers/text
        )
      
      #update file name
      ggsave(paste("Output/predictions/global/randomEffects/",TG,"_Tow-Survey-",date,".png",sep=""), plot = re_tow_plot_point_density,
             width = 7, height = 8, dpi = 300)
        
    }  
  }
      
  plot_TowSlopeAndIntercept(mdl_list)
    
#05 to plot the residuals with and without the random effects
    
  plot_residuals <- function(model_list){
  for(i in 1:length(model_list)){
    vc <- VarCorr(model_list[[i]])
    TG <- names(model_list[i])
    print(TG) 
    
    ifelse(TG == "Filter", 
           df_noNA <- df_modified %>% filter(!is.na(chla_sqrt)) %>% filter(!is.na(RFF_SVT_zib)) %>% filter(!is.na(tow_days)),
                      ifelse(TG == "Carni", 
                            df_noNA <- df %>% filter(!is.na(chla_sqrt)) %>% filter(!is.na(RCO_SVT_zib)) %>% filter(!is.na(tow_days)), 
                            df_noNA <- df %>% filter(!is.na(chla_sqrt)) %>% filter(!is.na(ROC_SVT_zib)) %>% filter(!is.na(tow_days))))
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
    
    ggsave(paste("Output/map/residuals_map/",TG,"_zoop_withRandom_",date,".png",sep=""), plot = residuals_withRandom,
           width = 9, height = 6, dpi = 300)
    print(paste("File saved: Output/map/residuals_map/",TG,"_zoop_withRandom_",date,".png",sep=""))
    
    # Model without random effects 
    ifelse(TG == "Filter", 
           Mdl_noRandom <- glmmTMB(RFF_SVT_zib ~ chla_sqrt + survey, 
                                          ziformula = ~1,
                                          data = df_noNA, family = beta_family(link = "logit")),
           ifelse(TG == "Carni", 
                  Mdl_noRandom <- glmmTMB(RCO_SVT_zib ~ chla_sqrt + survey, 
                                          ziformula = ~1,
                                          data = df_noNA, family = beta_family(link = "logit")), 
                  Mdl_noRandom <- glmmTMB(ROC_SVT_zib ~ chla_sqrt + survey, 
                                          ziformula = ~1,
                                          data = df_noNA, family = beta_family(link = "logit"))))
    
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
    
    #export plot of residuals without random effects
    ggsave(paste("Output/map/residuals_map/",TG,"_zoop_noRandom_",date,".png",sep=""), plot = residuals_noRandom,
           width = 9, height = 6, dpi = 300)
    print(paste("File saved: Output/map/residuals_map/",TG,"_zoop_noRandom_",date,".png",sep=""))
  
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
      ggsave(paste("Output/map/residuals_map/residualMap_",TG,"_",date,".png",sep=""), plot = residualMap,
             width = 8, height = 10, dpi = 300)
      print(paste("File saved: residualMap_",TG,"_",date,".png",sep=""))
      
    }
  }
  
  plot_residuals(mdl_list)

# ###Export selected models
#     write_rds(Filter_mdl_zib, "Output/previousModels/mdl_export/Filter_mdl_zib.rds")
#     write_rds(Carni_mdl_zib, "Output/previousModels/mdl_export/Carni_mdl_zib.rds")
#     write_rds(Omni_mdl_zib, "Output/previousModels/mdl_export/Omni_mdl_zib.rds")
  
## End ##
    