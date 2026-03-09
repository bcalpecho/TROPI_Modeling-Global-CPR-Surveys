# ---
# title: Wrangle ESMs using R package 'hotrstuff'
# author: Bryan Alpecho
# date: 2025-
# output: ESM projection for Chl-a 
# ---

#Workflow adapted from: 
# https://github.com/SnBuenafe/hotrstuff.
# Citation: Buenafe K, Schoeman D, Everett J (2025). hotrstuff: Facilitates the rapid download, wrangling and processing of Earth System Model (ESM) output from the Coupled Model Intercomparison Project (CMIP). R package version 0.0.2, https://github.com/SnBuenafe/hotrstuff.

# install 'hotrstuff' package
devtools::install_github("SnBuenafe/hotrstuff")

#load packages
library(hotrstuff)

# define variable
var <- "chlos"

#List of ESMs used for generating an ensemble
## ACCESS-ESM1-5
## CanESM5
## CanESM5-1
## IPSL-CM6A-LR
## CESM2-WACCM
## CMCC-ESM2
## MPI-ESM1-2-LR
## MPI-ESM-2-HR
## NorESM2-LM
## NorESM2-MM

#prior step: download wget scripts from ESGF server

#downloading ESM outputs
  base_dir <- "/home/bcalp/UQ" 
  
  htr_download_ESM(
    hpc = NA,
    indir = file.path(base_dir, "data", "wget"), #directory for wget scripts
    outdir = file.path(base_dir, "data", "raw", var), #directory for raw ESM output
    quiet = TRUE,
    security = FALSE)
  
#merge files
  htr_merge_files(
    indir = file.path(base_dir, "data", "raw", var), # input directory
    outdir = file.path(base_dir, "data", "proc", "merged", var), # output directory
    year_start = 2015, # earliest year across all the scenarios considered (e.g., historical, ssp126, ssp245, ssp585)
    year_end = 2100 # latest year across all the scenarios considered
  )

#Adjust and reframe time periods
  htr_slice_period(
    indir = file.path(base_dir, "data", "proc", "merged", var), # input directory
    outdir = file.path(base_dir, "data", "proc", "sliced", var), # output directory
    freq = "Omon", # ocean, daily
    scenario = "ssp",
    year_start = 2015,
    year_end = 2100,
    overwrite = FALSE
  )

#Fix calendar periods (if needed)
  htr_fix_calendar(indir = file.path(base_dir, "data", "proc", "sliced", var)) # will be rewritten

# #Changing frequency of climate data
  htr_change_freq(
    freq = "yearly",
    indir = file.path(base_dir, "data", "proc", "sliced", var), # input directory
    outdir = file.path(base_dir, "data", "proc", "yearly", var)
  )

#Regridding
  htr_regrid_esm(
    indir = file.path(base_dir, "data", "proc", "yearly", var),
    outdir = file.path(base_dir, "data", "proc", "regridded", "yearly", var),
    cell_res = 1.0,
    layer = "annual"
  )

#Generate ensemble
  #ensemble by median
  htr_create_ensemble(
    indir = file.path(base_dir, "data", "proc", "regridded", "yearly", var), # input directory
    outdir = file.path(base_dir, "data", "proc", "ensemble", "median", var), # output directory
    model_list = c("ACCESS-ESM1-5","CanESM5","CanESM5-1","CESM2-WACCM", "CMCC-ESM2","IPSL-CM6A-LR","MPI-ESM1-2-HR","MPI-ESM1-2-LR","NorESM2-LM","NorESM2-MM"), # list of models for ensemble
    variable = var, # variable name
    freq = "Omon", # original frequency of data
    scenario = "ssp370", # scenario
    mean = FALSE # if false, takes the median
  )
  ensemble_model <- list.files(file.path(base_dir, "data", "proc", "ensemble", "median", var), full.names = TRUE)
  ensemble <- rast(ensemble_model)
  plot(ensemble$chlos_86)  #to plot ensemble at year 2100 (value of 'time' = 86)
  
##
  