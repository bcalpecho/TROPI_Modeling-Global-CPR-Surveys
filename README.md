# Modeling Global CPR Surveys

Repository for MSc thesis project on global zooplankton community resolved to trophic groups (omnivores, carnivores, and gelatinous filter-feeders) across the satellite-derived chlorophyll-a gradient. The project used the combined data of CPR Surveys of the Australian CPR, the Atlantic CPR, the North Pacific CPR, and the Southern Ocean CPR. 

### Data availability
Our findings are supported by open-access datasets. The data of Continuous Plankton Recorder surveys are available in their respective online repositories. 
* ```Australian CPR Survey``` through the Australian Ocean Data Network (https://portal.aodn.org.au/search?uuid=bf287dfe-9ce4-4969-9c59-51c39ea4d011%20on%2029/06/202)
* ```Atlantic CPR Survey``` through the Archive for Marine Species and Habitats Data (https://doi.mba.ac.uk/data/2962)
* ```North Pacific CPR Survey``` through the Archive for Marine Species and Habitats Data (https://doi.mba.ac.uk/data/3086)
* ```Southern Ocean CPR Survey``` through the Australian Antarctic Data Center (https://data.aad.gov.au/aadc/cpr/)

The ESA OC-CCI chl-a data can be accessed at https://www.oceancolour.org/thredds/ncss/grid/CCI_ALL-v6.0-8DAY/dataset.html.  
The CMIP6 data are archived and distributed by the Earth System Grid Federation (ESGF; https://aims2.llnl.gov/search). 

### Scripts.
The data preparation and analysis can be performed in two ways. 
* Through a single script ```runCompleteAnalysis.R``` apart from the preparation of ensemble (Step 6a);
* Through modular scripts (Step-by-step scripts 1 - 7). 

  * ```data_input``` contains the CPR abundance tables, trait table, and combined data of the CPR Surveys.
  * ```1_generate_traits.R``` assign trophic groups to provided zooplankton taxon list based on trait tables.
  * ```2_extract_chla.R``` aggregates, extracts, and fill-up gaps of OC-CCI chlorophyll-a values.
  * ```3_generate_completeDF.R``` finalizes the complete data frame comprised of zooplankton (relative abundance and assigned traits) and chlorophyll data.
  * ```4_model_globalCPR.R``` fit the model and generates the model predictions.
  * ```5_predict_globalCPR.R``` generate model predictions given the ensemble data.
  * ```6a_generate_ensemble.R``` generate a median ensemble of chlos projections from ten CMIP6 ESM models.
  * ```6b_assess_models.R``` assess the quality of fit of GLMs and produce the supplementary figures for model assessment.
  * ```7_plot_modelsummary.R``` produce visual summary of selected models (Figures 3, 4, and 5).
  * ```helper_CPR.R``` additional functions for visualizing and preparing the zooplankton data. Includes creation of Figure 1 (Map of Global CPR).
  * ```helper_ensemble.R``` additional functions for visualizing and summarizing the ensemble data.
  
### Acknowledgements

The Global Alliance of Continuous Plankton Recorders have provided the data of the AusCPR, the Atlantic CPR, the North Pacific CPR, and the Southern Ocean CPR surveys. Data for AusCPR were sourced from Australia’s Integrated Marine Observing System (IMOS) – IMOS is enabled by the National Collaborative Research Infrastructure Strategy (NCRIS). Data for the Atlantic CPR were sourced from the Marine Biological Association of the United Kingdom through the UK Archive for Marine Species and Habitats Data (DASSH). Data for the North Pacific CPR were sourced from the North Pacific Marine Science Organization through the DASSH. Data for the Southern Ocean CPR were sourced from the Scientific Committee on Antarctic Research (SCAR) sponsored Southern Ocean CPR (SO-CPR) Survey Database hosted by the Australian Antarctic Data Centre. 

We hope that this work provides insights into how zooplankton is likely to respond to ongoing and future accelerating climate-driven changes in phytoplankton abundance and cell size. 
If you have any questions or comments, please send a mail at bryan.alpecho@ulb.be. 
