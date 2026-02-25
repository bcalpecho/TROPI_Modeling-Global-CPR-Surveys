# Modeling Global CPR Surveys

Repository for MSc thesis project on global zooplankton community reslved to trophic groups across the satellite-derived chlorophyll-a gradient.
The project used the combined data of CPR Surveys of the AusCPR, the Atlantic CPR, the North Pacific CPR, and the Southern Ocean CPR. 

### Step-by-step preparation of data and modeling.

```/data_input``` contains the CPR abundance tables, trait table, and combined data of the CPR Surveys.
For the chlorophyll-a data, the netcdf files of OC-CCI can be downloaded through [OC-CCI](https://www.oceancolour.org/thredds/ncss/grid/CCI_ALL-v6.0-8DAY/dataset.html). 

* ```1_generate_traits.R``` assign trophic groups to provided zooplankton taxon list based on trait tables.
* ```2_extract_chla.R``` aggregates, extracts, and fill-up gaps of OC-CCI chlorophyll-a values.
* ```3_generate_completeDF.R``` finalizes the complete data frame comprised of zooplankton (relative abundance and assigned traits) and chlorophyll data.
* ```4_model_globalCPR.R``` fit the model and generates the model predictions.
* ```5_predict_globalCPR.R``` generate model predictions given the ensemble data.
* ```6_assess_models.R``` assess the quality of fit of GLMs and produce the supplementary figures for model assessment.
* ```7_plot_modelsummary.R``` produce visual summary of selected models (Figure 3, 4, and 5).
* ```0_wrangling_CPR.R``` additional functions for visualizing and preparing the zooplankton data. Includes creation of Figure 1 (Map of Global CPR).
* ```0_wrangling_ensemble.R``` additional functions for visualizing and summarizing the ensemble data.
  
### Acknowledgements

The Global Alliance of Continuous Plankton Recorders have provided the data of the AusCPR, the Atlantic CPR, the North Pacific CPR, and the Southern Ocean CPR surveys. Data for AusCPR were sourced from Australia’s Integrated Marine Observing System (IMOS) – IMOS is enabled by the National Collaborative Research Infrastructure Strategy (NCRIS). Data for the Atlantic CPR were sourced from the Marine Biological Association of the United Kingdom through the UK Archive for Marine Species and Habitats Data (DASSH). Data for the North Pacific CPR were sourced from the North Pacific Marine Science Organization through the DASSH. Data for the Southern Ocean CPR were sourced from the Scientific Committee on Antarctic Research (SCAR) sponsored Southern Ocean CPR (SO-CPR) Survey Database hosted by the Australian Antarctic Data Centre. 

We hope that this work provides insights into how zooplankton is likely to respond to ongoing and future accelerating climate-driven changes in phytoplankton abundance and cell size. 
If you have any questions or comments, please send a mail at bryan.alpecho@ulb.be. 
