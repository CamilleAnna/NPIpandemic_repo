# NPIsocialDilemma_repo
 
 This repository is for data, codes, simulations output and raw figures for "Social dilemmas in non-pharmaceutical interventions to tackle epidemics" (Camille Simonet, Amy Sweeny & Luke McNally, 2021)

 
## scripts

All statistical analyses and simulations were run in R version 4.0.1

 - **Statistical analyses**
 	- **1_mobilityAnalysis.R**: mobility data analysis and codes producing manuscript figures
	- **CleanMCMC.R**: utility function sourced in above script
	
- **SIR models**
	- **2_mathematical_model.R**: main model simulations
	- **3_figures_main_heatmap.R**: make main text heatmap figures form simulations output
	- **4_figures_supplementary_heatmaps.R**: supplementary heatmap figures (run simulations + make figures)
	- **5_mathematic_model_and_figures_dynamicDecisionMaking.R**: simulations and figure script for dynamic models (main text and supplements)

- **Game theory toy model**
	- **6_gameTheory_regionPlots_BOX1.nb**: game theory toy model (Mathematic notebook, a pdf version is also provided)
	

## output

- **figures**: raw figures output, can be generated from scripts listed above, using output in ./output/model_runs

- **model_runs**:
	- **SIMULATION_SIR_SD_June01_nbIni_60000_alpha_4.RData**: simulations main model, social distancing intervention, with various initial frequency of compliance to generate (main heatmap, figure 3)
	- **SIMULATION_SIR_TTI_June01_nbIni_60000_tau_5.1_ptest_0.7.RData**: simulations main model, Test-Trace-Isolate intervention, with various initial frequency of compliance to generate (main heatmap, figure 4)
	- **SD_various_alphas_nbIni_60000_timeSim_800.RData**: runs of social distancing intervention with various C and alpha values (supp. figure heatmaps)
	- **TTI_various_tau_ptest_0.7_nbIni_60000_timeSim_800.RDat**a: runs of Test-Trace-Isolate intervention with various C and Tau values (supp. figure heatmaps)

- **mobility_analysis**: summary of statistical model for mobility analysis + R object file of MCMCglmm model output ( all outputs of script "1_mobilityAnalysis.R")

## data

data sources in mobility analysis script.

- **countryList_csv.csv**: countries codes list to use with "ggflags" in codes producing figure 1
- **MobilityInterventionsDF_MS.csv**:  mobility dataset used in mobility analysis
Interventions policies and times loaded from shared drive at: https://docs.google.com/spreadsheets/d/1LKTv_4-aC1I_Ld2HZAWs7yBkOb9fBELIV8bIEdDTykc/edit?usp=sharing

