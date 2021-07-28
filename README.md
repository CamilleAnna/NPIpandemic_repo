# NPIsocialDilemma_repo
 
 This repository is for data, codes, simulations output and raw figures for "Social dilemmas in non-pharmaceutical interventions to tackle epidemics" (Camille Simonet, Amy Sweeny & Luke McNally, 2021)

 
## scripts
 
 - **Statistical analyses**
 	- TBC
	
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


## data

TBC ...

