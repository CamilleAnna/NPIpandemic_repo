# NPIsocialDilemma_repo
 
 This repository is for data, codes, simulations output and raw figures for "Social dilemmas in non-pharmaceutical interventions to tackle epidemics" (Camille Simonet, Amy Sweeny & Luke McNally, 2021)

 
## scripts
 
 - Statistical analyses
 	- script 1
	- script 2
	
- SIR models
	- dd
	- dd
	
- Game theory toy model
	-dd

## output

- **figures**: raw figures output, can be generated from scripts listed above, using output in ./output/model_runs

- **model_runs**:
	- SIMULATION_SIR_SD_June01_nbIni_60000_alpha_4.RData: simulations main model, social distancing intervention, with various initial frequency of compliance to generate (main heatmap, figure 3)
	- SIMULATION_SIR_TTI_June01_nbIni_60000_tau_5.1_ptest_0.7.RData: simulations main model, Test-Trace-Isolate intervention, with various initial frequency of compliance to generate (main heatmap, figure 4)
	- SD_various_alphas_nbIni_60000_timeSim_800.RData: runs of social distancing intervention with various C and alpha values (supp. figure heatmaps)
	- TTI_various_tau_ptest_0.7_nbIni_60000_timeSim_800.RData: runs of Test-Trace-Isolate intervention with various C and Tau values (supp. figure heatmaps)


## data

TBC ...

