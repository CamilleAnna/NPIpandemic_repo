# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#   Social dilemmas in non-pharmaceutical interventions to tackle epidemics
#                     Simonet C, Sweeny A, McNally L (2021)
#                       Codes for mathematical models
#
#                       last edited: July, 27th, 2021
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Script summary
# [1] defines models for each intervention in a single code wrapper, define parameters
# [2] run SD model on range of initial frequency of compliance
# [3] run TTI model on range of initial frequency of compliance


setwd('path/to/repo') # setwd('~/Documents/PhD/Research/Covid19/NPIpandemic_repo/')
library(readxl); library(ggplot2); library(dplyr); library(tidyr);library(Rmisc); library(lubridate);library(plotly);library(yarrr); library(knitr);library(plyr);library(utils);library(tidyverse);library(gridExtra);library(deSolve); library(stringi); library(stringr);
library(ggthemes)


# [1] SETUP ----

# ... (a) Model  ----

# Code wrapper to run model

# Single wrapper for all possible scenarios:
#   npi = 'sd' runs social distancing scenario
#   npi = 'tti' runs Test=Trace-Isolate scenario

# null model (no intervention), can be run by setting appropriate parameter values:
#   if running under 'sd', setting a = 1 will run the null model
#   if running under 'tti', setting either tau to large value or p_test to 0 will run the null model

# Initial conditions are compatible with both sd and tti scenario, (so has all compartments, including 'test' and 'trace' compartments, but change rates for those is 0 if running 'sd'

sir_npi <- function(time, state, parameters, npi) {
  
  with(as.list(c(state, parameters)), {
    
    if(npi == 'sd'){
      
      if(x > 0){
        # BENEFITS: probability of not getting infected
        p_not_infected_c = (1 - (1-exp(-b*time_sim*(In + (1/a)*Ic)*(1/a)))) 
        p_not_infected_n = (1 - (1-exp(-b*time_sim*(In + (1/a)*Ic)*(1))))
        
        
        # COSTS: prop of next n days not recovered
        prop_time_not_recovered_c = ifelse(((1/(b*time_sim*(1/a)*(In + (1/a)*Ic)))+(1/g)) < time_sim,
                                           ((1/(b*time_sim*(1/a)*(In + (1/a)*Ic)))+(1/g))/time_sim,
                                           1)
        
        prop_time_not_recovered_n = ifelse(((1/(b*time_sim*(1)*(In + (1/a)*Ic)))+(1/g)) < time_sim,
                                           ((1/(b*time_sim*(1)*(In + (1/a)*Ic)))+(1/g))/time_sim,
                                           1)
        
        # "EFFECTIVE" COST: I use a parameter ip to "activate" the immunity passport
        # If ip = 0, no immunity passport, cost is C or E
        # if ip = 1, activate immunity passport, multiply by prop_time_not_recovered
        
        C_eff<- ifelse(ip == 0, C, C * prop_time_not_recovered_c)
        E_eff<- ifelse(ip == 0, E, E * prop_time_not_recovered_n)
        
        # PAYOFF = benefit - cost
        pi_n = p_not_infected_n - E_eff
        pi_c = p_not_infected_c - C_eff
      }else{
        pi_n = 0
        pi_c = 0
      }
      
      # Current fraction of compliers among the susceptible
      #fc = Sc/(Sc+Sn)  
      
      # ODE, with replicator dynamic for the S "activated" if x > 1
      
      dSn = -b*(In + (1/a)*Ic)*Sn + x*(Sn/(Sn+Sc)) * (pi_n - ((pi_n *(Sn/(Sn+Sc))) + (pi_c *(Sc/(Sn+Sc)))))
      dSc = -b*(In + (1/a)*Ic)*(1/a)*Sc + x*(Sc/(Sn+Sc)) * (pi_c - ((pi_n *(Sn/(Sn+Sc))) + (pi_c *(Sc/(Sn+Sc)))))
      
      dIn = b*(In + (1/a)*Ic)*Sn - g * In
      dIn_test = 0
      dIc = b*(In + (1/a)*Ic)*(1/a)*Sc - g * Ic
      dIc_test = 0
      dIc_trace = 0
      
      dRn = g * In
      dRc = g * Ic
      
    }else if(npi == 'tti'){
      
      qtest  = (g*tau)/(1 + (g*tau))
      qtrace = (g*tau)/(2 + (g*tau))
      
      if(x > 0){
        # BENEFIT - coded explicitly but in fact in TTI, it is the same for both compliers and non-compliers
        p_not_infected_c = (1 - (1 - exp(-b*time_sim*(In + Ic  + qtest*Ic_test + qtest*In_test + qtrace*Ic_trace))))
        p_not_infected_n = (1 - (1 - exp(-b*time_sim*(In + Ic  + qtest*Ic_test + qtest*In_test + qtrace*Ic_trace))))
        
        
        # COSTS - coded explicitly but in fact in TTI, it is the same for both compliers and non-compliers
        prop_time_not_recovered_c = ifelse(((1/(b*time_sim*(In + Ic  + qtest*Ic_test + qtest*In_test + qtrace*Ic_trace)))+(1/g)) < time_sim,
                                           ((1/(b*time_sim*(In + Ic  + qtest*Ic_test + qtest*In_test + qtrace*Ic_trace)))+(1/g))/time_sim,
                                           1)
        
        prop_time_not_recovered_n = ifelse(((1/(b*time_sim*(In + Ic  + qtest*Ic_test + qtest*In_test + qtrace*Ic_trace)))+(1/g)) < time_sim,
                                           ((1/(b*time_sim*(In + Ic  + qtest*Ic_test + qtest*In_test + qtrace*Ic_trace)))+(1/g))/time_sim,
                                           1)
        
        
        # "EFFECTIVE" COST: I use a parameter ip to "activate" the immunitypassport
        # If ip = 0, no immunity passport, cost is C or E
        # if ip = 1, activate immunity passport, multiply by prop_time_not_recovered
        
        C_eff<- ifelse(ip == 0, C, C * prop_time_not_recovered_c)
        E_eff<- ifelse(ip == 0, E, E * prop_time_not_recovered_n)
        
        
        # PAYOFF = benefit - cost
        pi_n = p_not_infected_c - E_eff
        pi_c = p_not_infected_n - C_eff
      }else{
        pi_n = 0
        pi_c = 0
      }
      
      # Sc/(Sc+Sn)=  Current fraction of compliers among the susceptible
      
      # ODE, with replicator dynamic for the S "activated" if x > 1              
      dSn = -b*(In + Ic  + qtest*Ic_test + qtest*In_test + qtrace*Ic_trace)*Sn + x*(Sn/(Sn+Sc)) * (pi_n - ((pi_n *(Sn/(Sn+Sc))) + (pi_c *(Sc/(Sn+Sc)))))
      dSc = -b*(In + Ic  + qtest*Ic_test + qtest*In_test + qtrace*Ic_trace)*Sc + x*(Sc/(Sn+Sc)) * (pi_c - ((pi_n *(Sn/(Sn+Sc))) + (pi_c *(Sc/(Sn+Sc)))))
      
      dIn = (1-p_test)*(b*(In + Ic  + qtest*Ic_test + qtest*In_test + qtrace*Ic_trace)*Sn)- g * In
      dIn_test = (p_test)*(b*(In + Ic  + qtest*Ic_test + qtest*In_test + qtrace*Ic_trace)*Sn)- g * In_test
      
      dIc = (1-p_test)*(b*(In + Ic  + qtest*In_test)*Sc)- g * Ic
      dIc_test = (p_test)*(b*(In + Ic  + qtest*In_test)*Sc)- g * Ic_test
      dIc_trace = b*(qtest*Ic_test + qtrace*Ic_trace)*Sc - g*Ic_trace
      
      dRn  = g * (In + In_test)
      dRc  = g * (Ic + Ic_test + Ic_trace)
      
    }else{
      
      print('npi argument must be either sd or tti')
    }
    
    #return(list(c(dSn, dSc, dEn, dEc, dEc_trace, dIn, dIn_test, dIc, dIc_test, dIc_trace, dRn, dRc)))})
    return(list(c(dSn, dSc, dIn, dIn_test, dIc, dIc_test, dIc_trace, dRn, dRc)))})
  
  
}

# code wrapper to plot model output.
# type = 'sir'  plots the sums of different classes (compliers/non-compliers, trace,test when releavant) to plot S, I, R
# type = 'sircn' plot each specific category
# set real_pop = TRUE plots actual number, FALSE plot proportion of pop
plot_epi_sir<- function(mod.out, type = 'sir', real_pop = FALSE, pop = 67886004){
  if(real_pop == FALSE){
    if(type == 'sir'){
      mod.out %>%
        as.data.frame() %>%
        mutate(S = Sn+Sc,
               #E = En+Ec+Ec_trace,
               I = In+Ic+In_test+Ic_test+Ic_trace,
               R = Rn + Rc) %>%
        gather(class, amount, 11:13) %>%
        ggplot(., aes(x = time, y = amount, col = class))+
        geom_line()
      
    }else if(type == 'sircn'){
      mod.out %>%
        as.data.frame() %>%
        gather(class, amount, 2:10) %>%
        ggplot(., aes(x = time, y = amount, col = class))+
        geom_line()
    }else{
      print('invalid type')
    }
  }else if(real_pop == TRUE){
    if(type == 'sir'){
      mod.out %>%
        as.data.frame() %>%
        mutate(Sn = Sn*pop,
               Sc = Sc*pop,
               #En = En*pop,
               #Ec = Ec*pop,
               #Ec_trace = Ec_trace*pop,
               In = In*pop,
               In_test = In_test*pop,
               Ic = Ic*pop,
               Ic_test = Ic_test*pop,
               Ic_trace = Ic_trace*pop,
               Rn = Rn*pop,
               Rc = Rc*pop) %>%
        mutate(S = Sn+Sc,
               #E = En+Ec+Ec_trace,
               I = In+Ic+In_test+Ic_test+Ic_trace,
               R = Rn + Rc) %>%
        gather(class, amount, 11:13) %>%
        ggplot(., aes(x = time, y = amount, col = class))+
        geom_line()
      
    }else if(type == 'sircn'){
      mod.out %>%
        as.data.frame() %>%
        mutate(Sn = Sn*pop,
               Sc = Sc*pop,
               #En = En*pop,
               #Ec = Ec*pop,
               #Ec_trace = Ec_trace*pop,
               In = In*pop,
               In_test = In_test*pop,
               Ic = Ic*pop,
               Ic_test = Ic_test*pop,
               Ic_trace = Ic_trace*pop,
               Rn = Rn*pop,
               Rc = Rc*pop) %>%
        gather(class, amount, 2:10) %>%
        ggplot(., aes(x = time, y = amount, col = class))+
        geom_line()
    }else{
      print('invalid type')
    }
    
  }else{
    print('real_pop must be TRUE/FALSE')
  }
}
integral<- function(x, strains){ 
  out<- as.data.frame(x)
  out_sum<- data.frame(time = out$time, sum = rowSums(as.data.frame(out[,strains])))
  sum(diff(out_sum$time) * (head(out_sum$sum,-1)+tail(out_sum$sum,-1)))/2
}


# ... (b) Parameters  ----

# Epidemic parameters
r0 = 3   # (reproduction number)
g = 1/14 # (recovery rate)
b = r0*g # (transmission rate estimated from R0)
N = 67886004 # used to define initial condition
nb_cases_ini = 60000#1000 # number of cases on March 16th when PM first announce to limit contacts = 58286


time_sim = 800
timestep = 0.1
times.run <- seq(0, time_sim, by = timestep)


# NO INTERVENTION

# b, g  are common to all simulations
# a only plays a role in SD model
# tau and p_test only play a role in TTI models
# x > 0 'activates' the replicator dynamic
# C, E, time_sim and ip only play a role when model is run with replicator dynamic, where ip = 1 means that immunity passport is activated (for static model, effect +/- immunity passport can be evaluate afterward --> see in figures script)

# ~~~~~~~~~~ Define intervention scenarios parameter sets ~~~~~~~~~~

# parameters for static models
# x, C, E, time_sim and ip do not play a role in this case, but defining them to dummy values otherwise function output NA
pars_null<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0, x = 0)
pars_sd<-   c(b = b, g = g, a = 4, tau = NA,   p_test = NA, x = 0)
pars_tti<-  c(b = b, g = g, a = NA,  tau = 5.1,  p_test = 0.7, x = 0)


# Quick test of effect of each intervention
prop_n = 0.9
ini<- c(Sn = prop_n - (0.5/N),
        Sc = (1 - prop_n) - (0.5/N),
        In = 0.5/N,
        In_test = 0,
        Ic = 0.5/N,
        Ic_test = 0,
        Ic_trace = 0,
        Rn = 0,
        Rc = 0)


# run static model
out_null<- ode(y = ini,  times = times.run, func = sir_npi, parms = pars_null, npi = 'sd') # does not matter which npi when running null model
out_sd<-   ode(y = ini,  times = times.run, func = sir_npi, parms = pars_sd, npi = 'sd')
out_tti<-  ode(y = ini,  times = times.run, func = sir_npi, parms = pars_tti, npi = 'tti')


grid.arrange(ncol = 3,
             plot_epi_sir(out_null, type = 'sir'),
             plot_epi_sir(out_sd, type = 'sir'),
             plot_epi_sir(out_tti, type = 'sir'))



# [2] SD STATIC FULL SIMULATION ----

# Run on a dataframe of various initial proportion of compliance
# record each time the bits needed to compute the payoff from each simulation afterward

min = 0.001 # was 0.001, test with 0.1
df.sd<- data.frame(Sn = seq(0+min, 1-min, min),
                   p_not_infected_c = NA,
                   p_not_infected_n = NA,
                   prop_time_not_recovered_c = NA,
                   prop_time_not_recovered_n = NA)



#pdf(paste0('./scratch/epidemic_output_static_SEIR_SD_', it, '.pdf'), width = 12, height = 4)
for(i in 1:nrow(df.sd)){
  
  print(i)
  
  # redefine ini at each loop iteration,
  # consider number of cases reported in early march in the UK (321)
  # split the initial number of cases among En and Ec according to the Sn and Sc proportions 
  
  ini<- c(Sn = df.sd$Sn[i] - ((df.sd$Sn[i]*nb_cases_ini)/N),
          Sc = (1 - df.sd$Sn[i]) - (((1 - df.sd$Sn[i])*nb_cases_ini)/N),
          In = (df.sd$Sn[i]*nb_cases_ini)/N,
          In_test = 0,
          Ic = ((1 - df.sd$Sn[i])*nb_cases_ini)/N,
          Ic_test = 0,
          Ic_trace = 0,
          Rn = 0,
          Rc = 0)
  
  
  out<- ode(y = ini, times = times.run, func = sir_npi, parms = pars_sd, npi = 'sd')
  
  #grid.arrange(ncol = 2,
  #             plot_epi_sir(out, type = 'seircn')+ylim(0, 1),
  #             plot_epi_sir(out, type = 'seir')+ylim(0,1))
  
  
  # EVALUATE BENEFIT AND COST AT EQUILIBRIUM (tmax), so compute PAYOFFS AFTERWARD
  
  # In SD model
  # proba of individual is still susceptible by the end of the epidemic
  df.sd$p_not_infected_c[i]<- tail(out, 1)[,c('Sc')]/sum(tail(out, 1)[,c('Sc', 'Ic', 'Rc')]) 
  df.sd$p_not_infected_n[i]<- tail(out, 1)[,c('Sn')]/sum(tail(out, 1)[,c('Sn', 'In', 'Rn')])
  
  # Using integral, average proportion of time individual spent not recovered
  df.sd$prop_time_not_recovered_c[i]<- integral(out, c('Sc', 'Ic'))/integral(out, c('Sc', 'Ic', 'Rc')) 
  df.sd$prop_time_not_recovered_n[i]<- integral(out, c('Sn', 'In'))/integral(out, c('Sn', 'In', 'Rn'))
  
}
#dev.off()


save.image('./output/model_runs/SIMULATION_SIR_SD_June01_nbIni_60000_alpha_4.RData')



# [3] TTI STATIC FULL SIMULATION ----

# Run on a dataframe of various initial proportion of compliance
# record each time the bits needed to compute the payoff from each simulation afterward

min = 0.001 # was 0.001, test with 0.1
df.tti<- data.frame(Sn = seq(0+min, 1-min, min),
                    p_not_infected_c = NA,
                    p_not_infected_n = NA,
                    prop_time_not_recovered_c = NA,
                    prop_time_not_recovered_n = NA)


#pdf(paste0('./scratch/epidemic_output_static_SEIR_TTI_', it, '.pdf'), width = 12, height = 4)
for(i in 1:nrow(df.tti)){
  print(i)
  
  # redefine ini at each loop iteration,
  # consider number of cases reported in early march in the UK (321)
  # split the initial number of cases among En and Ec according to the Sn and Sc proportions 
  
  ini<- c(Sn = df.tti$Sn[i] - ((df.tti$Sn[i]*nb_cases_ini)/N),
          Sc = (1 - df.tti$Sn[i]) - (((1 - df.tti$Sn[i])*nb_cases_ini)/N),
          In = (df.tti$Sn[i]*nb_cases_ini)/N,
          In_test = 0,
          Ic = ((1 - df.tti$Sn[i])*nb_cases_ini)/N,
          Ic_test = 0,
          Ic_trace = 0,
          Rn = 0,
          Rc = 0)
  
  
  out<- ode(y = ini, times = times.run, func = sir_npi, parms = pars_tti, npi = 'tti')
  
  #grid.arrange(ncol = 2,
  #             plot_epi_seir(out, type = 'seircn')+ylim(0, 1),
  #             plot_epi_seir(out, type = 'seir')+ylim(0,1))
  
  
  # EVALUATE BENEFIT AND COST AT EQUILIBRIUM (tmax), so compute PAYOFFS AFTERWARD
  # In TTI models
  # Evaluation at last time t, proba of individual is still susceptible by the end of the epidemic
  df.tti$p_not_infected_c[i]<- tail(out, 1)[,c('Sc')]/sum(tail(out, 1)[,c('Sc','Ic', 'Ic_test', 'Ic_trace', 'Rc')])
  df.tti$p_not_infected_n[i]<- tail(out, 1)[,c('Sn')]/sum(tail(out, 1)[,c('Sn','In', 'In_test','Rn')])
  
  # Using integral, average proportion of time individual spent not recovered
  df.tti$prop_time_not_recovered_c[i]<- integral(out, c('Sc', 'Ic', 'Ic_test', 'Ic_trace'))/integral(out, c('Sc', 'Ic',  'Ic_test', 'Ic_trace', 'Rc')) 
  df.tti$prop_time_not_recovered_n[i]<- integral(out, c('Sn', 'In', 'In_test'))/integral(out, c('Sn', 'In', 'In_test','Rn'))
}

#dev.off()

save.image('./output/model_runs/SIMULATION_SIR_TTI_June01_nbIni_60000_tau_5.1_ptest_0.7.RData')

