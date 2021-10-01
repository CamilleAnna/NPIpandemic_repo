# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#   Social dilemmas in non-pharmaceutical interventions to tackle epidemics
#                     Simonet C, Sweeny A, McNally L (2021)
#                 Codes for dynamic model simulations and figures
#
#                       last edited: July, 27th, 2021
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Script summary
# 

setwd('path/to/repo') # setwd('~/Documents/PhD/Research/Covid19/NPIsocialDilemma_repo/')
library(readxl); library(ggplot2); library(dplyr); library(tidyr);library(Rmisc); library(lubridate);library(plotly);library(yarrr); library(knitr);library(plyr);library(utils);library(tidyverse);library(gridExtra);library(deSolve); library(stringi); library(stringr);
library(ggthemes); library(deSolve)


# [1] SETUP ----

# ... (a) Model functions

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
      #dSn = -b*(In + (1/a)*Ic)*Sn + x*(Sn+Sc)*(1-fc) * (pi_n - (pi_n * (1 - fc) + pi_c * fc))
      #dSc = -b*(In + (1/a)*Ic)*(1/a)*Sc + x*(Sn+Sc)*fc * (pi_c - (pi_n * (1 - fc) + pi_c * fc))
      
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
      
      #fc = Sc/(Sc+Sn)  # Current fraction of compliers among the susceptible
      
      # ODE, with replicator dynamic for the S "activated" if x > 1              
      dSn = -b*(In + Ic  + qtest*Ic_test + qtest*In_test + qtrace*Ic_trace)*Sn + x*(Sn/(Sn+Sc)) * (pi_n - ((pi_n *(Sn/(Sn+Sc))) + (pi_c *(Sc/(Sn+Sc)))))
      dSc = -b*(In + Ic  + qtest*Ic_test + qtest*In_test + qtrace*Ic_trace)*Sc + x*(Sc/(Sn+Sc)) * (pi_c - ((pi_n *(Sn/(Sn+Sc))) + (pi_c *(Sc/(Sn+Sc)))))
      
      #dEn =  b*(In + Ic  + qtest*Ic_test + qtest*In_test + qtrace*Ic_trace)*Sn - s * En
      #dEc =  b*(In + Ic  + qtest*In_test)*Sc - s*Ec
      #dEc_trace = b*(qtest*Ic_test + qtrace*Ic_trace)*Sc - s*Ec_trace
      
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


# ... (b) Parameters
# Epidemic parameters
r0 = 3
g = 1/14
b = r0*g # (transmission rate estimated from R0)
N = 67886004 


# [2.1] SD DYNAMICS, main ----
#source('./scripts/2021_MODELSetup_alter.R')


# Initial conditions
nb_cases_ini = 2000
prop_n = 0.7
ini<- c(Sn = prop_n - ((prop_n*nb_cases_ini)/N),
        Sc = (1 - prop_n) - (((1 - prop_n)*nb_cases_ini)/N),
        In = (prop_n*nb_cases_ini)/N,
        In_test = 0,
        Ic = ((1 - prop_n)*nb_cases_ini)/N,
        Ic_test = 0,
        Ic_trace = 0,
        Rn = 0,
        Rc = 0)


# Time parameters

time_burnin = 30
time_sd = 400
startip = 60 # 
time_sd.ip = time_sd - startip
time_sim = 60 # t_eval in MS main text
timestep = 0.1
times.burnin <- seq(0, time_burnin, by = timestep)
times.run.sd <- seq(0, time_sd, by = timestep)
times.run.sd.ip <- seq(0, time_sd.ip, by = timestep)
times.noIntervention <- seq(0, time_burnin+time_sd, by = timestep)


# Plot dynamics
m<- rbind(
  c(0, 0.65, 0.555, 1),
  c(0, 0.65, 0, 0.444),
  c(0.7, 1, 0.555, 1)
)



pdf('./output/figures/Figure5_part1_SD_dynamic.pdf', width = (10+1+1)/2.54, height = (9+1+1)/2.54, pointsize = 7)


close.screen(all.screens = TRUE)
split.screen(m)
par(omi = rep(1/2.54,4))
for(i in 1:nrow(m)){
  screen(i)
  par(mar = c(0, 0, 0, 0))
  plot(1, axes = FALSE, type = "n")
  #box()
}


#cols = c('gold', 'red', 'dodgerblue', 'seagreen', 'darkorchid')
#darkcols = c('gold', 'red', 'dodgerblue', 'seagreen', 'darkorchid')

cols = rev(c('gold', 'red', 'dodgerblue', 'seagreen'))
darkcols = rev(c('gold', 'red', 'dodgerblue', 'seagreen'))

# Get those save colors version without having to used transparency settings after
colorRampPalette(c("white", "seagreen"))(10)[4]
cols = rev(c('#FFE871', '#FF7171', '#82C1FF', '#8ABEA1'))
darkcols = rev(c('#FFF1AA', '#FFAAAA', '#B4DAFF', '#B9D8C7'))



focus.x = 0.3
focus.C = 1.2
focus.a = 4

Es = c(1, 1.2, 1.3, 1.4)



# SCREEN 1: compliance dynamic ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
screen(1)

plot('', ylim = c(0, 1), xlim = c(0, time_burnin+time_sd), ylab = '', xlab = '', cex.axis = 0.9, xaxt = 'n')
#plot('', ylim = c(0, 1), xlim = c(0, 90), ylab = '', xlab = '', cex.axis = 0.9, xaxt = 'n')

for(i in 1:length(Es)){
  
  # parameters for dynamic model, no immunity passport
  pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
  #pars_null.dyn<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0,   x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_sd.dyn<-   c(b = b, g = g, a = focus.a, tau = NA,   p_test = NA,  x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_sd.dyn.ip<-c(b = b, g = g, a = focus.a, tau = NA,   p_test = NA,  x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 1)
  
  #out.burnin.background <- ode(y = ini,  times = times.burnin.background, func = sir_npi, parms = pars_null.dyn, npi = 'sd')  
  out.burnin <- ode(y = ini,  times = times.burnin, func = sir_npi, parms = pars_null, npi = 'sd')  
  out.run.sd<- ode(y = tail(out.burnin,1)[,-1], times = times.run.sd, func = sir_npi, parms = pars_sd.dyn, npi = 'sd')
  out.run.sd.ip<- ode(y = out.run.sd[which(out.run.sd[,'time'] == startip),-1], times = times.run.sd.ip, func = sir_npi, parms = pars_sd.dyn.ip, npi = 'sd')
  
  # Assemble output
  out.1<- rbind(out.burnin, out.run.sd[-1,])
  out.2<- rbind(out.burnin, out.run.sd[2:which(out.run.sd[,'time'] == startip),], out.run.sd.ip[-1,])
  out.1[,'time']<- seq(0, tail(time_burnin,1)+tail(time_sd, 1), timestep)
  out.2[,'time']<- seq(0, tail(time_burnin,1)+tail(time_sd, 1), timestep)
  
  # compliance level among non-recovered (in SD model, _test and _trace = 0)
  lines((Sc + Ic + Ic_test + Ic_trace)/(Sn + Sc + In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.1, lwd = 1.5, col = cols[i], lty = 1)
  out.test<- out.2[which(out.run.sd[,'time'] == I(startip+time_burnin)):nrow(out.2),]
  lines((Sc + Ic + Ic_test + Ic_trace)/(Sn + Sc + In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.test, lwd = 1.5, col = darkcols[i], lty = 2)
  
  
}
abline(v = time_burnin, lty = 2)
abline(v = time_burnin+startip, lty = 2)
abline(h = 1 - prop_n, col = 'grey')
mtext(side = 2, 'Compliance frequency', font = 2, cex = 1, line =  2.2)
mtext(side = 3, '+SD', font = 2, cex = 0.7, line =  0, at = time_burnin)#, adj = 0)
mtext(side = 3, '+IP', font = 2, cex = 0.7, line =  0, at = time_burnin+startip)#, adj = 0)
axis(1, at = seq(0, time_burnin+time_sd, 30))


# SCREEN 2: epidemic dynamic ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

screen(2) # alpha(cols[es], .3), lwd = 3,  xaxt = 'n'
plot('', ylim = c(0, 0.15), xlim = c(0, time_burnin+time_sd), ylab = '', xlab = '', cex.axis = 0.9, xaxt = 'n')
#plot('', ylim = c(0, 0.17), xlim = c(0, 90), ylab = '', xlab = '', cex.axis = 0.9, xaxt = 'n')

# plot epidemic in the absence of intervention
pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
out.Nointervention <- ode(y = ini,  times = times.noIntervention, func = sir_npi, parms = pars_null, npi = 'sd')  
#lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.Nointervention, lwd = 3, col = 'lightgrey', lty = 1)
x =  out.Nointervention[,'time']
y = rowSums(out.Nointervention[,c('In', 'In_test', 'Ic', 'Ic_test', 'Ic_trace')])
#polygon(c(x[x>=0], max(x), 0), c(y[x>=0], 0, 0), col = '#F0F0F0', border = NA)

polygon(c(x[x>=0], max(x), 0), c(ifelse(y[x>=0] > 0.155, 0.155, y[x>=0]), 0, 0), col = '#F0F0F0', border = NA)



for(i in 1:length(Es)){
  
  # parameters for dynamic model, no immunity passport
  pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
  #pars_null.dyn<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0,   x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_sd.dyn<-   c(b = b, g = g, a = focus.a, tau = NA,   p_test = NA,  x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_sd.dyn.ip<-c(b = b, g = g, a = focus.a, tau = NA,   p_test = NA,  x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 1)
  
  out.burnin <- ode(y = ini,  times = times.burnin, func = sir_npi, parms = pars_null, npi = 'sd')  
  out.run.sd<- ode(y = tail(out.burnin,1)[,-1], times = times.run.sd, func = sir_npi, parms = pars_sd.dyn, npi = 'sd')
  out.run.sd.ip<- ode(y = out.run.sd[which(out.run.sd[,'time'] == startip),-1], times = times.run.sd.ip, func = sir_npi, parms = pars_sd.dyn.ip, npi = 'sd')
  
  # Assemble output
  out.1<- rbind(out.burnin, out.run.sd[-1,])
  out.2<- rbind(out.burnin, out.run.sd[2:which(out.run.sd[,'time'] == startip),], out.run.sd.ip[-1,])
  out.1[,'time']<- seq(0, tail(time_burnin,1)+tail(time_sd, 1), timestep)
  out.2[,'time']<- seq(0, tail(time_burnin,1)+tail(time_sd, 1), timestep)
  
  # amount of infected (in SD model, _test and _trace = 0)
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.1, lwd = 1.5, col = cols[i], lty = 1)
  out.test<- out.2[which(out.run.sd[,'time'] == I(startip+time_burnin)):nrow(out.2),]
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.test, lwd = 1.5, col = darkcols[i], lty = 2)
  
}
abline(v = time_burnin, lty = 2)
abline(v = time_burnin+startip, lty = 2)
mtext(side = 2, 'Prevalence', line = 2.2, font = 2, cex = 1)
mtext(side = 1, 'Days', line = 2.2, font = 2, cex = 1)
axis(1, at = seq(0, time_burnin+time_sd, 30))
#axis(2, at = c(0, 0.005, 0.01),
#     labels = c(0, 0.005, 0.01))


# SCREEN 3: legend  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

screen(3)
plot(1, axes = FALSE, type = "n", ylim = c(0, 1), xlim = c(0, 1))

text(x = 0, y = 0.95, labels = 'Enforcement level:', adj = 0, cex = 0.8, font = 2)

lines(x = c(0, 0.1), y = rep(0.85, 2), lwd = 2, col = cols[1])
lines(x = c(0, 0.1), y = rep(0.8, 2),  lwd = 2, col = cols[2])
lines(x = c(0, 0.1), y = rep(0.75, 2), lwd = 2, col = cols[3])
lines(x = c(0, 0.1), y = rep(0.7, 2),  lwd = 2, col = cols[4])
#lines(x = c(0, 0.1), y = rep(0.65, 2),  lwd = 2, col = alpha(cols[5], .6))

text(x = 0.15, y = 0.85, labels = ': E < C', adj = 0, cex = 0.75)
text(x = 0.15, y = 0.8, labels = ': E = C', adj = 0, cex = 0.75)
text(x = 0.15, y = 0.75, labels = ': E > C', adj = 0, cex = 0.75)
text(x = 0.15, y = 0.7, labels = ': E >> C', adj = 0, cex = 0.75)
#text(x = 0.15, y = 0.65, labels = ': E = 5', adj = 0, cex = 0.75)

text(x = 0, y = 0.55, labels = 'From t = 90:', adj = 0, cex = 0.8, font = 2)

lines(x = c(0, 0.1), y = rep(0.45, 2),  lwd = 2, col = 'black')
lines(x = c(0, 0.1), y = rep(0.40, 2),  lwd = 2, col = 'black', lty = 3)

text(x = 0.15, y = 0.45, labels = ': without immunity passport', adj = 0, cex = 0.75)
text(x = 0.15, y = 0.40, labels = ': with immunity passport', adj = 0, cex = 0.75)


dev.off()



# [2.2] TTI DYNAMICS, main ----

# Initial conditions
nb_cases_ini = 2000
prop_n = 0.7
ini<- c(Sn = prop_n - ((prop_n*nb_cases_ini)/N),
        Sc = (1 - prop_n) - (((1 - prop_n)*nb_cases_ini)/N),
        In = (prop_n*nb_cases_ini)/N,
        In_test = 0,
        Ic = ((1 - prop_n)*nb_cases_ini)/N,
        Ic_test = 0,
        Ic_trace = 0,
        Rn = 0,
        Rc = 0)


# Time parameters

time_burnin = 30
time_NPI = 400
startip = 60 # IP introduced seven days after implementing Social distancing
time_NPI.ip = time_NPI - startip
time_sim = 60
timestep = 0.1
times.burnin <- seq(0, time_burnin, by = timestep)
times.run.NPI <- seq(0, time_NPI, by = timestep)
times.run.NPI.ip <- seq(0, time_NPI.ip, by = timestep)
times.noIntervention <- seq(0, time_burnin+time_NPI, by = timestep)


# Plot dynamics
m<- rbind(
  c(0, 0.65, 0.555, 1),
  c(0, 0.65, 0, 0.444),
  c(0.7, 1, 0.555, 1)
)


pdf('./output/figures/Figure5_part2_TTI_dynamic.pdf', width = (10+1+1)/2.54, height = (9+1+1)/2.54, pointsize = 7)

close.screen(all.screens = TRUE)
split.screen(m)
par(omi = rep(1/2.54,4))
for(i in 1:nrow(m)){
  screen(i)
  par(mar = c(0, 0, 0, 0))
  plot(1, axes = FALSE, type = "n")
  #box()
}


#cols = c('gold', 'red', 'dodgerblue', 'seagreen', 'darkorchid')
#darkcols = c('gold', 'red', 'dodgerblue', 'seagreen', 'darkorchid')

cols = rev(c('#FFE871', '#FF7171', '#82C1FF', '#8ABEA1'))
darkcols = rev(c('#FFF1AA', '#FFAAAA', '#B4DAFF', '#B9D8C7'))


focus.x = 0.3
focus.C = 1.2
#focus.a = 4

Es = c(1.1, 1.2, 1.25, 1.3)
#Es = c(1, 1.2, 1.3, 1.4)


focus.ptest = 0.7
focus.tau = 5.1

# SCREEN 1: compliance dynamic ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
screen(1)

plot('', ylim = c(0, 1), xlim = c(0, time_NPI+time_burnin), ylab = '', xlab = '', cex.axis = 0.9, xaxt = 'n')
for(i in 1:length(Es)){
  
  # parameters for dynamic model, no immunity passport
  pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
  #pars_null.dyn<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0,   x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_tti.dyn<-  c(b = b, g = g, a = NA,  tau = focus.tau,  p_test = focus.ptest, x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_tti.dyn.ip<-  c(b = b, g = g, a = NA,  tau = focus.tau,  p_test = focus.ptest, x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 1)
  
  out.burnin <- ode(y = ini,  times = times.burnin, func = sir_npi, parms = pars_null, npi = 'tti')  
  out.run.NPI<- ode(y = tail(out.burnin,1)[,-1], times = times.run.NPI, func = sir_npi, parms = pars_tti.dyn, npi = 'tti')
  out.run.NPI.ip<- ode(y = out.run.NPI[which(out.run.NPI[,'time'] == startip),-1], times = times.run.NPI.ip, func = sir_npi, parms = pars_tti.dyn.ip, npi = 'tti')
  
  # Assemble output
  out.1<- rbind(out.burnin, out.run.NPI[-1,])
  out.2<- rbind(out.burnin, out.run.NPI[2:which(out.run.NPI[,'time'] == startip),], out.run.NPI.ip[-1,])
  out.1[,'time']<- seq(0, tail(time_burnin,1)+tail(time_NPI, 1), timestep)
  out.2[,'time']<- seq(0, tail(time_burnin,1)+tail(time_NPI, 1), timestep)
  
  # compliance level among non-recovered (in SD model, _test and _trace = 0)
  lines((Sc + Ic + Ic_test + Ic_trace)/(Sn + Sc + In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.1, lwd = 1.5, col = cols[i], lty = 1)
  out.test<- out.2[which(out.run.NPI[,'time'] == I(startip+time_burnin)):nrow(out.2),]
  lines((Sc + Ic + Ic_test + Ic_trace)/(Sn + Sc + In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.test, lwd = 1.5, col = darkcols[i], lty = 2)
  
  
}
abline(v = time_burnin, lty = 2)
abline(v = time_burnin+startip, lty = 2)
abline(h = 1 - prop_n, col = 'grey')
mtext(side = 2, 'Compliance frequency', font = 2, cex = 1, line =  2.2)
mtext(side = 3, '+TTI', font = 2, cex = 0.7, line =  0, at = time_burnin)#, adj = 0)
mtext(side = 3, '+IP', font = 2, cex = 0.7, line =  0, at = time_burnin+startip)#, adj = 0)
axis(1, at = seq(0, time_NPI+time_burnin, 30))


# SCREEN 2: epidemic dynamic ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

screen(2) # alpha(cols[es], .3), lwd = 3,  xaxt = 'n'
plot('', ylim = c(0, 0.15), xlim = c(0, time_NPI+time_burnin), ylab = '', xlab = '', cex.axis = 0.9, xaxt = 'n')
# plot epidemic in the absence of intervention
pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
out.Nointervention <- ode(y = ini,  times = times.noIntervention, func = sir_npi, parms = pars_null, npi = 'tti')  
#lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.Nointervention, lwd = 3, col = 'lightgrey', lty = 1)
x =  out.Nointervention[,'time']
y = rowSums(out.Nointervention[,c('In', 'In_test', 'Ic', 'Ic_test', 'Ic_trace')])
#polygon(c(x[x>=0], max(x), 0), c(y[x>=0], 0, 0), col = '#F0F0F0', border = NA)

polygon(c(x[x>=0], max(x), 0), c(ifelse(y[x>=0] > 0.155, 0.155, y[x>=0]), 0, 0), col = '#F0F0F0', border = NA)


for(i in 1:length(Es)){
  
  # parameters for dynamic model, no immunity passport
  pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
  #pars_null.dyn<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0,   x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_tti.dyn<-  c(b = b, g = g, a = NA,  tau = focus.tau,  p_test = focus.ptest, x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_tti.dyn.ip<-  c(b = b, g = g, a = NA,  tau = focus.tau,  p_test = focus.ptest, x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 1)
  
  out.burnin <- ode(y = ini,  times = times.burnin, func = sir_npi, parms = pars_null, npi = 'tti')  
  out.run.NPI<- ode(y = tail(out.burnin,1)[,-1], times = times.run.NPI, func = sir_npi, parms = pars_tti.dyn, npi = 'tti')
  out.run.NPI.ip<- ode(y = out.run.NPI[which(out.run.NPI[,'time'] == startip),-1], times = times.run.NPI.ip, func = sir_npi, parms = pars_tti.dyn.ip, npi = 'tti')
  
  # Assemble output
  out.1<- rbind(out.burnin, out.run.NPI[-1,])
  out.2<- rbind(out.burnin, out.run.NPI[2:which(out.run.NPI[,'time'] == startip),], out.run.NPI.ip[-1,])
  out.1[,'time']<- seq(0, tail(time_burnin,1)+tail(time_NPI, 1), timestep)
  out.2[,'time']<- seq(0, tail(time_burnin,1)+tail(time_NPI, 1), timestep)
  
  # amount of infected (in SD model, _test and _trace = 0)
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.1, lwd = 1.5, col = cols[i], lty = 1)
  out.test<- out.2[which(out.run.NPI[,'time'] == I(startip+time_burnin)):nrow(out.2),]
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.test, lwd = 1.5, col = darkcols[i], lty = 2)
  
}

abline(v = time_burnin, lty = 2)
abline(v = time_burnin+startip, lty = 2)
mtext(side = 2, 'Prevalence', line = 2.2, font = 2, cex = 1)
mtext(side = 1, 'Days', line = 2.2, font = 2, cex = 1)
axis(1, at = seq(0, time_NPI+time_burnin, 30))


# SCREEN 3: legend  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

screen(3)
plot(1, axes = FALSE, type = "n", ylim = c(0, 1), xlim = c(0, 1))

text(x = 0, y = 0.95, labels = 'Enforcement level:', adj = 0, cex = 0.8, font = 2)

lines(x = c(0, 0.1), y = rep(0.85, 2), lwd = 2, col = cols[1])
lines(x = c(0, 0.1), y = rep(0.8, 2),  lwd = 2, col = cols[2])
lines(x = c(0, 0.1), y = rep(0.75, 2), lwd = 2, col = cols[3])
lines(x = c(0, 0.1), y = rep(0.7, 2),  lwd = 2, col = cols[4])
#lines(x = c(0, 0.1), y = rep(0.65, 2),  lwd = 2, col = alpha(cols[5], .6))

text(x = 0.15, y = 0.85, labels = ': E < C', adj = 0, cex = 0.75)
text(x = 0.15, y = 0.8, labels = ': E = C', adj = 0, cex = 0.75)
text(x = 0.15, y = 0.75, labels = ': E > C', adj = 0, cex = 0.75)
text(x = 0.15, y = 0.7, labels = ': E >> C', adj = 0, cex = 0.75)
#text(x = 0.15, y = 0.65, labels = ': E = 5', adj = 0, cex = 0.75)

text(x = 0, y = 0.55, labels = 'From t = 90:', adj = 0, cex = 0.8, font = 2)

lines(x = c(0, 0.1), y = rep(0.45, 2),  lwd = 1.5, col = 'black')
lines(x = c(0, 0.1), y = rep(0.40, 2),  lwd = 1.5, col = 'black', lty = 3)

text(x = 0.15, y = 0.45, labels = ': without immunity passport', adj = 0, cex = 0.75)
text(x = 0.15, y = 0.40, labels = ': with immunity passport', adj = 0, cex = 0.75)


dev.off()



# [3.1] SD DYNAMICS VARY TIMING ----

# need to empty enviroenment and re-run part [1] otherwise I think some params mess up
#source('./scripts/2021_MODELSetup_alter.R')

#cols = rev(c('gold', 'red', 'dodgerblue', 'seagreen'))
#darkcols = rev(c('gold', 'red', 'dodgerblue', 'seagreen'))

focus.x = 0.3
focus.C = 1.2
focus.a = 4

Es = c(1, 1.2, 1.3, 1.4)


pdf('./output/figures/FigureS5_SD_dynamic_vary.pdf', width = (11+1+1)/2.54, height = (11+1+1)/2.54, pointsize = 7)


par(mfrow = c(3,3),
    mar = rep(2.5,4))

# 1 ----

# Initial conditions
nb_cases_ini = 100
prop_n = 0.7
ini<- c(Sn = prop_n - ((prop_n*nb_cases_ini)/N),
        Sc = (1 - prop_n) - (((1 - prop_n)*nb_cases_ini)/N),
        In = (prop_n*nb_cases_ini)/N,
        In_test = 0,
        Ic = ((1 - prop_n)*nb_cases_ini)/N,
        Ic_test = 0,
        Ic_trace = 0,
        Rn = 0,
        Rc = 0)


# Time parameters

time_burnin = 30
time_sd = 400
startip = 30 # 
time_sd.ip = time_sd - startip
time_sim = 60
timestep = 0.1
times.burnin <- seq(0, time_burnin, by = timestep)
times.run.sd <- seq(0, time_sd, by = timestep)
times.run.sd.ip <- seq(0, time_sd.ip, by = timestep)
times.noIntervention <- seq(0, time_burnin+time_sd, by = timestep)


plot('', ylim = c(0, 0.17), xlim = c(0, time_burnin+time_sd), ylab = '', xlab = '', cex.axis = 1)

# plot epidemic in the absence of intervention
pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
out.Nointervention <- ode(y = ini,  times = times.noIntervention, func = sir_npi, parms = pars_null, npi = 'sd')  
#lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.Nointervention, lwd = 3, col = 'lightgrey', lty = 1)
x =  out.Nointervention[,'time']
y = rowSums(out.Nointervention[,c('In', 'In_test', 'Ic', 'Ic_test', 'Ic_trace')])
polygon(c(x[x>=0], max(x), 0), c(ifelse(y[x>=0] > 0.155, 0.155, y[x>=0]), 0, 0), col = '#F0F0F0', border = NA)



for(i in 1:length(Es)){
  
  # parameters for dynamic model, no immunity passport
  pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
  #pars_null.dyn<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0,   x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_sd.dyn<-   c(b = b, g = g, a = focus.a, tau = NA,   p_test = NA,  x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_sd.dyn.ip<-c(b = b, g = g, a = focus.a, tau = NA,   p_test = NA,  x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 1)
  
  out.burnin <- ode(y = ini,  times = times.burnin, func = sir_npi, parms = pars_null, npi = 'sd')  
  out.run.sd<- ode(y = tail(out.burnin,1)[,-1], times = times.run.sd, func = sir_npi, parms = pars_sd.dyn, npi = 'sd')
  out.run.sd.ip<- ode(y = out.run.sd[which(out.run.sd[,'time'] == startip),-1], times = times.run.sd.ip, func = sir_npi, parms = pars_sd.dyn.ip, npi = 'sd')
  
  # Assemble output
  out.1<- rbind(out.burnin, out.run.sd[-1,])
  out.2<- rbind(out.burnin, out.run.sd[2:which(out.run.sd[,'time'] == startip),], out.run.sd.ip[-1,])
  out.1[,'time']<- seq(0, tail(time_burnin,1)+tail(time_sd, 1), timestep)
  out.2[,'time']<- seq(0, tail(time_burnin,1)+tail(time_sd, 1), timestep)
  
  # amount of infected (in SD model, _test and _trace = 0)
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.1, lwd = 1.5, col = cols[i], lty = 1)
  out.test<- out.2[which(out.run.sd[,'time'] == I(startip+time_burnin)):nrow(out.2),]
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.test, lwd = 1.5, col = darkcols[i], lty = 2)
  
}
abline(v = time_burnin, lty = 2)
abline(v = time_burnin+startip, lty = 2)
#mtext(side = 2, 'Prevalence', line = 2.2, font = 2, cex = 1)
#mtext(side = 1, 'Days', line = 2.2, font = 2, cex = 1)
#axis(1, at = seq(0, time_burnin+time_sd, 30), labels = FALSE)
#axis(2, at = c(0, 0.005, 0.01),
#     labels = c(0, 0.005, 0.01))

mtext(side = 3, '+SD', font = 2, cex = 0.9, line =  0.3, at = time_burnin)#, adj = 0)
mtext(side = 3, '+IP 30 days later', font = 2, cex = 0.9, line =  0.3, at = time_burnin+startip, adj = 0)


# 2 ----

# Initial conditions
nb_cases_ini = 100
prop_n = 0.7
ini<- c(Sn = prop_n - ((prop_n*nb_cases_ini)/N),
        Sc = (1 - prop_n) - (((1 - prop_n)*nb_cases_ini)/N),
        In = (prop_n*nb_cases_ini)/N,
        In_test = 0,
        Ic = ((1 - prop_n)*nb_cases_ini)/N,
        Ic_test = 0,
        Ic_trace = 0,
        Rn = 0,
        Rc = 0)


# Time parameters

time_burnin = 30
time_sd = 400
startip = 60 # 
time_sd.ip = time_sd - startip
time_sim = 60
timestep = 0.1
times.burnin <- seq(0, time_burnin, by = timestep)
times.run.sd <- seq(0, time_sd, by = timestep)
times.run.sd.ip <- seq(0, time_sd.ip, by = timestep)
times.noIntervention <- seq(0, time_burnin+time_sd, by = timestep)


plot('', ylim = c(0, 0.17), xlim = c(0, time_burnin+time_sd), ylab = '', xlab = '', cex.axis = 1)

# plot epidemic in the absence of intervention
pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
out.Nointervention <- ode(y = ini,  times = times.noIntervention, func = sir_npi, parms = pars_null, npi = 'sd')  
#lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.Nointervention, lwd = 3, col = 'lightgrey', lty = 1)
x =  out.Nointervention[,'time']
y = rowSums(out.Nointervention[,c('In', 'In_test', 'Ic', 'Ic_test', 'Ic_trace')])
polygon(c(x[x>=0], max(x), 0), c(ifelse(y[x>=0] > 0.155, 0.155, y[x>=0]), 0, 0), col = '#F0F0F0', border = NA)

for(i in 1:length(Es)){
  
  # parameters for dynamic model, no immunity passport
  pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
  #pars_null.dyn<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0,   x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_sd.dyn<-   c(b = b, g = g, a = focus.a, tau = NA,   p_test = NA,  x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_sd.dyn.ip<-c(b = b, g = g, a = focus.a, tau = NA,   p_test = NA,  x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 1)
  
  out.burnin <- ode(y = ini,  times = times.burnin, func = sir_npi, parms = pars_null, npi = 'sd')  
  out.run.sd<- ode(y = tail(out.burnin,1)[,-1], times = times.run.sd, func = sir_npi, parms = pars_sd.dyn, npi = 'sd')
  out.run.sd.ip<- ode(y = out.run.sd[which(out.run.sd[,'time'] == startip),-1], times = times.run.sd.ip, func = sir_npi, parms = pars_sd.dyn.ip, npi = 'sd')
  
  # Assemble output
  out.1<- rbind(out.burnin, out.run.sd[-1,])
  out.2<- rbind(out.burnin, out.run.sd[2:which(out.run.sd[,'time'] == startip),], out.run.sd.ip[-1,])
  out.1[,'time']<- seq(0, tail(time_burnin,1)+tail(time_sd, 1), timestep)
  out.2[,'time']<- seq(0, tail(time_burnin,1)+tail(time_sd, 1), timestep)
  
  # amount of infected (in SD model, _test and _trace = 0)
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.1, lwd = 1.5, col = cols[i], lty = 1)
  out.test<- out.2[which(out.run.sd[,'time'] == I(startip+time_burnin)):nrow(out.2),]
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.test, lwd = 1.5, col = darkcols[i], lty = 2)
  
}
abline(v = time_burnin, lty = 2)
abline(v = time_burnin+startip, lty = 2)
#mtext(side = 2, 'Prevalence', line = 2.2, font = 2, cex = 1)
#mtext(side = 1, 'Days', line = 2.2, font = 2, cex = 1)
#axis(1, at = seq(0, time_burnin+time_sd, 30))
#axis(2, at = c(0, 0.005, 0.01),
#     labels = c(0, 0.005, 0.01))

mtext(side = 3, '+SD', font = 2, cex = 0.9, line =  0.3, at = time_burnin)#, adj = 0)
mtext(side = 3, '+IP 60 days later', font = 2, cex = 0.9, line =  0.3, at = time_burnin+startip, adj = 0)


# 3 ----

# Initial conditions
nb_cases_ini = 100
prop_n = 0.7
ini<- c(Sn = prop_n - ((prop_n*nb_cases_ini)/N),
        Sc = (1 - prop_n) - (((1 - prop_n)*nb_cases_ini)/N),
        In = (prop_n*nb_cases_ini)/N,
        In_test = 0,
        Ic = ((1 - prop_n)*nb_cases_ini)/N,
        Ic_test = 0,
        Ic_trace = 0,
        Rn = 0,
        Rc = 0)


# Time parameters

time_burnin = 30
time_sd = 400
startip = 90 # 
time_sd.ip = time_sd - startip
time_sim = 60
timestep = 0.1
times.burnin <- seq(0, time_burnin, by = timestep)
times.run.sd <- seq(0, time_sd, by = timestep)
times.run.sd.ip <- seq(0, time_sd.ip, by = timestep)
times.noIntervention <- seq(0, time_burnin+time_sd, by = timestep)


plot('', ylim = c(0, 0.17), xlim = c(0, time_burnin+time_sd), ylab = '', xlab = '', cex.axis = 1)

# plot epidemic in the absence of intervention
pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
out.Nointervention <- ode(y = ini,  times = times.noIntervention, func = sir_npi, parms = pars_null, npi = 'sd')  
#lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.Nointervention, lwd = 3, col = 'lightgrey', lty = 1)
x =  out.Nointervention[,'time']
y = rowSums(out.Nointervention[,c('In', 'In_test', 'Ic', 'Ic_test', 'Ic_trace')])
polygon(c(x[x>=0], max(x), 0), c(ifelse(y[x>=0] > 0.155, 0.155, y[x>=0]), 0, 0), col = '#F0F0F0', border = NA)

for(i in 1:length(Es)){
  
  # parameters for dynamic model, no immunity passport
  pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
  #pars_null.dyn<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0,   x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_sd.dyn<-   c(b = b, g = g, a = focus.a, tau = NA,   p_test = NA,  x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_sd.dyn.ip<-c(b = b, g = g, a = focus.a, tau = NA,   p_test = NA,  x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 1)
  
  out.burnin <- ode(y = ini,  times = times.burnin, func = sir_npi, parms = pars_null, npi = 'sd')  
  out.run.sd<- ode(y = tail(out.burnin,1)[,-1], times = times.run.sd, func = sir_npi, parms = pars_sd.dyn, npi = 'sd')
  out.run.sd.ip<- ode(y = out.run.sd[which(out.run.sd[,'time'] == startip),-1], times = times.run.sd.ip, func = sir_npi, parms = pars_sd.dyn.ip, npi = 'sd')
  
  # Assemble output
  out.1<- rbind(out.burnin, out.run.sd[-1,])
  out.2<- rbind(out.burnin, out.run.sd[2:which(out.run.sd[,'time'] == startip),], out.run.sd.ip[-1,])
  out.1[,'time']<- seq(0, tail(time_burnin,1)+tail(time_sd, 1), timestep)
  out.2[,'time']<- seq(0, tail(time_burnin,1)+tail(time_sd, 1), timestep)
  
  # amount of infected (in SD model, _test and _trace = 0)
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.1, lwd = 1.5, col = cols[i], lty = 1)
  out.test<- out.2[which(out.run.sd[,'time'] == I(startip+time_burnin)):nrow(out.2),]
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.test, lwd = 1.5, col = darkcols[i], lty = 2)
  
}
abline(v = time_burnin, lty = 2)
abline(v = time_burnin+startip, lty = 2)
#mtext(side = 2, 'Prevalence', line = 2.2, font = 2, cex = 1)
#mtext(side = 1, 'Days', line = 2.2, font = 2, cex = 1)
#axis(1, at = seq(0, time_burnin+time_sd, 30))
#axis(2, at = c(0, 0.005, 0.01),
#     labels = c(0, 0.005, 0.01))

mtext(side = 3, '+SD', font = 2, cex = 0.9, line =  0.3, at = time_burnin)#, adj = 0)
mtext(side = 3, '+IP 90 days later', font = 2, cex = 0.9, line =  0.3, at = time_burnin+startip, adj = 0)

mtext(side = 4, 'Ini = 100', font = 2, cex = 1, line = 1)




# 4 ----

# Initial conditions
nb_cases_ini = 1000
prop_n = 0.7
ini<- c(Sn = prop_n - ((prop_n*nb_cases_ini)/N),
        Sc = (1 - prop_n) - (((1 - prop_n)*nb_cases_ini)/N),
        In = (prop_n*nb_cases_ini)/N,
        In_test = 0,
        Ic = ((1 - prop_n)*nb_cases_ini)/N,
        Ic_test = 0,
        Ic_trace = 0,
        Rn = 0,
        Rc = 0)


# Time parameters

time_burnin = 30
time_sd = 400
startip = 30 # 
time_sd.ip = time_sd - startip
time_sim = 60
timestep = 0.1
times.burnin <- seq(0, time_burnin, by = timestep)
times.run.sd <- seq(0, time_sd, by = timestep)
times.run.sd.ip <- seq(0, time_sd.ip, by = timestep)
times.noIntervention <- seq(0, time_burnin+time_sd, by = timestep)


plot('', ylim = c(0, 0.17), xlim = c(0, time_burnin+time_sd), ylab = '', xlab = '', cex.axis = 1)

# plot epidemic in the absence of intervention
pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
out.Nointervention <- ode(y = ini,  times = times.noIntervention, func = sir_npi, parms = pars_null, npi = 'sd')  
#lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.Nointervention, lwd = 3, col = 'lightgrey', lty = 1)
x =  out.Nointervention[,'time']
y = rowSums(out.Nointervention[,c('In', 'In_test', 'Ic', 'Ic_test', 'Ic_trace')])
polygon(c(x[x>=0], max(x), 0), c(ifelse(y[x>=0] > 0.155, 0.155, y[x>=0]), 0, 0), col = '#F0F0F0', border = NA)

for(i in 1:length(Es)){
  
  # parameters for dynamic model, no immunity passport
  pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
  #pars_null.dyn<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0,   x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_sd.dyn<-   c(b = b, g = g, a = focus.a, tau = NA,   p_test = NA,  x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_sd.dyn.ip<-c(b = b, g = g, a = focus.a, tau = NA,   p_test = NA,  x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 1)
  
  out.burnin <- ode(y = ini,  times = times.burnin, func = sir_npi, parms = pars_null, npi = 'sd')  
  out.run.sd<- ode(y = tail(out.burnin,1)[,-1], times = times.run.sd, func = sir_npi, parms = pars_sd.dyn, npi = 'sd')
  out.run.sd.ip<- ode(y = out.run.sd[which(out.run.sd[,'time'] == startip),-1], times = times.run.sd.ip, func = sir_npi, parms = pars_sd.dyn.ip, npi = 'sd')
  
  # Assemble output
  out.1<- rbind(out.burnin, out.run.sd[-1,])
  out.2<- rbind(out.burnin, out.run.sd[2:which(out.run.sd[,'time'] == startip),], out.run.sd.ip[-1,])
  out.1[,'time']<- seq(0, tail(time_burnin,1)+tail(time_sd, 1), timestep)
  out.2[,'time']<- seq(0, tail(time_burnin,1)+tail(time_sd, 1), timestep)
  
  # amount of infected (in SD model, _test and _trace = 0)
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.1, lwd = 1.5, col = cols[i], lty = 1)
  out.test<- out.2[which(out.run.sd[,'time'] == I(startip+time_burnin)):nrow(out.2),]
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.test, lwd = 1.5, col = darkcols[i], lty = 2)
  
}
abline(v = time_burnin, lty = 2)
abline(v = time_burnin+startip, lty = 2)
#mtext(side = 2, 'Prevalence', line = 3, font = 2, cex = 1)
#mtext(side = 1, 'Days', line = 2.2, font = 2, cex = 1)
#axis(1, at = seq(0, time_burnin+time_sd, 30))
#axis(2, at = c(0, 0.005, 0.01),
#     labels = c(0, 0.005, 0.01))







# 5 ----

# Initial conditions
nb_cases_ini = 1000
prop_n = 0.7
ini<- c(Sn = prop_n - ((prop_n*nb_cases_ini)/N),
        Sc = (1 - prop_n) - (((1 - prop_n)*nb_cases_ini)/N),
        In = (prop_n*nb_cases_ini)/N,
        In_test = 0,
        Ic = ((1 - prop_n)*nb_cases_ini)/N,
        Ic_test = 0,
        Ic_trace = 0,
        Rn = 0,
        Rc = 0)


# Time parameters

time_burnin = 30
time_sd = 400
startip = 60 # 
time_sd.ip = time_sd - startip
time_sim = 60
timestep = 0.1
times.burnin <- seq(0, time_burnin, by = timestep)
times.run.sd <- seq(0, time_sd, by = timestep)
times.run.sd.ip <- seq(0, time_sd.ip, by = timestep)
times.noIntervention <- seq(0, time_burnin+time_sd, by = timestep)


plot('', ylim = c(0, 0.17), xlim = c(0, time_burnin+time_sd), ylab = '', xlab = '', cex.axis = 1)

# plot epidemic in the absence of intervention
pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
out.Nointervention <- ode(y = ini,  times = times.noIntervention, func = sir_npi, parms = pars_null, npi = 'sd')  
#lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.Nointervention, lwd = 3, col = 'lightgrey', lty = 1)
x =  out.Nointervention[,'time']
y = rowSums(out.Nointervention[,c('In', 'In_test', 'Ic', 'Ic_test', 'Ic_trace')])
polygon(c(x[x>=0], max(x), 0), c(ifelse(y[x>=0] > 0.155, 0.155, y[x>=0]), 0, 0), col = '#F0F0F0', border = NA)

for(i in 1:length(Es)){
  
  # parameters for dynamic model, no immunity passport
  pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
  #pars_null.dyn<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0,   x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_sd.dyn<-   c(b = b, g = g, a = focus.a, tau = NA,   p_test = NA,  x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_sd.dyn.ip<-c(b = b, g = g, a = focus.a, tau = NA,   p_test = NA,  x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 1)
  
  out.burnin <- ode(y = ini,  times = times.burnin, func = sir_npi, parms = pars_null, npi = 'sd')  
  out.run.sd<- ode(y = tail(out.burnin,1)[,-1], times = times.run.sd, func = sir_npi, parms = pars_sd.dyn, npi = 'sd')
  out.run.sd.ip<- ode(y = out.run.sd[which(out.run.sd[,'time'] == startip),-1], times = times.run.sd.ip, func = sir_npi, parms = pars_sd.dyn.ip, npi = 'sd')
  
  # Assemble output
  out.1<- rbind(out.burnin, out.run.sd[-1,])
  out.2<- rbind(out.burnin, out.run.sd[2:which(out.run.sd[,'time'] == startip),], out.run.sd.ip[-1,])
  out.1[,'time']<- seq(0, tail(time_burnin,1)+tail(time_sd, 1), timestep)
  out.2[,'time']<- seq(0, tail(time_burnin,1)+tail(time_sd, 1), timestep)
  
  # amount of infected (in SD model, _test and _trace = 0)
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.1, lwd = 1.5, col = cols[i], lty = 1)
  out.test<- out.2[which(out.run.sd[,'time'] == I(startip+time_burnin)):nrow(out.2),]
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.test, lwd = 1.5, col = darkcols[i], lty = 2)
  
}
abline(v = time_burnin, lty = 2)
abline(v = time_burnin+startip, lty = 2)
#mtext(side = 2, 'Prevalence', line = 2.2, font = 2, cex = 1)
#mtext(side = 1, 'Days', line = 2.2, font = 2, cex = 1)
#axis(1, at = seq(0, time_burnin+time_sd, 30))
#axis(2, at = c(0, 0.005, 0.01),
#     labels = c(0, 0.005, 0.01))



# 6 ----

# Initial conditions
nb_cases_ini = 1000
prop_n = 0.7
ini<- c(Sn = prop_n - ((prop_n*nb_cases_ini)/N),
        Sc = (1 - prop_n) - (((1 - prop_n)*nb_cases_ini)/N),
        In = (prop_n*nb_cases_ini)/N,
        In_test = 0,
        Ic = ((1 - prop_n)*nb_cases_ini)/N,
        Ic_test = 0,
        Ic_trace = 0,
        Rn = 0,
        Rc = 0)


# Time parameters

time_burnin = 30
time_sd = 400
startip = 90 # 
time_sd.ip = time_sd - startip
time_sim = 60
timestep = 0.1
times.burnin <- seq(0, time_burnin, by = timestep)
times.run.sd <- seq(0, time_sd, by = timestep)
times.run.sd.ip <- seq(0, time_sd.ip, by = timestep)
times.noIntervention <- seq(0, time_burnin+time_sd, by = timestep)


plot('', ylim = c(0, 0.17), xlim = c(0, time_burnin+time_sd), ylab = '', xlab = '', cex.axis = 1)

# plot epidemic in the absence of intervention
pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
out.Nointervention <- ode(y = ini,  times = times.noIntervention, func = sir_npi, parms = pars_null, npi = 'sd')  
#lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.Nointervention, lwd = 3, col = 'lightgrey', lty = 1)
x =  out.Nointervention[,'time']
y = rowSums(out.Nointervention[,c('In', 'In_test', 'Ic', 'Ic_test', 'Ic_trace')])
polygon(c(x[x>=0], max(x), 0), c(ifelse(y[x>=0] > 0.155, 0.155, y[x>=0]), 0, 0), col = '#F0F0F0', border = NA)

for(i in 1:length(Es)){
  
  # parameters for dynamic model, no immunity passport
  pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
  #pars_null.dyn<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0,   x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_sd.dyn<-   c(b = b, g = g, a = focus.a, tau = NA,   p_test = NA,  x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_sd.dyn.ip<-c(b = b, g = g, a = focus.a, tau = NA,   p_test = NA,  x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 1)
  
  out.burnin <- ode(y = ini,  times = times.burnin, func = sir_npi, parms = pars_null, npi = 'sd')  
  out.run.sd<- ode(y = tail(out.burnin,1)[,-1], times = times.run.sd, func = sir_npi, parms = pars_sd.dyn, npi = 'sd')
  out.run.sd.ip<- ode(y = out.run.sd[which(out.run.sd[,'time'] == startip),-1], times = times.run.sd.ip, func = sir_npi, parms = pars_sd.dyn.ip, npi = 'sd')
  
  # Assemble output
  out.1<- rbind(out.burnin, out.run.sd[-1,])
  out.2<- rbind(out.burnin, out.run.sd[2:which(out.run.sd[,'time'] == startip),], out.run.sd.ip[-1,])
  out.1[,'time']<- seq(0, tail(time_burnin,1)+tail(time_sd, 1), timestep)
  out.2[,'time']<- seq(0, tail(time_burnin,1)+tail(time_sd, 1), timestep)
  
  # amount of infected (in SD model, _test and _trace = 0)
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.1, lwd = 1.5, col = cols[i], lty = 1)
  out.test<- out.2[which(out.run.sd[,'time'] == I(startip+time_burnin)):nrow(out.2),]
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.test, lwd = 1.5, col = darkcols[i], lty = 2)
  
}
abline(v = time_burnin, lty = 2)
abline(v = time_burnin+startip, lty = 2)
#mtext(side = 2, 'Prevalence', line = 2.2, font = 2, cex = 1)
#mtext(side = 1, 'Days', line = 2.2, font = 2, cex = 1)
#axis(1, at = seq(0, time_burnin+time_sd, 30))
#axis(2, at = c(0, 0.005, 0.01),
#     labels = c(0, 0.005, 0.01))

mtext(side = 4, 'Ini = 1000', font = 2, cex = 1, line = 1)






# 7 ----

# Initial conditions
nb_cases_ini = 10000
prop_n = 0.7
ini<- c(Sn = prop_n - ((prop_n*nb_cases_ini)/N),
        Sc = (1 - prop_n) - (((1 - prop_n)*nb_cases_ini)/N),
        In = (prop_n*nb_cases_ini)/N,
        In_test = 0,
        Ic = ((1 - prop_n)*nb_cases_ini)/N,
        Ic_test = 0,
        Ic_trace = 0,
        Rn = 0,
        Rc = 0)


# Time parameters

time_burnin = 30
time_sd = 400
startip = 30 # 
time_sd.ip = time_sd - startip
time_sim = 60
timestep = 0.1
times.burnin <- seq(0, time_burnin, by = timestep)
times.run.sd <- seq(0, time_sd, by = timestep)
times.run.sd.ip <- seq(0, time_sd.ip, by = timestep)
times.noIntervention <- seq(0, time_burnin+time_sd, by = timestep)


plot('', ylim = c(0, 0.17), xlim = c(0, time_burnin+time_sd), ylab = '', xlab = '', cex.axis = 1)

# plot epidemic in the absence of intervention
pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
out.Nointervention <- ode(y = ini,  times = times.noIntervention, func = sir_npi, parms = pars_null, npi = 'sd')  
#lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.Nointervention, lwd = 3, col = 'lightgrey', lty = 1)
x =  out.Nointervention[,'time']
y = rowSums(out.Nointervention[,c('In', 'In_test', 'Ic', 'Ic_test', 'Ic_trace')])
polygon(c(x[x>=0], max(x), 0), c(ifelse(y[x>=0] > 0.155, 0.155, y[x>=0]), 0, 0), col = '#F0F0F0', border = NA)

for(i in 1:length(Es)){
  
  # parameters for dynamic model, no immunity passport
  pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
  #pars_null.dyn<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0,   x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_sd.dyn<-   c(b = b, g = g, a = focus.a, tau = NA,   p_test = NA,  x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_sd.dyn.ip<-c(b = b, g = g, a = focus.a, tau = NA,   p_test = NA,  x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 1)
  
  out.burnin <- ode(y = ini,  times = times.burnin, func = sir_npi, parms = pars_null, npi = 'sd')  
  out.run.sd<- ode(y = tail(out.burnin,1)[,-1], times = times.run.sd, func = sir_npi, parms = pars_sd.dyn, npi = 'sd')
  out.run.sd.ip<- ode(y = out.run.sd[which(out.run.sd[,'time'] == startip),-1], times = times.run.sd.ip, func = sir_npi, parms = pars_sd.dyn.ip, npi = 'sd')
  
  # Assemble output
  out.1<- rbind(out.burnin, out.run.sd[-1,])
  out.2<- rbind(out.burnin, out.run.sd[2:which(out.run.sd[,'time'] == startip),], out.run.sd.ip[-1,])
  out.1[,'time']<- seq(0, tail(time_burnin,1)+tail(time_sd, 1), timestep)
  out.2[,'time']<- seq(0, tail(time_burnin,1)+tail(time_sd, 1), timestep)
  
  # amount of infected (in SD model, _test and _trace = 0)
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.1, lwd = 1.5, col = cols[i], lty = 1)
  out.test<- out.2[which(out.run.sd[,'time'] == I(startip+time_burnin)):nrow(out.2),]
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.test, lwd = 1.5, col = darkcols[i], lty = 2)
  
}
abline(v = time_burnin, lty = 2)
abline(v = time_burnin+startip, lty = 2)
#mtext(side = 2, 'Prevalence', line = 2.2, font = 2, cex = 1)
#mtext(side = 1, 'Days', line = 2.2, font = 2, cex = 1)
#axis(1, at = seq(0, time_burnin+time_sd, 30))
#axis(2, at = c(0, 0.005, 0.01),
#     labels = c(0, 0.005, 0.01))







# 8 ----

# Initial conditions
nb_cases_ini = 10000
prop_n = 0.7
ini<- c(Sn = prop_n - ((prop_n*nb_cases_ini)/N),
        Sc = (1 - prop_n) - (((1 - prop_n)*nb_cases_ini)/N),
        In = (prop_n*nb_cases_ini)/N,
        In_test = 0,
        Ic = ((1 - prop_n)*nb_cases_ini)/N,
        Ic_test = 0,
        Ic_trace = 0,
        Rn = 0,
        Rc = 0)


# Time parameters

time_burnin = 30
time_sd = 400
startip = 60 # 
time_sd.ip = time_sd - startip
time_sim = 60
timestep = 0.1
times.burnin <- seq(0, time_burnin, by = timestep)
times.run.sd <- seq(0, time_sd, by = timestep)
times.run.sd.ip <- seq(0, time_sd.ip, by = timestep)
times.noIntervention <- seq(0, time_burnin+time_sd, by = timestep)


plot('', ylim = c(0, 0.17), xlim = c(0, time_burnin+time_sd), ylab = '', xlab = '', cex.axis = 1)

# plot epidemic in the absence of intervention
pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
out.Nointervention <- ode(y = ini,  times = times.noIntervention, func = sir_npi, parms = pars_null, npi = 'sd')  
#lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.Nointervention, lwd = 3, col = 'lightgrey', lty = 1)
x =  out.Nointervention[,'time']
y = rowSums(out.Nointervention[,c('In', 'In_test', 'Ic', 'Ic_test', 'Ic_trace')])
polygon(c(x[x>=0], max(x), 0), c(ifelse(y[x>=0] > 0.155, 0.155, y[x>=0]), 0, 0), col = '#F0F0F0', border = NA)

for(i in 1:length(Es)){
  
  # parameters for dynamic model, no immunity passport
  pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
  #pars_null.dyn<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0,   x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_sd.dyn<-   c(b = b, g = g, a = focus.a, tau = NA,   p_test = NA,  x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_sd.dyn.ip<-c(b = b, g = g, a = focus.a, tau = NA,   p_test = NA,  x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 1)
  
  out.burnin <- ode(y = ini,  times = times.burnin, func = sir_npi, parms = pars_null, npi = 'sd')  
  out.run.sd<- ode(y = tail(out.burnin,1)[,-1], times = times.run.sd, func = sir_npi, parms = pars_sd.dyn, npi = 'sd')
  out.run.sd.ip<- ode(y = out.run.sd[which(out.run.sd[,'time'] == startip),-1], times = times.run.sd.ip, func = sir_npi, parms = pars_sd.dyn.ip, npi = 'sd')
  
  # Assemble output
  out.1<- rbind(out.burnin, out.run.sd[-1,])
  out.2<- rbind(out.burnin, out.run.sd[2:which(out.run.sd[,'time'] == startip),], out.run.sd.ip[-1,])
  out.1[,'time']<- seq(0, tail(time_burnin,1)+tail(time_sd, 1), timestep)
  out.2[,'time']<- seq(0, tail(time_burnin,1)+tail(time_sd, 1), timestep)
  
  # amount of infected (in SD model, _test and _trace = 0)
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.1, lwd = 1.5, col = cols[i], lty = 1)
  out.test<- out.2[which(out.run.sd[,'time'] == I(startip+time_burnin)):nrow(out.2),]
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.test, lwd = 1.5, col = darkcols[i], lty = 2)
  
}
abline(v = time_burnin, lty = 2)
abline(v = time_burnin+startip, lty = 2)
#mtext(side = 2, 'Prevalence', line = 2.2, font = 2, cex = 1)
#mtext(side = 1, 'Days', line = 2, font = 2, cex = 1)
#axis(1, at = seq(0, time_burnin+time_sd, 30))
#axis(2, at = c(0, 0.005, 0.01),
#     labels = c(0, 0.005, 0.01))



# 9 ----

# Initial conditions
nb_cases_ini = 10000
prop_n = 0.7
ini<- c(Sn = prop_n - ((prop_n*nb_cases_ini)/N),
        Sc = (1 - prop_n) - (((1 - prop_n)*nb_cases_ini)/N),
        In = (prop_n*nb_cases_ini)/N,
        In_test = 0,
        Ic = ((1 - prop_n)*nb_cases_ini)/N,
        Ic_test = 0,
        Ic_trace = 0,
        Rn = 0,
        Rc = 0)


# Time parameters

time_burnin = 30
time_sd = 400
startip = 90 # 
time_sd.ip = time_sd - startip
time_sim = 60
timestep = 0.1
times.burnin <- seq(0, time_burnin, by = timestep)
times.run.sd <- seq(0, time_sd, by = timestep)
times.run.sd.ip <- seq(0, time_sd.ip, by = timestep)
times.noIntervention <- seq(0, time_burnin+time_sd, by = timestep)


plot('', ylim = c(0, 0.17), xlim = c(0, time_burnin+time_sd), ylab = '', xlab = '', cex.axis = 1)

# plot epidemic in the absence of intervention
pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
out.Nointervention <- ode(y = ini,  times = times.noIntervention, func = sir_npi, parms = pars_null, npi = 'sd')  
#lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.Nointervention, lwd = 3, col = 'lightgrey', lty = 1)
x =  out.Nointervention[,'time']
y = rowSums(out.Nointervention[,c('In', 'In_test', 'Ic', 'Ic_test', 'Ic_trace')])
polygon(c(x[x>=0], max(x), 0), c(ifelse(y[x>=0] > 0.155, 0.155, y[x>=0]), 0, 0), col = '#F0F0F0', border = NA)

for(i in 1:length(Es)){
  
  # parameters for dynamic model, no immunity passport
  pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
  #pars_null.dyn<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0,   x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_sd.dyn<-   c(b = b, g = g, a = focus.a, tau = NA,   p_test = NA,  x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_sd.dyn.ip<-c(b = b, g = g, a = focus.a, tau = NA,   p_test = NA,  x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 1)
  
  out.burnin <- ode(y = ini,  times = times.burnin, func = sir_npi, parms = pars_null, npi = 'sd')  
  out.run.sd<- ode(y = tail(out.burnin,1)[,-1], times = times.run.sd, func = sir_npi, parms = pars_sd.dyn, npi = 'sd')
  out.run.sd.ip<- ode(y = out.run.sd[which(out.run.sd[,'time'] == startip),-1], times = times.run.sd.ip, func = sir_npi, parms = pars_sd.dyn.ip, npi = 'sd')
  
  # Assemble output
  out.1<- rbind(out.burnin, out.run.sd[-1,])
  out.2<- rbind(out.burnin, out.run.sd[2:which(out.run.sd[,'time'] == startip),], out.run.sd.ip[-1,])
  out.1[,'time']<- seq(0, tail(time_burnin,1)+tail(time_sd, 1), timestep)
  out.2[,'time']<- seq(0, tail(time_burnin,1)+tail(time_sd, 1), timestep)
  
  # amount of infected (in SD model, _test and _trace = 0)
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.1, lwd = 1.5, col = cols[i], lty = 1)
  out.test<- out.2[which(out.run.sd[,'time'] == I(startip+time_burnin)):nrow(out.2),]
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.test, lwd = 1.5, col = darkcols[i], lty = 2)
  
}
abline(v = time_burnin, lty = 2)
abline(v = time_burnin+startip, lty = 2)
#mtext(side = 2, 'Prevalence', line = 2.2, font = 2, cex = 1)
#mtext(side = 1, 'Days', line = 2.2, font = 2, cex = 1)
#axis(1, at = seq(0, time_burnin+time_sd, 30))
#axis(2, at = c(0, 0.005, 0.01),
#     labels = c(0, 0.005, 0.01))

mtext(side = 4, 'Ini = 10 000', font = 2, cex = 1, line = 1)






dev.off()



# [3.2] TTI DYNAMICS VARY TIMING ----

# Same here, empty environeemnt and re-run part [1] first to make sure no mess up between parts of codes

cols = rev(c('#FFE871', '#FF7171', '#82C1FF', '#8ABEA1'))
darkcols = rev(c('#FFF1AA', '#FFAAAA', '#B4DAFF', '#B9D8C7'))


focus.x = 0.3
focus.C = 1.2
#focus.a = 4
focus.ptest = 0.7
focus.tau = 5.1

#Es = c(1, 1.2, 1.3, 1.4)
Es = c(1.1, 1.2, 1.25, 1.3)

pdf('./output/figures/FigureS6_TTI_dynamic_vary.pdf', width = (11+1+1)/2.54, height = (11+1+1)/2.54, pointsize = 7)


par(mfrow = c(3,3),
    mar = rep(2.5,4))


# 1 ----
# Initial conditions
nb_cases_ini = 100
prop_n = 0.7
ini<- c(Sn = prop_n - ((prop_n*nb_cases_ini)/N),
        Sc = (1 - prop_n) - (((1 - prop_n)*nb_cases_ini)/N),
        In = (prop_n*nb_cases_ini)/N,
        In_test = 0,
        Ic = ((1 - prop_n)*nb_cases_ini)/N,
        Ic_test = 0,
        Ic_trace = 0,
        Rn = 0,
        Rc = 0)


# Time parameters

time_burnin = 30
time_NPI = 400
startip = 30 # IP introduced seven days after implementing Social distancing
time_NPI.ip = time_NPI - startip
time_sim = 60
timestep = 0.1
times.burnin <- seq(0, time_burnin, by = timestep)
times.run.NPI <- seq(0, time_NPI, by = timestep)
times.run.NPI.ip <- seq(0, time_NPI.ip, by = timestep)
times.noIntervention <- seq(0, time_burnin+time_NPI, by = timestep)


plot('', ylim = c(0, 0.15), xlim = c(0, time_NPI+time_burnin), ylab = '', xlab = '', cex.axis = 1)
# plot epidemic in the absence of intervention
pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
out.Nointervention <- ode(y = ini,  times = times.noIntervention, func = sir_npi, parms = pars_null, npi = 'tti')  
#lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.Nointervention, lwd = 3, col = 'lightgrey', lty = 1)
x =  out.Nointervention[,'time']
y = rowSums(out.Nointervention[,c('In', 'In_test', 'Ic', 'Ic_test', 'Ic_trace')])
polygon(c(x[x>=0], max(x), 0), c(ifelse(y[x>=0] > 0.155, 0.155, y[x>=0]), 0, 0), col = '#F0F0F0', border = NA)


for(i in 1:length(Es)){
  
  # parameters for dynamic model, no immunity passport
  pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
  #pars_null.dyn<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0,   x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_tti.dyn<-  c(b = b, g = g, a = NA,  tau = focus.tau,  p_test = focus.ptest, x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_tti.dyn.ip<-  c(b = b, g = g, a = NA,  tau = focus.tau,  p_test = focus.ptest, x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 1)
  
  out.burnin <- ode(y = ini,  times = times.burnin, func = sir_npi, parms = pars_null, npi = 'tti')  
  out.run.NPI<- ode(y = tail(out.burnin,1)[,-1], times = times.run.NPI, func = sir_npi, parms = pars_tti.dyn, npi = 'tti')
  out.run.NPI.ip<- ode(y = out.run.NPI[which(out.run.NPI[,'time'] == startip),-1], times = times.run.NPI.ip, func = sir_npi, parms = pars_tti.dyn.ip, npi = 'tti')
  
  # Assemble output
  out.1<- rbind(out.burnin, out.run.NPI[-1,])
  out.2<- rbind(out.burnin, out.run.NPI[2:which(out.run.NPI[,'time'] == startip),], out.run.NPI.ip[-1,])
  out.1[,'time']<- seq(0, tail(time_burnin,1)+tail(time_NPI, 1), timestep)
  out.2[,'time']<- seq(0, tail(time_burnin,1)+tail(time_NPI, 1), timestep)
  
  # amount of infected (in SD model, _test and _trace = 0)
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.1, lwd = 1.5, col = cols[i], lty = 1)
  out.test<- out.2[which(out.run.NPI[,'time'] == I(startip+time_burnin)):nrow(out.2),]
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.test, lwd = 1.5, col = darkcols[i], lty = 2)
  
}

abline(v = time_burnin, lty = 2)
abline(v = time_burnin+startip, lty = 2)
#mtext(side = 2, 'Prevalence', line = 2.2, font = 2, cex = 1)
#mtext(side = 1, 'Days', line = 2.2, font = 2, cex = 1)
#axis(1, at = seq(0, time_NPI+time_burnin, 30))

mtext(side = 3, '+TTI', font = 2, cex = 0.9, line =  0.3, at = time_burnin)#, adj = 0)
mtext(side = 3, '+IP 30 days later', font = 2, cex = 0.9, line =  0.3, at = time_burnin+startip, adj = 0)




# 2 ----
# Initial conditions
nb_cases_ini = 100
prop_n = 0.7
ini<- c(Sn = prop_n - ((prop_n*nb_cases_ini)/N),
        Sc = (1 - prop_n) - (((1 - prop_n)*nb_cases_ini)/N),
        In = (prop_n*nb_cases_ini)/N,
        In_test = 0,
        Ic = ((1 - prop_n)*nb_cases_ini)/N,
        Ic_test = 0,
        Ic_trace = 0,
        Rn = 0,
        Rc = 0)


# Time parameters

time_burnin = 30
time_NPI = 400
startip = 60 # IP introduced seven days after implementing Social distancing
time_NPI.ip = time_NPI - startip
time_sim = 60
timestep = 0.1
times.burnin <- seq(0, time_burnin, by = timestep)
times.run.NPI <- seq(0, time_NPI, by = timestep)
times.run.NPI.ip <- seq(0, time_NPI.ip, by = timestep)
times.noIntervention <- seq(0, time_burnin+time_NPI, by = timestep)


plot('', ylim = c(0, 0.15), xlim = c(0, time_NPI+time_burnin), ylab = '', xlab = '', cex.axis = 1)
# plot epidemic in the absence of intervention
pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
out.Nointervention <- ode(y = ini,  times = times.noIntervention, func = sir_npi, parms = pars_null, npi = 'tti')  
#lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.Nointervention, lwd = 3, col = 'lightgrey', lty = 1)
x =  out.Nointervention[,'time']
y = rowSums(out.Nointervention[,c('In', 'In_test', 'Ic', 'Ic_test', 'Ic_trace')])
polygon(c(x[x>=0], max(x), 0), c(ifelse(y[x>=0] > 0.155, 0.155, y[x>=0]), 0, 0), col = '#F0F0F0', border = NA)


for(i in 1:length(Es)){
  
  # parameters for dynamic model, no immunity passport
  pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
  #pars_null.dyn<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0,   x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_tti.dyn<-  c(b = b, g = g, a = NA,  tau = focus.tau,  p_test = focus.ptest, x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_tti.dyn.ip<-  c(b = b, g = g, a = NA,  tau = focus.tau,  p_test = focus.ptest, x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 1)
  
  out.burnin <- ode(y = ini,  times = times.burnin, func = sir_npi, parms = pars_null, npi = 'tti')  
  out.run.NPI<- ode(y = tail(out.burnin,1)[,-1], times = times.run.NPI, func = sir_npi, parms = pars_tti.dyn, npi = 'tti')
  out.run.NPI.ip<- ode(y = out.run.NPI[which(out.run.NPI[,'time'] == startip),-1], times = times.run.NPI.ip, func = sir_npi, parms = pars_tti.dyn.ip, npi = 'tti')
  
  # Assemble output
  out.1<- rbind(out.burnin, out.run.NPI[-1,])
  out.2<- rbind(out.burnin, out.run.NPI[2:which(out.run.NPI[,'time'] == startip),], out.run.NPI.ip[-1,])
  out.1[,'time']<- seq(0, tail(time_burnin,1)+tail(time_NPI, 1), timestep)
  out.2[,'time']<- seq(0, tail(time_burnin,1)+tail(time_NPI, 1), timestep)
  
  # amount of infected (in SD model, _test and _trace = 0)
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.1, lwd = 1.5, col = cols[i], lty = 1)
  out.test<- out.2[which(out.run.NPI[,'time'] == I(startip+time_burnin)):nrow(out.2),]
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.test, lwd = 1.5, col = darkcols[i], lty = 2)
  
}

abline(v = time_burnin, lty = 2)
abline(v = time_burnin+startip, lty = 2)
#mtext(side = 2, 'Prevalence', line = 2.2, font = 2, cex = 1)
#mtext(side = 1, 'Days', line = 2.2, font = 2, cex = 1)
#axis(1, at = seq(0, time_NPI+time_burnin, 30))

mtext(side = 3, '+TTI', font = 2, cex = 0.9, line =  0.3, at = time_burnin)#, adj = 0)
mtext(side = 3, '+IP 60 days later', font = 2, cex = 0.9, line =  0.3, at = time_burnin+startip, adj = 0)


# 3 ----
# Initial conditions
nb_cases_ini = 100
prop_n = 0.7
ini<- c(Sn = prop_n - ((prop_n*nb_cases_ini)/N),
        Sc = (1 - prop_n) - (((1 - prop_n)*nb_cases_ini)/N),
        In = (prop_n*nb_cases_ini)/N,
        In_test = 0,
        Ic = ((1 - prop_n)*nb_cases_ini)/N,
        Ic_test = 0,
        Ic_trace = 0,
        Rn = 0,
        Rc = 0)


# Time parameters

time_burnin = 30
time_NPI = 400
startip = 90 # IP introduced seven days after implementing Social distancing
time_NPI.ip = time_NPI - startip
time_sim = 60
timestep = 0.1
times.burnin <- seq(0, time_burnin, by = timestep)
times.run.NPI <- seq(0, time_NPI, by = timestep)
times.run.NPI.ip <- seq(0, time_NPI.ip, by = timestep)
times.noIntervention <- seq(0, time_burnin+time_NPI, by = timestep)


plot('', ylim = c(0, 0.15), xlim = c(0, time_NPI+time_burnin), ylab = '', xlab = '', cex.axis = 1)
# plot epidemic in the absence of intervention
pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
out.Nointervention <- ode(y = ini,  times = times.noIntervention, func = sir_npi, parms = pars_null, npi = 'tti')  
#lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.Nointervention, lwd = 3, col = 'lightgrey', lty = 1)
x =  out.Nointervention[,'time']
y = rowSums(out.Nointervention[,c('In', 'In_test', 'Ic', 'Ic_test', 'Ic_trace')])
polygon(c(x[x>=0], max(x), 0), c(ifelse(y[x>=0] > 0.155, 0.155, y[x>=0]), 0, 0), col = '#F0F0F0', border = NA)


for(i in 1:length(Es)){
  
  # parameters for dynamic model, no immunity passport
  pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
  #pars_null.dyn<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0,   x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_tti.dyn<-  c(b = b, g = g, a = NA,  tau = focus.tau,  p_test = focus.ptest, x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_tti.dyn.ip<-  c(b = b, g = g, a = NA,  tau = focus.tau,  p_test = focus.ptest, x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 1)
  
  out.burnin <- ode(y = ini,  times = times.burnin, func = sir_npi, parms = pars_null, npi = 'tti')  
  out.run.NPI<- ode(y = tail(out.burnin,1)[,-1], times = times.run.NPI, func = sir_npi, parms = pars_tti.dyn, npi = 'tti')
  out.run.NPI.ip<- ode(y = out.run.NPI[which(out.run.NPI[,'time'] == startip),-1], times = times.run.NPI.ip, func = sir_npi, parms = pars_tti.dyn.ip, npi = 'tti')
  
  # Assemble output
  out.1<- rbind(out.burnin, out.run.NPI[-1,])
  out.2<- rbind(out.burnin, out.run.NPI[2:which(out.run.NPI[,'time'] == startip),], out.run.NPI.ip[-1,])
  out.1[,'time']<- seq(0, tail(time_burnin,1)+tail(time_NPI, 1), timestep)
  out.2[,'time']<- seq(0, tail(time_burnin,1)+tail(time_NPI, 1), timestep)
  
  # amount of infected (in SD model, _test and _trace = 0)
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.1, lwd = 1.5, col = cols[i], lty = 1)
  out.test<- out.2[which(out.run.NPI[,'time'] == I(startip+time_burnin)):nrow(out.2),]
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.test, lwd = 1.5, col = darkcols[i], lty = 2)
  
}

abline(v = time_burnin, lty = 2)
abline(v = time_burnin+startip, lty = 2)
#mtext(side = 2, 'Prevalence', line = 2.2, font = 2, cex = 1)
#mtext(side = 1, 'Days', line = 2.2, font = 2, cex = 1)
#axis(1, at = seq(0, time_NPI+time_burnin, 30))


mtext(side = 3, '+TTI', font = 2, cex = 0.9, line =  0.3, at = time_burnin)#, adj = 0)
mtext(side = 3, '+IP 90 days later', font = 2, cex = 0.9, line =  0.3, at = time_burnin+startip, adj = 0)

mtext(side = 4, 'Ini = 100', font = 2, cex = 1, line = 1)





# 4 ----
# Initial conditions
nb_cases_ini = 1000
prop_n = 0.7
ini<- c(Sn = prop_n - ((prop_n*nb_cases_ini)/N),
        Sc = (1 - prop_n) - (((1 - prop_n)*nb_cases_ini)/N),
        In = (prop_n*nb_cases_ini)/N,
        In_test = 0,
        Ic = ((1 - prop_n)*nb_cases_ini)/N,
        Ic_test = 0,
        Ic_trace = 0,
        Rn = 0,
        Rc = 0)


# Time parameters

time_burnin = 30
time_NPI = 400
startip = 30 # IP introduced seven days after implementing Social distancing
time_NPI.ip = time_NPI - startip
time_sim = 60
timestep = 0.1
times.burnin <- seq(0, time_burnin, by = timestep)
times.run.NPI <- seq(0, time_NPI, by = timestep)
times.run.NPI.ip <- seq(0, time_NPI.ip, by = timestep)
times.noIntervention <- seq(0, time_burnin+time_NPI, by = timestep)


plot('', ylim = c(0, 0.15), xlim = c(0, time_NPI+time_burnin), ylab = '', xlab = '', cex.axis = 1)
# plot epidemic in the absence of intervention
pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
out.Nointervention <- ode(y = ini,  times = times.noIntervention, func = sir_npi, parms = pars_null, npi = 'tti')  
#lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.Nointervention, lwd = 3, col = 'lightgrey', lty = 1)
x =  out.Nointervention[,'time']
y = rowSums(out.Nointervention[,c('In', 'In_test', 'Ic', 'Ic_test', 'Ic_trace')])
polygon(c(x[x>=0], max(x), 0), c(ifelse(y[x>=0] > 0.155, 0.155, y[x>=0]), 0, 0), col = '#F0F0F0', border = NA)


for(i in 1:length(Es)){
  
  # parameters for dynamic model, no immunity passport
  pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
  #pars_null.dyn<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0,   x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_tti.dyn<-  c(b = b, g = g, a = NA,  tau = focus.tau,  p_test = focus.ptest, x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_tti.dyn.ip<-  c(b = b, g = g, a = NA,  tau = focus.tau,  p_test = focus.ptest, x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 1)
  
  out.burnin <- ode(y = ini,  times = times.burnin, func = sir_npi, parms = pars_null, npi = 'tti')  
  out.run.NPI<- ode(y = tail(out.burnin,1)[,-1], times = times.run.NPI, func = sir_npi, parms = pars_tti.dyn, npi = 'tti')
  out.run.NPI.ip<- ode(y = out.run.NPI[which(out.run.NPI[,'time'] == startip),-1], times = times.run.NPI.ip, func = sir_npi, parms = pars_tti.dyn.ip, npi = 'tti')
  
  # Assemble output
  out.1<- rbind(out.burnin, out.run.NPI[-1,])
  out.2<- rbind(out.burnin, out.run.NPI[2:which(out.run.NPI[,'time'] == startip),], out.run.NPI.ip[-1,])
  out.1[,'time']<- seq(0, tail(time_burnin,1)+tail(time_NPI, 1), timestep)
  out.2[,'time']<- seq(0, tail(time_burnin,1)+tail(time_NPI, 1), timestep)
  
  # amount of infected (in SD model, _test and _trace = 0)
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.1, lwd = 1.5, col = cols[i], lty = 1)
  out.test<- out.2[which(out.run.NPI[,'time'] == I(startip+time_burnin)):nrow(out.2),]
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.test, lwd = 1.5, col = darkcols[i], lty = 2)
  
}

abline(v = time_burnin, lty = 2)
abline(v = time_burnin+startip, lty = 2)
#mtext(side = 2, 'Prevalence', line = 3, font = 2, cex = 1)
#mtext(side = 1, 'Days', line = 2.2, font = 2, cex = 1)
#axis(1, at = seq(0, time_NPI+time_burnin, 30))





# 5 ----
# Initial conditions
nb_cases_ini = 1000
prop_n = 0.7
ini<- c(Sn = prop_n - ((prop_n*nb_cases_ini)/N),
        Sc = (1 - prop_n) - (((1 - prop_n)*nb_cases_ini)/N),
        In = (prop_n*nb_cases_ini)/N,
        In_test = 0,
        Ic = ((1 - prop_n)*nb_cases_ini)/N,
        Ic_test = 0,
        Ic_trace = 0,
        Rn = 0,
        Rc = 0)


# Time parameters

time_burnin = 30
time_NPI = 400
startip = 60 # IP introduced seven days after implementing Social distancing
time_NPI.ip = time_NPI - startip
time_sim = 60
timestep = 0.1
times.burnin <- seq(0, time_burnin, by = timestep)
times.run.NPI <- seq(0, time_NPI, by = timestep)
times.run.NPI.ip <- seq(0, time_NPI.ip, by = timestep)
times.noIntervention <- seq(0, time_burnin+time_NPI, by = timestep)


plot('', ylim = c(0, 0.15), xlim = c(0, time_NPI+time_burnin), ylab = '', xlab = '', cex.axis = 1)
# plot epidemic in the absence of intervention
pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
out.Nointervention <- ode(y = ini,  times = times.noIntervention, func = sir_npi, parms = pars_null, npi = 'tti')  
#lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.Nointervention, lwd = 3, col = 'lightgrey', lty = 1)
x =  out.Nointervention[,'time']
y = rowSums(out.Nointervention[,c('In', 'In_test', 'Ic', 'Ic_test', 'Ic_trace')])
polygon(c(x[x>=0], max(x), 0), c(ifelse(y[x>=0] > 0.155, 0.155, y[x>=0]), 0, 0), col = '#F0F0F0', border = NA)


for(i in 1:length(Es)){
  
  # parameters for dynamic model, no immunity passport
  pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
  #pars_null.dyn<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0,   x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_tti.dyn<-  c(b = b, g = g, a = NA,  tau = focus.tau,  p_test = focus.ptest, x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_tti.dyn.ip<-  c(b = b, g = g, a = NA,  tau = focus.tau,  p_test = focus.ptest, x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 1)
  
  out.burnin <- ode(y = ini,  times = times.burnin, func = sir_npi, parms = pars_null, npi = 'tti')  
  out.run.NPI<- ode(y = tail(out.burnin,1)[,-1], times = times.run.NPI, func = sir_npi, parms = pars_tti.dyn, npi = 'tti')
  out.run.NPI.ip<- ode(y = out.run.NPI[which(out.run.NPI[,'time'] == startip),-1], times = times.run.NPI.ip, func = sir_npi, parms = pars_tti.dyn.ip, npi = 'tti')
  
  # Assemble output
  out.1<- rbind(out.burnin, out.run.NPI[-1,])
  out.2<- rbind(out.burnin, out.run.NPI[2:which(out.run.NPI[,'time'] == startip),], out.run.NPI.ip[-1,])
  out.1[,'time']<- seq(0, tail(time_burnin,1)+tail(time_NPI, 1), timestep)
  out.2[,'time']<- seq(0, tail(time_burnin,1)+tail(time_NPI, 1), timestep)
  
  # amount of infected (in SD model, _test and _trace = 0)
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.1, lwd = 1.5, col = cols[i], lty = 1)
  out.test<- out.2[which(out.run.NPI[,'time'] == I(startip+time_burnin)):nrow(out.2),]
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.test, lwd = 1.5, col = darkcols[i], lty = 2)
  
}

abline(v = time_burnin, lty = 2)
abline(v = time_burnin+startip, lty = 2)
#mtext(side = 2, 'Prevalence', line = 2.2, font = 2, cex = 1)
#mtext(side = 1, 'Days', line = 2.2, font = 2, cex = 1)
#axis(1, at = seq(0, time_NPI+time_burnin, 30))



# 6 ----
# Initial conditions
nb_cases_ini = 1000
prop_n = 0.7
ini<- c(Sn = prop_n - ((prop_n*nb_cases_ini)/N),
        Sc = (1 - prop_n) - (((1 - prop_n)*nb_cases_ini)/N),
        In = (prop_n*nb_cases_ini)/N,
        In_test = 0,
        Ic = ((1 - prop_n)*nb_cases_ini)/N,
        Ic_test = 0,
        Ic_trace = 0,
        Rn = 0,
        Rc = 0)


# Time parameters

time_burnin = 30
time_NPI = 400
startip = 90 # IP introduced seven days after implementing Social distancing
time_NPI.ip = time_NPI - startip
time_sim = 60
timestep = 0.1
times.burnin <- seq(0, time_burnin, by = timestep)
times.run.NPI <- seq(0, time_NPI, by = timestep)
times.run.NPI.ip <- seq(0, time_NPI.ip, by = timestep)
times.noIntervention <- seq(0, time_burnin+time_NPI, by = timestep)


plot('', ylim = c(0, 0.15), xlim = c(0, time_NPI+time_burnin), ylab = '', xlab = '', cex.axis = 1)
# plot epidemic in the absence of intervention
pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
out.Nointervention <- ode(y = ini,  times = times.noIntervention, func = sir_npi, parms = pars_null, npi = 'tti')  
#lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.Nointervention, lwd = 3, col = 'lightgrey', lty = 1)
x =  out.Nointervention[,'time']
y = rowSums(out.Nointervention[,c('In', 'In_test', 'Ic', 'Ic_test', 'Ic_trace')])
polygon(c(x[x>=0], max(x), 0), c(ifelse(y[x>=0] > 0.155, 0.155, y[x>=0]), 0, 0), col = '#F0F0F0', border = NA)


for(i in 1:length(Es)){
  
  # parameters for dynamic model, no immunity passport
  pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
  #pars_null.dyn<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0,   x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_tti.dyn<-  c(b = b, g = g, a = NA,  tau = focus.tau,  p_test = focus.ptest, x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_tti.dyn.ip<-  c(b = b, g = g, a = NA,  tau = focus.tau,  p_test = focus.ptest, x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 1)
  
  out.burnin <- ode(y = ini,  times = times.burnin, func = sir_npi, parms = pars_null, npi = 'tti')  
  out.run.NPI<- ode(y = tail(out.burnin,1)[,-1], times = times.run.NPI, func = sir_npi, parms = pars_tti.dyn, npi = 'tti')
  out.run.NPI.ip<- ode(y = out.run.NPI[which(out.run.NPI[,'time'] == startip),-1], times = times.run.NPI.ip, func = sir_npi, parms = pars_tti.dyn.ip, npi = 'tti')
  
  # Assemble output
  out.1<- rbind(out.burnin, out.run.NPI[-1,])
  out.2<- rbind(out.burnin, out.run.NPI[2:which(out.run.NPI[,'time'] == startip),], out.run.NPI.ip[-1,])
  out.1[,'time']<- seq(0, tail(time_burnin,1)+tail(time_NPI, 1), timestep)
  out.2[,'time']<- seq(0, tail(time_burnin,1)+tail(time_NPI, 1), timestep)
  
  # amount of infected (in SD model, _test and _trace = 0)
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.1, lwd = 1.5, col = cols[i], lty = 1)
  out.test<- out.2[which(out.run.NPI[,'time'] == I(startip+time_burnin)):nrow(out.2),]
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.test, lwd = 1.5, col = darkcols[i], lty = 2)
  
}

abline(v = time_burnin, lty = 2)
abline(v = time_burnin+startip, lty = 2)
#mtext(side = 2, 'Prevalence', line = 2.2, font = 2, cex = 1)
#mtext(side = 1, 'Days', line = 2.2, font = 2, cex = 1)
#axis(1, at = seq(0, time_NPI+time_burnin, 30))


mtext(side = 4, 'Ini = 1000', font = 2, cex = 1, line = 1)



# 7 ----
# Initial conditions
nb_cases_ini = 10000
prop_n = 0.7
ini<- c(Sn = prop_n - ((prop_n*nb_cases_ini)/N),
        Sc = (1 - prop_n) - (((1 - prop_n)*nb_cases_ini)/N),
        In = (prop_n*nb_cases_ini)/N,
        In_test = 0,
        Ic = ((1 - prop_n)*nb_cases_ini)/N,
        Ic_test = 0,
        Ic_trace = 0,
        Rn = 0,
        Rc = 0)


# Time parameters

time_burnin = 30
time_NPI = 400
startip = 30 # IP introduced seven days after implementing Social distancing
time_NPI.ip = time_NPI - startip
time_sim = 60
timestep = 0.1
times.burnin <- seq(0, time_burnin, by = timestep)
times.run.NPI <- seq(0, time_NPI, by = timestep)
times.run.NPI.ip <- seq(0, time_NPI.ip, by = timestep)
times.noIntervention <- seq(0, time_burnin+time_NPI, by = timestep)


plot('', ylim = c(0, 0.15), xlim = c(0, time_NPI+time_burnin), ylab = '', xlab = '', cex.axis = 1)
# plot epidemic in the absence of intervention
pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
out.Nointervention <- ode(y = ini,  times = times.noIntervention, func = sir_npi, parms = pars_null, npi = 'tti')  
#lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.Nointervention, lwd = 3, col = 'lightgrey', lty = 1)
x =  out.Nointervention[,'time']
y = rowSums(out.Nointervention[,c('In', 'In_test', 'Ic', 'Ic_test', 'Ic_trace')])
polygon(c(x[x>=0], max(x), 0), c(ifelse(y[x>=0] > 0.155, 0.155, y[x>=0]), 0, 0), col = '#F0F0F0', border = NA)


for(i in 1:length(Es)){
  
  # parameters for dynamic model, no immunity passport
  pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
  #pars_null.dyn<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0,   x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_tti.dyn<-  c(b = b, g = g, a = NA,  tau = focus.tau,  p_test = focus.ptest, x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_tti.dyn.ip<-  c(b = b, g = g, a = NA,  tau = focus.tau,  p_test = focus.ptest, x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 1)
  
  out.burnin <- ode(y = ini,  times = times.burnin, func = sir_npi, parms = pars_null, npi = 'tti')  
  out.run.NPI<- ode(y = tail(out.burnin,1)[,-1], times = times.run.NPI, func = sir_npi, parms = pars_tti.dyn, npi = 'tti')
  out.run.NPI.ip<- ode(y = out.run.NPI[which(out.run.NPI[,'time'] == startip),-1], times = times.run.NPI.ip, func = sir_npi, parms = pars_tti.dyn.ip, npi = 'tti')
  
  # Assemble output
  out.1<- rbind(out.burnin, out.run.NPI[-1,])
  out.2<- rbind(out.burnin, out.run.NPI[2:which(out.run.NPI[,'time'] == startip),], out.run.NPI.ip[-1,])
  out.1[,'time']<- seq(0, tail(time_burnin,1)+tail(time_NPI, 1), timestep)
  out.2[,'time']<- seq(0, tail(time_burnin,1)+tail(time_NPI, 1), timestep)
  
  # amount of infected (in SD model, _test and _trace = 0)
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.1, lwd = 1.5, col = cols[i], lty = 1)
  out.test<- out.2[which(out.run.NPI[,'time'] == I(startip+time_burnin)):nrow(out.2),]
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.test, lwd = 1.5, col = darkcols[i], lty = 2)
  
}

abline(v = time_burnin, lty = 2)
abline(v = time_burnin+startip, lty = 2)
#mtext(side = 2, 'Prevalence', line = 2.2, font = 2, cex = 1)
#mtext(side = 1, 'Days', line = 2.2, font = 2, cex = 1)
#axis(1, at = seq(0, time_NPI+time_burnin, 30))





# 8 ----
# Initial conditions
nb_cases_ini = 10000
prop_n = 0.7
ini<- c(Sn = prop_n - ((prop_n*nb_cases_ini)/N),
        Sc = (1 - prop_n) - (((1 - prop_n)*nb_cases_ini)/N),
        In = (prop_n*nb_cases_ini)/N,
        In_test = 0,
        Ic = ((1 - prop_n)*nb_cases_ini)/N,
        Ic_test = 0,
        Ic_trace = 0,
        Rn = 0,
        Rc = 0)


# Time parameters

time_burnin = 30
time_NPI = 400
startip = 60 # IP introduced seven days after implementing Social distancing
time_NPI.ip = time_NPI - startip
time_sim = 60
timestep = 0.1
times.burnin <- seq(0, time_burnin, by = timestep)
times.run.NPI <- seq(0, time_NPI, by = timestep)
times.run.NPI.ip <- seq(0, time_NPI.ip, by = timestep)
times.noIntervention <- seq(0, time_burnin+time_NPI, by = timestep)


plot('', ylim = c(0, 0.15), xlim = c(0, time_NPI+time_burnin), ylab = '', xlab = '', cex.axis = 1)
# plot epidemic in the absence of intervention
pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
out.Nointervention <- ode(y = ini,  times = times.noIntervention, func = sir_npi, parms = pars_null, npi = 'tti')  
#lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.Nointervention, lwd = 3, col = 'lightgrey', lty = 1)
x =  out.Nointervention[,'time']
y = rowSums(out.Nointervention[,c('In', 'In_test', 'Ic', 'Ic_test', 'Ic_trace')])
polygon(c(x[x>=0], max(x), 0), c(ifelse(y[x>=0] > 0.155, 0.155, y[x>=0]), 0, 0), col = '#F0F0F0', border = NA)


for(i in 1:length(Es)){
  
  # parameters for dynamic model, no immunity passport
  pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
  #pars_null.dyn<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0,   x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_tti.dyn<-  c(b = b, g = g, a = NA,  tau = focus.tau,  p_test = focus.ptest, x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_tti.dyn.ip<-  c(b = b, g = g, a = NA,  tau = focus.tau,  p_test = focus.ptest, x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 1)
  
  out.burnin <- ode(y = ini,  times = times.burnin, func = sir_npi, parms = pars_null, npi = 'tti')  
  out.run.NPI<- ode(y = tail(out.burnin,1)[,-1], times = times.run.NPI, func = sir_npi, parms = pars_tti.dyn, npi = 'tti')
  out.run.NPI.ip<- ode(y = out.run.NPI[which(out.run.NPI[,'time'] == startip),-1], times = times.run.NPI.ip, func = sir_npi, parms = pars_tti.dyn.ip, npi = 'tti')
  
  # Assemble output
  out.1<- rbind(out.burnin, out.run.NPI[-1,])
  out.2<- rbind(out.burnin, out.run.NPI[2:which(out.run.NPI[,'time'] == startip),], out.run.NPI.ip[-1,])
  out.1[,'time']<- seq(0, tail(time_burnin,1)+tail(time_NPI, 1), timestep)
  out.2[,'time']<- seq(0, tail(time_burnin,1)+tail(time_NPI, 1), timestep)
  
  # amount of infected (in SD model, _test and _trace = 0)
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.1, lwd = 1.5, col = cols[i], lty = 1)
  out.test<- out.2[which(out.run.NPI[,'time'] == I(startip+time_burnin)):nrow(out.2),]
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.test, lwd = 1.5, col = darkcols[i], lty = 2)
  
}

abline(v = time_burnin, lty = 2)
abline(v = time_burnin+startip, lty = 2)
#mtext(side = 2, 'Prevalence', line = 2.2, font = 2, cex = 1)
#mtext(side = 1, 'Days', line = 2, font = 2, cex = 1)
#axis(1, at = seq(0, time_NPI+time_burnin, 30))



# 9 ----
# Initial conditions
nb_cases_ini = 10000
prop_n = 0.7
ini<- c(Sn = prop_n - ((prop_n*nb_cases_ini)/N),
        Sc = (1 - prop_n) - (((1 - prop_n)*nb_cases_ini)/N),
        In = (prop_n*nb_cases_ini)/N,
        In_test = 0,
        Ic = ((1 - prop_n)*nb_cases_ini)/N,
        Ic_test = 0,
        Ic_trace = 0,
        Rn = 0,
        Rc = 0)


# Time parameters

time_burnin = 30
time_NPI = 400
startip = 90 # IP introduced seven days after implementing Social distancing
time_NPI.ip = time_NPI - startip
time_sim = 60
timestep = 0.1
times.burnin <- seq(0, time_burnin, by = timestep)
times.run.NPI <- seq(0, time_NPI, by = timestep)
times.run.NPI.ip <- seq(0, time_NPI.ip, by = timestep)
times.noIntervention <- seq(0, time_burnin+time_NPI, by = timestep)


plot('', ylim = c(0, 0.15), xlim = c(0, time_NPI+time_burnin), ylab = '', xlab = '', cex.axis = 1)
# plot epidemic in the absence of intervention
pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
out.Nointervention <- ode(y = ini,  times = times.noIntervention, func = sir_npi, parms = pars_null, npi = 'tti')  
#lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.Nointervention, lwd = 3, col = 'lightgrey', lty = 1)
x =  out.Nointervention[,'time']
y = rowSums(out.Nointervention[,c('In', 'In_test', 'Ic', 'Ic_test', 'Ic_trace')])
polygon(c(x[x>=0], max(x), 0), c(ifelse(y[x>=0] > 0.155, 0.155, y[x>=0]), 0, 0), col = '#F0F0F0', border = NA)


for(i in 1:length(Es)){
  
  # parameters for dynamic model, no immunity passport
  pars_null<- c(b = b, g = g, a = 1, tau = 1000, p_test = 0, x = 0)
  #pars_null.dyn<- c(b = b, g = g, a = 1,   tau = 1000, p_test = 0,   x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_tti.dyn<-  c(b = b, g = g, a = NA,  tau = focus.tau,  p_test = focus.ptest, x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 0)
  pars_tti.dyn.ip<-  c(b = b, g = g, a = NA,  tau = focus.tau,  p_test = focus.ptest, x = focus.x, C = focus.C, E = Es[i], time_sim = time_sim, ip = 1)
  
  out.burnin <- ode(y = ini,  times = times.burnin, func = sir_npi, parms = pars_null, npi = 'tti')  
  out.run.NPI<- ode(y = tail(out.burnin,1)[,-1], times = times.run.NPI, func = sir_npi, parms = pars_tti.dyn, npi = 'tti')
  out.run.NPI.ip<- ode(y = out.run.NPI[which(out.run.NPI[,'time'] == startip),-1], times = times.run.NPI.ip, func = sir_npi, parms = pars_tti.dyn.ip, npi = 'tti')
  
  # Assemble output
  out.1<- rbind(out.burnin, out.run.NPI[-1,])
  out.2<- rbind(out.burnin, out.run.NPI[2:which(out.run.NPI[,'time'] == startip),], out.run.NPI.ip[-1,])
  out.1[,'time']<- seq(0, tail(time_burnin,1)+tail(time_NPI, 1), timestep)
  out.2[,'time']<- seq(0, tail(time_burnin,1)+tail(time_NPI, 1), timestep)
  
  # amount of infected (in SD model, _test and _trace = 0)
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.1, lwd = 1.5, col = cols[i], lty = 1)
  out.test<- out.2[which(out.run.NPI[,'time'] == I(startip+time_burnin)):nrow(out.2),]
  lines((In + In_test + Ic + Ic_test + Ic_trace) ~ time, data = out.test, lwd = 1.5, col = darkcols[i], lty = 2)
  
}

abline(v = time_burnin, lty = 2)
abline(v = time_burnin+startip, lty = 2)
#mtext(side = 2, 'Prevalence', line = 2.2, font = 2, cex = 1)
#mtext(side = 1, 'Days', line = 2.2, font = 2, cex = 1)
#axis(1, at = seq(0, time_NPI+time_burnin, 30))


mtext(side = 4, 'Ini = 10000', font = 2, cex = 1, line = 1)


dev.off()
