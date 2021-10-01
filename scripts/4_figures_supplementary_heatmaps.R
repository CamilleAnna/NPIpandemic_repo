# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#   Social dilemmas in non-pharmaceutical interventions to tackle epidemics
#                     Simonet C, Sweeny A, McNally L (2021)
#                     Codes for supplementary figures output
#
#                       last edited: July, 27th, 2021
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Script summary
# Runs simulation of each intervention for a range of C and intervention efficacy parameters values
# efficacy in SD varies with alpha, efficacy in TTI varies with tau
#   (notes:)
#   run same model as main model script, but simplified code for faster run (two distinct functions for SD and TTI rather than single big code wrapper)
#   again, payoff with/without immunity passport evaluate afterwards so does not change anything within the model functions

# Then load output and plot heatmaps for all combination of C and efficacy parameter values


setwd('path/to/repo') # setwd('~/Documents/PhD/Research/Covid19/NPIpandemic_repo/')
library(readxl); library(ggplot2); library(dplyr); library(tidyr);library(Rmisc); library(lubridate);library(plotly);library(yarrr); library(knitr);library(plyr);library(utils);library(tidyverse);library(gridExtra);library(deSolve); library(stringi); library(stringr);
library(ggthemes); library(deSolve)


# [1] Model functions and parameters ----
# same a in main script model but simplified code for quicker run
sir_npi_quick_sd <- function(time, state, parameters, npi) {
  
  with(as.list(c(state, parameters)), {
    
    
    dSn = -b*(In + (1/a)*Ic)*Sn
    dSc = -b*(In + (1/a)*Ic)*(1/a)*Sc
    
    dIn = b*(In + (1/a)*Ic)*Sn - g * In
    dIc = b*(In + (1/a)*Ic)*(1/a)*Sc - g * Ic
    
    dRn = g * In
    dRc = g * Ic
    
    return(list(c(dSn, dSc, dIn, dIc, dRn, dRc)))})
  
}
sir_npi_quick_tti <- function(time, state, parameters, npi) {
  
  with(as.list(c(state, parameters)), {
    
    qtest  = (g*tau)/(1 + (g*tau))
    qtrace = (g*tau)/(2 + (g*tau))
    
    
    # ODE, with replicator dynamic for the S "activated" if x > 1              
    dSn = -b*(In + Ic  + qtest*Ic_test + qtest*In_test + qtrace*Ic_trace)*Sn
    dSc = -b*(In + Ic  + qtest*Ic_test + qtest*In_test + qtrace*Ic_trace)*Sc
    
    dIn = (1-p_test)*(b*(In + Ic  + qtest*Ic_test + qtest*In_test + qtrace*Ic_trace)*Sn)- g * In
    dIn_test = (p_test)*(b*(In + Ic  + qtest*Ic_test + qtest*In_test + qtrace*Ic_trace)*Sn)- g * In_test
    
    dIc = (1-p_test)*(b*(In + Ic  + qtest*In_test)*Sc)- g * Ic
    dIc_test = (p_test)*(b*(In + Ic  + qtest*In_test)*Sc)- g * Ic_test
    dIc_trace = b*(qtest*Ic_test + qtrace*Ic_trace)*Sc - g*Ic_trace
    
    dRn  = g * (In + In_test)
    dRc  = g * (Ic + Ic_test + Ic_trace)
    
    
    return(list(c(dSn, dSc, dIn, dIn_test, dIc, dIc_test, dIc_trace, dRn, dRc)))})
  
  
}
integral<- function(x, strains){ 
  out<- as.data.frame(x)
  out_sum<- data.frame(time = out$time, sum = rowSums(as.data.frame(out[,strains])))
  sum(diff(out_sum$time) * (head(out_sum$sum,-1)+tail(out_sum$sum,-1)))/2
}


# define general

# Epidemic parameters
r0 = 3   # (reproduction number)
g = 1/14 # (recovery rate)
b = r0*g # (transmission rate estimated from R0)
N = 67886004 # used to define initial condition
nb_cases_ini = 60000 # number of cases on March 16th when PM first announce to limit contacts = 58286

# Simulation parameters
time_sim = 800
timestep = 0.1
times.run <- seq(0, time_sim, by = timestep)


# [2.1] Run - SOCIAL DISTANCING ----
# Simulation dataframe to fill
reductions<- c(0.05, 0.1, 0.25,0.5, 0.75, 0.9, 0.95)
alphas<- 1/(1-reductions) # corresponds to reduction of beta of 5%, 10%, 25%, 50%, 75%, 90%, 95%
df.sd.list<- vector('list', length = length(alphas))
names(df.sd.list)<- paste0('alpha_', round(alphas, 2))


# Simulation
for(j in 1:length(alphas)){
  
  print(j)
  
  # redefine a new parameter set with a different alpha and re-initialise an entry of the df.sd.list()
  pars_sd.focal<- c(b = b, g = g, a = alphas[j], tau = NA,   p_test = NA, x = 0)
  
  
  min = 0.001 # was 0.001, test with 0.1
  df.sd.list[[j]]<- data.frame(Sn = seq(0+min, 1-min, min),
                             p_not_infected_c = NA,
                             p_not_infected_n = NA,
                             prop_time_not_recovered_c = NA,
                             prop_time_not_recovered_n = NA)
  
  
  # run through the df.sd.list[[j]] for different initial compliance
  for(i in 1:nrow(df.sd.list[[j]])){
    
    #print(i)
    
    # redefine ini at each loop iteration, with different frequencies of compliers
    ini<- c(Sn = df.sd.list[[j]]$Sn[i] - ((df.sd.list[[j]]$Sn[i]*nb_cases_ini)/N),
            Sc = (1 - df.sd.list[[j]]$Sn[i]) - (((1 - df.sd.list[[j]]$Sn[i])*nb_cases_ini)/N),
            In = (df.sd.list[[j]]$Sn[i]*nb_cases_ini)/N,
            Ic = ((1 - df.sd.list[[j]]$Sn[i])*nb_cases_ini)/N,
            Rn = 0,
            Rc = 0)
    
    # Run
    out<- ode(y = ini, times = times.run, func = sir_npi_quick_sd, parms = pars_sd.focal)

    # EVALUATE BENEFIT AND COST AT EQUILIBRIUM (tmax), to compute PAYOFFS AFTERWARD
    # proba that an individual is still susceptible by the end of the epidemic
    df.sd.list[[j]]$p_not_infected_c[i]<- tail(out, 1)[,c('Sc')]/sum(tail(out, 1)[,c('Sc', 'Ic', 'Rc')]) 
    df.sd.list[[j]]$p_not_infected_n[i]<- tail(out, 1)[,c('Sn')]/sum(tail(out, 1)[,c('Sn', 'In', 'Rn')])
    
    # Using integral, average proportion of time individual spent not recovered
    df.sd.list[[j]]$prop_time_not_recovered_c[i]<- integral(out, c('Sc', 'Ic'))/integral(out, c('Sc', 'Ic', 'Rc')) 
    df.sd.list[[j]]$prop_time_not_recovered_n[i]<- integral(out, c('Sn', 'In'))/integral(out, c('Sn', 'In', 'Rn'))
    
  }
  
  
}

#save.image('./output/model_runs/SD_various_alphas_nbIni_60000_timeSim_800.RData')


# [2.2] Run - TTI ----

# Simulation dataframe to fill
taus<- c(1, 2, 5, 7, 8, 10, 12)
df.tti.list<- vector('list', length = length(taus))
names(df.tti.list)<- paste0('tau_', taus)


# Simulation
for(j in 1:length(taus)){
  
  print(j)
  
  # redefine a new parameter set with a different alpha and re-initialise an entry of the df.tti.list()
  
  pars_tti.focal<- c(b = b, g = g, a = NA, tau = taus[j], p_test = 0.7, x = 0)
  
  
  min = 0.001 # was 0.001, test with 0.1
  df.tti.list[[j]]<- data.frame(Sn = seq(0+min, 1-min, min),
                               p_not_infected_c = NA,
                               p_not_infected_n = NA,
                               prop_time_not_recovered_c = NA,
                               prop_time_not_recovered_n = NA)
  
  
  # run through the df.tti.list[[j]] for different initial compliance
  for(i in 1:nrow(df.tti.list[[j]])){
    
    #print(i)
    
    # redefine ini at each loop iteration, with different frequencies of compliers
    ini<- c(Sn = df.tti.list[[j]]$Sn[i] - ((df.tti.list[[j]]$Sn[i]*nb_cases_ini)/N),
            Sc = (1 - df.tti.list[[j]]$Sn[i]) - (((1 - df.tti.list[[j]]$Sn[i])*nb_cases_ini)/N),
            In = (df.tti.list[[j]]$Sn[i]*nb_cases_ini)/N,
            In_test = 0,
            Ic = ((1 - df.tti.list[[j]]$Sn[i])*nb_cases_ini)/N,
            Ic_test = 0,
            Ic_trace = 0,
            Rn = 0,
            Rc = 0)
    
    # Run (no burnin phase)
    out<- ode(y = ini, times = times.run, func = sir_npi_quick_tti, parms = pars_tti.focal)
    
    # EVALUATE BENEFIT AND COST AT EQUILIBRIUM (tmax), to compute PAYOFFS AFTERWARD
    # proba that an individual is still susceptible by the end of the epidemic
    df.tti.list[[j]]$p_not_infected_c[i]<- tail(out, 1)[,c('Sc')]/sum(tail(out, 1)[,c('Sc', 'Ic', 'Rc')]) 
    df.tti.list[[j]]$p_not_infected_n[i]<- tail(out, 1)[,c('Sn')]/sum(tail(out, 1)[,c('Sn', 'In', 'Rn')])
    
    # Using integral, average proportion of time individual spent not recovered
    df.tti.list[[j]]$prop_time_not_recovered_c[i]<- integral(out, c('Sc', 'Ic'))/integral(out, c('Sc', 'Ic', 'Rc')) 
    df.tti.list[[j]]$prop_time_not_recovered_n[i]<- integral(out, c('Sn', 'In'))/integral(out, c('Sn', 'In', 'Rn'))
    
  }
  
  
}

#save.image('./output/model_runs/TTI_various_tau_ptest_0.7_nbIni_60000_timeSim_800.RData')



# [3] PLOTTING ----

# same function wrappers as in main figures script, just with slightly different specs for size & layout to fit
plot.heatmap<- function(payoff.df, focal.payoff.comp, zlim.fixed = FALSE, zlims = NA, draw.y.axis){
  
  # Plot the heatmap
  payoffs.comp.df<- payoff.df[,c('Sc', 'E', focal.payoff.comp)]
  colnames(payoffs.comp.df)<- c('Sc', 'E', 'focal.payoff.comp')
  
  payoff.matrix<- payoffs.comp.df %>%
    spread(key = Sc, value = focal.payoff.comp)
  row.names(payoff.matrix)<- payoff.matrix$E
  payoff.matrix<- payoff.matrix %>% select(-E) %>% as.matrix()
  
  #cols<-rev(colorRampPalette(c("navy", 'blue', "white", "red", "darkred"))(300))
  cols<-rev(colorRampPalette(c("blue", "white", "red"))(300))
  rotate <- function(x) t(apply(x, 2, rev))
  
  
  if(zlim.fixed == FALSE){
    
    image(payoff.matrix,col=cols, useRaster=T,
          zlim=max(abs(payoff.matrix)) * c(-1, 1),
          #zlim=c(-2,2),
          xaxt = 'n', yaxt = 'n')
  }else{
    
    image(payoff.matrix,col=cols, useRaster=T,
          #zlim=max(abs(payoff.matrix)) * c(-1, 1),
          zlim=c(zlims[1],zlims[2]),
          xaxt = 'n', yaxt = 'n')
  }
  
  
  # Add the axes, with re-scaling appropriate to match 0-1 range to whatever E range
 
  if(draw.y.axis == TRUE){
  axis(2, at = c(0, 0.25, 0.5, 0.75, 1), font = 1, lwd = 1,  cex.axis = 0.7, tck = -0.05, padj = 0.95, col = 'black')
  }
  axis(1, at = c(0, 0.25, 0.5, 0.75, 1), labels = quantile(unique(payoffs.comp.df$E)), cex.axis = 0.7, font = 1, lwd = 1, tck = -0.05, padj = -2, col = 'black')
  axis(3, at = c(0, 1), labels = c('', ''), lwd = 1, lwd.ticks = 0)
  axis(4, at = c(0, 1), labels = c('', ''), lwd = 1, lwd.ticks = 0)
  
  # Add the 0 cline line
  # find all points where there is a transition from negative to positive values of payoff difference
  # select payoff we want
  
  foo<- payoffs.comp.df
  
  # find negative to positive transition points
  foo.neg<- foo %>%
    filter(focal.payoff.comp < 0) %>%
    group_by(Sc) %>%
    filter(focal.payoff.comp == max(focal.payoff.comp)) %>%
    ungroup()
  
  foo.pos<- foo %>%
    filter(focal.payoff.comp > 0) %>%
    group_by (Sc) %>%
    filter(focal.payoff.comp == min(focal.payoff.comp)) %>%
    ungroup()
  
  
  d<- data.frame(plot.scale = quantile(seq(0,1, length.out = 10000), p = seq(0,1,0.0001))) %>%
    mutate(my.scale = quantile(seq(min(payoffs.comp.df$E), max(payoffs.comp.df$E), length.out = 10000), p = seq(0,1,0.0001)))
  
  # compute isocline value and rescale it on a 0-1 scale
  foo2<- inner_join(foo.neg, foo.pos, by = 'Sc') %>%
    mutate(cline = E.x, #(E.x + E.y)/2,
           cline.scaled = NA) %>%
    arrange(Sc)
  
  for(i in 1:nrow(foo2)){
    
    foo2$cline.scaled[i] = d[which.min(abs(d$my.scale - foo2$cline[i])),'plot.scale']
    
  } 
  
  foo2<- foo2 %>%
    mutate(keep = seq(1,nrow(.)) %% 40) %>%
    filter(keep == 1)
  
  # Draw the line
  lines(Sc ~ cline.scaled, data = foo2, lwd = 1, lty = 2)
  
}
legend<- function(payoff.df.focal, focal.payoff.comp){
  
  cols<-rev(colorRampPalette(c("navy","blue", "white", "red", "darkred"))(300))
  clbr<-matrix(ncol=length(cols),nrow=2)
  clbr[1,]<-seq(0,1,length.out=length(cols))
  clbr[2,]<-seq(0,1,length.out=length(cols))
  
  
  plot(NA,NA,xlim=c(0,1),ylim=c(0,1),yaxt="n", xaxt="n", xlab="",ylab="",cex.lab=1.5,bty="n")
  
  image(z=clbr,y=seq(0, 1,length.out=length(cols)),
        x=c(0,0.17),
        col=cols,
        zlim=c(0,1),
        yaxt="n",xlab="",ylab="",useRaster=T,cex.lab=1.5,add=T)
  
  axis(side = 4, pos = 0.3, at=seq(0, 1, by=0.25), srt = 45, labels = FALSE,tck = -.09, lwd = 1)
  
  
  text(y = c(0, 0.25, 0.75, 1),
       x = 0.5,
       labels = substr(formatC(
         round(max(abs(payoff.df.focal[,which(colnames(payoff.df.focal) == focal.payoff.comp)])) * c(-1, -0.5, 0.5, 1), 2), digit = 2, format = 's'),
         1, 4),
       font = 1, cex = 0.4, srt = 0, adj = 0)
  
  text(y = 0.5,
       x = 0.5,
       labels = '0',
       font = 1, cex = 0.4, srt = 0, adj = 0)
  
}
prep_df<- function(model.df.output, test.C, Es.range){
  
es<- data.frame(Sn = model.df.output$Sn,
                E = rep(seq(Es.range[1], Es.range[2], 0.001), each = nrow(model.df.output)))

cs<- data.frame(Sn = model.df.output$Sn,
                C = rep(c(test.C), each = nrow(model.df.output)))

cb<- left_join(cs, es, by = 'Sn')

payoff.df<- model.df.output %>%
  left_join(cb, by = 'Sn') %>%
  arrange(E) %>%
  mutate(pi_c = p_not_infected_c - C,
         pi_n = p_not_infected_n - E,
         pi_c.ip = p_not_infected_c - C * prop_time_not_recovered_c,
         pi_n.ip = p_not_infected_n - E * prop_time_not_recovered_n,
         Sc = 1 - Sn) %>%
  mutate(payoff.comp = pi_c - pi_n,
         payoff.comp.ip = pi_c.ip - pi_n.ip)


payoff.df.focal.C<- payoff.df[payoff.df$C == test.C,]

return(payoff.df.focal.C)

}

# define figure layout
m  <- rbind(
  c(0, 0.2, 0.75, 0.95),     # 1
  c(0.205, 0.24, 0.85, 0.95), # 2
  c(0.25, 0.45, 0.75, 0.95), # 3
  c(0.455, 0.49,0.85, 0.95),  # 4
  c(0.5, 0.7, 0.75, 0.95),   # 5
  c(0.705, 0.74, 0.85, 0.95), # 6
  c(0.75, 0.95, 0.75, 0.95), # 7
  c(0.955, 0.99, 0.85, 0.95), # 8
  
  c(0, 0.2, 0.5, 0.7),     # 1
  c(0.205, 0.24, 0.6, 0.7), # 2
  c(0.25, 0.45, 0.5, 0.7), # 3
  c(0.455, 0.49,0.6, 0.7),  # 4
  c(0.5, 0.7, 0.5, 0.7),   # 5
  c(0.705, 0.74, 0.6, 0.7), # 6
  c(0.75, 0.95, 0.5, 0.7), # 7
  c(0.955, 0.99, 0.6, 0.7), # 8
  
  c(0, 0.2, 0.25, 0.45),     # 1
  c(0.205, 0.24, 0.35, 0.45), # 2
  c(0.25, 0.45, 0.25, 0.45), # 3
  c(0.455, 0.49,0.35, 0.45),  # 4
  c(0.5, 0.7, 0.25, 0.45),   # 5
  c(0.705, 0.74, 0.35, 0.45), # 6
  c(0.75, 0.95, 0.25, 0.45), # 7
  c(0.955, 0.99, 0.35, 0.45), # 8
  
  c(0, 0.2, 0, 0.2),     # 1
  c(0.205, 0.24, 0.1, 0.2), # 2
  c(0.25, 0.45, 0, 0.2), # 3
  c(0.455, 0.49,0.1, 0.2),  # 4
  c(0.5, 0.7, 0, 0.2),   # 5
  c(0.705, 0.74, 0.1, 0.2), # 6
  c(0.75, 0.95, 0, 0.2), # 7
  c(0.955, 0.99, 0.1, 0.2)) # 8


# Plot output SD ----

# load the model output from above run
load('./output/model_runs/SD_various_alphas_nbIni_60000_timeSim_800.RData')

close.screen(all.screens = TRUE)
fty = 2
cex.labels = 1


# plots a 4x4 heatmap grid for increase values of cost of compliance (C) and increasing intervention efficacy (alpha)

pdf('./output/figures/FigureS1_SD_heatmap_vary_alpha.pdf', width = (13+1+1)/2.54, height = (13+1+1)/2.54, pointsize = 7)
close.screen(all.screens = TRUE)
split.screen(m)
par(omi = rep(1/2.54,4))
for(i in 1:nrow(m)){
  screen(i)
  par(mar = c(0, 0, 0, 0))
  plot(1, axes = FALSE, type = "n")
  #box()
}

#dev.off()

# ALPHA = 1.33 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LEGENDS, screens 2,4,6,8
screen(2); legend(prep_df(df.sd.list$alpha_1.33, test.C = 0.5, Es.range = c(0,1)), 'payoff.comp')
screen(4); legend(prep_df(df.sd.list$alpha_1.33, test.C = 1,   Es.range = c(0,2)), 'payoff.comp')
mtext('SOCIAL DISTANCING', side = 3, line = 3, font = 2, padj = 0, col = 'black')
screen(6); legend(prep_df(df.sd.list$alpha_1.33, test.C = 1.5, Es.range = c(0,2)), 'payoff.comp')
screen(8); legend(prep_df(df.sd.list$alpha_1.33, test.C = 2,   Es.range = c(0,3)), 'payoff.comp')

# HEATMAPS, screens 1,3,5,7

screen(1)
plot.heatmap(prep_df(df.sd.list$alpha_1.33, test.C = 0.5, Es.range = c(0,1)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = TRUE)

mtext('C = 0.5', side = 3, line = 1, font = 2)

screen(3)
plot.heatmap(prep_df(df.sd.list$alpha_1.33, test.C = 1, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)


mtext('C = 1', side = 3, line = 1, font = 2)


screen(5)
plot.heatmap(prep_df(df.sd.list$alpha_1.33, test.C = 1.5, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)

mtext('C = 1.5', side = 3, line = 1, font = 2)

screen(7)
plot.heatmap(prep_df(df.sd.list$alpha_1.33, test.C = 2, Es.range = c(0,3)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose), 
             draw.y.axis = FALSE)

mtext('C = 2', side = 3, line = 1, font = 2)

mtext('a = 1.33', side = 4, line = 2.5, font = fty, cex = cex.labels)


# ALPHA = 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LEGENDS, screens 10,12,14,16
screen(10); legend(prep_df(df.sd.list$alpha_2, test.C = 0.5, Es.range = c(0,1)), 'payoff.comp')
screen(12); legend(prep_df(df.sd.list$alpha_2, test.C = 1,   Es.range = c(0,2)), 'payoff.comp')
screen(14); legend(prep_df(df.sd.list$alpha_2, test.C = 1.5, Es.range = c(0,2)), 'payoff.comp')
screen(16); legend(prep_df(df.sd.list$alpha_2, test.C = 2,   Es.range = c(0,3)), 'payoff.comp')


# HEATMAPS, screens 9,11,13,15

screen(9)
plot.heatmap(prep_df(df.sd.list$alpha_2, test.C = 0.5, Es.range = c(0,1)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = TRUE)

screen(11)
plot.heatmap(prep_df(df.sd.list$alpha_2, test.C = 1, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)

screen(13)
plot.heatmap(prep_df(df.sd.list$alpha_2, test.C = 1.5, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)

screen(15)
plot.heatmap(prep_df(df.sd.list$alpha_2, test.C = 2, Es.range = c(0,3)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose), 
             draw.y.axis = FALSE)

mtext('a = 2', side = 4, line = 2.5, font = fty, cex = cex.labels)


# ALPHA = 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LEGENDS, screens 18,20,22,24
screen(18); legend(prep_df(df.sd.list$alpha_4, test.C = 0.5, Es.range = c(0,1)), 'payoff.comp')
screen(20); legend(prep_df(df.sd.list$alpha_4, test.C = 1,   Es.range = c(0,2)), 'payoff.comp')
screen(22); legend(prep_df(df.sd.list$alpha_4, test.C = 1.5, Es.range = c(0,2)), 'payoff.comp')
screen(24); legend(prep_df(df.sd.list$alpha_4, test.C = 2,   Es.range = c(0,3)), 'payoff.comp')

# HEATMAPS, screens 17,19,21,23

screen(17)
plot.heatmap(prep_df(df.sd.list$alpha_4, test.C = 0.5, Es.range = c(0,1)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = TRUE)

mtext('Initial compliance frequency', side = 2, line = 2.2, font = fty, cex = cex.labels, adj = -1)

screen(19)
plot.heatmap(prep_df(df.sd.list$alpha_4, test.C = 1, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)

screen(21)
plot.heatmap(prep_df(df.sd.list$alpha_4, test.C = 1.5, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)

screen(23)
plot.heatmap(prep_df(df.sd.list$alpha_4, test.C = 2, Es.range = c(0,3)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose), 
             draw.y.axis = FALSE)


mtext('a = 4', side = 4, line = 2.5, font = fty, cex = cex.labels)


# ALPHA = 10 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LEGENDS, screens 26, 28, 30, 32
screen(26); legend(prep_df(df.sd.list$alpha_10, test.C = 0.5, Es.range = c(0,1)), 'payoff.comp')
screen(28); legend(prep_df(df.sd.list$alpha_10, test.C = 1,   Es.range = c(0,2)), 'payoff.comp')
mtext('Enforcement level (E)', side = 1, line = 6, font = fty, cex = cex.labels, padj = 0)
screen(30); legend(prep_df(df.sd.list$alpha_10, test.C = 1.5, Es.range = c(0,2)), 'payoff.comp')
screen(32); legend(prep_df(df.sd.list$alpha_10, test.C = 2,   Es.range = c(0,3)), 'payoff.comp')

# HEATMAPS, screens 25, 27, 29, 31

screen(25)
plot.heatmap(prep_df(df.sd.list$alpha_10, test.C = 0.5, Es.range = c(0,1)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = TRUE)

screen(27)
plot.heatmap(prep_df(df.sd.list$alpha_10, test.C = 1, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)


screen(29)
plot.heatmap(prep_df(df.sd.list$alpha_10, test.C = 1.5, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)

screen(31)
plot.heatmap(prep_df(df.sd.list$alpha_10, test.C = 2, Es.range = c(0,3)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose), 
             draw.y.axis = FALSE)


mtext('a = 10', side = 4, line = 2.5, font = fty, cex = cex.labels)


dev.off()


# Plot output SD+IP ----

load('./output/model_runs/SD_various_alphas_nbIni_60000_timeSim_800.RData')

close.screen(all.screens = TRUE)
fty = 2
cex.labels = 1


pdf('./output/figures/FigureS2_SD_IP_heatmap_vary_alpha.pdf', width = (13+1+1)/2.54, height = (13+1+1)/2.54, pointsize = 7)
close.screen(all.screens = TRUE)
split.screen(m)
par(omi = rep(1/2.54,4))
for(i in 1:nrow(m)){
  screen(i)
  par(mar = c(0, 0, 0, 0))
  plot(1, axes = FALSE, type = "n")
  #box()
}

#dev.off()

# ALPHA = 1.33 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LEGENDS, screens 2,4,6,8
screen(2); legend(prep_df(df.sd.list$alpha_1.33, test.C = 0.5, Es.range = c(0,1)), 'payoff.comp.ip')
screen(4); legend(prep_df(df.sd.list$alpha_1.33, test.C = 1,   Es.range = c(0,2)), 'payoff.comp.ip')
mtext('SOCIAL DISTANCING + IMMUNITY PASSPORT', side = 3, line = 3, font = 2, padj = 0, col = 'black')
screen(6); legend(prep_df(df.sd.list$alpha_1.33, test.C = 1.5, Es.range = c(0,2)), 'payoff.comp.ip')
screen(8); legend(prep_df(df.sd.list$alpha_1.33, test.C = 2,   Es.range = c(0,3)), 'payoff.comp.ip')

# HEATMAPS, screens 1,3,5,7

screen(1)
plot.heatmap(prep_df(df.sd.list$alpha_1.33, test.C = 0.5, Es.range = c(0,1)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = TRUE)

mtext('C = 0.5', side = 3, line = 1, font = 2)

screen(3)
plot.heatmap(prep_df(df.sd.list$alpha_1.33, test.C = 1, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)


mtext('C = 1', side = 3, line = 1, font = 2)


screen(5)
plot.heatmap(prep_df(df.sd.list$alpha_1.33, test.C = 1.5, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)

mtext('C = 1.5', side = 3, line = 1, font = 2)

screen(7)
plot.heatmap(prep_df(df.sd.list$alpha_1.33, test.C = 2, Es.range = c(0,3)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose), 
             draw.y.axis = FALSE)

mtext('C = 2', side = 3, line = 1, font = 2)

mtext('a = 1.33', side = 4, line = 2.5, font = fty, cex = cex.labels)


# ALPHA = 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LEGENDS, screens 10,12,14,16
screen(10); legend(prep_df(df.sd.list$alpha_2, test.C = 0.5, Es.range = c(0,1)), 'payoff.comp.ip')
screen(12); legend(prep_df(df.sd.list$alpha_2, test.C = 1,   Es.range = c(0,2)), 'payoff.comp.ip')
screen(14); legend(prep_df(df.sd.list$alpha_2, test.C = 1.5, Es.range = c(0,2)), 'payoff.comp.ip')
screen(16); legend(prep_df(df.sd.list$alpha_2, test.C = 2,   Es.range = c(0,3)), 'payoff.comp.ip')

# HEATMAPS, screens 9,11,13,15

screen(9)
plot.heatmap(prep_df(df.sd.list$alpha_2, test.C = 0.5, Es.range = c(0,1)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = TRUE)

screen(11)
plot.heatmap(prep_df(df.sd.list$alpha_2, test.C = 1, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)

screen(13)
plot.heatmap(prep_df(df.sd.list$alpha_2, test.C = 1.5, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)

screen(15)
plot.heatmap(prep_df(df.sd.list$alpha_2, test.C = 2, Es.range = c(0,3)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose), 
             draw.y.axis = FALSE)

mtext('a = 2', side = 4, line = 2.5, font = fty, cex = cex.labels)


# ALPHA = 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LEGENDS, screens 18,20,22,24
screen(18); legend(prep_df(df.sd.list$alpha_4, test.C = 0.5, Es.range = c(0,1)), 'payoff.comp.ip')
screen(20); legend(prep_df(df.sd.list$alpha_4, test.C = 1,   Es.range = c(0,2)), 'payoff.comp.ip')
screen(22); legend(prep_df(df.sd.list$alpha_4, test.C = 1.5, Es.range = c(0,2)), 'payoff.comp.ip')
screen(24); legend(prep_df(df.sd.list$alpha_4, test.C = 2,   Es.range = c(0,3)), 'payoff.comp.ip')

# HEATMAPS, screens 17,19,21,23

screen(17)
plot.heatmap(prep_df(df.sd.list$alpha_4, test.C = 0.5, Es.range = c(0,1)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = TRUE)

mtext('Initial compliance frequency', side = 2, line = 2.2, font = fty, cex = cex.labels, adj = -1)

screen(19)
plot.heatmap(prep_df(df.sd.list$alpha_4, test.C = 1, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)

screen(21)
plot.heatmap(prep_df(df.sd.list$alpha_4, test.C = 1.5, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)

screen(23)
plot.heatmap(prep_df(df.sd.list$alpha_4, test.C = 2, Es.range = c(0,3)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose), 
             draw.y.axis = FALSE)


mtext('a = 4', side = 4, line = 2.5, font = fty, cex = cex.labels)


# ALPHA = 10 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LEGENDS, screens 26, 28, 30, 32
screen(26); legend(prep_df(df.sd.list$alpha_10, test.C = 0.5, Es.range = c(0,1)), 'payoff.comp.ip')
screen(28); legend(prep_df(df.sd.list$alpha_10, test.C = 1,   Es.range = c(0,2)), 'payoff.comp.ip')
mtext('Enforcement level (E)', side = 1, line = 6, font = fty, cex = cex.labels, padj = 0)
screen(30); legend(prep_df(df.sd.list$alpha_10, test.C = 1.5, Es.range = c(0,2)), 'payoff.comp.ip')
screen(32); legend(prep_df(df.sd.list$alpha_10, test.C = 2,   Es.range = c(0,3)), 'payoff.comp.ip')

# HEATMAPS, screens 25, 27, 29, 31

screen(25)
plot.heatmap(prep_df(df.sd.list$alpha_10, test.C = 0.5, Es.range = c(0,1)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = TRUE)

screen(27)
plot.heatmap(prep_df(df.sd.list$alpha_10, test.C = 1, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)


screen(29)
plot.heatmap(prep_df(df.sd.list$alpha_10, test.C = 1.5, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)

screen(31)
plot.heatmap(prep_df(df.sd.list$alpha_10, test.C = 2, Es.range = c(0,3)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose), 
             draw.y.axis = FALSE)


mtext('a = 10', side = 4, line = 2.5, font = fty, cex = cex.labels)


dev.off()


# Plot output TTI ----

load('./output/model_runs/TTI_various_tau_ptest_0.7_nbIni_60000_timeSim_800.RData')

close.screen(all.screens = TRUE)
fty = 2
cex.labels = 1


pdf('./output/figures/FigureS3_TTI_heatmap_vary_tau.pdf', width = (13+1+1)/2.54, height = (13+1+1)/2.54, pointsize = 7)
close.screen(all.screens = TRUE)
split.screen(m)
par(omi = rep(1/2.54,4))
for(i in 1:nrow(m)){
  screen(i)
  par(mar = c(0, 0, 0, 0))
  plot(1, axes = FALSE, type = "n")
  #box()
}

#dev.off()



# TAU = 10 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LEGENDS, screens 2 4 6 8
screen(2); legend(prep_df(df.tti.list$tau_10, test.C = 0.5, Es.range = c(0,1)), 'payoff.comp')
screen(4); legend(prep_df(df.tti.list$tau_10, test.C = 1,   Es.range = c(0,2)), 'payoff.comp')
mtext('TEST-TRACE-ISOLATE', side = 3, line = 3, font = 2, padj = 0, col = 'black')
screen(6); legend(prep_df(df.tti.list$tau_10, test.C = 1.5, Es.range = c(0,2)), 'payoff.comp')
screen(8); legend(prep_df(df.tti.list$tau_10, test.C = 2,   Es.range = c(0,3)), 'payoff.comp')

# HEATMAPS, screens 1, 3, 5, 7

screen(1)
plot.heatmap(prep_df(df.tti.list$tau_10, test.C = 0.5, Es.range = c(0,1)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = TRUE)

mtext('C = 0.5', side = 3, line = 1, font = 2)

screen(3)
plot.heatmap(prep_df(df.tti.list$tau_10, test.C = 1, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)

mtext('C = 1', side = 3, line = 1, font = 2)

screen(5)
plot.heatmap(prep_df(df.tti.list$tau_10, test.C = 1.5, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)

mtext('C = 1.5', side = 3, line = 1, font = 2)

screen(7)
plot.heatmap(prep_df(df.tti.list$tau_10, test.C = 2, Es.range = c(0,3)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose), 
             draw.y.axis = FALSE)

mtext('C = 2', side = 3, line = 1, font = 2)

mtext('tau = 10', side = 4, line = 2.5, font = fty, cex = cex.labels)



# TAU = 8 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LEGENDS, screens 10, 12, 14, 16
screen(10); legend(prep_df(df.tti.list$tau_8, test.C = 0.5, Es.range = c(0,1)), 'payoff.comp')
screen(12); legend(prep_df(df.tti.list$tau_8, test.C = 1,   Es.range = c(0,2)), 'payoff.comp')
screen(14); legend(prep_df(df.tti.list$tau_8, test.C = 1.5, Es.range = c(0,2)), 'payoff.comp')
screen(16); legend(prep_df(df.tti.list$tau_8, test.C = 2,   Es.range = c(0,3)), 'payoff.comp')

# HEATMAPS, screens 9, 11, 13, 15

screen(9)
plot.heatmap(prep_df(df.tti.list$tau_8, test.C = 0.5, Es.range = c(0,1)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = TRUE)


screen(11)
plot.heatmap(prep_df(df.tti.list$tau_8, test.C = 1, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)

screen(13)
plot.heatmap(prep_df(df.tti.list$tau_8, test.C = 1.5, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)

screen(15)
plot.heatmap(prep_df(df.tti.list$tau_8, test.C = 2, Es.range = c(0,3)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose), 
             draw.y.axis = FALSE)


mtext('tau = 8', side = 4, line = 2.5, font = fty, cex = cex.labels)



# TAU = 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LEGENDS, screens 18, 20, 22, 24
screen(10); legend(prep_df(df.tti.list$tau_5, test.C = 0.5, Es.range = c(0,1)), 'payoff.comp')
screen(12); legend(prep_df(df.tti.list$tau_5, test.C = 1,   Es.range = c(0,2)), 'payoff.comp')
screen(14); legend(prep_df(df.tti.list$tau_5, test.C = 1.5, Es.range = c(0,2)), 'payoff.comp')
screen(16); legend(prep_df(df.tti.list$tau_5, test.C = 2,   Es.range = c(0,3)), 'payoff.comp')

# HEATMAPS, screens 17, 19, 21, 23

screen(17)
plot.heatmap(prep_df(df.tti.list$tau_5, test.C = 0.5, Es.range = c(0,1)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = TRUE)

mtext('Initial compliance frequency', side = 2, line = 2.2, font = fty, cex = cex.labels, adj = -1)

screen(19)
plot.heatmap(prep_df(df.tti.list$tau_5, test.C = 1, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)

screen(21)
plot.heatmap(prep_df(df.tti.list$tau_5, test.C = 1.5, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)

screen(23)
plot.heatmap(prep_df(df.tti.list$tau_5, test.C = 2, Es.range = c(0,3)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose), 
             draw.y.axis = FALSE)

mtext('tau = 5', side = 4, line = 2.5, font = fty, cex = cex.labels)




# TAU = 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LEGENDS, screens 26, 28, 30, 32
screen(26); legend(prep_df(df.tti.list$tau_2, test.C = 0.5, Es.range = c(0,1)), 'payoff.comp')
screen(28); legend(prep_df(df.tti.list$tau_2, test.C = 1,   Es.range = c(0,2)), 'payoff.comp')
mtext('Enforcement level (E)', side = 1, line = 6, font = fty, cex = cex.labels, padj = 0)
screen(30); legend(prep_df(df.tti.list$tau_2, test.C = 1.5, Es.range = c(0,2)), 'payoff.comp')
screen(32); legend(prep_df(df.tti.list$tau_2, test.C = 2,   Es.range = c(0,3)), 'payoff.comp')

# HEATMAPS, screens 25, 27, 29, 31

screen(25)
plot.heatmap(prep_df(df.tti.list$tau_2, test.C = 0.5, Es.range = c(0,1)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = TRUE)

screen(27)
plot.heatmap(prep_df(df.tti.list$tau_2, test.C = 1, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)



screen(29)
plot.heatmap(prep_df(df.tti.list$tau_2, test.C = 1.5, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)


screen(31)
plot.heatmap(prep_df(df.tti.list$tau_2, test.C = 2, Es.range = c(0,3)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose), 
             draw.y.axis = FALSE)


mtext('tau = 2', side = 4, line = 2.5, font = fty, cex = cex.labels)


dev.off()


# Plot output TTI+IP ----

load('./output/model_runs/TTI_various_tau_ptest_0.7_nbIni_60000_timeSim_800.RData')

close.screen(all.screens = TRUE)
fty = 2
cex.labels = 1


pdf('./output/figures/FigureS4_TTI_IP_heatmap_vary_tau.pdf', width = (13+1+1)/2.54, height = (13+1+1)/2.54, pointsize = 7)
close.screen(all.screens = TRUE)
split.screen(m)
par(omi = rep(1/2.54,4))
for(i in 1:nrow(m)){
  screen(i)
  par(mar = c(0, 0, 0, 0))
  plot(1, axes = FALSE, type = "n")
  #box()
}

#dev.off()


# TAU = 10 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LEGENDS, screens 2,4,6,8
screen(2); legend(prep_df(df.tti.list$tau_10, test.C = 0.5, Es.range = c(0,1)), 'payoff.comp.ip')
screen(4); legend(prep_df(df.tti.list$tau_10, test.C = 1,   Es.range = c(0,2)), 'payoff.comp.ip')
mtext('TEST-TRACE-ISOLATE + IMMUNITY PASSPORT', side = 3, line = 3, font = 2, padj = 0, col = 'black')
screen(6); legend(prep_df(df.tti.list$tau_10, test.C = 1.5, Es.range = c(0,2)), 'payoff.comp.ip')
screen(8); legend(prep_df(df.tti.list$tau_10, test.C = 2,   Es.range = c(0,3)), 'payoff.comp.ip')

# HEATMAPS, screens 1,3,5,7

screen(1)
plot.heatmap(prep_df(df.tti.list$tau_10, test.C = 0.5, Es.range = c(0,1)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = TRUE)
mtext('C = 0.5', side = 3, line = 1, font = 2)

screen(3)
plot.heatmap(prep_df(df.tti.list$tau_10, test.C = 1, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)
mtext('C = 1', side = 3, line = 1, font = 2)


screen(5)
plot.heatmap(prep_df(df.tti.list$tau_10, test.C = 1.5, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)
mtext('C = 1.5', side = 3, line = 1, font = 2)

screen(7)
plot.heatmap(prep_df(df.tti.list$tau_10, test.C = 2, Es.range = c(0,3)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose), 
             draw.y.axis = FALSE)

mtext('C = 2', side = 3, line = 1, font = 2)

mtext('tau = 10', side = 4, line = 2.5, font = fty, cex = cex.labels)

# TAU = 8 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LEGENDS, screens 10, 12, 14, 16
screen(10); legend(prep_df(df.tti.list$tau_8, test.C = 0.5, Es.range = c(0,1)), 'payoff.comp.ip')
screen(12); legend(prep_df(df.tti.list$tau_8, test.C = 1,   Es.range = c(0,2)), 'payoff.comp.ip')
screen(14); legend(prep_df(df.tti.list$tau_8, test.C = 1.5, Es.range = c(0,2)), 'payoff.comp.ip')
screen(16); legend(prep_df(df.tti.list$tau_8, test.C = 2,   Es.range = c(0,3)), 'payoff.comp.ip')

# HEATMAPS, screens 9, 11, 13, 15

screen(9)
plot.heatmap(prep_df(df.tti.list$tau_8, test.C = 0.5, Es.range = c(0,1)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = TRUE)

screen(11)
plot.heatmap(prep_df(df.tti.list$tau_8, test.C = 1, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)

screen(13)
plot.heatmap(prep_df(df.tti.list$tau_8, test.C = 1.5, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)

screen(15)
plot.heatmap(prep_df(df.tti.list$tau_8, test.C = 2, Es.range = c(0,3)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose), 
             draw.y.axis = FALSE)


mtext('tau = 8', side = 4, line = 2.5, font = fty, cex = cex.labels)


# TAU = 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LEGENDS, screens 18, 20, 22, 24
screen(18); legend(prep_df(df.tti.list$tau_5, test.C = 0.5, Es.range = c(0,1)), 'payoff.comp.ip')
screen(20); legend(prep_df(df.tti.list$tau_5, test.C = 1,   Es.range = c(0,2)), 'payoff.comp.ip')
screen(22); legend(prep_df(df.tti.list$tau_5, test.C = 1.5, Es.range = c(0,2)), 'payoff.comp.ip')
screen(24); legend(prep_df(df.tti.list$tau_5, test.C = 2,   Es.range = c(0,3)), 'payoff.comp.ip')

# HEATMAPS, screens 17, 19, 21, 23

screen(17)
plot.heatmap(prep_df(df.tti.list$tau_5, test.C = 0.5, Es.range = c(0,1)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = TRUE)
mtext('Initial compliance frequency', side = 2, line = 2.2, font = fty, cex = cex.labels, adj = -1)

screen(19)
plot.heatmap(prep_df(df.tti.list$tau_5, test.C = 1, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)

screen(21)
plot.heatmap(prep_df(df.tti.list$tau_5, test.C = 1.5, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)

screen(23)
plot.heatmap(prep_df(df.tti.list$tau_5, test.C = 2, Es.range = c(0,3)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose), 
             draw.y.axis = FALSE)

mtext('tau = 5', side = 4, line = 2.5, font = fty, cex = cex.labels)


# TAU = 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LEGENDS, screens 26, 28, 30, 32
screen(26); legend(prep_df(df.tti.list$tau_2, test.C = 0.5, Es.range = c(0,1)), 'payoff.comp.ip')
screen(28); legend(prep_df(df.tti.list$tau_2, test.C = 1,   Es.range = c(0,2)), 'payoff.comp.ip')
mtext('Enforcement level (E)', side = 1, line = 6, font = fty, cex = cex.labels, padj = 0)
screen(30); legend(prep_df(df.tti.list$tau_2, test.C = 1.5, Es.range = c(0,2)), 'payoff.comp.ip')
screen(32); legend(prep_df(df.tti.list$tau_2, test.C = 2,   Es.range = c(0,3)), 'payoff.comp.ip')

# HEATMAPS, screens 25, 27, 29, 31

screen(25)
plot.heatmap(prep_df(df.tti.list$tau_2, test.C = 0.5, Es.range = c(0,1)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = TRUE)


screen(27)
plot.heatmap(prep_df(df.tti.list$tau_2, test.C = 1, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)


screen(29)
plot.heatmap(prep_df(df.tti.list$tau_2, test.C = 1.5, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = FALSE)


screen(31)
plot.heatmap(prep_df(df.tti.list$tau_2, test.C = 2, Es.range = c(0,3)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose), 
             draw.y.axis = FALSE)

mtext('tau = 2', side = 4, line = 2.5, font = fty, cex = cex.labels)



dev.off()


