# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#   Social dilemmas in non-pharmaceutical interventions to tackle epidemics
#                     Simonet C, Sweeny A, McNally L (2021)
#                       Codes for figures output
#
#                       last edited: July, 27th, 2021
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Script summary
# [1] defines code wrappers to compute payoff from model run output, and plot the heatmap, and define figure layout
# [2] plot main SD heatmap figure
# [3] plot main TTI heatmap figure


setwd('path/to/repo') # setwd('~/Documents/PhD/Research/Covid19/NPIpandemic_repo/')
library(readxl); library(ggplot2); library(dplyr); library(tidyr);library(Rmisc); library(lubridate);library(plotly);library(yarrr); library(knitr);library(plyr);library(utils);library(tidyverse);library(gridExtra);library(deSolve); library(stringi); library(stringr);
library(ggthemes)


# [1] setup  ----

# Code wrapper to take model output (from '2_mathematical_model.R') and prep it for plotting
#   defines a range of values of E and C to compute payoff over
#   compute payoff with and withou immunity passport

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
  
  if(is.na(test.C) == TRUE){
    payoff.df.focal.C<- payoff.df
    
  }else{
    payoff.df.focal.C<- payoff.df[payoff.df$C == test.C,]
    
  }
  
  return(payoff.df.focal.C)
  
}


# code wrapper to make main text heatmap figure
# payoff.df =  output of 'prep_df' code wrapper
# focal.payoff.comp = 
# zlim.fixed = weather range for zlims (heatmap color min-max range) should be fixed
# zlims = vector of two zlims of fixed = TRUE

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
  
  axis(2, at = c(0, 0.25, 0.5, 0.75, 1), cex.axis = 0.9, font = 1, lwd = 1)
  axis(1, at = c(0, 0.25, 0.5, 0.75, 1), labels = quantile(unique(payoffs.comp.df$E)), cex.axis = 0.9, font = 1, lwd = 1)
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
legend<- function(payoff.df.focal, focal.payoff.comp, zlim.fixed = FALSE, zlims = NA){
  
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
  
  if(zlim.fixed == TRUE){
    
    text(y = seq(0, 1, by=0.25),
         x = 0.5,
         labels = round(zlims[2]* c(-1, -0.5, -0, 0.5, 1), 2),
         font = 1, cex = 0.7, srt = 0, adj = 0)
    
  }else{
  text(y = c(0, 0.25, 0.75, 1),
       x = 0.5,
       labels = substr(formatC(
         round(max(abs(payoff.df.focal[,which(colnames(payoff.df.focal) == focal.payoff.comp)])) * c(-1, -0.5, 0.5, 1), 2), digit = 2, format = 's'),
         1, 4),
       font = 1, cex = 0.7, srt = 0, adj = 0)
  
  text(y = 0.5,
       x = 0.5,
       labels = '0',
       font = 1, cex = 0.7, srt = 0, adj = 0)
  
  }
  
  
}


# Some stuff for complicated figure layout
w = 13
h = 12
s = (5-0.6)/3

horizontal<- c(0,
               s,
               s+0.3,
               s+0.3+s,
               s+0.3+s+0.3,
               s+0.3+s+0.3+s,
               s+0.3+s+0.3+s+0.2,
               s+0.3+s+0.3+s+0.2+0.8,
               s+0.3+s+0.3+s+0.2+0.8+1,
               s+0.3+s+0.3+s+0.2+0.8+1+s,
               s+0.3+s+0.3+s+0.2+0.8+1+s+0.3,
               s+0.3+s+0.3+s+0.2+0.8+1+s+0.3+s,
               s+0.3+s+0.3+s+0.2+0.8+1+s+0.3+s+0.3,
               s+0.3+s+0.3+s+0.2+0.8+1+s+0.3+s+0.3+s,
               s+0.3+s+0.3+s+0.2+0.8+1+s+0.3+s+0.3+s+0.2,
               s+0.3+s+0.3+s+0.2+0.8+1+s+0.3+s+0.3+s+0.2+0.8)

vertical<-
  c(0,
    s, 
    s+0.3,
    s+0.3+s,
    s+0.3+s+0.3,
    s+0.3+s+0.3+s,
    s+0.3+s+0.3+s+2,
    s+0.3+s+0.3+s+2+5)


round(horizontal/w, 3)
round(vertical/h, 3)


m  <- rbind(
  c(0.000, 0.385, 0.583, 1.000), # 1
  c(0.400, 0.462, 0.792, 1.000), # 2
  
  c(0.000, 0.113, 0.294, 0.417),# 3
  c(0.000, 0.113, 0.147, 0.269),# 4
  c(0.000, 0.113, 0.000, 0.122),# 5
  
  c(0.136, 0.249, 0.294, 0.417),# 6
  c(0.136, 0.249, 0.147, 0.269),# 7
  c(0.136, 0.249, 0.000, 0.122),# 8
  
  c(0.272, 0.385, 0.294, 0.417),# 9
  c(0.272, 0.385, 0.147, 0.269),# 10
  c(0.272, 0.385, 0.000, 0.122),# 11
  
  c(0.538, 0.923, 0.583, 1.000),# 12
  c(0.938, 1.000, 0.792, 1.000),# 13
  
  c(0.538, 0.651, 0.294, 0.417),# 14
  c(0.538, 0.651, 0.147, 0.269),# 15
  c(0.538, 0.651, 0.000, 0.122),# 16
  
  c(0.674, 0.787, 0.294, 0.417),# 17
  c(0.674, 0.787, 0.147, 0.269),# 18
  c(0.674, 0.787, 0.000, 0.122),# 19
  
  c(0.810, 0.923, 0.294, 0.417),# 20
  c(0.810, 0.923, 0.147, 0.269),# 21
  c(0.810, 0.923, 0.000, 0.122)) # 22


close.screen(all.screens = TRUE)
fty = 2
cex.labels = 1



# [2] Figure 3 - SD main heatmap ----
load('./output/model_runs/SIMULATION_SIR_SD_June01_nbIni_60000_alpha_4.RData')

zlim.choose = 2


pdf('./output/figures/Figure3_SD_MAIN.pdf', width = (13+1+1)/2.54, height = (12+1+1)/2.54, pointsize = 7)

close.screen(all.screens = TRUE)
split.screen(m)
par(omi = rep(1/2.54,4))
for(i in 1:nrow(m)){
  screen(i)
  par(mar = c(0, 0, 0, 0))
  plot(1, axes = FALSE, type = "n")
  #box()
  }


# main pannel SD ----
screen(1)

plot.heatmap(prep_df(df.sd, test.C = 1.2, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = TRUE)

mtext('Initial compliance frequency', side = 2, line = 2.2, font = fty, cex = cex.labels)
mtext('Enforcement level (E)', side = 1, line = 2.2, font = fty, cex = cex.labels)
mtext('SOCIAL DISTANCING', side = 3, line = 1, font = 2)

screen(2); legend(prep_df(df.sd, test.C = 1.2, Es.range = c(0,2)), 'payoff.comp', zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose))


# main pannel SD+IP ----

screen(12)

plot.heatmap(prep_df(df.sd, test.C = 1.2, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = TRUE)

mtext('Initial compliance frequency', side = 2, line = 2.2, font = fty, cex = cex.labels)
mtext('Enforcement level (E)', side = 1, line = 2.2, font = fty, cex = cex.labels)
mtext('SOCIAL DISTANCING + IMMUNITY PASSPORT', side = 3, line = 1, font = 2)

screen(13); legend(prep_df(df.sd, test.C = 1.2, Es.range = c(0,2)), 'payoff.comp.ip', zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose))



# insets: SD, low E ----

payoff.df<- prep_df(df.sd, test.C = 1.2, Es.range = c(0,2))
test.C = 1.2

# screen 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
payoff.df2<- payoff.df[payoff.df$E == 0.4,]; focal.E = 0.4

screen(3)
plot(p_not_infected_n ~ Sc, data = payoff.df2, col = 'red', type = 'l',
     ylab = '', xlab = '', ylim = c(0, 1),
     xaxt = 'n', yaxt = 'n')
lines(p_not_infected_c ~ Sc, data = payoff.df2, col = 'dodgerblue')
axis(2,payoff.df2[payoff.df2$Sc ==  min(payoff.df2$Sc),c('p_not_infected_c', 'p_not_infected_n')], labels = c('a', 'b'), lwd.ticks = 0, line = -0.8, lwd = 0,las = 2, cex.axis = 0.8)
mtext('Benefit', side = 2, line = 2, font = fty, cex = cex.labels)
mtext('Low E', side = 3, line = 0.25, font = fty, cex = cex.labels)


# screen 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

screen(4)
plot(C ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '', ylim = c(-0.25, 2.1),
     xaxt = 'n', yaxt = 'n')
lines(E ~ Sc, data = payoff.df2, col = 'red')
axis(2, c(test.C, focal.E), labels = c('c', 'd'), lwd.ticks = 0, line = -0.8, lwd = 0,las = 2, cex.axis = 0.8)
mtext('Cost', side = 2, line = 2, font = fty, cex = cex.labels)

# screen 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

payoff.lims<-c(floor(min(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n', 'pi_c')])),
ceiling(max(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n', 'pi_c')])))

screen(5)
plot(pi_c ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '',
     #ylim = payoff.lims,
     ylim = c(-1.8, 0.6),
     xaxt = 'n', yaxt = 'n')
lines(pi_n ~ Sc, data = payoff.df2, col = 'red')

axis(2,
     at = (payoff.df2[payoff.df2$Sc ==  min(payoff.df2$Sc),c('p_not_infected_c', 'p_not_infected_n')] - c(test.C, focal.E)),
     labels = c('a-c', 'b-d'),
     lwd.ticks = 0, line = -0.8, lwd = 0,las = 2, cex.axis = 0.8)

mtext('Payoff', side = 2, line = 2, font = fty, cex = cex.labels)


# insets: SD, mid E ----

# screen 6 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
payoff.df2<- payoff.df[payoff.df$E == 1,]; focal.E = 1
screen(6)
plot(p_not_infected_n ~ Sc, data = payoff.df2, col = 'red', type = 'l',
     ylab = '', xlab = '', ylim = c(0, 1),
     xaxt = 'n', yaxt = 'n')
lines(p_not_infected_c ~ Sc, data = payoff.df2, col = 'dodgerblue')
mtext('Intermediate E', side = 3, line = 0.25, font = fty, cex = cex.labels)

# screen 7 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
screen(7)
plot(C ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '', ylim = c(-0.25, 2.1),
     xaxt = 'n', yaxt = 'n')
lines(E ~ Sc, data = payoff.df2, col = 'red')

# screen 8 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

payoff.lims<-c(floor(min(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n', 'pi_c')])),
               ceiling(max(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n', 'pi_c')])))

screen(8)
plot(pi_c ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '',
     #ylim = payoff.lims,
     ylim = c(-1.8, 0.6),
     xaxt = 'n', yaxt = 'n')
lines(pi_n ~ Sc, data = payoff.df2, col = 'red')


focalpoint = which.min(abs(payoff.df2$pi_c - payoff.df2$pi_n))


points(x = payoff.df2[focalpoint,'Sc'],
       y = payoff.df2[focalpoint,'pi_c'],
       pch = 16, cex = 1.1)
points(x = payoff.df2[focalpoint,'Sc'],
       y = payoff.df2[focalpoint,'pi_c'],
       pch = 1, cex = 1.1)


mtext('Compliance frequency', side = 1, line = 1.5, font = fty, cex = cex.labels)


# insets: SD, high E ----

# screen 9 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
payoff.df2<- payoff.df[payoff.df$E == 1.8,]; focal.E = 1.8
screen(9)
plot(p_not_infected_n ~ Sc, data = payoff.df2, col = 'red', type = 'l',
     ylab = '', xlab = '', ylim = c(0, 1),
     xaxt = 'n', yaxt = 'n')
lines(p_not_infected_c ~ Sc, data = payoff.df2, col = 'dodgerblue')
mtext('High E', side = 3, line = 0.25, font = fty, cex = cex.labels)

# screen 10 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
screen(10)
plot(C ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '', ylim = c(-0.25, 2.1),
     xaxt = 'n', yaxt = 'n')
lines(E ~ Sc, data = payoff.df2, col = 'red')

# screen 11 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

payoff.lims<-c(floor(min(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n', 'pi_c')])),
               ceiling(max(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n', 'pi_c')])))

screen(11)
plot(pi_c ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '',
     #ylim = payoff.lims,
     ylim = c(-1.8, 0.6),
     xaxt = 'n', yaxt = 'n')
lines(pi_n ~ Sc, data = payoff.df2, col = 'red')




# insets: SD+IP, low E ----
# screen 14 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
payoff.df2<- payoff.df[payoff.df$E == 0.4,]; focal.E = 0.4

screen(14)
plot(p_not_infected_n ~ Sc, data = payoff.df2, col = 'red', type = 'l',
     ylab = '', xlab = '', ylim = c(0, 1),
     xaxt = 'n', yaxt = 'n')
lines(p_not_infected_c ~ Sc, data = payoff.df2, col = 'dodgerblue')

axis(2,payoff.df2[payoff.df2$Sc ==  min(payoff.df2$Sc),c('p_not_infected_c', 'p_not_infected_n')], labels = c('a', 'b'), lwd.ticks = 0, line = -0.8, lwd = 0,las = 2, cex.axis = 0.8)

mtext('Benefit', side = 2, line = 2, font = fty, cex = cex.labels)
mtext('Low E', side = 3, line = 0.25, font = fty, cex = cex.labels)

# screen 15 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
screen(15)
plot(prop_time_not_recovered_c * C ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '', ylim = c(-0.25, 2.1),
     xaxt = 'n', yaxt = 'n')
lines(prop_time_not_recovered_n * E ~ Sc, data = payoff.df2, col = 'red')
axis(2,
     c(test.C * tail(payoff.df2$prop_time_not_recovered_c, 1), focal.E * tail(payoff.df2$prop_time_not_recovered_n, 1)),
     labels = c('c', 'd'), lwd.ticks = 0, line = -0.8, lwd = 0,las = 2, cex.axis = 0.8)
mtext('Cost', side = 2, line = 2, font = fty, cex = cex.labels)


# screen 16 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

payoff.lims<-c(floor(min(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n.ip', 'pi_c.ip')])),
               ceiling(max(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n.ip', 'pi_c.ip')])))

screen(16)
plot(pi_c.ip ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '',
     #ylim = payoff.lims,
     ylim = c(-1, 0.6),
     xaxt = 'n', yaxt = 'n')
lines(pi_n.ip ~ Sc, data = payoff.df2, col = 'red')

axis(2,
     payoff.df2[payoff.df2$Sc ==  min(payoff.df2$Sc),c('p_not_infected_c', 'p_not_infected_n')] - c(test.C * tail(payoff.df2$prop_time_not_recovered_c, 1), focal.E * tail(payoff.df2$prop_time_not_recovered_n, 1)),
     labels = c('a-c', 'b-d'),
     lwd.ticks = 0, line = -0.8, lwd = 0,las = 2, cex.axis = 0.8)

axis(2,
     at = payoff.df2[payoff.df2$Sc ==  min(payoff.df2$Sc),c('p_not_infected_n')] - c(focal.E * tail(payoff.df2$prop_time_not_recovered_n, 1)),
     #at = 0.006,
     labels = 'b-d',
     lwd.ticks = 0, line = -0.8, lwd = 0,las = 2, cex.axis = 0.8)


mtext('Payoff', side = 2, line = 2, font = fty, cex = cex.labels)


# insets: SD+IP, mid E ----
# screen 17 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
payoff.df2<- payoff.df[payoff.df$E == 1.25,]
screen(17)
plot(p_not_infected_n ~ Sc, data = payoff.df2, col = 'red', type = 'l',
     ylab = '', xlab = '', ylim = c(0, 1),
     xaxt = 'n', yaxt = 'n')
lines(p_not_infected_c ~ Sc, data = payoff.df2, col = 'dodgerblue')
mtext('Intermediate E', side = 3, line = 0.25, font = fty, cex = cex.labels)

# screen 18 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
screen(18)
plot(prop_time_not_recovered_c * C ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '', ylim = c(-0.25, 2.1),
     xaxt = 'n', yaxt = 'n')
lines(prop_time_not_recovered_n * E ~ Sc, data = payoff.df2, col = 'red')
#axis(2,
#     c(1.5 * tail(payoff.df2$prop_time_not_recovered_c, 1), 0.25 * tail(payoff.df2$prop_time_not_recovered_n, 1)),
#     labels = c('c', 'd'), lwd.ticks = 0, line = -0.8, lwd = 0,las = 2, cex.axis = 0.8)
#mtext('Cost', side = 2, line = 2, font = fty, cex = cex.labels)


# screen 19 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

payoff.lims<-c(floor(min(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n.ip', 'pi_c.ip')])),
               ceiling(max(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n.ip', 'pi_c.ip')])))

screen(19)
plot(pi_c.ip ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '',
     #ylim = payoff.lims,
     ylim = c(-0.5, 0.1),
     xaxt = 'n', yaxt = 'n')
lines(pi_n.ip ~ Sc, data = payoff.df2, col = 'red')

#axis(2,
#     payoff.df2[payoff.df2$Sc ==  min(payoff.df2$Sc),c('p_not_infected_c', 'p_not_infected_n')] - c(1.5 * tail(payoff.df2$prop_time_not_recovered_c, 1), 1.6 * tail(payoff.df2$prop_time_not_recovered_n, 1)),
#     labels = c('a-c', 'b-d'),
#     lwd.ticks = 0, line = -0.8, lwd = 0,las = 2, cex.axis = 0.8)

focalpoint = which.min(abs(payoff.df2$pi_c.ip - payoff.df2$pi_n.ip))

points(x = payoff.df2[focalpoint,'Sc'],
       y = payoff.df2[focalpoint,'pi_c.ip'],
       pch = 1, cex = 1.1)
mtext('Compliance frequency', side = 1, line = 1.5, font = fty, cex = cex.labels)

# insets: SD+IP, high E ----
# screen 20 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

payoff.df2<- payoff.df[payoff.df$E == 1.999,]; focal.E = 1.999
screen(20)
plot(p_not_infected_n ~ Sc, data = payoff.df2, col = 'red', type = 'l',
     ylab = '', xlab = '', ylim = c(0, 1),
     xaxt = 'n', yaxt = 'n')
lines(p_not_infected_c ~ Sc, data = payoff.df2, col = 'dodgerblue')
mtext('High E', side = 3, line = 0.25, font = fty, cex = cex.labels)

# screen 21 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
screen(21)
plot(prop_time_not_recovered_c * C ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '', ylim = c(-0.25, 2.1),
     xaxt = 'n', yaxt = 'n')
lines(prop_time_not_recovered_n * E ~ Sc, data = payoff.df2, col = 'red')
#axis(2,
#     c(1.5 * tail(payoff.df2$prop_time_not_recovered_c, 1), 3 * tail(payoff.df2$prop_time_not_recovered_n, 1)),
#     labels = c('c', 'd'), lwd.ticks = 0, line = -0.8, lwd = 0,las = 2, cex.axis = 0.8)
#mtext('Cost', side = 2, line = 2, font = fty, cex = cex.labels)


# screen 22 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

payoff.lims<-c(floor(min(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n.ip', 'pi_c.ip')])),
               ceiling(max(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n.ip', 'pi_c.ip')])))

screen(22)
plot(pi_c.ip ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '',
     #ylim = payoff.lims,
     ylim = c(-1, 0.6),
     xaxt = 'n', yaxt = 'n')
lines(pi_n.ip ~ Sc, data = payoff.df2, col = 'red')

#axis(2,
#     payoff.df2[payoff.df2$Sc ==  min(payoff.df2$Sc),c('p_not_infected_c', 'p_not_infected_n')] - c(1.5 * tail(payoff.df2$prop_time_not_recovere#d_c, 1), 3 * tail(payoff.df2$prop_time_not_recovered_n, 1)),
#     labels = c('a-c', 'b-d'),
#     lwd.ticks = 0, line = -0.8, lwd = 0,las = 2, cex.axis = 0.8)


# close pdf ----
dev.off()



# TTI ----

load('./output/model_runs/SIMULATION_SIR_TTI_June01_nbIni_60000_tau_5.1_ptest_0.7.RData')

zlim.choose = 1.2

pdf('./output/figures/Figure4_TTI_MAIN.pdf', width = (13+1+1)/2.54, height = (12+1+1)/2.54, pointsize = 7)

close.screen(all.screens = TRUE)
split.screen(m)
par(omi = rep(1/2.54,4))
for(i in 1:nrow(m)){
  screen(i)
  par(mar = c(0, 0, 0, 0))
  plot(1, axes = FALSE, type = "n")
  #box()
}

# main pannel TTI ----

screen(1)

plot.heatmap(prep_df(df.tti, test.C = 1.2, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp',
             zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = TRUE)

mtext('Initial compliance frequency', side = 2, line = 2.2, font = fty, cex = cex.labels)
mtext('Enforcement level (E)', side = 1, line = 2.2, font = fty, cex = cex.labels)
mtext('TEST-TRACE-ISOLATE', side = 3, line = 1, font = 2)

screen(2); legend(prep_df(df.tti, test.C = 1.2, Es.range = c(0,2)), 'payoff.comp', zlim.fixed = FALSE, zlims = c(-zlim.choose, zlim.choose))

# main pannel TTI+IP ----

screen(12) # force zlim to 1.2 for this one, to compare better to left map

plot.heatmap(prep_df(df.tti, test.C = 1.2, Es.range = c(0,2)),
             focal.payoff.comp = 'payoff.comp.ip',
             zlim.fixed = TRUE, zlims = c(-zlim.choose, zlim.choose),
             draw.y.axis = TRUE)

mtext('Initial compliance frequency', side = 2, line = 2.2, font = fty, cex = cex.labels)
mtext('Enforcement level (E)', side = 1, line = 2.2, font = fty, cex = cex.labels)
mtext('TEST-TRACE-ISOLATE + IMMUNITY PASSPORT', side = 3, line = 1, font = 2)

screen(13); legend(prep_df(df.tti, test.C = 1.2, Es.range = c(0,2)), 'payoff.comp.ip', zlim.fixed = TRUE, zlims = c(-zlim.choose, zlim.choose))




# insets: TTI, low E ----

payoff.df<- prep_df(df.tti, test.C = 1.2, Es.range = c(0,2))
test.C = 1.2

# screen 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

payoff.df2<- payoff.df[payoff.df$E == 0.5,]; focal.E = 0.5
screen(3)
plot(p_not_infected_n ~ Sc, data = payoff.df2, type = 'l',
     ylab = '', xlab = '', ylim = c(0, 1),
     xaxt = 'n', yaxt = 'n',
     lwd = 1, lty = 1, col = 'red')
lines(p_not_infected_c ~ Sc, data = payoff.df2, lwd = 1, lty = 2, col = 'dodgerblue')

axis(2,
     payoff.df2[payoff.df2$Sc ==  min(payoff.df2$Sc),c('p_not_infected_c', 'p_not_infected_n')], labels = rep('a = b', 2), lwd.ticks = 0, line = -0.8, lwd = 0,las = 2, cex.axis = 0.8)

mtext('Benefit', side = 2, line = 2, font = fty, cex = cex.labels)
mtext('Low E', side = 3, line = 0.25, font = fty, cex = cex.labels)


# screen 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
screen(4)
plot(C ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '', ylim = c(-0.25, 2.1),
     xaxt = 'n', yaxt = 'n')
lines(E ~ Sc, data = payoff.df2, col = 'red')
axis(2, c(test.C, focal.E), labels = c('c', 'd'), lwd.ticks = 0, line = -0.8, lwd = 0,las = 2, cex.axis = 0.8)
mtext('Cost', side = 2, line = 2, font = fty, cex = cex.labels)


# screen 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

payoff.lims.screen5<-c(min(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n', 'pi_c')]),
               (max(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n', 'pi_c')]))) * c((1+0.1),(1+0.1))

screen(5)
plot(pi_c ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '',
     #ylim = payoff.lims.screen5,
     ylim = c(-2, 0.5),
     xaxt = 'n', yaxt = 'n')
lines(pi_n ~ Sc, data = payoff.df2, col = 'red')

axis(2,
     payoff.df2[payoff.df2$Sc ==  min(payoff.df2$Sc),c('p_not_infected_c', 'p_not_infected_n')] - c(test.C, focal.E),
     labels = c('a-c', 'b-d'),
     lwd.ticks = 0, line = -0.8, lwd = 0,las = 2, cex.axis = 0.8)

mtext('Payoff', side = 2, line = 2, font = fty, cex = cex.labels)


# insets: TTI, mid E ----
# screen 6 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
payoff.df2<- payoff.df[payoff.df$E == 1.2,]; focal.E = 1.2
screen(6)
plot(p_not_infected_n ~ Sc, data = payoff.df2, col = 'red', type = 'l',
     ylab = '', xlab = '', ylim = c(0, 1),
     xaxt = 'n', yaxt = 'n')
lines(p_not_infected_c ~ Sc, data = payoff.df2, col = 'dodgerblue', lty = 2)
mtext('Intermediate E', side = 3, line = 0.25, font = fty, cex = cex.labels)

# screen 7 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
screen(7)
plot(C ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '', ylim = c(-0.25, 2.1),
     xaxt = 'n', yaxt = 'n')
lines(E ~ Sc, data = payoff.df2, col = 'red', lty = 2)

# screen 8 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# payoff.lims<-c(min(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n', 'pi_c')]),
#                max(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n', 'pi_c')]))

screen(8)
plot(pi_c ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '',
     #ylim = payoff.lims.screen5,
     ylim = c(-2, 0.5),
     xaxt = 'n', yaxt = 'n')
lines(pi_n ~ Sc, data = payoff.df2, col = 'red', lty = 2)

mtext('Compliance frequency', side = 1, line = 1.5, font = fty, cex = cex.labels)


# insets: TTI, high E ----
# screen 9 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
payoff.df2<- payoff.df[payoff.df$E == 2,]; focal.E = 1.9
screen(9)
plot(p_not_infected_n ~ Sc, data = payoff.df2, col = 'red', type = 'l',
     ylab = '', xlab = '', ylim =  c(0, 1),
     xaxt = 'n', yaxt = 'n')
lines(p_not_infected_c ~ Sc, data = payoff.df2, col = 'dodgerblue', lty = 2)
mtext('High E', side = 3, line = 0.25, font = fty, cex = cex.labels)

# screen 10 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
screen(10)
plot(C ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '', ylim = c(-0.25, 2.1),
     xaxt = 'n', yaxt = 'n')
lines(E ~ Sc, data = payoff.df2, col = 'red')

# screen 11 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

payoff.lims.screen11<-c(min(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n', 'pi_c')]),
                       (max(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n', 'pi_c')]))) * c((1+0.1),(1-0.1))

screen(11)
plot(pi_c ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '',
     #ylim = payoff.lims.screen11,
     ylim = c(-2, 0.5),
     xaxt = 'n', yaxt = 'n')
lines(pi_n ~ Sc, data = payoff.df2, col = 'red')



# insets: TTI+IP, low E ----
# screen 14 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
payoff.df2<- payoff.df[payoff.df$E == 0.5,]; focal.E = 0.5
screen(14)
plot(p_not_infected_n ~ Sc, data = payoff.df2, col = 'red', type = 'l',
     ylab = '', xlab = '', ylim =  c(0, 1),
     xaxt = 'n', yaxt = 'n')
lines(p_not_infected_c ~ Sc, data = payoff.df2, col = 'dodgerblue', lty = 2)

axis(2,payoff.df2[payoff.df2$Sc ==  min(payoff.df2$Sc),c('p_not_infected_c', 'p_not_infected_n')], labels = rep('a = b', 2), lwd.ticks = 0, line = -0.8, lwd = 0,las = 2, cex.axis = 0.8)

mtext('Benefit', side = 2, line = 2, font = fty, cex = cex.labels)
mtext('Low E', side = 3, line = 0.25, font = fty, cex = cex.labels)

# screen 15 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
screen(15)
plot(prop_time_not_recovered_c * C ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '', ylim = c(-0.25, 2.1),
     xaxt = 'n', yaxt = 'n')
lines(prop_time_not_recovered_n * E ~ Sc, data = payoff.df2, col = 'red')
axis(2,
     c(test.C * tail(payoff.df2$prop_time_not_recovered_c, 1), focal.E * tail(payoff.df2$prop_time_not_recovered_n, 1)),
     labels = c('c', 'd'), lwd.ticks = 0, line = -0.8, lwd = 0,las = 2, cex.axis = 0.8)
mtext('Cost', side = 2, line = 2, font = fty, cex = cex.labels)


# screen 16 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

payoff.lims.screen16<-c(min(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n.ip', 'pi_c.ip')]),
               max(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n.ip', 'pi_c.ip')])) * c(1+0.1, 1+2)

#payoff.lims.screen16<- c(-2, 2)

screen(16)
plot(pi_c.ip ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '',
     #ylim = payoff.lims.screen16,
     ylim = c(-1, 0.5),
     xaxt = 'n', yaxt = 'n')
lines(pi_n.ip ~ Sc, data = payoff.df2, col = 'red')

axis(2,
     payoff.df2[payoff.df2$Sc ==  min(payoff.df2$Sc),c('p_not_infected_c', 'p_not_infected_n')] - c(test.C * tail(payoff.df2$prop_time_not_recovered_c, 1), focal.E * tail(payoff.df2$prop_time_not_recovered_n, 1)),
     labels = c('a-c', 'b-d'),
     lwd.ticks = 0, line = -0.8, lwd = 0,las = 2, cex.axis = 0.8)

mtext('Payoff', side = 2, line = 2, font = fty, cex = cex.labels)


# insets: TTI+IP, mid E ----
# screen 17 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
payoff.df2<- payoff.df[payoff.df$E == 1.2,]

screen(17)
plot(p_not_infected_n ~ Sc, data = payoff.df2, col = 'red', type = 'l',
     ylab = '', xlab = '', ylim =  c(0, 1),
     xaxt = 'n', yaxt = 'n')
lines(p_not_infected_c ~ Sc, data = payoff.df2, col = 'dodgerblue', lty = 2)
mtext('Intermediate E', side = 3, line = 0.25, font = fty, cex = cex.labels)

# screen 18 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
screen(18)
plot(prop_time_not_recovered_c * C ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '', ylim = c(-0.25, 2.1),
     xaxt = 'n', yaxt = 'n')
lines(prop_time_not_recovered_n * E ~ Sc, data = payoff.df2, col = 'red', lty = 2)
#axis(2,
#     c(1.5 * tail(payoff.df2$prop_time_not_recovered_c, 1), 0.25 * tail(payoff.df2$prop_time_not_recovered_n, 1)),
#     labels = c('c', 'd'), lwd.ticks = 0, line = -0.8, lwd = 0,las = 2, cex.axis = 0.8)
#mtext('Cost', side = 2, line = 2, font = fty, cex = cex.labels)


# screen 19 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#payoff.lims.screen19<- c(-0.5, -0.15)

payoff.lims.screen19<-c(min(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n.ip', 'pi_c.ip')]),
                        max(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n.ip', 'pi_c.ip')])) * c(1+0.1, 1+2)

screen(19)
plot(pi_c.ip ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '',
     #ylim = payoff.lims.screen16,
     ylim = c(-1, 0.5),
     xaxt = 'n', yaxt = 'n')
lines(pi_n.ip ~ Sc, data = payoff.df2, col = 'red', lty = 2)

#axis(2,
#     payoff.df2[payoff.df2$Sc ==  min(payoff.df2$Sc),c('p_not_infected_c', 'p_not_infected_n')] - c(1.5 * tail(payoff.df2$prop_time_not_recovered_c, 1), 1.6 * tail(payoff.df2$prop_time_not_recovered_n, 1)),
#     labels = c('a-c', 'b-d'),
#     lwd.ticks = 0, line = -0.8, lwd = 0,las = 2, cex.axis = 0.8)

mtext('Compliance frequency', side = 1, line = 1.5, font = fty, cex = cex.labels)

# insets: TTI+IP, high E ----
# screen 20 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
payoff.df2<- payoff.df[payoff.df$E == 1.999,]
screen(20)
plot(p_not_infected_n ~ Sc, data = payoff.df2, col = 'red', type = 'l',
     ylab = '', xlab = '', ylim =  c(0, 1),
     xaxt = 'n', yaxt = 'n')
lines(p_not_infected_c ~ Sc, data = payoff.df2, col = 'dodgerblue', lty = 2)
mtext('High E', side = 3, line = 0.25, font = fty, cex = cex.labels)

# screen 21 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
screen(21)
plot(prop_time_not_recovered_c * C ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '', ylim = c(-0.25, 2.1),
     xaxt = 'n', yaxt = 'n')
lines(prop_time_not_recovered_n * E ~ Sc, data = payoff.df2, col = 'red')
#axis(2,
#     c(1.5 * tail(payoff.df2$prop_time_not_recovered_c, 1), 3 * tail(payoff.df2$prop_time_not_recovered_n, 1)),
#     labels = c('c', 'd'), lwd.ticks = 0, line = -0.8, lwd = 0,las = 2, cex.axis = 0.8)
#mtext('Cost', side = 2, line = 2, font = fty, cex = cex.labels)


# screen 22 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

payoff.lims.screen22<-c(floor(min(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n.ip', 'pi_c.ip')])),
                ceiling(max(payoff.df2[payoff.df2$Sc %in% c(min(payoff.df2$Sc), max(payoff.df2$Sc)),c('pi_n.ip', 'pi_c.ip')])))
#payoff.lims.screen22<- c(-0.65, -0.4)

#payoff.lims.screen22<- c(-0.77, -0.33)
screen(22)
plot(pi_c.ip ~ Sc, data = payoff.df2, col = 'dodgerblue', type = 'l',
     ylab = '', xlab = '',
     #ylim = payoff.lims.screen22,
     ylim = c(-1, 0.5),
     xaxt = 'n', yaxt = 'n')
lines(pi_n.ip ~ Sc, data = payoff.df2, col = 'red')

#axis(2,
#     payoff.df2[payoff.df2$Sc ==  min(payoff.df2$Sc),c('p_not_infected_c', 'p_not_infected_n')] - c(1.5 * tail(payoff.df2$prop_time_not_recovere#d_c, 1), 3 * tail(payoff.df2$prop_time_not_recovered_n, 1)),
#     labels = c('a-c', 'b-d'),
#     lwd.ticks = 0, line = -0.8, lwd = 0,las = 2, cex.axis = 0.8)





# close pdf ----
dev.off()















