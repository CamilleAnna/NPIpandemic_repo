
## Analysing effects of mobility on stringency and vice versa 

# packages ----------------------------------------------------------------

library(tidyverse)
library(janitor)
library(lubridate)
library(MCMCglmm)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(patchwork)
library(DescTools)


# plot set-up ------------------------------------------------------------
require(nationalparkcolors)
require(ggthemes)
require(ggflags)
require(ggpubr)

dir.create("ManuscriptFigures")
dir.create("ManuscriptOutput")


sparkPal <- park_palette("GeneralGrant", 8)[-5]

addSmallLegend <- function(myPlot, pointSize = 1, lineSize=1.5, 
                           textSize = 10, spaceLegend = 0.8) {
  myPlot +
    guides(#shape = guide_legend(override.aes = list(size = pointSize)),
      #line= guide_legend(override.aes = list(size=lineSize)), 
      color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize))
  #legend.key.size = unit(spaceLegend, "lines"))
}


# functions ---------------------------------------------------------------

source("Functions/CleanMCMC.R")

shift <- function(x, n){
  c(x[-(seq(n))], rep(NA, n))
}

# read in data  -----------------------------------------------------------

mobility<- read.csv("ManuscriptData/MobilityInterventionDF_MS.csv") %>% 
  clean_names(case="upper_camel") %>% 
  mutate(Date=as.Date(Date)) %>% 
  filter(Date < "2020-06-01") # 13696 - filter to first few months 

# data prep ---------------------------------------------------------------

mobilityCats <- mobility %>% select(4:9, StringencyIndexForDisplay) %>% colnames() %>% str_replace(., "PercentChangeFromBaseline", "")
  
mobilityLocalities <-mobility %>% pull(CountryName) %>% unique() %>%  as.character()
CountryCodes <- mobility %>% select(CountryName, CountryCode2) %>% filter(!duplicated(CountryName))

CountryList <- read.csv("Data/countryList_csv.csv", stringsAsFactors=FALSE,
                        colClasses = c("character"),  na.strings="") %>% 
  rename(CountryName=Name, CountryCode2=Code) %>% mutate(CountryCode2=str_to_lower(CountryCode2)) # just in case 

setdiff(mobilityLocalities, CountryList$CountryName)

mobility %>% 
  select(-CountryCode2) %>% 
  left_join(., CountryCodes, by="CountryName") %>% 
  rename_at(vars(ends_with("PercentChangeFromBaseline")),.funs = funs(sub("PercentChangeFromBaseline", "", .))) %>% 
  group_by(CountryName) %>% 
  select(matches("Country"), StringencyIndexForDisplay, Date, mobilityCats) %>% #select(-StringencyIndex) %>% 
  filter(!is.na(GroceryAndPharmacy), !is.na(Residential), !is.na(StringencyIndexForDisplay))  %>% 
  mutate_at(mobilityCats, .funs=list(First= ~ .x-lag(.x)))  -> mobilityDF # 13482 

mobilityLocalities <-mobilityDF %>% pull(CountryName) %>% unique() %>%  as.character() #127, redo post clean 
 
mobilityCatsFirst <- paste(mobilityCats, "First", sep="_")

mobilityDF[, mobilityCatsFirst] <- lapply(mobilityDF[, mobilityCatsFirst], shift, 1)   

mobilityDF %>% filter(!is.na(GroceryAndPharmacy_First)) %>% 
  group_by(CountryName) %>% 
  mutate_at(mobilityCatsFirst, as.numeric) %>% 
  mutate_at(mobilityCats, as.numeric) %>% 
  rowwise() %>% 
  mutate(meanActivityFirst=mean(GroceryAndPharmacy_First, Parks_First, 
                                Residential_First, RetailAndRecreation_First, 
                                TransitStations_First, Workplaces_First),
         meanActivityFirstNonRes=mean(GroceryAndPharmacy_First, Parks_First, 
                                      RetailAndRecreation_First, 
                                      TransitStations_First, Workplaces_First), 
         meanActivityRaw=mean(GroceryAndPharmacy, Parks, 
                              Residential, RetailAndRecreation, 
                              TransitStations, Workplaces),
         meanActivityRawNonRes=mean(GroceryAndPharmacy, Parks, 
                                    RetailAndRecreation, 
                                    TransitStations, Workplaces),
         StringencyIndex_FirstBin=as.factor(ifelse(StringencyIndexForDisplay_First>0, 1, 
                                                   0))) %>% 
  mutate(Month=month(Date) %>% as.factor) %>% 
  mutate(CountryName = factor(CountryName)) %>% 
  as.data.frame()-> MobilityModDF  #13355 obs 



## MCMCGLMM models 

library(MCMCglmm)
library(plotMCMC)

hist(MobilityModDF$meanActivityFirst) # gaussian 
hist(MobilityModDF$meanActivityRaw) # not gaussian but dunno how to transform 
hist(MobilityModDF$StringencyIndexForDisplay_First) #way zero-inflated 
hist(MobilityModDF$StringencyIndexForDisplay) #p awful 


Prior1 <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100)))


Prior2 <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100),
                        G2 = list(V = diag(1), nu = 0.002, alpha.mu = rep(0,1), alpha.V = diag(1)*100)))


mcmcMobility<- as.formula(paste("meanActivityFirstNonRes", "~", paste(c("StringencyIndexForDisplay", "Month"), collapse="+")))
mcmcStringency<- as.formula(paste("StringencyIndexForDisplay_First", "~", paste(c("meanActivityRaw", "Month"), collapse="+")))

mf<-50
MobModsFin <- list()
MobModsFin[[1]] <- MCMCglmm(fixed= mcmcMobility, 
                         random= ~  CountryName, 
                         prior=Prior1,
                         data=MobilityModDF,
                         family = "gaussian",
                         #pr=TRUE, 
                         nitt = 13000*mf,#REMEMBER YOU'VE DONE THIS
                         thin = 10*mf,burnin=3000*mf)

summary(MobModsFin[[1]])

MobModsFin[[2]] <- MCMCglmm(fixed= mcmcStringency, 
                         random= ~  CountryName, #+ meanActivityRaw:CountryName, 
                         prior=Prior1,
                         data=MobilityModDF, 
                         family = "gaussian",
                         #trunc=TRUE, 
                         #pr=TRUE, 
                         nitt = 13000*mf,
                         thin = 10*mf,burnin=3000*mf)

summary(MobModsFin[[2]])

saveRDS(MobModsFin, "ManuscriptOutput/MobilityRegressionsMS.rds")


## moodel outputs 

MobModOut <- list()
MobModOut[[1]] <- clean.MCMC.GLMM(MobModsFin[[1]]) %>% 
  rename("lowerCI95"=3, "upperCI95"=4) %>% 
  mutate(model="Mobility First Difference")
MobModOut[[2]] <- clean.MCMC.GLMM(MobModsFin[[2]]) %>% 
  rename("lowerCI95"=3, "upperCI95"=4) %>% 
  mutate(model="Stringency First Difference")

MobModOutDF <- bind_rows(MobModOut[1:2])
write.csv(MobModDF, "ManuscriptOutput/MobilityGLMMOutput.csv", row.names = F)

ModEfx<- MobModOutDF %>% 
  mutate(post.mean=as.numeric(post.mean), 
         lowerCI95=as.numeric(lowerCI95),
        upperCI95=as.numeric(upperCI95)) %>% 
  filter(variable %in% c("StringencyIndexForDisplay", "meanActivityRawNonRes")) %>% 
  ggplot(., aes(x=variable, y=post.mean, colour=model)) + 
  geom_point(size=4, position = position_dodge(width=0.5)) + 
  geom_hline(linetype="dashed", yintercept=0)+
  geom_errorbar(aes(ymax=upperCI95, ymin=lowerCI95),
                width=0.3, size=1.2, position = position_dodge(width=0.5)) + 
  scale_x_discrete(limits=c("StringencyIndexForDisplay", "meanActivityRawNonRes"), 
                   labels=c("Stringency Index", "Mean Activity")) + 
  scale_colour_manual(values=sparkPal[7:6], name="Model Response") + 
  theme_tufte(base_size = 14) + 
  labs(y="estimate", x="Fixed Effect") + 
  theme(text=element_text(family= "Helvetica"))  
  #ggtitle("Mixed Model Effects")

ModEfx <-addSmallLegend(ModEfx)


# Correlations  -----------------------------------------------------------

CorrDF <- 
  mobilityDF %>% 
  group_by(CountryName, CountryCode2) %>% 
  summarise(GroceryAndPharmacy=mean(GroceryAndPharmacy),
            Parks=mean(Parks),
            Residential=mean(Residential),
            RetailAndRecreation=mean(RetailAndRecreation),
            TransitStations=mean(TransitStations),
            StringencyIndex=mean(StringencyIndexForDisplay),
            Workplaces=mean(Workplaces)) %>%  
  rename("Grocery & Pharmacy"=GroceryAndPharmacy, "Retail & Recreation"=RetailAndRecreation, 
         "Transit Stations"= TransitStations, "Stringency Index"=StringencyIndex) %>% 
  mutate(CountryCode2=str_to_lower(CountryCode2))

mobilityCatClean<-  c("Grocery & Pharmacy", "Parks", "Retail & Recreation", 
                      "Transit Stations", "Workplaces", "Residential", "Stringency Index") 


CorrPlots <- list()

PredMod<- list()
PredPlots <- list()
PredOut<-list()

for(i in 1:6){
  PredMod[[i]] <- lm(CorrDF %>% pull(mobilityCatClean[i]) ~ CorrDF %>% pull(mobilityCatClean[7]))
  PredPlots[[i]] <- data.frame(pred=predict(PredMod[[i]], CorrDF), Stringency=CorrDF$`Stringency Index`)
  PredOut[[i]] <- data.frame(r2=round(summary(PredMod[[i]])$r.squared,2), 
                             estimate=round(summary(PredMod[[i]])$coefficients[2], 2),
                             p=round(anova(PredMod[[i]])$`Pr(>F)`[1],3),
                             Category=mobilityCatClean[i])
}

PredOutAll <- bind_rows(PredOut)
write.csv(PredOutAll, "ManuscriptOutput/AggregateMobilityRegressionsUpdatedMS.csv", row.names = FALSE)

CorrPlots <- lapply(1:6,
                    function(i){
                      require(DescTools)
                      yVal<-CorrDF %>% pull(mobilityCatClean[i]) %>% as.numeric() %>% RoundTo(5, floor)
                      xVal<-CorrDF %>% pull(mobilityCatClean[7]) %>% as.numeric() %>% RoundTo(10, floor)
                      CorrPlots[[i]] <- 
                        CorrDF %>% 
                        ggplot(aes(y=CorrDF %>% pull(mobilityCatClean[i]) %>% as.numeric(), 
                                   x=CorrDF %>% pull(mobilityCatClean[7]) %>% as.numeric())) +
                        geom_smooth(method="lm", se=FALSE, colour="darkgrey")+
                        #geom_abline(color="darkgrey", intercept=PredMod[[i]]$coefficients[1], 
                        #           slope=PredMod[[i]]$coefficients[2])+ 
                        geom_flag(aes(y=CorrDF %>% pull(mobilityCatClean[i]) %>% as.numeric(), 
                                      x=CorrDF %>% pull(mobilityCatClean[7]) %>% as.numeric(),  
                                      country=CorrDF %>% pull(CountryCode2)), size=3.7) + 
                        theme_classic(base_size = 14) + 
                        theme(text=element_text(family= "Helvetica")) + 
                        ggpubr::stat_cor(label.x = 15, label.y = min(yVal)) + 
                        scale_y_continuous(breaks = seq(min(yVal), max(yVal), 
                                                        by = ((max(yVal)-min(yVal))/length(unique(yVal))*2) %>% RoundTo(., 5))) +
                        scale_x_continuous(breaks = seq(min(xVal), max(xVal), by = 10)) + 
                        labs(y=mobilityCatClean[i], x=mobilityCatClean[7])
                      #ggtitle(mobilityLocalities[i])
                    }) 


library(patchwork)

CountryCorr <- (CorrPlots[[1]] + CorrPlots[[2]] + CorrPlots[[3]] )/ (CorrPlots[[4]] + CorrPlots[[5]]+CorrPlots[[6]])

CountryCorr + ggsave("ManuscriptFigures/CountryMeanCorr.png", units="mm", height=230, width=300, dpi=400)


# sparklines --------------------------------------------------------------

mobilityDF %>% #na.omit() %>% 
  select(Date, CountryName, CountryCode, mobilityCats) %>% 
  gather(key="MobilityCategory", value="DiffBaseline", mobilityCats) %>% 
  mutate(MobilityCategory=as.factor(MobilityCategory),
         MobilityCategory=fct_relevel(MobilityCategory,
                                      "GroceryAndPharmacy", "Parks", "Residential", 
                                      "RetailAndRecreation", "TransitStations", 
                                      "Workplaces", "StringencyIndexForDisplay")) %>% 
  filter(!is.na(MobilityCategory)) -> mobilityLongRaw



levels(mobilityLongRaw$MobilityCategory)
mobilityCatClean<-  c("Grocery & Pharmacy", "Parks", "Residential", "Retail & Recreation", 
                      "Transit Stations", "Workplaces", "Stringency Index") 

datesXminor<- mobilityLongRaw %>% filter(CountryName==mobilityLocalities[120]) %>% 
  pull(Date) %>% factor() %>% unique() %>% as.character() 

datesXmajor <- mobilityLongRaw %>% filter(CountryName==mobilityLocalities[120]) %>% 
  slice(1, 107) %>% # first and last of one country and category 
  pull(Date) %>% factor() %>% unique() %>% as.character() 

# sparkline all 
starts <- mobilityLongRaw %>% 
  mutate(Date=as.Date(Date), 
         DiffBaseline=round(DiffBaseline, 2)) %>% 
  group_by(CountryName, MobilityCategory) %>% filter(Date== min(Date))

ends <- mobilityLongRaw %>% 
  mutate(Date=as.Date(Date), 
         DiffBaseline=round(DiffBaseline, 2)) %>% 
  group_by(CountryName, MobilityCategory) %>% filter(Date== max(Date))

SparkPlots <- list()
SparkPlots<-lapply(1:length(mobilityLocalities), 
                   function(i){
                     require(ggplot2)
                     require(ggthemes)
                     require(grid)
                     SparkPlots[[i]] <-
                       mobilityLongRaw %>% filter(CountryName==mobilityLocalities[i]) %>% 
                       ggplot()  + 
                       geom_line(aes(x=as.Date(Date), y=DiffBaseline, colour=MobilityCategory, group=MobilityCategory), 
                                 size=0.8) + 
                       geom_point(data=starts %>% filter(CountryName==mobilityLocalities[i]), 
                                  aes(x=as.Date(Date), y=DiffBaseline, colour=MobilityCategory), size=2) + 
                       geom_point(data=ends %>% filter(CountryName==mobilityLocalities[i]), 
                                  aes(x=as.Date(Date), y=DiffBaseline, colour=MobilityCategory), size=2) + 
                       geom_text(data=starts %>% filter(CountryName==mobilityLocalities[i]), 
                                 aes(x=as.Date(Date), y=DiffBaseline, #colour=MobilityCategory, 
                                     label=DiffBaseline), hjust=1.2, vjust=0.8, size=3) + 
                       geom_text(data=ends %>% filter(CountryName==mobilityLocalities[i]), 
                                 aes(x=as.Date(Date), y=DiffBaseline, #colour=MobilityCategory,
                                     label=DiffBaseline), hjust=-0.3, size=3) + 
                       facet_grid(MobilityCategory ~ ., scales="free_y") + 
                       #scale_x_date(breaks=as.Date(datesXmajor), name="Date") + 
                       scale_y_continuous(expand = c(0.4, 0)) + 
                       scale_x_date(expand=c(0,12)) + 
                       scale_colour_manual(values=sparkPal, name="Category", labels=mobilityCatClean)+ 
                       theme_tufte(base_size = 12) + 
                       theme(text=element_text(family="Helvetica"), 
                             legend.position = "none", 
                             axis.title = element_blank(), 
                             axis.text.y = element_blank(),
                             strip.text = element_blank(), 
                             axis.ticks.y=element_blank(), 
                             axis.text.x = element_blank(), 
                             axis.ticks.x=element_blank()) + 
                       ggtitle(mobilityLocalities[i]) 
                   }) 


legend <- get_legend(SparkPlots[[1]]) # have to grab before removing legend above run with then without legends 
SparkPlots[[128]] <-legend

require(cowplot)
plot_grid(plotlist=SparkPlots[1:128], nrow=12, ncol=12) + 
  ggsave("ManuscriptFigures/SparkPlotsAll.png", units="mm", height=1000, width=1000, dpi=400)


#sparkline UK 
require(ggthemes)

SparkPlotUK  <-
  mobilityLongRaw %>% filter(CountryName==mobilityLocalities[120]) %>% 
  ggplot()  + 
  geom_line(aes(x=as.Date(Date), y=DiffBaseline, colour=MobilityCategory, group=MobilityCategory), 
            size=0.8) + 
  geom_point(data=starts %>% filter(CountryName==mobilityLocalities[120]), 
             aes(x=as.Date(Date), y=DiffBaseline, colour=MobilityCategory), size=2) + 
  geom_point(data=ends %>% filter(CountryName==mobilityLocalities[120]), 
             aes(x=as.Date(Date), y=DiffBaseline, colour=MobilityCategory), size=2) + 
  geom_text(data=starts %>% filter(CountryName==mobilityLocalities[120]), 
            aes(x=as.Date(Date), y=DiffBaseline, #colour=MobilityCategory, 
                label=DiffBaseline), hjust=1.5, vjust=0.8, size=3) + 
  geom_text(data=ends %>% filter(CountryName==mobilityLocalities[120]), 
            aes(x=as.Date(Date), y=DiffBaseline, #colour=MobilityCategory,
                label=DiffBaseline), hjust=-0.4, size=3) + 
  facet_grid(MobilityCategory ~ ., scales="free_y") + 
  scale_x_date(breaks=as.Date(datesXmajor), name="Date", expand=c(0,10)) + 
  scale_y_continuous(expand = c(0.4, 0)) + 
  scale_colour_manual(values=sparkPal, name="Category", labels=mobilityCatClean)+ 
  theme_tufte(base_size = 12) + 
  #ggtitle(paste(mobilityLocalities[120], "Mobility Trends", sep=" ")) + 
  theme(text=element_text(family="Helvetica"), 
        #legend.position = "none", 
        axis.title = element_blank(), 
        axis.text.y = element_blank(),
        strip.text = element_blank(), 
        axis.ticks.y=element_blank(), 
        title = element_text(hjust=5)) 

SparkPlotUK <-addSmallLegend(SparkPlotUK) 
SparkPlotUK + ggsave("ManuscriptFigures/SparkPlotUK.png", units="mm", height=300, width=250, dpi=300)


## ComboFig 
right_col <- SparkPlotUK/ModEfx + plot_layout(heights = c(1.2, 0.6)) 
ComboFig <- (CountryCorr + plot_layout(tag_level = 'new') | (right_col)) + plot_layout(widths = c(2, 1)) #+ plot_layout(tag_level = 'new')  
ComboFig + plot_annotation(tag_levels = c("A", "1")) & theme(plot.tag = element_text(size = 12)) + 
  ggsave("ManuscriptFigures/PanelPlot.png", units="mm", height=300, width=360, dpi=300)



# Narrative Plots  --------------------------------------------------------
library(gsheet)
library(ggrepel)
InterventionTimes <- gsheet2tbl("https://docs.google.com/spreadsheets/d/1LKTv_4-aC1I_Ld2HZAWs7yBkOb9fBELIV8bIEdDTykc/edit?usp=sharing") %>% 
  as.data.frame() %>% filter(Type!="") %>% 
  mutate(Date=as.Date(Date, format="%d/%m/%Y"),
         zoom=TRUE, 
         DateDescription=paste(EventDescription, " (", Date, ")", sep=""))

NarrativeDF <- mobilityLongRaw %>% 
  filter(CountryName %in% InterventionTimes$Country, MobilityCategory!="StringencyIndexForDisplay") %>% as.data.frame() %>% 
  mutate(CountryName=factor(CountryName), 
         MobilityCategory=fct_recode(MobilityCategory, 
                    `Grocery & Pharmacy`="GroceryAndPharmacy", 
                    `Retail & Recreation`="RetailAndRecreation", 
                    `Transit Stations`="TransitStations"))


starts <- mobilityLongRaw %>% 
  mutate(Date=as.Date(Date), 
         DiffBaseline=round(DiffBaseline, 2)) %>% 
  group_by(CountryName, MobilityCategory) %>% filter(Date== min(Date)) %>% as.data.frame() %>% 
  mutate(MobilityCategory=fct_recode(MobilityCategory, 
                                     `Grocery & Pharmacy`="GroceryAndPharmacy", 
                                     `Retail & Recreation`="RetailAndRecreation", 
                                     `Transit Stations`="TransitStations"))

ends <- mobilityLongRaw %>% 
  mutate(Date=as.Date(Date), 
         DiffBaseline=round(DiffBaseline, 2)) %>% 
  group_by(CountryName, MobilityCategory) %>% filter(Date== max(Date)) %>% as.data.frame() %>% 
  mutate(MobilityCategory=fct_recode(MobilityCategory, 
                                     `Grocery & Pharmacy`="GroceryAndPharmacy", 
                                     `Retail & Recreation`="RetailAndRecreation", 
                                     `Transit Stations`="TransitStations"))

MobilityCats <- levels(NarrativeDF$MobilityCategory)[1:6]
InterventionCats <- InterventionTimes$Type %>% unique()
NarrativeLocalities <- InterventionTimes$Country %>% unique()
shadesOfGrey <- colorRampPalette(c("grey0", "grey60"))
NarrativePal <- c(sparkPal[1:6], shadesOfGrey(5))
Colours <- c(MobilityCats, InterventionCats)


datesXminor<- mobilityLongRaw %>% filter(CountryName==mobilityLocalities[120]) %>% 
  pull(Date) %>% factor() %>% unique() %>% as.character() 

datesXmajor <- mobilityLongRaw %>% filter(CountryName==mobilityLocalities[120]) %>% 
  slice(1, 107) %>% 
  pull(Date) %>% factor() %>% unique() %>% as.character() 

datesZoomMajor <- InterventionTimes %>% filter(Country==NarrativeLocalities[11]) %>% 
  pull(Date) %>% factor() %>% unique() %>% as.character() 
datesZoomMinor <- datesXminor[16:47]
datesZoomMajor <- datesZoomMinor[c(2,30)]

NarrativeSparkUK <- 
  NarrativeDF %>% filter(CountryName==NarrativeLocalities[11])  %>% 
  ggplot()  + 
  geom_line(aes(x=as.Date(Date), y=DiffBaseline, group=MobilityCategory), 
            size=0.8, colour="grey40") + 
  geom_vline(data=InterventionTimes %>% filter(Country==NarrativeLocalities[11]), 
             aes(xintercept=as.Date(Date), colour=Type), size=0.6) + 
  geom_point(data=starts %>% filter(CountryName==NarrativeLocalities[11], MobilityCategory %in% MobilityCats), 
             aes(x=as.Date(Date), y=DiffBaseline), size=2) + 
  geom_point(data=ends %>% filter(CountryName==NarrativeLocalities[11], MobilityCategory %in% MobilityCats), 
             aes(x=as.Date(Date), y=DiffBaseline), size=2) + 
  geom_text(data=starts %>% filter(CountryName==NarrativeLocalities[11], MobilityCategory %in% MobilityCats), 
            aes(x=as.Date(Date), y=DiffBaseline, #colour=MobilityCategory, 
                label=DiffBaseline), hjust=1.5, vjust=0.8, size=3) + 
  geom_text(data=ends %>% filter(CountryName==NarrativeLocalities[11], MobilityCategory %in% MobilityCats), 
            aes(x=as.Date(Date), y=DiffBaseline, #colour=MobilityCategory,
                label=DiffBaseline), hjust=-0.4, size=3) + 
  annotate("rect", xmin = as.Date(datesZoomMajor[1]), xmax = as.Date(datesZoomMajor[2]), ymin = -Inf, ymax = Inf, alpha = 0.1) +
  scale_x_date(breaks=as.Date(datesXmajor), name="Date", expand=c(0,10)) + 
  #scale_y_continuous(expand = c(0.4, 0)) + 
  facet_grid(MobilityCategory ~ ., scales="free_y", switch = "y", labeller=label_wrap_gen(width=10)) + 
  scale_colour_manual(values=NarrativePal, name="Intervention Type", limits=Colours[7:11], 
                      labels=c(InterventionCats)) + 
  theme_tufte(base_size = 12) + 
  #ggtitle(paste(NarrativeLocalities[11], "Event Narrative", sep=" ")) + 
  theme(text=element_text(family="Helvetica"), 
        legend.position = "none", 
        legend.box="vertical", legend.margin=margin(), 
        strip.text.y = element_text(size=9),
        #axis.title = element_blank(), 
        axis.text.y = element_blank(),
        #strip.text = element_blank(), 
        axis.ticks.y=element_blank(), 
        title = element_text(hjust=5)) +
  guides(colour=guide_legend(nrow=2))
  

#NarrativeSparkUK <- addSmallLegend(NarrativeSparkUK) 

legendNarrative <- get_legend(NarrativeSparkUK) # have to grab before removing legend above 


yvals <- MobilityModDF  %>% 
  filter(CountryName==NarrativeLocalities[[11]]) %>% 
  mutate(Date=as.character(Date)) %>% 
  filter(Date %in% datesZoomMinor) %>% 
  pull(meanActivityRawNonRes)
xvals <- InterventionTimes %>% filter(Country==NarrativeLocalities[11]) %>% arrange(Date) %>%  pull(Date)

InterventionTimesWr <- InterventionTimes %>% 
  mutate(DateDescription= strwrap(DateDescription, width=150, simplify = FALSE), 
         DateDescription=sapply(DateDescription, paste, collapse="\n  "))

NarrativeZoomUK <- MobilityModDF %>% 
  mutate(Date=as.character(Date)) %>% 
  filter(CountryName==NarrativeLocalities[[11]], Date %in% datesZoomMinor) %>% 
  mutate(Date=as.Date(Date)) %>% 
  ggplot() + 
  geom_line(aes(x=as.Date(Date), y=meanActivityRawNonRes), 
            size=0.8, colour="grey40") + 
  geom_vline(data=InterventionTimes %>% filter(Country==NarrativeLocalities[11]), 
             aes(xintercept=as.Date(Date), colour=Type, label=DateDescription), size=1, alpha=0.7) + 
  scale_x_date(breaks=as.Date(datesZoomMajor), minor_breaks = as.Date(datesZoomMinor), name="Date", expand=c(0,0.5)) + 
  scale_y_continuous(expand = c(0.2, 0)) + 
  scale_colour_manual(values=NarrativePal, name="Intervention Type", limits=Colours[7:11], 
                      labels=c(InterventionCats))+ 
  theme_tufte(base_size = 12) + 
  annotate("text", x=as.Date(xvals[1:4]), y=c(rep(RoundTo(min(yvals), 5, floor), 4)), 
                label=c("1", "2", "3", "4,5"), hjust=1.2, vjust=0.8, size=4) + 
  #labs(y="Mean Activity, % Baseline") + 
  theme(text=element_text(family="Helvetica"), 
        legend.position = "none", 
        axis.title = element_blank(), 
        #axis.text.y = element_blank(),
        #strip.text = element_blank(), 
        #axis.ticks.y=element_blank(), 
        title = element_text(hjust=5), 
        plot.margin = unit(c(20,50,20,20), "pt")) # t, r, b, l  

AnnotationDF <- list()
for(x in 1:length(NarrativeLocalities)) { 
  df <- InterventionTimesWr %>% filter(Country == NarrativeLocalities[[x]]) %>% select(Type, Date, DateDescription)
  AnnotationDF[[x]] <- df %>% arrange(Date)
  }

AnnotationsUK <- list()
UKgrobs<- list()
for(x in 1:length(unique(InterventionTimes$Type))) {
  df <- AnnotationDF[[11]]
  AnnotationsUK[[x]] <- df %>% slice(x) %>% pull(DateDescription)
  UKgrobs[[x]] = textGrob(paste(x, AnnotationsUK[[x]], sep="-"), x=0.01, y=0.5, 
                          just = c("left", "top"), 
                          gp = gpar(fontsize = 10, col =  "black", fontfamily = "Helvetica"))
}

labsFootUK <- list()
for(x in 1:length(UKgrobs)){
labsFootUK[[x]] = gTree(paste("labsFoot",x), children = gList(UKgrobs[[x]]))
} 

footer = matrix(list(labsFootUK[[1]],labsFootUK[[2]], labsFootUK[[3]], labsFootUK[[4]], labsFootUK[[5]]), nrow=5)
footer_fin = gtable::gtable_matrix(name="Footer", 
                            grobs=footer, 
                            widths=unit(1, "null"), 
                            heights=
                              unit.c(unit(1.1, "grobheight", labsFootUK[[1]]) + unit(1.5, "lines"), 
                            unit(1.1, "grobheight", labsFootUK[[2]]) + unit(1.5, "lines"), 
                            unit(1.1, "grobheight", labsFootUK[[3]]) +  unit(1.5, "lines"),
                            unit(1.1, "grobheight", labsFootUK[[4]]) +  unit(1.5, "lines"),
                            unit(1.1, "grobheight", labsFootUK[[5]]) +  unit(1.5, "lines")))

#NarrativeZoomUKGrob <- as_grob(NarrativeZoomUK)
#grid.newpage()
#grid.draw(NarrativeZoomUKGrob)
#NarrativeZoomUKText <- gtable_add_grob(NarrativeZoomUKGrob, footer_fin, t=14, l=4, r=7)
#grid.draw(NarrativeZoomUKText)  

require(patchwork)
require(cowplot)

right_col <- grid.arrange(NarrativeZoomUK, footer_fin, legendNarrative, ncol=1, 
                          heights=c(2,1, 1)) %>% as_ggplot()
NarrativeSparkUK + right_col + 
  ggsave("ManuscriptFigures/SparkPlotsNarrativeUK.pdf", units="mm", height=300, width=390, dpi=400)


