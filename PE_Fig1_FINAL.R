require(tidyr)
require(scales)
require(ggplot2)
require(dplyr)
require(stringr)
require(nlme)
require(grid)
require(broom)
require(agricolae)
library(devtools)
require(DiagrammeR)


# ggplot theme----
theme_figs <- theme_bw() +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        axis.line.x = element_line(colour = 'grey60', size = 0.35), axis.line.y = element_line(colour = 'grey60', size = 0.35), 
        axis.ticks = element_line(colour = 'grey60', size = 0.35), 
        strip.background = element_blank(), strip.text = element_text(hjust = 0, vjust = 1, size = 11)) 
theme_set(theme_figs)

update_geom_defaults('point', c(size = 1.5))
update_geom_defaults('errorbar', c(size = 0.5))

# set up data----
#Metadata
Plant.endo.dat <- read.csv("~/Downloads/Endo_Data_Cleaned.csv")
Plant.endo.dat <-Plant.endo.dat[!is.na(Plant.endo.dat$Weight), ] #11 taken out, need to consider this

Plant.endo.dat.pool<-Plant.endo.dat %>%
  group_by(Pool.number) %>%
  summarise(Pool.Endo_Prev = mean(Endo_Prev), Pool.Weight=sum(Weight))
names(Plant.endo.dat.pool)<-c("Pool", "Endo_Prev", "Pool.Weight")

poolprev.dat <- read.csv("~/Downloads/FullPREV_2018_1.csv")  %>%
  filter(Dam == "D" | Dam == "ND") %>%
  filter(Distance == "12" | Distance == "24" | Distance == "36")

poolprev.dat
summary(poolprev.dat)
poolprev.dat$Pool<-as.factor(poolprev.dat$Pool)
poolprev.dat$fDistance<-as.factor(poolprev.dat$Distance)
poolprev.dat$Distance<-as.numeric(poolprev.dat$Distance)
poolprev.dat2<-poolprev.dat[poolprev.dat$Group!="Control",]
poolprev.dat3<-poolprev.dat2[poolprev.dat2$Pool!="4",] 
poolprev.dat3$DAE<-as.numeric(poolprev.dat3$DAE)

#AUDPS Analysis
poolprev.dat4a<-group_by(poolprev.dat3,Pool, Pool_Endo, DAE)%>%
  arrange(Pool, DAE)%>%
  summarize(Rhiz_Prev=as.numeric(mean(Rhiz)))
poolprev.dat4b<-group_by(poolprev.dat4a, Pool)%>%
  arrange(Pool)%>%
  summarize(MAXRhiz_Prev=max(Rhiz_Prev))
poolprev.dat4a<-as.data.frame(poolprev.dat4a)

poolprev.dat4b<-as.data.frame(poolprev.dat4b)

require(plyr)
poolprev.dat5<-join(poolprev.dat4a, poolprev.dat4b, by = "Pool")
Peakday<-poolprev.dat5[poolprev.dat5$Rhiz_Prev==poolprev.dat5$MAXRhiz_Prev,]
Peakday$PEAK<-Peakday$DAE
Peakday$DAE<-NULL
Peakday$Pool_Endo<-NULL
Peakday$Rhiz_Prev<-NULL
Peakday$MAXRhiz_Prev<-NULL
poolprev.dat5b<-join(poolprev.dat5, Peakday, by = "Pool")
poolprev.dat5c<-poolprev.dat5b[poolprev.dat5b$DAE<=poolprev.dat5b$PEAK,]

detach(package:plyr)

GCdat.audpcL <- group_by(poolprev.dat5c, Pool, Pool_Endo) %>% 
  arrange(Pool, DAE) %>% 
  summarize(audps = as.numeric(paste(as.numeric(audps(Rhiz_Prev, DAE, type='absolute')), collapse='-')), #added as.numeric to get rid of the following code
            rel.audps = as.numeric(paste(as.numeric(audps(Rhiz_Prev, DAE, type='relative')), collapse='-'))) %>% 
  ungroup() %>% na.omit()

GCdat.audpcL$audps = as.numeric(GCdat.audpcL$audps)
GCdat.audpcL$Pool_Endo = as.factor(GCdat.audpcL$Pool_Endo)
GCdat.audpcL$logaudps<-log(GCdat.audpcL$audps)

audps.modL<-lm(audps ~ Pool_Endo, data=GCdat.audpcL)
audps.modLlog<-lm(logaudps ~ Pool_Endo,  data=GCdat.audpcL)

resids.fig(audps.modL, GCdat.audpcL) #not normal
resids.fig(audps.modLlog, GCdat.audpcL) #better

anova(audps.modLlog) 
#p=0.09799

#VIOLIN PLOT OF AUDPC data by timing group
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
my_y_title <- expression(paste(italic("bacteria X"), " isolates with corresponding types"))
.... + labs(y=my_y_title)
pdf(file = "~/Desktop/PE_Fig1.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4,
    useDingbats = FALSE) # The height of the plot in inches
ggplot(GCdat.audpcL, aes(factor(Pool_Endo), logaudps))+theme_bw()+
  ylab("AUDPS (log-transformed)")+
  scale_fill_manual(values=cbPalette)+
  geom_dotplot(aes(fill=as.factor(Pool_Endo)), dotsize=0.75, binaxis='y', stackdir='center')+
  geom_boxplot(width=0.3, alpha=0)+
  scale_x_discrete(labels=c("Endo_Free" = "Endophyte \n Free", "Endo_Inoc" = "Endophyte \nInoculated"))+
  theme(axis.text=element_text(size=16),
        axis.title.y= element_text(size=16),
        axis.ticks=element_line(size=1),
        legend.title = element_text(size=16), 
        legend.text=element_text(size=16),
        legend.background = element_blank(),
        legend.position="none", axis.title.x=element_blank())
dev.off()
require(cowplot)
plot_grid(severity, prev, labels = c('A', 'B'))
