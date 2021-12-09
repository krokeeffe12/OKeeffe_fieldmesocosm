require(tidyr)
require(scales)
require(ggplot2)
require(dplyr)
require(stringr)
require(nlme)
require(grid)
require(broom)

# ggplot theme----
theme_figs <- theme_bw() +
  theme(axis.title = element_text(size = 11), axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.position = c(1, 1), legend.justification = c(1, 0.5),
        legend.key.size = unit(0.35, 'cm'),
         legend.text = element_text(size = 11),
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

#Prevalence Data
setwd("~/Documents/Mitchell Lab/Brandon/Summer 2018/Pool Experiment/Data")

poolprev.dat <- read.csv("~/Downloads/FullPREV_2018_1.csv")  %>%
  filter(Dam == "D" | Dam == "ND") %>%
  filter(Distance == "12" | Distance == "24" | Distance == "36")

poolprev.dat
summary(poolprev.dat)
poolprev.dat$Pool<-as.factor(poolprev.dat$Pool)
poolprev.dat$fDistance<-as.factor(poolprev.dat$Distance)
poolprev.dat$Distance<-as.numeric(poolprev.dat$Distance)
poolprev.dat2<-poolprev.dat[poolprev.dat$Group=="Exp",]
poolprev.dat3<-poolprev.dat2[poolprev.dat2$Pool!="4",] #contaminated

# restructure data
poolprev.long <- gather(poolprev.dat3, key = Parasite, value = Infected, Rhiz, Coll, Pucc, Sept) %>% 
  mutate(Parasite = factor(case_when(.$Parasite == "Coll" ~ "Colletotrichum",
                                     .$Parasite == "Rhiz"~ "Rhizoctonia",
                                     .$Parasite == "Pucc" ~ "Puccinia",
                                     .$Parasite == "Sept" ~ "Septoria"), levels = c("Colletotrichum", "Rhizoctonia", "Puccinia", "Septoria")))

poolprev.long2<-poolprev.long[poolprev.long$Parasite=="Rhizoctonia",]

#Pool_Endo only
PP.sum<-poolprev.long2 %>%
  group_by(Pool, DAE, Pool_Endo) %>%
  summarize(Rhiz.Prev = mean(Infected, na.rm = T)) %>%
  ungroup()

new2<-merge(x = Plant.endo.dat.pool, y = PP.sum, by="Pool", all.y=TRUE)
new2$Endo_Prev<-as.numeric(new2$Endo_Prev)

#Peak Prevalence (model-based)
new3<-new2[new2$DAE==25,]
library(ggplot2)

#Exploratory graphs
ggplot(new3, aes(Pool_Endo, Rhiz.Prev))+
  geom_boxplot()

ggplot(new3, aes(Pool.Weight, Rhiz.Prev))+
  geom_point()+
  theme(legend.position = "right")+
  geom_smooth(method="lm")

#Linear models
M0<-lm(Rhiz.Prev~1, data=new3) #null
M1<-lm(Rhiz.Prev~Pool_Endo*Pool.Weight, data=new3) #interactive
M1a<-lm(Rhiz.Prev~Pool.Weight+Pool_Endo, data=new3) #additive
M2<-lm(Rhiz.Prev~Pool_Endo, data=new3) #indiivdual
M3<-lm(Rhiz.Prev~Pool.Weight, data=new3) #indiivdual

AIC(M0, M1, M1a, M2, M3)
anova(M1, M1a, M2, M3) 

anova(M1) #interaction is insignificant

M1a.predict <- cbind(new3, predict(M1, interval = 'confidence'))

# Step 1: Call the pdf command to start the plot
pdf(file = "~/Desktop/PE_Fig5.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5,
    useDingbats = FALSE) # The height of the plot in inches

ggplot(M1a.predict, aes(Pool.Weight,Rhiz.Prev, group=Pool_Endo, color=Pool_Endo))+
  geom_point(size=3)+
  geom_line(aes(Pool.Weight, fit))+
  geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.2)+
  scale_color_manual(name="Endophyte \nInoculation", 
                    labels = c("Absent", 
                               "Present"), 
                    values = c("Endo_Free"="#E69F00", 
                               "Endo_Inoc"="#56B4E9"))+
  ylab("Peak Rhizoctonia Prevalence") +
  xlab("Aboveground Biomass (g)")+
  theme(axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16),
        axis.title.y= element_text(size=16),
        axis.ticks=element_line(size=1),
        legend.position = c(0.2, 0.85),
        legend.title = element_text(size=16), 
        legend.text=element_text(size=16),
        legend.background = element_blank(),
        axis.title.x=element_text(size=16))
  
dev.off()



