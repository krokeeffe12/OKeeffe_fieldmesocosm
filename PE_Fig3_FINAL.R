library(ggplot2)
library(beeswarm)
library(dplyr)
library(lme4)
library(nlme)
require(ggpubr)

#Set theme
theme_oeco <- theme_bw() +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), 
        axis.line.x = element_line(size = 0.35, colour = 'grey50'), axis.line.y = element_line(size = 0.35, colour = 'black'),
        axis.ticks = element_line(size = 0.5, colour = 'grey50'), legend.text=element_text(size=12),
        axis.title.y = element_text(margin = margin(0, 20, 0, 0)))
theme_set(theme_oeco)

#Upload weight data
pool.weights <- read.csv("~/Downloads/Endo_Data_Cleaned2.csv")
pool.weights<-pool.weights[!is.na(pool.weights$Weight), ]

#Upload prevalence data
poolprev.dat <- read.csv("~/Downloads/PoolPREV_1_All.csv")  %>%
  filter(Dam == "D" | Dam == "ND") 

new2<-merge(x = pool.weights, y = poolprev.dat, by="Plant", all.x=TRUE)
new_sub<-new2[1:23]
new_sub$DAE<-NULL
new_sub$Survey_Date<-NULL
new_sub2<-unique(new_sub) %>%
  filter(Control=="Experimental") #Remove control pools

head(new_sub2)
pool.dat<-group_by(new_sub2, Pool.number, Endo_Prev) %>% 
  arrange(Pool.number) %>% 
  summarize(Pool.Weight = sum(Weight))%>% 
ungroup() %>% na.omit()

fit0<-gls(Pool.Weight~1, method="ML", data=pool.dat)
fit1<-lm(Pool.Weight~Endo_Prev,
          data=pool.dat)

AIC(fit0, fit1)
anova(fit0, fit1) #Endophyte group is a bit better than null model
anova(fit1) #Endophyte group significant predictor
plot(fit1) #looks good

# Step 1: Call the pdf command to start the plot
pdf(file = "~/Desktop/PE_Fig3.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5,
    useDingbats = FALSE) # The height of the plot in inches

# Step 2: Create the plot with R code
ggplot(pool.dat, aes(x=Endo_Prev, y=Pool.Weight))+
  ylim(0,175)+
  stat_smooth(method="lm", color="black")+
  geom_point(size=2)+
  xlab("Confirmed Endophyte Infection Prevalence")+
  ylab("Population-level\nAboveground Biomass (g)")+
  theme(axis.title=element_text(size=14, face="bold"),
        axis.text=element_text(size=14))

# Step 3: Run dev.off() to create the file!
dev.off()





