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
Plant.endo.dat <-Plant.endo.dat[!is.na(Plant.endo.dat$Weight), ] 
Plant.endo.dat$fPool.number<-as.factor(Plant.endo.dat$Pool.number)

Plant.endo.dat.pool<-Plant.endo.dat %>% dplyr::group_by(fPool.number) %>%
  dplyr::summarise(Pool.Endo_Prev = mean(Endo_Prev), Pool.Weight=sum(Weight))
names(Plant.endo.dat.pool)<-c("Pool", "Endo_Prev", "Pool.Weight")
#Prevalence Data

poolprev.datA <- read.csv("~/Downloads/FullPREV_2018_1.csv")  %>%
  filter(Dam == "D" | Dam == "ND") %>%
  filter(Distance == "12" | Distance == "24" | Distance == "36")

poolprev.datA
summary(poolprev.datA)
poolprev.datA$Pool<-as.factor(poolprev.datA$Pool)
poolprev.datA$fDistance<-as.factor(poolprev.datA$Distance)
poolprev.datA$Distance<-as.numeric(poolprev.datA$Distance)
poolprev.datA2<-poolprev.datA[poolprev.datA$Group=="Exp",]
poolprev.datA3<-poolprev.datA2[poolprev.datA2$Pool!="4",] 

# restructure data
poolprev.long <- gather(poolprev.datA3, key = Parasite, value = Infected, Rhiz, Coll, Pucc, Sept) %>% 
  mutate(Parasite = factor(case_when(.$Parasite == "Coll" ~ "Colletotrichum",
                                     .$Parasite == "Rhiz"~ "Rhizoctonia",
                                     .$Parasite == "Pucc" ~ "Puccinia",
                                     .$Parasite == "Sept" ~ "Septoria"), levels = c("Colletotrichum", "Rhizoctonia", "Puccinia", "Septoria")))

#Restructure again to summarize the pool data 
poolprev.long2<-poolprev.long[poolprev.long$Parasite=="Rhizoctonia",]

#Pool_Endo only
PP.sum<-poolprev.long2 %>%
  dplyr::group_by(Pool, DAE, Pool_Endo) %>%
  dplyr::summarize(Rhiz.Prev = mean(Infected, na.rm = T)) %>%
  ungroup()

#Combine prevalence and endophyte data
new2<-merge(x = Plant.endo.dat.pool, y = PP.sum, by="Pool", all.y=TRUE)
new2$Endo_Prev<-as.numeric(new2$Endo_Prev)
new2<-new2[new2$DAE<28,]
#write.csv(new2, "2018.csv")
new3<-read.csv("~/Desktop/2018.csv", header=TRUE)
#Figure out random structure first
#null model
Rmod.null <-gls(Rhiz.Prev~1, method="REML", data=new2)
#Add in fixed effects
Rmod.gls <-gls(Rhiz.Prev~Pool_Endo*DAE, method="REML", data=new2)
#Random intercepts
Rmod.ri<-lme(Rhiz.Prev~Pool_Endo*DAE, 
             random=~1|Pool, method="REML", 
             data=new2)
#Random slopes
Rmod.slopes<-lme(Rhiz.Prev~Pool_Endo*DAE, 
                 random=~DAE|Pool, method="REML", 
                 control=list(opt="optim",maxIter=5000000, msMaxIter=500000), #let it run longer than the defaults)
                 data=new2)
#Random slopes with interaction term in the random structure
Rmod.slopesb<-lme(Rhiz.Prev~Pool_Endo*DAE, 
                 random=~Pool_Endo:DAE|Pool, method="REML", 
                 control=list(opt="optim",maxIter=5000000, msMaxIter=500000), #let it run longer than the defaults)
                 data=new2)

AIC(Rmod.null, Rmod.gls, Rmod.ri, Rmod.slopes, Rmod.slopesb)
anova(Rmod.null, Rmod.gls, Rmod.ri, Rmod.slopes, Rmod.slopesb)
#random intercepts has lowest AIC, but not significant with random slopes
#random slopes takes into account repeated measures so keeping
#keeping interaction because best random structure
resids.fig(Rmod.slopesb, new2) 

ctrl <- lmeControl(opt='optim')

#Now go to ML
Rmod.slopes.ML<-lme(Rhiz.Prev~Pool_Endo*DAE, 
                  random=~Pool_Endo:DAE|Pool, method="ML", 
                  control=list(opt="optim",maxIter=5000000, msMaxIter=500000), #let it run longer than the defaults)
                  data=new2)

#Add in temporal autocorrelation structure
Rmod.slopes.ML.corr<-lme(Rhiz.Prev~Pool_Endo*DAE, 
             random=~Pool_Endo:DAE|Pool, method="ML", 
             correlation=corCAR1(form=~DAE|Pool),
             control=list(opt="optim",maxIter=5000000, msMaxIter=500000),
             data=new2)

anova(Rmod.slopes.ML, Rmod.slopes.ML.corr) #autocorrelation structure is better

#What polynomial for time
Rmod.slopes1<-lme(Rhiz.Prev~Pool_Endo*poly(DAE, 1, raw=TRUE), 
                   random=~Pool_Endo:DAE|Pool, method="ML", 
                   correlation=corCAR1(form=~DAE|Pool),control=list(opt="optim",maxIter=5000000, msMaxIter=500000),
                   data=new2)

Rmod.slopes2b<-lme(Rhiz.Prev~Pool_Endo*poly(DAE, 2, raw=TRUE), 
                  random=~Pool_Endo:DAE|Pool, method="ML", 
                  correlation=corCAR1(form=~DAE|Pool),control=list(opt="optim",maxIter=5000000, msMaxIter=500000),
                  data=new2)

Rmod.slopes3<-lme(Rhiz.Prev~Pool_Endo*poly(DAE, 3, raw=TRUE), 
              random=~DAE|Pool, method="ML", 
              correlation=corCAR1(form=~DAE|Pool),control=list(opt="optim",maxIter=5000000, msMaxIter=500000),
              data=new3)

Rmod.slopes4<-lme(Rhiz.Prev~Pool_Endo*poly(DAE, 4, raw=TRUE), 
                  random=~DAE|Pool, method="ML", 
                  correlation=corCAR1(form=~DAE|Pool),control=list(opt="optim",maxIter=5000000, msMaxIter=500000),
                  data=new2)

AIC(Rmod.slopes1, Rmod.slopes2b, Rmod.slopes3, Rmod.slopes4)
anova(Rmod.slopes1, Rmod.slopes2b, Rmod.slopes3, Rmod.slopes4) #correlation structure is significantly better

Rmod.slopes3L<-lme(Rhiz.Prev~Pool_Endo*poly(DAE, 3, raw=TRUE), 
                  random=~Pool_Endo:DAE|Pool, method="ML", 
                  control=list(opt="optim",maxIter=5000000, msMaxIter=500000),
                  data=new3)


AIC(Rmod.slopes3L, Rmod.slopes3)
anova(Rmod.slopes3L, Rmod.slopes3) #correlation structure is significantly better

cbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
Prev<-Rmod.slopes3 %>% broom.mixed::augment() %>% 
  ggplot(aes(x=DAE,y=Rhiz.Prev))+
  theme_bw() + 
  #geom_line(aes(color=Pool_Endo, group=Pool),size=1, alpha=0.2)+
  #geom_line(data =  PP.sum, aes(x=DAE, y=Rhiz.Prev, group = Plant, color=Pool_Endo, alpha=0.00001))+
  #geom_smooth(aes(y=.fixed),size=2, method = "lm", se = F)+ #this is what you want to change.
  geom_point(aes(color=Pool_Endo),size=3, position=position_jitter(width=0.5), alpha=0.3)+
  geom_line(aes(y=.fixed, color=Pool_Endo),size=3)+
  #geom_line(data=df2017A, aes(y=.fixed, color=Pool_Endo),size=3, lty="dashed")+
  scale_color_manual(values=cbPalette)+
  ylab("Proportion of Leaves \nExhibiting Rhizoctonia Lesions") +
  xlab("Days after Rhizoctonia Inoculation")+
  guides(color=guide_legend(title="Endophyte Inoculation"))+
  theme(axis.text=element_text(size=16),
        axis.title= element_text(size=16),
        axis.ticks=element_line(size=2),
        legend.title = element_text(size=16),
        legend.position=c(0.27, 0.83),
        legend.text=element_text(size=16),
        legend.background = element_blank())+
  scale_color_manual(name="Endophyte Inoculation", 
                          labels = c("Absent", 
                                     "Present"), 
                          values = c("Endo_Free"="#E69F00", 
                                     "Endo_Inoc"="#56B4E9"))

#Maximum prevalence
new4<-new2 %>%
  dplyr::group_by(Pool, Pool_Endo) %>%
  dplyr::summarize(MaxRhiz.Prev = max(Rhiz.Prev, na.rm = T)) %>%
  ungroup()

t.test(MaxRhiz.Prev~Pool_Endo, data=new4)

require(plotrix)
myData <- aggregate(new4$MaxRhiz.Prev,
                    by = list(cyl = new4$Pool_Endo),
                    FUN = function(x) c(mean = mean(x), sd = sd(x),
                                        n = length(x)))
myData <- do.call(data.frame, myData)
myData$se <- myData$x.sd / sqrt(myData$x.n)
colnames(myData) <- c("Pool_Endo", "mean", "sd", "n", "se")
myData$names <- c(paste(myData$Pool_Endo, "Endophyte Inoculation"))

Max<-ggplot(myData)+
  geom_bar(aes(x=Pool_Endo, y=mean, fill=Pool_Endo), width=0.6, stat="identity", alpha=0.7)+
  geom_errorbar(aes(x=Pool_Endo, ymin=mean-se, ymax=mean+se), alpha=0.9, size=1.1, width=0.3)+
  scale_fill_manual(name="Endophyte Inoculation", 
                     labels = c("Absent", 
                                "Present"), 
                     values = c("Endo_Free"="#E69F00", 
                                "Endo_Inoc"="#56B4E9"))+
  ylab("Peak Rhizoctonia Prevalence") +
  xlab("Endophyte Inoculation ")+
  theme_bw()+
  theme(legend.position="none")+
  theme(axis.text=element_text(size=16),
        axis.title= element_text(size=16),
        axis.ticks=element_line(size=2),
        #legend.title = element_text(size=16),
        legend.position="none")+
        #legend.text=element_text(size=16),
        #legend.background = element_blank())+
  scale_x_discrete(breaks=c("Endo_Free","Endo_Inoc"),
                   labels=c("Absent", "Present"))


require(cowplot)
plot_grid(Prev, Max, labels = c('A', 'B'))

resids.fig <- function(mod, df) {
  residdf <- dplyr::mutate(df, resids = residuals(mod, type = 'normalized'),
                           fits = fitted(mod))
  fig2 <-ggplot(residdf, aes(x = fits, y = resids)) + geom_point() +
    labs(x = 'Fitted values', y = '')
  
  fig3 <- ggplot(residdf) + stat_qq(aes(sample = resids)) +
    labs(x = 'Theoretical Quantiles', y = 'Sample Quantiles')
  
  # qqline plot = FALSE, according to James should work
  
  fig4 <- ggplot(residdf, aes(x = resids)) + geom_histogram(aes(y=..density..), colour = 'grey50') +
    labs(x = 'Residuals', y = 'Frequency') + scale_y_continuous(expand = c(0, 0)) +
    stat_function(fun = dnorm, color = "red", args = list(mean = mean(residdf$resids),
                                                          sd = sd(residdf$resids)))
  grid.draw(rbind(ggplotGrob(fig2), ggplotGrob(fig3), ggplotGrob(fig4), size = 'first'))
  
  return(summary(mod))
}

