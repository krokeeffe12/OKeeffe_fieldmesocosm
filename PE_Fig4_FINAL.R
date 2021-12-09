require(tidyr)
require(scales)
require(ggplot2)
require(dplyr)
require(stringr)
library(agricolae)
require(nlme)
require(ggplot2)
require(grid)
require(broom)

# ggplot theme----
theme_figs <- theme_classic() +
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

## Load in data
Leaves<-read.csv(file="~/Downloads/Pool_Experiment_BW_Cum.csv", sep=",", header=T)
names(Leaves)
Leaves2<-Leaves[,1:10]
Leaves2<-Leaves2[Leaves2$Days_After_Infection==65,]

Weights<-read.csv(file="~/Downloads/Workbook1.csv", sep=",", header=T)
df_merged<-merge(x = Leaves2, y = Weights, by = "Plant_ID", all.x = TRUE)
df_merged$Total_Leaves<-as.numeric(df_merged$Total_Leaves)
df_merged<-df_merged[df_merged$Mass<0.3,]
df_test<-df_merged[!duplicated(df_merged),]
df_test<-df_test[!is.na(df_test$Mass), ]
df_test$Mass_g<-df_test$Mass*1000
df_test<-df_test 
df_test$fPool<-as.factor(df_test$Pool_Number.x)
df_pool<-df_test %>%
  dplyr::group_by(fPool)%>%
  summarize(pBiomass=sum(Mass_g), pLeaves=sum(Total_Leaves))

#df_pool$mean_PoolBiomass<-df_pool$pBiomass*13
df_pool$Pool_Number_num<-as.numeric(df_pool$fPool)

require(lmodel2)
Ex1.res <- lmodel2(log(pLeaves) ~ log(pBiomass), range.y="interval", range.x = "interval", data=df_pool, nperm=999)
print.lmodel2(Ex1.res)
plot(Ex1.res,  "RMA", xlab="Population-level aboveground biomass", main="",ylab="Population-level number of leaves")

reg <- Ex1.res$regression.results
names(reg) <- c("method", "intercept", "slope", "angle", "p-value")
reg

pdf(file = "~/Desktop/PE_Fig4.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4,
    useDingbats = FALSE) # The height of the plot in inches

ggplot(df_pool,aes(x=log(pBiomass),y=log(pLeaves))) + geom_point() + 
  geom_abline(intercept=Ex1.res$confidence.intervals[4,2],
              slope=Ex1.res$confidence.intervals[4,5],colour="grey", lty="dashed", size=1) + 
  geom_abline(intercept=Ex1.res$confidence.intervals[4,3],
              slope=Ex1.res$confidence.intervals[4,4],colour="grey", lty="dashed", size=1) +
  geom_abline(intercept=Ex1.res$regression.results[4,2],
              slope=Ex1.res$regression.results[4,3],colour="black", size=1) +
  theme_bw()+
  ylim(5.7,6.6)+
  xlab("Population-level Aboveground Biomass\n(log-transformed)")+
  ylab("Population-level Number of Leaves\n(log-transformed)")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))

dev.off()


