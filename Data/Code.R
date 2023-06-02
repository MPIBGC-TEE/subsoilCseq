rm(list=ls())

if (dir.exists("F:/OneDrive - mails.ucas.edu.cn/AA_Research/AB_Wang_Luo_Age")){
  setwd("F:/OneDrive - mails.ucas.edu.cn/AA_Research/AB_Wang_Luo_Age")
} else {
  setwd("C:/Users/86136/OneDrive - mails.ucas.edu.cn/AA_Research/AB_Wang_Luo_Age")
}

library(tidyr)
library(ggplot2)
library(plotly)
library(dplyr)
library(RColorBrewer)
library("scales")


###
### Fig.a: Mean age and transit times
###
df <- read.csv("./Data_CS/SoilRadioCarbon_Profile_Mean_Age_Tau_Fivelayers.cs.csv")

MeanAge = c(df$MeanAge20cm,df$MeanAge40cm,df$MeanAge60cm,df$MeanAge80cm,df$MeanAge100cm)

MeanTau = c(df$TransitTime_20,df$TransitTime_40,df$TransitTime_60,df$TransitTime_80,df$TransitTime_100)

Depth = c(rep(20,910),rep(40,910),rep(60,910),rep(80,910),rep(100,910))

df.m <- data.frame(Value = c(MeanAge,MeanTau),
                   Variable = c(rep("Age",length(MeanAge)),rep("Tau",length(MeanTau))),
                   Depth = c(Depth,Depth),
                   Biome = c(rep(df$Biome,10),rep(df$Biome,10)))
library(plyr)

adata <- ddply(df.m,c("Depth","Variable"),
               summarise,
               mean=mean(Value,na.rm=T),
               sd = sd(Value,na.rm=T),
               se=sd/sqrt(sum(!is.na(Value)))) %>%
  mutate(Variable = factor(Variable, levels=c("Tau","Age")))

write.csv(adata,"./Data_CS/Data_figa.csv",row.names = F)
################################################

Fig.a <- ggplot()+
  geom_point(data=adata,aes(x=Depth, y=mean, color=Variable,group=Variable),shape=1,size=2)+
  geom_smooth(data=adata,aes(x=Depth, y=mean, color=Variable),show.legend = F, size=0.1)+
  geom_errorbar(data=adata,aes(ymin = mean - se, ymax = mean + se,x=Depth, color=Variable,group=Variable), width=2,
                position=position_dodge(0.05)) +
  scale_y_continuous(position = "right",limits = c(0,5000)) +
  scale_x_reverse(labels=c("20" = "0-20", "40" = "20-40","60" = "40-60", "80"="60-80","100"="80-100")) +
  theme_classic()+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = c(0.8, 0.8),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        axis.ticks.length=unit(.25, "cm"),
        axis.text.y = element_text(angle = 45))+ 
  labs(y = "Mean soil carbon age and transit time (yr)")+
  labs(x = "Depth (cm)") +
  coord_flip() +
  scale_fill_manual(values=c("aquamarine4","chocolate2"))+
  scale_color_manual(values=c("aquamarine4","chocolate2")) 


###
### Fig.b: Mean transit time by biomes
###

bdata <- df.m[df.m$Variable=="Tau",]

bdata <- ddply(bdata, c("Depth","Biome"),
               summarise,
               mean=mean(Value,na.rm=T),
               sd = sd(Value,na.rm=T),
               se=sd/sqrt(sum(!is.na(Value)))) %>%
  mutate(Biome = factor(Biome, levels=c("Tropical/Subtropical forests","Tropical/Subtropical grasslands/savannas","Temperate forests",
                                        "Temperate grasslands","Mediterranean/Montane shrublands","Boreal forests","Tundra","Deserts",
                                        "Crops")))

write.csv(bdata,"./Data_CS/Data_figb.csv",row.names = F)

col.figb <- c("darkgreen","darkolivegreen1","olivedrab","yellowgreen","palegreen",
              "darkslategrey","cadetblue","moccasin","gold")


Fig.b <- ggplot()+
  geom_point(data=bdata,aes(x=Depth, y=mean, color=Biome,group=Biome),shape=1,size=2)+
  geom_smooth(data=bdata,aes(x=Depth, y=mean, color=Biome),show.legend = F,size=0.1)+
  geom_errorbar(data=bdata,aes(ymin = mean - se, ymax = mean + se,x=Depth, color=Biome,group=Biome), width=10,
                position=position_dodge(0.05)) +
  scale_y_continuous(position = "right",trans='log10',breaks=c(1,10,100,1000,10000)) +
  scale_x_reverse(labels=c("20" = "0-20", "40" = "20-40","60" = "40-60", "80"="60-80","100"="80-100")) +
  theme_classic()+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = c(0.8, 0.8),
        legend.key.size = unit(0.2, 'cm'),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        axis.ticks.length=unit(.25, "cm"),
        axis.text.y = element_text(angle = 45))+ 
  labs(y = "Mean transit time of soil carbon input (yr)")+
  labs(x = "Depth (cm)") +
  coord_flip() +
  scale_fill_manual(values=col.figb)+
  scale_color_manual(values=col.figb) 

#plot, further modified manually in PDF
library(cowplot)
windows(h=4,w=9)
par(mgp = c(2.5,0.8,0),mar=c(4,4.5,2,1))
plot_grid(Fig.a, Fig.b, ncol=2,align = "vh")

