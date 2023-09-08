rm(list=ls())

# load packages
library(tidyr)
library(dplyr)
library(plyr)
library(ggplot2)
library(viridis)
library(cowplot)
library(reshape2)
library(randomForest)
library(edarf)


#### initial processing ####
load("dfNa.RData") # no management 4 degree increase simulation
load("dfNb.RData") # no management 4 degree increase simulation
load("dfNc.RData") # no management 4 degree increase simulation
dfN <- rbind(dfNa,dfNb,dfNc) # full dataframe
rm(dfNa,dfNb,dfNc) # remove the pieces

load("dfHa.RData") # heterogeneity management 4 degree increase simulation
load("dfHb.RData") # heterogeneity management 4 degree increase simulation
load("dfHc.RData") # heterogeneity management 4 degree increase simulation
dfH <- rbind(dfHa,dfHb,dfHc) # full dataframe
rm(dfHa,dfHb,dfHc) # remove the pieces

load("dfEa.RData") # stochasticity management 4 degree increase simulation
load("dfEb.RData") # stochasticity management 4 degree increase simulation
load("dfEc.RData") # stochasticity management 4 degree increase simulation
dfE <- rbind(dfEa,dfEb,dfEc) # full dataframe
rm(dfEa,dfEb,dfEc) # remove the pieces

# find the log10 of the dispersal because it is a little easier for some plots
dfN$ldisp <- log10(dfN$disp)
dfH$ldisp <- log10(dfH$disp)
dfE$ldisp <- log10(dfE$disp)

# simplify dfH and dfE
dfH$NfNA<-rep(dfH$Nf[dfH$h1==dfH$h2],each=2) # adds another column that is the final population size for simulation with no action
dfH <- dfH[dfH$h1!=dfH$h2,] # remove the rows with no management
dfE$NfNA<-rep(dfE$Nf[dfE$e1==dfE$e2],each=2) # adds another column that is the final population size for simulation with no action
dfE <- dfE[dfE$e1!=dfE$e2,] # remove the rows with no management


#### random forests plots (figures 2, S2-S4) ####
#### abundance regression ####
# new dataframes for rf

# create new dataframe for relative change regression
dfN99 <- dfN %>%
  mutate(NRatio = Nf/Ni)

df99 <- df99 %>%
  subset(Ni>0 & NRatio<quantile(NRatio,0.99))   # only keep populations that begin positive
  dfN2[dfN2$NRatio<quantile(dfN2$NRatio,0.99),] # remove outlier NRatios for regression


# set the seed and sample for the random forest
set.seed(1)
sampN <- 1:100000

# random forest regression
rfN<-randomForest(NRatio~ldisp+toler+r+E+H,dfN99[sampN,],importance=T,do.trace=T)

# find the importance and put in dataframe
rfImpsN <- as.data.frame(importance(rfN))
rfImpsN$variableName <- c('gamma','sigma','rho','S','H')

# reorganize it into a better dataframe
dfVarImpN<-data.frame(variable=rownames(rfImpsN),
                      variableName=rfImpsN[,3],
                      IncMSE=rfImpsN[,1],
                      IncNodePurity=rfImpsN[,2],
                      rank=rank(-rfImpsN[,1]),
                      order=c(1,2,3,5,4))

# set seed for partial dependence
set.seed(2)
# calculate partial dependence with sample+100000
partDepN <- partial_dependence(rfN,c('ldisp','toler','r','E','H'),data=dfN99[100000+sampN,])

# final dataframe for plotting
pdpN <- partDepN %>%
  pivot_longer(cols = -NRatio, names_to = "variable", values_to = "value") %>%
  drop_na(value) %>%
  dplyr::select(variable, value, NRatio)

# Plot the data
mseN <- ggplot(dfVarImpN,aes(x=order,y=IncMSE))+
  geom_point(size=3)+
  geom_text(aes(label=variableName),parse=T,vjust = 0, nudge_y = 25)+
  geom_segment(linetype=2,aes(xend=order,y=0,yend=IncMSE))+
  scale_x_continuous(name="")+
  scale_y_continuous(name="Variable importance\n(percent increase in\nmean squared error)",
                     expand=c(0,0),limits=c(0,max(dfVarImpN$IncMSE)*1.2))+
  theme_classic()+
  theme(legend.position = "none",axis.text.x = element_blank())

purityN <- ggplot(dfVarImpN,aes(x=order,y=IncNodePurity))+
  geom_point(size=3)+
  geom_text(aes(label=variableName),parse=T,vjust = 0, nudge_y = 50)+
  geom_segment(linetype=2,aes(xend=order,y=0,yend=IncNodePurity))+
  scale_x_continuous(name="")+
  scale_y_continuous(name="Variable importance\n(increase in node purity)",
                     expand=c(0,0),limits=c(0,max(dfVarImpN$IncNodePurity)*1.2))+
  theme_classic()+
  theme(legend.position = "none",axis.text.x = element_blank())

dispN <- ggplot(subset(pdpN,variable=='ldisp'),aes(x=value,y=NRatio))+
  geom_line(linewidth=1.5)+
  scale_x_continuous(name=(expression(gamma: " Mean dispersal distance")),
                     breaks=c(-4,-2,0),labels=c("0.0001","0.01","1"),expand=c(0,0),
                     limits=c(-4,0))+
  scale_y_continuous(name='Relative final\npopulation size',limits=c(0,1),expand=c(0,0))+
  theme_classic()

tolerN <- ggplot(subset(pdpN,variable=='toler'),aes(x=value,y=NRatio))+
  geom_line(linewidth=1.5)+
  scale_x_continuous(name=(expression(sigma: " Thermal tolerance breadth")),expand=c(0,0),
                     limits=c(2,5),breaks=c(2,3,4,5))+
  scale_y_continuous(name='Relative final\npopulation size',limits=c(0,1),expand=c(0,0))+
  theme_classic()

fecundN <- ggplot(subset(pdpN,variable=='r'),aes(x=value,y=NRatio))+
  geom_line(linewidth=1.5)+
  scale_x_continuous(name=(expression(rho: " Fecundity")),expand=c(0,0),
                     limits=c(1.5,6),breaks=c(2,4,6))+
  scale_y_continuous(name='Relative final\npopulation size',limits=c(0,1),expand=c(0,0))+
  theme_classic()

heterogN <- ggplot(subset(pdpN,variable=='H'),aes(x=value,y=NRatio))+
  geom_line(linewidth=1.5)+
  scale_x_continuous(name=(expression(H: " Local temperature heterogeneity")),
                     limits=c(0,2),breaks=c(0,1,2),expand=c(0,0))+
  scale_y_continuous(name='Relative final\npopulation size',limits=c(0,1),expand=c(0,0))+
  theme_classic()

stochastN <- ggplot(subset(pdpN,variable=='E'),aes(x=value,y=NRatio))+
  geom_line(linewidth=1.5)+
  scale_x_continuous(name=(expression(S: " Environmental stochasticity")),expand=c(0,0),
                     limits=c(0,1),breaks=c(0,0.5,1))+
  scale_y_continuous(name='Relative final\npopulation size',limits=c(0,1),expand=c(0,0))+
  theme_classic()

hetNL <- ggplot(sdfH,aes(x= as.numeric(dispCat2),y= as.numeric(hetCat2),fill=NfM/NiM))+
  geom_tile()+
  scale_fill_viridis_c(name="Relative final\npopulation size",limits=c(0,1),
                       guide = guide_colorbar(direction = "horizontal"))+
  scale_x_continuous(name="Mean dispersal distance",breaks=c(0.5,8.5,16.5),labels=c("0.0001","0.01","1"),expand=c(0,0))+
  scale_y_continuous(name=(expression(H: " Local temperature heterogeneity")),
                     breaks=c(0.5,8.5,16.5),labels=c(0,1,2),expand=c(0,0))+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5,direction='horizontal')) +
  coord_fixed()

legendN<-get_legend(hetNL)

hetN <- ggplot(sdfH,aes(x= as.numeric(dispCat2),y= as.numeric(hetCat2),fill=NfM/NiM))+
  geom_tile()+
  scale_fill_viridis_c(name="Average\npercent\npopulation\nremaining",limits=c(0,1),guide=F)+
  scale_x_continuous(name="Mean dispersal distance",breaks=c(0.5,8.5,16.5),labels=c("0.0001","0.01","1"),expand=c(0,0))+
  scale_y_continuous(name=(expression(H: " Local temperature heterogeneity")),
                     breaks=c(0.5,8.5,16.5),labels=c(0,1,2),expand=c(0,0))+
  coord_fixed()+ theme(plot.margin = unit(c(0,1,1,1), "lines"))
stoN <- ggplot(sdfE,aes(x= as.numeric(dispCat2),y= as.numeric(stoCat2),fill=NfM/NiM))+
  geom_tile()+
  scale_fill_viridis_c(name="Average\npercent\npopulation\nremaining",limits=c(0,1),guide=F)+
  scale_x_continuous(name="Mean dispersal distance",breaks=c(0.5,8.5,16.5),labels=c("0.0001","0.01","1"),expand=c(0,0))+
  scale_y_continuous(name=(expression(S: " Environmental stochasticity")),
                     breaks=c(0.5,8.5,16.5),labels=c(0,0.5,1),expand=c(0,0))+
  coord_fixed()+ theme(plot.margin = unit(c(0,1,1,1), "lines"))


# final version of figure 2 (or S2)
f2N <- plot_grid(plot_grid(mseN,dispN,tolerN,fecundN,heterogN,stochastN,ncol=2,
                           labels=c('a','b','c','d','e','f')),
                 plot_grid(legendN,
                           plot_grid(hetN,stoN,labels=c('h','g'),align='v',ncol=1),
                           ncol=1,rel_heights = c(1,4)),
                 rel_widths = c(2.2,1),ncol=2)

#### extinction classification ####

# sample and seed for persistence classification
set.seed(3)
samps<-10000

# random forest classification
rfX<-randomForest(factor(ext)~ldisp+toler+r+E+H,dfN,importance=T,sampsize=c(samps,samps),do.trace=T)

# find the importance and put in dataframe
rfImpsX <- as.data.frame(rfX$importance)
rfImpsX$variableName <- c('gamma','sigma','rho','S','H')

dfVarImpX<-data.frame(variable=rownames(rfImpsX),
                      variableName=rfImpsX[,5],
                      decreaseAccuracy=rfImpsX[,3],
                      decreaseGini=rfImpsX[,4],
                      rank=rank(-rfImpsX[,3]),
                      order=c(1,2,3,5,4))

# partial dependence
set.seed(4)
partDepX <- partial_dependence(rfX,c('ldisp','toler','r','E','H'),data=dfN2[sampN+100000,])

# dataframe for partial dependence
pdpX <- partDepX %>%
  pivot_longer(cols = -c("0","1"), names_to = "variable", values_to = "value") %>%
  drop_na(value) %>%
  dplyr::select(variable, value, "1")

colnames(pdpX) <- c('variable','value','p')

# Plot the data

s1X<-ggplot(dfVarImpX,aes(x=order,y=decreaseAccuracy))+
  geom_point(size=3)+
  geom_text(aes(label=variableName),parse=T,vjust = 0, nudge_y = 0.0025)+
  geom_segment(linetype=2,aes(xend=order,y=0,yend=decreaseAccuracy))+
  scale_x_continuous(name="")+
  scale_y_continuous(name="Variable importance\n(mean decrease in\naccuracy)",
                     expand=c(0,0),limits=c(0,max(dfVarImpX$decreaseAccuracy)*1.2))+
  theme_classic()+
  theme(legend.position = "none",axis.text.x = element_blank())

s2X <- ggplot(dfVarImpX,aes(x=rank,y=decreaseGini))+
  geom_point(size=3)+
  geom_text(aes(label=variableName),parse=T,vjust = 0, nudge_y = 250)+
  geom_segment(linetype=2,aes(xend=rank,y=0,yend=decreaseGini))+
  #geom_hline(yintercept=0,linetype=3)+
  scale_x_continuous(name="")+
  scale_y_continuous(name="Variable importance\n(mean decrease in Gini)",
                     expand=c(0,0),limits=c(0,max(dfVarImpX$decreaseGini)*1.2))+
  theme_classic()+
  theme(legend.position = "none",axis.text.x = element_blank())

s3X <- ggplot(subset(pdpX,variable=='ldisp'),aes(x=value,y=1-p))+
  geom_line(linewidth=1.5)+
  scale_x_continuous(name=(expression(gamma: " Mean dispersal distance")),
                     breaks=c(-4,-2,0),labels=c("0.0001","0.01","1"),expand=c(0,0),
                     limits=c(-4,0))+
  scale_y_continuous(name="Probability of\npersisting",limits=c(0,1),expand=c(0,0))+
  theme_classic()

s4X <-  ggplot(subset(pdpX,variable=='toler'),aes(x=value,y=1-p))+
  geom_line(linewidth=1.5)+
  scale_x_continuous(name=(expression(sigma: " Thermal tolerance breadth")),
                     limits=c(2,5),breaks=c(2,3,4,5),expand=c(0,0))+
  scale_y_continuous(name="Probability of\npersisting",limits=c(0,1),expand=c(0,0))+
  theme_classic()

s5X <- ggplot(subset(pdpX,variable=='r'),aes(x=value,y=1-p))+
  geom_line(linewidth=1.5)+
  scale_x_continuous(name=(expression(rho: " Fecundity")),expand=c(0,0),
                     limits=c(1.5,6),breaks=c(2,4,6))+
  scale_y_continuous(name="Probability of\npersisting",limits=c(0,1),expand=c(0,0))+
  theme_classic()

s6X <- ggplot(subset(pdpX,variable=='H'),aes(x=value,y=1-p))+
  geom_line(linewidth=1.5)+
  scale_x_continuous(name=(expression(H: " Local temperature heterogeneity")),
                     expand=c(0,0),limits=c(0,2),breaks=c(0,1,2))+
  scale_y_continuous(name="Probability of\npersisting",limits=c(0,1),expand=c(0,0))+
  theme_classic()

s7X <- ggplot(subset(pdpX,variable=='E'),aes(x=value,y=1-p))+
  geom_line(linewidth=1.5)+
  scale_x_continuous(name=(expression(S: " Environmental stochasticity")),
                     expand=c(0,0),limits=c(0,1),breaks=c(0,0.5,1))+
  scale_y_continuous(name="Probability of\npersisting",limits=c(0,1),expand=c(0,0))+
  theme_classic()


hetPersL <- ggplot(sdfH,aes(x= as.numeric(dispCat2),y= as.numeric(hetCat2),fill=1-extM))+
  geom_tile()+
  scale_fill_viridis_c(name="Probability of\npersisting",limits=c(0,1),option="magma")+
  scale_x_continuous(name="Mean dispersal distance",breaks=c(0.5,8.5,16.5),labels=c("0.0001","0.01","1"),expand=c(0,0))+
  scale_y_continuous(name="Local climate heterogeneity",breaks=c(0.5,8.5,16.5),labels=c(0,1,2),expand=c(0,0))+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5,direction='horizontal')) +
  coord_fixed()

legendP<-get_legend(hetPersL)

hetPers <- ggplot(sdfH,aes(x= as.numeric(dispCat2),y= as.numeric(hetCat2),fill=1-extM))+
  geom_tile()+
  scale_fill_viridis_c(name="Probability of\npersisting",limits=c(0,1),option="magma",guide=F)+
  scale_x_continuous(name="Mean dispersal distance",breaks=c(0.5,8.5,16.5),labels=c("0.0001","0.01","1"),expand=c(0,0))+
  scale_y_continuous(name="Local climate heterogeneity",breaks=c(0.5,8.5,16.5),labels=c(0,1,2),expand=c(0,0))+
  coord_fixed()

stoPers <- ggplot(sdfE,aes(x= as.numeric(dispCat2),y= as.numeric(stoCat2),fill=1-extM))+
  geom_tile()+
  scale_fill_viridis_c(name="Persistence\nlikelihood",limits=c(0,1),guide=F,option="magma")+
  scale_x_continuous(name=(expression(gamma: " Mean dispersal distance")),
                     breaks=c(0.5,8.5,16.5),labels=c("0.0001","0.01","1"),expand=c(0,0))+
  scale_y_continuous(name=(expression(S: " Environmental stochasticity")),
                     breaks=c(0.5,8.5,16.5),labels=c(0,0.5,1),expand=c(0,0))+
  coord_fixed()

# final version of figure S3 (or S4)
f2X <- plot_grid(plot_grid(s1X,s3X,s4X,s5X,s6X,s7X,ncol=2,
                           labels=c('a','b','c','d','e','f')),
                 plot_grid(legendP,
                           plot_grid(hetPers,stoPers,labels=c('h','g'),align='v',ncol=1),
                           ncol=1,rel_heights = c(1,4)),
                 rel_widths = c(2.2,1),ncol=2)


#### figure 4 (and S6) ####
# We want to categorize dfH by whether management was positive or negative and how strong that was
dfH$percImp<-dfH$Nf/dfH$NfNA # how much does management improve the final population? final population divided by final population with no action
# sorry, this is a bit weird, but it works
# we want to categorize the outcomes by percent changed
# first for management that modifies heterogeneity
dfH$impCat <- NaN # make a vector to store categorizations
dfH$impCat[dfH$percImp==0]<-0 # if exctinct with management but not without, category 0
dfH$impCat[dfH$percImp>0 & dfH$percImp<0.1]<-1 # if managed is <10% of managed, category 1
dfH$impCat[dfH$percImp>=0.1 & dfH$percImp<0.5]<-2 # 10-50%, category 2
dfH$impCat[dfH$percImp>=0.5 & dfH$percImp<0.95]<-3 # 50-95%, category 3
dfH$impCat[dfH$percImp>=0.95 & dfH$percImp<=1.05]<-4 # 95-105% (roughly neutral), category 4
dfH$impCat[dfH$percImp>1.05 & dfH$percImp<=2]<-5 # 105-200%, category 5
dfH$impCat[dfH$percImp>2 & dfH$percImp<=10]<-6 # 200-1000%, category 6
dfH$impCat[dfH$percImp>10 &  dfH$percImp<=10000000000000]<-7 # >1000%, category 7
dfH$impCat[is.infinite(dfH$percImp) & dfH$Ni>0]<-8 # managed persists when unmanaged went extinct, category 8
# anything left are populations that go extinct with or without management

# now do the same for management that modifies stochasticity
dfE$percImp<-dfE$Nf/dfE$NfNA # how much does management improve the final population? final population divided by final population with no action
# we want to categorize the outcomes by percent changed
dfE$impCat <- NaN # make a vector to store categorizations
dfE$impCat[dfE$percImp==0]<-0 # if exctinct with management but not without, category 0
dfE$impCat[dfE$percImp>0 & dfE$percImp<0.1]<-1 # if managed is <10% of managed, category 1
dfE$impCat[dfE$percImp>=0.1 & dfE$percImp<0.5]<-2 # 10-50%, category 2
dfE$impCat[dfE$percImp>=0.5 & dfE$percImp<0.95]<-3 # 50-95%, category 3
dfE$impCat[dfE$percImp>=0.95 & dfE$percImp<=1.05]<-4 # 95-105% (roughly neutral), category 4
dfE$impCat[dfE$percImp>1.05 & dfE$percImp<=2]<-5 # 105-200%, category 5
dfE$impCat[dfE$percImp>2 & dfE$percImp<=10]<-6 # 200-1000%, category 6
dfE$impCat[dfE$percImp>10 &  dfE$percImp<=10000000000000]<-7 # >1000%, category 7
dfE$impCat[is.infinite(dfE$percImp) & dfE$Ni>0]<-8 # managed persists when unmanaged went extinct, category 8
# anything left are populations that go extinct with or without management

# now let's plot them
# increased heterogeneity over dispersal
dc1<-ggplot(dfH[dfH$h1==1 & dfH$h2==1 & dfH$Ni>0,],aes(x=(disp),stat(count),fill=factor(-cat2)))+
  geom_histogram(bins=100,position="fill")+
  scale_fill_manual(name="Relative change in\npopulation size\nbetween management\nand no action",
                    values = c("#000080","#0000FF","#0080FF","#00FFFF","#808080","#FFFF00","#FF8000","#FF0000","#800000","#000000"),
                    na.value="white",
                    labels=c("Extinction prevented",">10","2-10","1.05-2 (increased)","0.95-1.05 (neutral)","0.5-0.95 (decreased)","0.1-0.5","<0.1","Extinction caused","Always extinct"))+
  scale_y_continuous(expand=c(0,0),name="Proportion\nof simulations",breaks=c(0,0.25,0.5,0.75,1),labels=c(0,"",0.5,"",1))+
  scale_x_log10(name="Mean dispersal",limit=c(0.0001,1),breaks=c(0.0001,0.001,0.01,0.1,1),label=c("0.0001","","0.01","","1"),expand=c(0,0))+
  guides(fill = guide_legend(override.aes = list(colour = "black")),fill=guide_legend(ncol=2))+
  ggtitle("Increased heterogeneity")
# increased heterogeneity over thermal tolerance
tc1<-ggplot(dfH[dfH$h1==1 & dfH$h2==1 & dfH$Ni>0,],aes(x=(toler),stat(count),fill=factor(-cat2)))+
  geom_histogram(bins=100,position="fill")+
  scale_fill_manual(name="Effect of increasing\nlocal heterogeneity",
                    values = c("#000080","#0000FF","#0080FF","#00FFFF","#808080","#FFFF00","#FF8000","#FF0000","#800000","#000000"),
                    na.value = "white",
                    guide=F)+
  scale_y_continuous(expand=c(0,0),name="Proportion\nof simulations",breaks=c(0,0.25,0.5,0.75,1),labels=c(0,"",0.5,"",1))+
  scale_x_continuous(name="Thermal tolerance",limit=c(2,5),breaks=c(2,3,4,5),label=c(2,3,4,5),expand=c(0,0))
# increased heterogeneity over fecundity
rc1<-ggplot(dfH[dfH$h1==1 & dfH$h2==1 & dfH$Ni>0,],aes(x=(r),stat(count),fill=factor(-cat2)))+
  geom_histogram(bins=100,position="fill")+
  scale_fill_manual(name="Relative change in\npopulation size\nbetween management\nand no action",
                    values = c("#000080","#0000FF","#0080FF","#00FFFF","#808080","#FFFF00","#FF8000","#FF0000","#800000","#000000"),
                    na.value="white",
                    guide=F)+
  scale_y_continuous(expand=c(0,0),name="Proportion\nof simulations",breaks=c(0,0.25,0.5,0.75,1),labels=c(0,"",0.5,"",1))+
  scale_x_continuous(name="Fecundity",limit=c(1.5,6),breaks=c(2,3,4,5,6),label=c(2,"",4,"",6),expand=c(0,0))
# increased heterogeneity over stochasticity
sc1<-ggplot(dfH[dfH$h1==1 & dfH$h2==1 & dfH$Ni>0,],aes(x=(E),stat(count),fill=factor(-cat2)))+
  geom_histogram(bins=100,position="fill")+
  scale_fill_manual(name="Relative change in\npopulation size\nbetween management\nand no action",
                    values = c("#000080","#0000FF","#0080FF","#00FFFF","#808080","#FFFF00","#FF8000","#FF0000","#800000","#000000"),
                    na.value="white",
                    guide=F)+
  scale_y_continuous(expand=c(0,0),name="Proportion\nof simulations",breaks=c(0,0.25,0.5,0.75,1),labels=c(0,"",0.5,"",1))+
  scale_x_continuous(name="Environmental stochasticity",limit=c(0,1),breaks=c(0,0.25,0.5,0.75,1),label=c(0,"",0.5,"",1),expand=c(0,0))

# decreased stochasticity over dispersal
dc2<-ggplot(dfE[dfE$e1==0.5 & dfE$e2==0.25 & dfE$Ni>0,],aes(x=(disp),stat(count),fill=factor(-cat2)))+
  geom_histogram(bins=100,position="fill")+
  scale_fill_manual(name="Relative change in population size\nfrom reducing environmental stochasticity",
                    values = c("#000080","#0000FF","#0080FF","#00FFFF","#808080","#FFFF00","#FF8000","#FF0000","#800000","#000000"),
                    na.value="white",
                    labels=c("Extinction prevented",">10","2-10","1.05-2 (increased)","0.95-1.05 (neutral)","0.5-0.95 (decreased)","0.1-0.5","<0.1","Extinction caused","Always extinct"))+
  scale_y_continuous(expand=c(0,0),name="Proportion\nof simulations",breaks=c(0,0.25,0.5,0.75,1),labels=c(0,"",0.5,"",1))+
  scale_y_continuous(expand=c(0,0),name="Proportion\nof simulations",breaks=c(0,0.25,0.5,0.75,1),labels=c(0,"",0.5,"",1))+
  scale_x_log10(name="Mean dispersal",limit=c(0.0001,1),breaks=c(0.0001,0.001,0.01,0.1,1),label=c("0.0001","","0.01","","1"),expand=c(0,0))+
  guides(fill = guide_legend(override.aes = list(colour = "black")))+
  ggtitle("Decreased stochasticity")
# decreased stochasticity over thermal tolerance
tc2<-ggplot(dfE[dfE$e1==0.5 & dfE$e2==0.25 & dfE$Ni>0,],aes(x=(toler),stat(count),fill=factor(-cat2)))+
  geom_histogram(bins=100,position="fill")+
  scale_fill_manual(name="Relative change in\npopulation size\nbetween management\nand no action",
                    values = c("#000080","#0000FF","#0080FF","#00FFFF","#808080","#FFFF00","#FF8000","#FF0000","#800000","#000000"),
                    na.value="white",
                    guide=F)+
  scale_y_continuous(expand=c(0,0),name="Proportion\nof simulations",breaks=c(0,0.25,0.5,0.75,1),labels=c(0,"",0.5,"",1))+
  scale_x_continuous(name="Thermal tolerance",limit=c(2,5),breaks=c(2,3,4,5),label=c(2,3,4,5),expand=c(0,0))
# decreased stochasticity over fecundity
rc2<-ggplot(dfE[dfE$e1==0.5 & dfE$e2==0.25 & dfE$Ni>0,],aes(x=(r),stat(count),fill=factor(-cat2)))+
  geom_histogram(bins=100,position="fill")+
  scale_fill_manual(name="Relative change in\npopulation size\nbetween management\nand no action",
                    values = c("#000080","#0000FF","#0080FF","#00FFFF","#808080","#FFFF00","#FF8000","#FF0000","#800000","#000000"),
                    na.value="white",
                    guide=F)+
  scale_y_continuous(expand=c(0,0),name="Proportion\nof simulations",breaks=c(0,0.25,0.5,0.75,1),labels=c(0,"",0.5,"",1))+
  scale_x_continuous(name="Fecundity",limit=c(1.5,6),breaks=c(2,3,4,5,6),label=c(2,"",4,"",6),expand=c(0,0))
# decreased stochasticity over heterogeneity
ec1<-ggplot(dfE[dfE$e1==0.5 & dfE$e2==0.25 & dfE$Ni>0,],aes(x=(H),stat(count),fill=factor(-cat2)))+
  geom_histogram(bins=100,position="fill")+
  scale_fill_manual(name="Relative change in\npopulation size\nbetween management\nand no action",
                    values = c("#000080","#0000FF","#0080FF","#00FFFF","#808080","#FFFF00","#FF8000","#FF0000","#800000","#000000"),
                    na.value="white",
                    guide=F)+
  scale_y_continuous(expand=c(0,0),name="Proportion\nof simulations",breaks=c(0,0.25,0.5,0.75,1),labels=c(0,"",0.5,"",1))+
  scale_x_continuous(name="Local heterogeneity",limit=c(0,2),breaks=c(0,0.5,1,1.5,2),label=c(0,"",1,"",2),expand=c(0,0))


# now find the legends
legend1<-get_legend(dc1+theme(legend.title = element_text(size = 12),legend.text = element_text(size = 8),legend.box.margin = ggplot2::margin(0, 30, 30, 30)))
legend2<-get_legend(dc2+theme(legend.title = element_text(size = 12),legend.text = element_text(size = 8),legend.box.margin = ggplot2::margin(0, 30, 30, 30)))

# plot everything without the legends
noleg <- plot_grid(dc1+theme(legend.position = "none"),dc2+theme(legend.position = "none"),tc1,tc2,rc1,rc2,nrow=3,labels=c("a","b","c","d","e","f"),align="h")
# final figure 4 (or S6)
pg3 <- plot_grid(noleg,legend1,nrow=1,rel_widths = c(2, 0.6))


#### Figure 3 ####
# make sure you run figure 4 first
# Separate dfN into categories again, but this time depending on whether they are in the lowest or highest 25% quantiles
# for dispersal
dfN$dispCat2 <- 0-(dfN$disp<quantile(dfN$disp,0.25))
dfN$dispCat2 <- dfN$dispCat2+(dfN$disp>quantile(dfN$disp,0.75))
# for thermal tolerance
dfN$tolerCat <- 0-(dfN$toler<quantile(dfN$toler,0.25))
dfN$tolerCat <- dfN$tolerCat+(dfN$toler>quantile(dfN$toler,0.75))
# and for fecundity
dfN$rCat <- 0-(dfN$r<quantile(dfN$r,0.25))+0
dfN$rCat <- dfN$rCat+(dfN$r>quantile(dfN$r,0.75))+0
# categorize by whether the population benefits from each type of management
dfN$benH <- 0+(dfH$impCat>4)
dfN$benE <- 0+(dfE$impCat>4)
# categorize by whether the population is hurt by the management
dfN$badH <- 0+(dfH$impCat<4)
dfN$badE <- 0+(dfE$impCat<4)

# this summarizes the percent benefit and percent harmed by management 

sdf <- ddply(dfN[dfN$Ni>0 & dfN$dispCat %in% c(-1,1) & dfN$tolerCat %in% c(-1,1) & dfN$rCat %in% c(-1,1),],c("dispCat","tolerCat","rCat"),summarize,N=length(i),
             persP=1-sum(ext)/length(ext),
             meanN=mean(Nf/Ni,na.rm=T),
             medianN=median(Nf/Ni,na.rm=T),
             benH=sum(benH,na.rm=T)/length(benH),
             benE=sum(benE,na.rm=T)/length(benE),
             badH=sum(badH,na.rm=T)/length(badH),
             badE=sum(badE,na.rm=T)/length(badE))

# this summary data was converted into figure 3 in the paper
