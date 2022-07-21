rm(list=ls())

require(plyr)
require(ggplot2)
require(viridis)
require(cowplot)
require(reshape2)

#### initial processing
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


#### Figure 2a,b,c,d,e,f ####
require(randomForest)
samps<-100000 # sample size in random forest
rf<-randomForest(factor(ext)~ldisp+toler+r+E+H,dfN,importance=T,sampsize=c(samps,samps),scale=F,do.trace=T) # random forest calculation

# Function to calculate partial dependence
# Largely based on the partialPlot funciton in the randomforest package
partialDep <- function(rf,pred.data,x.var,log=F,L=201){
  x <- pred.data[,x.var]
  n <- nrow(pred.data)
  if(log){
    x.pt <- exp(seq(log(min(x)),log(max(x)),length=L))
  } else{
    x.pt <- seq(min(x),max(x),length=L)
  }
  y.pt<-numeric(L)
  for(i in 1:L){
    x.data<-pred.data
    x.data[,x.var] <- x.pt[i]
    pr<-predict(rf,x.data,type="prob")
    y.pt[i]<-mean(log(pr[,1])-rowMeans(log(pr)),na.rm=T)
    y.pt[i]<-mean(log(ifelse(pr[,1] == 0,
                             .Machine$double.eps, pr[,1]))
                  - rowMeans(log(ifelse(pr == 0, .Machine$double.eps, pr))),na.rm=TRUE)
  }
  retFrame<-data.frame(matrix(c(x.pt,y.pt),L))
  colnames(retFrame)<-c(x.var,"y")
  return(retFrame)
}

# Turn random forest variable importance results into a dataframe
dfVarImp<-data.frame(varImp=rf$importance[,3],rank=rank(-rf$importance[,3]),vars=rownames(rf$importance))
dfVarImp$vars <- rownames(dfVarImp) # this is the name for each of the independent variables, but you might want to change it to the symbols used in the paper
# plot the relative variable importance for each of the independent variables
s1<-ggplot(dfVarImp,aes(x=rank,y=varImp))+
  geom_point(size=3)+
  geom_text(aes(label=vars),parse=T,vjust = 0, nudge_y = 0.002)+
  geom_segment(linetype=2,aes(xend=rank,y=0,yend=varImp))+
  geom_hline(yintercept=0,linetype=3)+
  scale_x_continuous(name="Importance rank")+
  scale_y_continuous(name="Variable importance")+
  theme_classic()+
  theme(legend.position = "none")

# Now calculate the partial dependence for each of the independent variables
samp <- 100001:110000
pdDisp<-partialDep(rf,dfN[samp,],'ldisp',L=201)
pdToler<-partialDep(rf,dfN[samp,],'toler',L=201)
pdr<-partialDep(rf,dfN[samp,],'r',L=201)
pdH<-partialDep(rf,dfN[samp,],'H',L=201)
pdE<-partialDep(rf,dfN[samp,],'E',L=201)

# And plot each of them
# dispersal distance (log transformed)
s2 <- ggplot(pdDisp,aes(x=ldisp,y=y))+
  geom_line(size=1.5)+
  geom_hline(yintercept = 0,linetype=2)+
  scale_x_continuous(name=(expression(gamma: " Mean dispersal distance")))+
  scale_y_continuous(name='Partial dependence',limits=c(-12,12))+
  theme_classic()
# thermal tolerance breadth
s3 <- ggplot(pdToler,aes(x=toler,y=y))+
  geom_line(size=1.5)+
  geom_hline(yintercept = 0,linetype=2)+
  scale_x_continuous(name=(expression(sigma: " Thermal tolerance breadth")))+
  scale_y_continuous(name='Partial dependence',limits=c(-12,12))+
  theme_classic()
# fecundity
s4 <- ggplot(pdr,aes(x=r,y=y))+
  geom_line(size=1.5)+
  geom_hline(yintercept = 0,linetype=2)+
  scale_x_continuous(name=(expression(rho: " Fecundity")))+
  scale_y_continuous(name='Partial dependence',limits=c(-12,12))+
  theme_classic()
# local temperature heterogeneity
s5 <- ggplot(pdH,aes(x=H,y=y))+
  geom_line(size=1.5)+
  geom_hline(yintercept = 0,linetype=2)+
  scale_x_continuous(name=(expression(H: " Local temperature heterogeneity")))+
  scale_y_continuous(name='Partial dependence',limits=c(-12,12))+
  theme_classic()
# environmental stochasticity
s6 <- ggplot(pdE,aes(x=E,y=y))+
  geom_line(size=1.5)+
  geom_hline(yintercept = 0,linetype=2)+
  scale_x_continuous(name=(expression(S: " Environmental stochasticity")))+
  scale_y_continuous(name='Partial dependence',limits=c(-12,12))+
  theme_classic()

#plot them all together
plot_grid(s1,s2,s3,s4,s5,s6,nrow=3,labels=c('a','b','c','d','e','f'))


#### Figure 2g,h ####
# discretize the model results into categories for average dispersal, local temperature heterogeneity, and environmental stochasticity
dfN$dispCat<-cut(dfN$ldisp, breaks=seq(-4,0,by=0.25),labels=seq(1,16),include.lowest = T)
dfN$hetCat<-cut(dfN$H, breaks=seq(0,2,by=0.125),labels=seq(1,16),include.lowest = T)
dfN$stoCat<-cut(dfN$E, breaks=seq(0,1,by=0.0625),labels=seq(1,16),include.lowest = T)

# Now summarize the results for each combination of dispersal and heterogeneity category
sdfH<-ddply(dfN,c("dispCat","hetCat"),summarize,N=length(i),
            NiM=mean(Ni,na.rm=T),NfM=mean(Nf,na.rm=T), # mean initial population and mean final population
            RiM=mean(Ri,na.rm=T),RfM=mean(Rf,na.rm=T), # mean initial range and mean final range
            extM=mean(ext,na.rm=T),dispM=mean(disp,na.rm=T)) # proportion extinct and mean dispersal

# Same for environmental stochasticity
sdfE<-ddply(dfN,c("dispCat","stoCat"),summarize,N=length(i),
            NiM=mean(Ni,na.rm=T),NfM=mean(Nf,na.rm=T),
            RiM=mean(Ri,na.rm=T),RfM=mean(Rf,na.rm=T),
            extM=mean(ext,na.rm=T),dispM=mean(disp,na.rm=T))

# Plot the persistence likelihood depending on dispersal and heterogeneity
hetPers <- ggplot(sdfH,aes(x= as.numeric(dispCat),y= as.numeric(hetCat),fill=1-extM))+
  geom_tile()+
  scale_fill_viridis_c(name="Persistence\nlikelihood",limits=c(0,1),option="magma")+
  scale_x_continuous(name="Dispersal ability",breaks=c(0.5,8.5,16.5),labels=c(0.0001,0.01,1),expand=c(0,0))+
  scale_y_continuous(name="Local climate heterogeneity",breaks=c(0.5,8.5,16.5),labels=c(0,1,2),expand=c(0,0))+
  coord_fixed()
# Grab the legend for a multiplot
legendP<-get_legend(hetPers)
# but remove the legend so we don't have it later
hetPers <- hetPers+theme(legend.position="none")

# Plot the average proportion of remaining population depending on dispersal and heterogeneity
hetN <- ggplot(sdfH,aes(x= as.numeric(dispCat),y= as.numeric(hetCat),fill=NfM/NiM))+
  geom_tile()+
  scale_fill_viridis_c(name="Average\npercent\npopulation\nremaining",limits=c(0,1))+
  scale_x_continuous(name="Dispersal ability",breaks=c(0.5,8.5,16.5),labels=c(0.0001,0.01,1),expand=c(0,0))+
  scale_y_continuous(name="Local climate heterogeneity",breaks=c(0.5,8.5,16.5),labels=c(0,1,2),expand=c(0,0))+
  coord_fixed()
# Grab the legend for a multiplot
legendN<-get_legend(hetN)
# but remove the legend so we don't have it later
hetN <- hetN+theme(legend.position="none")

# Plot the persistence likelihood depending on dispersal and stochasticity
stoPers <- ggplot(sdfE,aes(x= as.numeric(dispCat),y= as.numeric(stoCat),fill=1-extM))+
  geom_tile()+
  scale_fill_viridis_c(name="Persistence\nlikelihood",limits=c(0,1),option="magma")+
  scale_x_continuous(name="Dispersal ability",breaks=c(0.5,8.5,16.5),labels=c(0.0001,0.01,1),expand=c(0,0))+
  scale_y_continuous(name="Envrionmental stochasticity",breaks=c(0.5,8.5,16.5),labels=c(0,0.5,1),expand=c(0,0))+
  coord_fixed()+
  theme(legend.position="none")

# Plot the average proportion of remaining population depending on dispersal and stochasticity
stoN <- ggplot(sdfE,aes(x= as.numeric(dispCat),y= as.numeric(stoCat),fill=NfM/NiM))+
  geom_tile()+
  scale_fill_viridis_c(name="Average\npercent\npopulation\nremaining",limits=c(0,1))+
  scale_x_continuous(name="Dispersal ability",breaks=c(0.5,8.5,16.5),labels=c(0.0001,0.01,1),expand=c(0,0))+
  scale_y_continuous(name="Envrionmental stochasticity",breaks=c(0.5,8.5,16.5),labels=c(0,0.5,1),expand=c(0,0))+
  coord_fixed()+
  theme(legend.position="none")

# Plot them all together
plot_grid(hetPers,stoPers,legendP, hetN,stoN,legendN,nrow=2,labels=c("a","b","","c","d",""),rel_widths = c(3,3,2))



#### Figure 4 ####
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
dc1<-ggplot(dfH[dfH$Ni>0,],aes(x=(disp),stat(count),fill=factor(-impCat)))+
  geom_histogram(bins=100,position="fill")+
  scale_fill_manual(name="Relative change in\npopulation size\nbetween management\nand no action",
                    values = c("#800000","#FF0000","#FF8000","#FFFF00","#808080","#00FFFF","#0080FF","#0000FF","#000080","#000000"),
                    na.value="white",
                    labels=c("Extinction prevented",">10","2-10","1.05-2 (increased)","0.95-1.05 (neutral)","0.5-0.95 (decreased)","0.1-0.5","<0.1","Extinction caused","Always extinct"))+
  scale_y_continuous(expand=c(0,0),name="Proportion\nof simulations",breaks=c(0,0.25,0.5,0.75,1),labels=c(0,"",0.5,"",1))+
  scale_x_log10(name="Mean dispersal",limit=c(0.0001,1),breaks=c(0.0001,0.001,0.01,0.1,1),label=c("0.0001","","0.01","","1"),expand=c(0,0))+
  guides(fill = guide_legend(override.aes = list(colour = "black")),fill=guide_legend(ncol=2))+
  ggtitle("Increased heterogeneity")
# increased heterogeneity over thermal tolerance
tc1<-ggplot(dfH[dfH$Ni>0,],aes(x=(toler),stat(count),fill=factor(-impCat)))+
  geom_histogram(bins=100,position="fill")+
  scale_fill_manual(name="Effect of increasing\nlocal heterogeneity",
                    values = c("#800000","#FF0000","#FF8000","#FFFF00","#808080","#00FFFF","#0080FF","#0000FF","#000080","#000000"),
                    na.value = "white",
                    guide=F)+
  scale_y_continuous(expand=c(0,0),name="Proportion\nof simulations",breaks=c(0,0.25,0.5,0.75,1),labels=c(0,"",0.5,"",1))+
  scale_x_continuous(name="Thermal tolerance",limit=c(2,5),breaks=c(2,3,4,5),label=c(2,3,4,5),expand=c(0,0))
# increased heterogeneity over fecundity
rc1<-ggplot(dfH[dfH$Ni>0,],aes(x=(r),stat(count),fill=factor(-impCat)))+
  geom_histogram(bins=100,position="fill")+
  scale_fill_manual(name="Relative change in\npopulation size\nbetween management\nand no action",
                    values = c("#800000","#FF0000","#FF8000","#FFFF00","#808080","#00FFFF","#0080FF","#0000FF","#000080","#000000"),
                    na.value="white",
                    guide=F)+
  scale_y_continuous(expand=c(0,0),name="Proportion\nof simulations",breaks=c(0,0.25,0.5,0.75,1),labels=c(0,"",0.5,"",1))+
  scale_x_continuous(name="Fecundity",limit=c(1.5,6),breaks=c(2,3,4,5,6),label=c(2,"",4,"",6),expand=c(0,0))
# increased heterogeneity over stochasticity
sc1<-ggplot(dfH[dfH$Ni>0,],aes(x=(E),stat(count),fill=factor(-impCat)))+
  geom_histogram(bins=100,position="fill")+
  scale_fill_manual(name="Relative change in\npopulation size\nbetween management\nand no action",
                    values = c("#800000","#FF0000","#FF8000","#FFFF00","#808080","#00FFFF","#0080FF","#0000FF","#000080","#000000"),
                    na.value="white",
                    guide=F)+
  scale_y_continuous(expand=c(0,0),name="Proportion\nof simulations",breaks=c(0,0.25,0.5,0.75,1),labels=c(0,"",0.5,"",1))+
  scale_x_continuous(name="Environmental stochasticity",limit=c(0,1),breaks=c(0,0.25,0.5,0.75,1),label=c(0,"",0.5,"",1),expand=c(0,0))
# decreased stochasticity over dispersal
dc2<-ggplot(dfE[dfE$Ni>0,],aes(x=(disp),stat(count),fill=factor(-impCat)))+
  geom_histogram(bins=100,position="fill")+
  scale_fill_manual(name="Relative change in population size\nfrom reducing environmental stochasticity",
                    values = c("#800000","#FF0000","#FF8000","#FFFF00","#808080","#00FFFF","#0080FF","#0000FF","#000080","#000000"),
                    na.value="white",
                    labels=c("Extinction prevented",">10","2-10","1.05-2 (increased)","0.95-1.05 (neutral)","0.5-0.95 (decreased)","0.1-0.5","<0.1","Extinction caused","Always extinct"))+
  scale_y_continuous(expand=c(0,0),name="Proportion\nof simulations",breaks=c(0,0.25,0.5,0.75,1),labels=c(0,"",0.5,"",1))+
  scale_x_log10(name="Mean dispersal",limit=c(0.0001,1),breaks=c(0.0001,0.001,0.01,0.1,1),label=c("0.0001","","0.01","","1"),expand=c(0,0))+
  guides(fill = guide_legend(override.aes = list(colour = "black")))+
  ggtitle("Decreased stochasticity")
# decreased stochasticity over thermal tolerance
tc2<-ggplot(dfE[dfE$Ni>0,],aes(x=(toler),stat(count),fill=factor(-impCat)))+
  geom_histogram(bins=100,position="fill")+
  scale_fill_manual(name="Relative change in\npopulation size\nbetween management\nand no action",
                    values = c("#800000","#FF0000","#FF8000","#FFFF00","#808080","#00FFFF","#0080FF","#0000FF","#000080","#000000"),
                    na.value="white",
                    guide=F)+
  scale_y_continuous(expand=c(0,0),name="Proportion\nof simulations",breaks=c(0,0.25,0.5,0.75,1),labels=c(0,"",0.5,"",1))+
  scale_x_continuous(name="Thermal tolerance",limit=c(2,5),breaks=c(2,3,4,5),label=c(2,3,4,5),expand=c(0,0))
# decreased stochasticity over fecundity
rc2<-ggplot(dfE[dfE$Ni>0,],aes(x=(r),stat(count),fill=factor(-impCat)))+
  geom_histogram(bins=100,position="fill")+
  scale_fill_manual(name="Relative change in\npopulation size\nbetween management\nand no action",
                    values = c("#800000","#FF0000","#FF8000","#FFFF00","#808080","#00FFFF","#0080FF","#0000FF","#000080","#000000"),
                    na.value="white",
                    guide=F)+
  scale_y_continuous(expand=c(0,0),name="Proportion\nof simulations",breaks=c(0,0.25,0.5,0.75,1),labels=c(0,"",0.5,"",1))+
  scale_x_continuous(name="Fecundity",limit=c(1.5,6),breaks=c(2,3,4,5,6),label=c(2,"",4,"",6),expand=c(0,0))
# decreased stochasticity over heterogeneity
ec2<-ggplot(dfE[dfE$Ni>0,],aes(x=(H),stat(count),fill=factor(-impCat)))+
  geom_histogram(bins=100,position="fill")+
  scale_fill_manual(name="Relative change in\npopulation size\nbetween management\nand no action",
                    values = c("#800000","#FF0000","#FF8000","#FFFF00","#808080","#00FFFF","#0080FF","#0000FF","#000080","#000000"),
                    na.value="white",
                    guide=F)+
  scale_y_continuous(expand=c(0,0),name="Proportion\nof simulations",breaks=c(0,0.25,0.5,0.75,1),labels=c(0,"",0.5,"",1))+
  scale_x_continuous(name="Local heterogeneity",limit=c(0,2),breaks=c(0,0.5,1,1.5,2),label=c(0,"",1,"",2),expand=c(0,0))

# now find the legends
legend1<-get_legend(dc1+theme(legend.title = element_text(size = 12),legend.text = element_text(size = 8),legend.box.margin = ggplot2::margin(0, 30, 30, 30)))
legend2<-get_legend(dc2+theme(legend.title = element_text(size = 12),legend.text = element_text(size = 8),legend.box.margin = ggplot2::margin(0, 30, 30, 30)))

# plot everything without the legends
noleg <- plot_grid(dc1+theme(legend.position = "none"),dc2+theme(legend.position = "none"),tc1,tc2,rc1,rc2,nrow=3,labels=c("a","b","c","d","e","f"),align="h")
plot_grid(noleg,legend1,nrow=1,rel_widths = c(2, 0.6)) # add the legends


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


#### SOMETHING IS VERY WRONG AND NEEDS FXING?!
# this summarizes the percent benefit and percent harmed by management 
sdf <- ddply(dfN[dfN$Ni>0 & dfN$dispCat2 %in% c(-1,1) & dfN$tolerCat %in% c(-1,1) & dfN$rCat %in% c(-1,1),],c("dispCat2","tolerCat","rCat"),summarize,N=length(i),
             persP=1-sum(ext)/length(ext),
             benH=sum(benH,na.rm=T)/length(benH),
             badH=sum(badH,na.rm=T)/length(badH),
             benE=sum(benE,na.rm=T)/length(benE),
             badE=sum(badE,na.rm=T)/length(badE))

# this summary data was converted into figure 3 in the paper

