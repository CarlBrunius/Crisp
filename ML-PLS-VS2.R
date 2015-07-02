rm(list=ls())
setwd('C:/R/Crisp')

### Load data (See commented code below)
load(file='trtData.RData')
load(file='ML_BigModels.RData')
load(file='permutation.RData')

## Maka data tables for Husman vs Celle analysis
featMin=round(HusmanBig$nFeat$minModel[1])
featMid=round(HusmanBig$nFeat$midModel[1])
featMax=round(HusmanBig$nFeat$maxModel[1])
featsMin=HusmanBig$VIPRank$minModel[,1][order(HusmanBig$VIPRank$minModel[,1])[1:featMin]]
featsMid=HusmanBig$VIPRank$midModel[,1][order(HusmanBig$VIPRank$midModel[,1])[1:featMid]]
featsMax=HusmanBig$VIPRank$maxModel[,1][order(HusmanBig$VIPRank$maxModel[,1])[1:featMax]]
VIPMin=HusmanBig$VIPRank$minModel[,1]
VIPMid=HusmanBig$VIPRank$midModel[,1]
VIPMax=HusmanBig$VIPRank$maxModel[,1]
VIPMean=(VIPMin*VIPMid*VIPMax)^(1/3)
VIPs=data.frame(Min=VIPMin,Mid=VIPMid,Max=VIPMax,Geom=VIPMean)
write.xlsx(VIPs,file='HusmanBig.xlsx','VIPs')
n=10
VIPOrdMin=VIPMin[order(VIPMin)][1:n]
VIPOrdMid=VIPMid[order(VIPMid)][1:n]
VIPOrdMax=VIPMax[order(VIPMax)][1:n]
library(xlsx)
write.xlsx(featsMin,file='HusmanBig.xlsx','featsMin',append=TRUE)
write.xlsx(featsMid,file='HusmanBig.xlsx','featsMid',append=TRUE)
write.xlsx(featsMax,file='HusmanBig.xlsx','featsMax',append=TRUE)
write.xlsx(VIPOrdMin,file='HusmanBig.xlsx','VIPMin10',append=TRUE)
write.xlsx(VIPOrdMid,file='HusmanBig.xlsx','VIPMid10',append=TRUE)
write.xlsx(VIPOrdMax,file='HusmanBig.xlsx','VIPMax10',append=TRUE)
unique10=sort(unique(c(names(VIPOrdMin),names(VIPOrdMid),names(VIPOrdMax))))
VIPOrdMin=VIPMin[names(VIPMin)%in%unique10]
VIPOrdMid=VIPMid[names(VIPMid)%in%unique10]
VIPOrdMax=VIPMax[names(VIPMax)%in%unique10]
VIPOrdMean=(VIPOrdMin*VIPOrdMid*VIPOrdMax)^(1/3)
VIP10=data.frame(Min=VIPOrdMin,Mid=VIPOrdMid,Max=VIPOrdMax,Geom=VIPOrdMean)
write.xlsx(VIP10,file='HusmanBig.xlsx','VIP10',append=TRUE)


## Maka data tables for Deli vs Celle analysis
featMin=round(DeliBig$nFeat$minModel[1])
featMid=round(DeliBig$nFeat$midModel[1])
featMax=round(DeliBig$nFeat$maxModel[1])
featsMin=DeliBig$VIPRank$minModel[,1][order(DeliBig$VIPRank$minModel[,1])[1:featMin]]
featsMid=DeliBig$VIPRank$midModel[,1][order(DeliBig$VIPRank$midModel[,1])[1:featMid]]
featsMax=DeliBig$VIPRank$maxModel[,1][order(DeliBig$VIPRank$maxModel[,1])[1:featMax]]
VIPMin=DeliBig$VIPRank$minModel[,1]
VIPMid=DeliBig$VIPRank$midModel[,1]
VIPMax=DeliBig$VIPRank$maxModel[,1]
VIPMean=(VIPMin*VIPMid*VIPMax)^(1/3)
VIPs=data.frame(Min=VIPMin,Mid=VIPMid,Max=VIPMax,Geom=VIPMean)
write.xlsx(VIPs,file='DeliBig.xlsx','VIPs')
n=10
VIPOrdMin=VIPMin[order(VIPMin)][1:n]
VIPOrdMid=VIPMid[order(VIPMid)][1:n]
VIPOrdMax=VIPMax[order(VIPMax)][1:n]
library(xlsx)
write.xlsx(featsMin,file='DeliBig.xlsx','featsMin',append=TRUE)
write.xlsx(featsMid,file='DeliBig.xlsx','featsMid',append=TRUE)
write.xlsx(featsMax,file='DeliBig.xlsx','featsMax',append=TRUE)
write.xlsx(VIPOrdMin,file='DeliBig.xlsx','VIPMin10',append=TRUE)
write.xlsx(VIPOrdMid,file='DeliBig.xlsx','VIPMid10',append=TRUE)
write.xlsx(VIPOrdMax,file='DeliBig.xlsx','VIPMax10',append=TRUE)
unique10=sort(unique(c(names(VIPOrdMin),names(VIPOrdMid),names(VIPOrdMax))))
VIPOrdMin=VIPMin[names(VIPMin)%in%unique10]
VIPOrdMid=VIPMid[names(VIPMid)%in%unique10]
VIPOrdMax=VIPMax[names(VIPMax)%in%unique10]
VIPOrdMean=(VIPOrdMin*VIPOrdMid*VIPOrdMax)^(1/3)
VIP10=data.frame(Min=VIPOrdMin,Mid=VIPOrdMid,Max=VIPOrdMax,Geom=VIPOrdMean)
write.xlsx(VIP10,file='DeliBig.xlsx','VIP10',append=TRUE)

### Take out H0 from Deli and Husman permutations
deliMiss=deliMissMin=deliMissMid=deliMissMax=husmanMiss=husmanMissMin=husmanMissMid=husmanMissMax=numeric()
deliCompMin=deliCompMid=deliCompMax=husmanCompMin=husmanCompMid=husmanCompMax=numeric()
for (p in 1:length(deliPerm)) {
  # deliMiss=c(deliMiss,deliPerm[[p]]$Y
  deliMissMin=c(deliMissMin,deliPerm[[p]]$misClass$minModel)
  deliMissMid=c(deliMissMid,deliPerm[[p]]$misClass$midModel)
  deliMissMax=c(deliMissMax,deliPerm[[p]]$misClass$maxModel)
  deliCompMin=c(deliCompMin,deliPerm[[p]]$nComp$minModel[1])
  deliCompMid=c(deliCompMid,deliPerm[[p]]$nComp$midModel[1])
  deliCompMax=c(deliCompMax,deliPerm[[p]]$nComp$maxModel[1])
  husmanMissMin=c(husmanMissMin,husmanPerm[[p]]$misClass$minModel)
  husmanMissMid=c(husmanMissMid,husmanPerm[[p]]$misClass$midModel)
  husmanMissMax=c(husmanMissMax,husmanPerm[[p]]$misClass$maxModel)
  husmanCompMin=c(husmanCompMin,husmanPerm[[p]]$nComp$minModel[1])
  husmanCompMid=c(husmanCompMid,husmanPerm[[p]]$nComp$midModel[1])
  husmanCompMax=c(husmanCompMax,husmanPerm[[p]]$nComp$maxModel[1])
}
png(file='DeliMiss.png',width=1024,height=1024,pointsize=25)
hist(deliMissMin,10,col=rgb(0,0,0,0.8),xlim=c(0,20),ylim=c(0,15),axes=FALSE,main='',xlab='Permutation misclassifications (Deli)') # Min in black
hist(deliMissMid,10,col=rgb(1,0,0,0.5),add=TRUE) # Mid in red
hist(deliMissMax,10,col=rgb(1,1,1,0.3),add=TRUE) # Max in white
axis(1)
axis(2,las=1)
box(bty='l')
legend('topright',legend=c('MinModel','MidModel','MaxModel'),fill=c(rgb(0,0,0,0.8),rgb(1,0,0,0.5),rgb(1,1,1,0.3)))
dev.off()

png(file='HusmanMiss.png',width=1024,height=1024,pointsize=25)
hist(husmanMissMin,10,col=rgb(0,0,0,0.8),xlim=c(0,20),ylim=c(0,18),axes=FALSE,main='',xlab='Permutation misclassifications (Husman)') # Min in black
hist(husmanMissMid,10,col=rgb(1,0,0,0.5),add=TRUE) # Mid in red
hist(husmanMissMax,10,col=rgb(1,1,1,0.3),add=TRUE) # Max in white
axis(1)
axis(2,las=1)
box(bty='l')
legend('topright',legend=c('MinModel','MidModel','MaxModel'),fill=c(rgb(0,0,0,0.8),rgb(1,0,0,0.5),rgb(1,1,1,0.3)))
dev.off()

### calculate permutation p-values for multilevel models
pDeli=sum(deliMissMin<=4)/length(deliMissMin)
pHusman=sum(husmanMissMin<=4)/length(husmanMissMin)

