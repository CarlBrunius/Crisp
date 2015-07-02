rm(list=ls())
setwd('C:/R/Crisp')

## Continued from Drift2.r

## Load packages required for parallell processing and multivariate analysis
library(doParallel)
library(foreach)
library(mixOmics)
library(pROC)

### Load data 
load(file='trtData.RData')  # Data file corresponding to LCMSDataCorrected, decoded with respect to individual and treatment

###Prepare data for PCA-analysis
idData=matrix(nrow=length(unique(trtData$ID)),ncol=ncol(trtData)-2)
ID=numeric(length(unique(trtData$ID)))
### Make PCAs to check data quality
n=0
for (i in unique(trtData$ID)) {
  n=n+1
  ID[n]=i
  idData[n,]=colMeans(trtData[trtData$ID==i,-1:-2])
}
colnames(idData)=colnames(trtData)[-1:-2]
idData=as.data.frame(idData)
idData=cbind(ID,idData)
library(xlsx)
write.xlsx(trtData,'pcaData.xlsx','trtData')
write.xlsx(idData,'pcaData.xlsx','idData',append=T)
write.csv2(trtData,'trtData.csv')
write.csv2(idData,'idData.csv')
## PCA in SIMCA
#  id:  308 excluded from model (too many NAs)
#       728 outside t2
#       708 high DModX
#  trt: 308 excluded from model (too many NAs)
#       728 celle slightly outside t2
#       554 celle outside t2
#       469 husman high DModX (1 imputed value at 65 mins)
#       415 husman high DModX
#       728 husman high DModX
#       824 deli high DModX

### Exclude ID 308 - too many imputed values!
trtData=trtData[trtData$ID!=308,]

### Divide according to treatments
celle=trtData[trtData$treatment=='celle',]
celle=celle[order(celle$ID),]
Xcelle=celle[,-1:-2]
deli=trtData[trtData$treatment=='deli',]
# > deli[,1]-celle[,1]
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
deli=deli[order(deli$ID),]
Xdeli=deli[,-1:-2]
husman=trtData[trtData$treatment=='husman',]
# > husman[,1]-celle[,1]
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
husman=husman[order(husman$ID),]
Xhusman=husman[,-1:-2]
str(celle)

### Make contrasts for multilevel analyses
diffDeli=Xdeli-Xcelle
diffHusman=Xhusman-Xcelle
diffDH=Xhusman-Xdeli
dim(diffDeli)
dim(diffHusman)
dim(diffDH)

### Setup parallel processing
# cl=makeCluster(detectCores())  # For using all cores
cl=makeCluster(3)  # For using 3 cores specifically - so you can use your computer for other stuff when you're running statistics.
registerDoParallel(cl)

### Perform 'Big' Multilevel PLS-VS analysis
nRep=60
comps=2
featRatio=0.9
resampInner=TRUE
DeliBig=mplsDAPar(diffDeli,nRep=nRep,comps=comps,featRatio=featRatio,metric='miss',resampInner=resampInner)
HusmanBig=mplsDAPar(diffHusman,nRep=nRep,comps=comps,featRatio=featRatio,metric='miss',resampInner=resampInner)
stopCluster(cl)
save(list=ls()[grep('Big',ls())],file='ML_BigModels.RData')

### Permutation analysis
nRep=8
comps=2
featRatio=0.75
resampInner=TRUE
pStart=1
pMax=80
deliPerm=husmanPerm=list()
cl=makeCluster(4)
registerDoParallel(cl)
for (p in pStart:pMax) {
  cat(paste('\n Permutation ',p,' (from ',pStart,'-',pMax,'):\n',sep=''))
  deliPerm[[p]]=mplsDAPar(diffDeli,nRep=nRep,comps=comps,featRatio=featRatio,metric='miss',resampInner=resampInner,perm=T)
  husmanPerm[[p]]=mplsDAPar(diffHusman,nRep=nRep,comps=comps,featRatio=featRatio,metric='miss',resampInner=resampInner,perm=T)
}
stopCluster(cl)
save(deliPerm,husmanPerm,file='permutation.RData')
load(file='permutation.RData')

# Prepare and plot permutation H0 distributions
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
range=range(c(deliMissMin,deliMissMid,deliMissMax))
range[1]=floor(range[1]/2)*2
range[2]=ceiling(range[2]/2)*2
hist(deliMissMin,10,col=rgb(0,0,0,0.8),xlim=range,ylim=c(0,15),axes=FALSE,main='',xlab='Permutation misclassifications (Deli)') # Min in black
hist(deliMissMid,10,col=rgb(1,0,0,0.5),add=TRUE) # Mid in red
hist(deliMissMax,10,col=rgb(1,1,1,0.3),add=TRUE) # Max in white
axis(1)
axis(2,las=1)
box(bty='l')
legend('topright',legend=c('MinModel','MidModel','MaxModel'),fill=c(rgb(0,0,0,0.8),rgb(1,0,0,0.5),rgb(1,1,1,0.3)))
dev.off()

png(file='HusmanMiss.png',width=1024,height=1024,pointsize=25)
range=range(c(husmanMissMin,husmanMissMid,husmanMissMax))
range[1]=floor(range[1]/2)*2
range[2]=ceiling(range[2]/2)*2
hist(husmanMissMin,10,col=rgb(0,0,0,0.8),xlim=range,ylim=c(0,12),axes=FALSE,main='',xlab='Permutation misclassifications (Husman)') # Min in black
hist(husmanMissMid,10,col=rgb(1,0,0,0.5),add=TRUE) # Mid in red
hist(husmanMissMax,10,col=rgb(1,1,1,0.3),add=TRUE) # Max in white
axis(1)
axis(2,las=1)
box(bty='l')
legend('topright',legend=c('MinModel','MidModel','MaxModel'),fill=c(rgb(0,0,0,0.8),rgb(1,0,0,0.5),rgb(1,1,1,0.3)))
dev.off()
