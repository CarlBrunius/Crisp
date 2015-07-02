rm(list=ls())
setwd('C:/R/Crisp')

# Continued from Drift1.r

library(xcms)
library(CAMERA)

# Function to calculate coefficients of variation per column (ie feature) in a feature matrix
cv=function(mat) {
  mean=apply(mat,2,mean)
  sd=apply(mat,2,sd)
  return(cv=sd/mean)
}

load(file='Data5.RData')
load(file='BatchData.Rdata')
venn.diagram(list(B1a=B1aFeats$finalVars,B1b=B1bFeats$finalVars,B2=B2Feats$finalVars,B3=B3Feats$finalVars,B4=B4Feats$finalVars),file='batches.png')

### Perform between batch correction
# Final feature matrix -> AFeats
batches=c('B1a','B1b','B2','B3','B4')
BatchA=get(paste(batches[1],'Feats',sep=''))
AFeats=BatchA$TestFeatsFinal
for (b in batches[2:length(batches)]) {
  BatchB=get(paste(b,'Feats',sep=''))
  Avars=colnames(AFeats)
  Bvars=BatchB$finalVars
  ABvars=Avars[Avars%in%Bvars]
  AFeats=subset(AFeats,select=ABvars)
  BFeats=subset(BatchB$TestFeatsFinal,select=ABvars)
  Asamps=rownames(AFeats)
  Bsamps=rownames(BatchB$TestFeatsFinal)
  AQCInd=grep('Q',Asamps)
  BQCInd=grep('Q',Bsamps)
  AQCFeatMeans=colMeans(AFeats[AQCInd,])
  BQCFeatMeans=colMeans(BFeats[BQCInd,])
  corrVect=AQCFeatMeans/BQCFeatMeans
  BFeats=t(t(BFeats)*corrVect)
  AFeats=rbind(AFeats,BFeats)
}

sampNames=rownames(AFeats)
featNames=colnames(AFeats)

## Take out QC variables and calculate CVs
QCInd=grep('Q',sampNames)
QCFeats=AFeats[QCInd,]
QCCV=cv(QCFeats)

## Take out sample data for multivariate analysis
CorrData=AFeats[-QCInd,]
save(CorrData,file='LCMSDataCorrected.RData')
