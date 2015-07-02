rm(list=ls())
setwd('C:/R/Crisp')

## Code to install some necessary packages
# source("http://bioconductor.org/biocLite.R")
# biocLite("xcms")
# biocLite("CAMERA")

## Load packages
library(xcms)
library(CAMERA)
require("XML")

## Load data
load(file='Data5.RData')     # LC-MS original data (XCMS-set object)
load(file='injTable.RData')  # Data on injection tables
# str(injTable)
# 'data.frame':   493 obs. of  13 variables:
 # $ Line              : num  1 2 6 7 8 9 10 11 12 13 ...
 # $ CheckToRun        : num  0 0 0 0 0 0 0 0 0 0 ...
 # $ Injections        : num  8 1 1 1 1 1 1 1 1 1 ...
 # $ Volume            : num  4 4 4 4 4 4 4 4 4 4 ...
 # $ Position          : chr  "BD2" "BD1" "BD1" "RA1" ...
 # $ SampleID          : chr  "QC" "QC" "QC" "A30" ...
 # $ DataPath          : chr  "rppos1" "rppos1" "rppos1" "rppos1" ...
 # $ ResultDatafile    : chr  "QC_BD2_08_6042.d" "QC_BD1_01_6043.d" "QC_BD1_01_6047.d" "A30_RA1_01_6048.d" ...
 # $ NumberInSequence  : num  8 9 13 14 15 16 17 18 19 20 ...
 # $ Method            : chr  "CaBr_RP_2.m" "CaBr_RP_2.m" "CaBr_RP_2.m" "CaBr_RP_2.m" ...
 # $ MS_Method         : chr  "CaBr_Pos_50to1200_1_SME.m" "CaBr_Pos_50to1200_1_SME.m" "CaBr_Pos_50to1200_1_SME.m" "CaBr_Pos_50to1200_1_SME.m" ...
 # $ MSProcessingMethod: chr  "CaBr_DA_ESPos_1.m" "CaBr_DA_ESPos_1.m" "CaBr_DA_ESPos_1.m" "CaBr_DA_ESPos_1.m" ...

## Function to take out peaktable from data
peakTab=function(XS) {
	pTab=peakTable(XS)
	id=paste('RP',pTab$mz,'@',pTab$rt,sep='')
	X=t(pTab[,-1:-12])
	colnames(X)=id
	return(X)
}

## Make peak table from XCMS object
PT=peakTab(Data5)
# Sort according to injection table order
PT=PT[order(match(rownames(PT),injTable$ResultDatafile)),]
# > identical(rownames(PT),injTable$ResultDatafile) [1] TRUE

### BATCH 1A
## Bring out batch 1a
B1a=PT[injTable$DataPath=='rppos1a',]
B1aMeta=injTable[injTable$DataPath=='rppos1a',]
## Bring out indices for QCs
QRInd=grep('Q',B1aMeta$SampleID)
table(B1aMeta[QRInd,5:6])
qcInd=grep('BD1',B1aMeta$Position)
minIndex=min(qcInd)
maxIndex=max(qcInd)
# refInd=QRInd[!QRInd%in%qcInd]
## Grab QCs
inj=B1aMeta$NumberInSequence[qcInd]
Feats=B1a[qcInd,]
FeatsScaled=scale(Feats,center=FALSE)
NAs=colSums(is.na(FeatsScaled))>0
# sum(NAs)  =0 -> no missing data or zero variance!
QCObj=list(inj=inj,Feats=FeatsScaled,FeatsRaw=Feats)
## Grab entire batch
B1aCorrInd=CorrInd=minIndex:nrow(B1a)
CorrObj=list(inj=B1aMeta$NumberInSequence[CorrInd],Feats=B1a[CorrInd,])
## Perform Within Batch Correction
## Drift functions from clustFuncBare2.r
report=TRUE
A=clust(QCObj$inj,QCObj$Feats,report=report)
B=driftCalc(A,spar=0.2,smoothFunc='spline',report=report)
C=cleanClust(B,report=report)
CM=C$corMat
CMadd=matrix(rep(CM[47,],11),nrow=11,byrow=T)
CM2=rbind(CM,CMadd)
C$corMat=CM2
D2=driftCorr2(C,CorrObj=CorrObj,report=report)
E2=cleanVar2(D2,report=report)
B1aFeats=E2
save(B1a,B1aMeta,B1aFeats,file='B1a.Rdata')

### BATCH 1B
## Bring out batch 1b
B1b=PT[injTable$DataPath=='rppos1b',]
B1b[2:3,]=B1b[3:2,]
B1bMeta=injTable[injTable$DataPath=='rppos1b',]
B1bMeta[2:3,-9]=B1bMeta[3:2,-9]
## Bring out indices for QCs and Refs
QRInd=grep('Q',B1bMeta$SampleID)
table(B1bMeta[QRInd,5:6])
qcInd=grep('BC7',B1bMeta$Position)
minIndex=min(qcInd)
maxIndex=max(qcInd)
## Grab QCs
B1bQCinj=inj=B1bMeta$NumberInSequence[qcInd]
Feats=B1b[qcInd,]
FeatsScaled=scale(Feats,center=FALSE)
NAs=colSums(is.na(FeatsScaled))>0
# sum(NAs)=2 -> missing data or zero variance!
QCObj=list(inj=inj,Feats=FeatsScaled[,NAs==0],FeatsRaw=Feats[,NAs==0])
## Grab entire batch
CorrInd=minIndex:maxIndex
CorrObj=list(inj=B1bMeta$NumberInSequence[CorrInd],Feats=B1b[CorrInd,NAs==0])
## Perform Within Batch Correction
## Drift functions from clustFuncBare2.r 
report=TRUE
A=clust(QCObj$inj,QCObj$Feats,report=report)
B=driftCalc(A,spar=0.2,smoothFunc='spline',report=report)
C=cleanClust(B,report=report)
D=driftCorr2(C,CorrObj=CorrObj,report=report)
E=cleanVar2(D,report=report)
B1bFeats=E
save(B1b,B1bMeta,B1bFeats,file='B1b.Rdata')

### BATCH 2
## Bring out batch 2
B2=PT[injTable$DataPath=='rppos2',]
B2Meta=injTable[injTable$DataPath=='rppos2',]
## Bring out indices for QCs and Refs
QRInd=grep('Q',B2Meta$SampleID)
## qcInd=grep('QCE',B2Meta$SampleID)
qcInd=grep('BE7',B2Meta$Position)
minIndex=min(qcInd)
maxIndex=max(qcInd)
## Grab QCs
inj=B2Meta$NumberInSequence[qcInd]
Feats=B2[qcInd,]
FeatsScaled=scale(Feats,center=FALSE)
NAs=colSums(is.na(FeatsScaled))>0
# sum(NAs)=1 -> no missing data or zero variance!
QCObj=list(inj=inj,Feats=FeatsScaled[,NAs==0],FeatsRaw=Feats[,NAs==0])
## Grab entire batch
CorrInd=minIndex:maxIndex
CorrObj=list(inj=B2Meta$NumberInSequence[CorrInd],Feats=B2[CorrInd,NAs==0])
## Perform Within Batch Correction
## Drift functions from clustFuncBare2.r 
report=FALSE
A=clust(QCObj$inj,QCObj$Feats,report=report)
B=driftCalc(A,spar=0.2,smoothFunc='spline',report=report)
C=cleanClust(B,report=report)
D=driftCorr2(C,CorrObj=CorrObj,report=report)
E=cleanVar2(D,report=report)
B2Feats=E
save(B2,B2Meta,B2Feats,file='B2.Rdata')

### BATCH 3
## Bring out batch 3
B3=PT[injTable$DataPath=='rppos3',]
B3Meta=injTable[injTable$DataPath=='rppos3',]
## Bring out indices for QCs and Refs
QRInd=grep('Q',B3Meta$SampleID)
table(B3Meta[QRInd,5:6])
## qcInd=grep('QCE',B3Meta$SampleID)
qcInd=grep('BE7',B3Meta$Position)
minIndex=min(qcInd)
maxIndex=max(qcInd)
## Grab QCs
inj=B3Meta$NumberInSequence[qcInd]
Feats=B3[qcInd,]
FeatsScaled=scale(Feats,center=FALSE)
NAs=colSums(is.na(FeatsScaled))>0
# sum(NAs)=1 -> no missing data or zero variance!
QCObj=list(inj=inj,Feats=FeatsScaled[,NAs==0],FeatsRaw=Feats[,NAs==0])
## Grab entire batch
CorrInd=minIndex:maxIndex
CorrObj=list(inj=B3Meta$NumberInSequence[CorrInd],Feats=B3[CorrInd,NAs==0])
## Perform Within Batch Correction
## Drift functions from clustFuncBare2.r in QCData folder
report=FALSE
A=clust(QCObj$inj,QCObj$Feats,report=report)
B=driftCalc(A,spar=0.2,smoothFunc='spline',report=report)
C=cleanClust(B,report=report)
D=driftCorr2(C,CorrObj=CorrObj,report=report)
E=cleanVar2(D,report=report)
B3Feats=E
save(B3,B3Meta,B3Feats,file='B3.Rdata')

### BATCH 4
## Bring out batch 4
B4=PT[injTable$DataPath=='rppos4',]
B4Meta=injTable[injTable$DataPath=='rppos4',]
## Bring out indices for QCs and Refs
QRInd=grep('Q',B4Meta$SampleID)
table(B4Meta[QRInd,5:6])
## qcInd=grep('QCE',B4Meta$SampleID)
qcInd=grep('BE7',B4Meta$Position)
minIndex=min(qcInd)
maxIndex=max(qcInd)
## Grab QCs
inj=B4Meta$NumberInSequence[qcInd]
Feats=B4[qcInd,]
FeatsScaled=scale(Feats,center=FALSE)
NAs=colSums(is.na(FeatsScaled))>0
# sum(NAs)=1 -> no missing data or zero variance!
QCObj=list(inj=inj,Feats=FeatsScaled[,NAs==0],FeatsRaw=Feats[,NAs==0])
## Grab entire batch
B4CorrInd=CorrInd=minIndex:maxIndex
CorrObj=list(inj=B4Meta$NumberInSequence[CorrInd],Feats=B4[CorrInd,NAs==0])
## Perform Within Batch Correction
## Drift functions from clustFuncBare2.r 
report=FALSE
A=clust(QCObj$inj,QCObj$Feats,report=report)
B=driftCalc(A,spar=0.2,smoothFunc='spline',report=report)
C=cleanClust(B,report=report)
D=driftCorr2(C,CorrObj=CorrObj,report=report)
E=cleanVar2(D,report=report)
B4Feats=E
save(B4,B4Meta,B4Feats,file='B4.Rdata')

save(list=ls(pattern='B.')[substr(ls(pattern='B.'),1,1)=='B'],file='BatchData.Rdata')
load(file='BatchData.Rdata')
