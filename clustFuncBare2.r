

library(xcms)
library(reshape)
library(mclust)
library(rgl)

## Make peaklist
peakTab=function(XS) {
	pTab=peakTable(XS)
	id=paste('RP',pTab$mz,'@',pTab$rt,sep='')
	X=t(pTab[,-1:-13])
	colnames(X)=id
	return(X)
}

## Bring out 'QC' group: Scaled and NAs removed
grabQC=function(XS,batch,grp='QC') {
	incl=(XS@phenoData[,1]==batch & XS@phenoData[,2]==grp)
	peakTab=peakTab(XS)
	QC=peakTab[incl,]
	QCCV=cv(QC)
	QCscale=scale(QC,center=FALSE)
	NAs=colSums(is.na(QCscale))>0
	QCRawNaRm=QC[,!NAs]
	QCFeats=QCscale[,!NAs]
	inj=as.numeric(substr(matrix(unlist(strsplit(rownames(QC),'_')),ncol=6,byrow=T)[,6],1,3))
	return(list(batch=batch,inj=inj,Feats=QCFeats,RawFeats=QC,RawFeatsNaRm=QCRawNaRm,NAs=NAs))
}

## Bring out Reference samples from same batch as QCs
grabRef=function(XS,QC,grp='Ref') {
	Feats=peakTab(XS)[XS@phenoData[,1]==QC$batch & XS@phenoData[,2]==grp,!QC$NAs]
	inj=as.numeric(substr(matrix(unlist(strsplit(rownames(Feats),'_')),ncol=6,byrow=T)[,6],1,3))
	return(list(inj=inj,Feats=Feats))
}

## Bring out entire batch same as QCs
grabBatch=function(XS,QC) {
	Feats=peakTab(XS)[XS@phenoData[,1]==QC$batch,!QC$NAs]
	inj=as.numeric(substr(matrix(unlist(strsplit(rownames(Feats),'_')),ncol=6,byrow=T)[,6],1,3))
	return(list(inj=inj,Feats=Feats))
}

## Simple function for calculating cv per column (ie variable)
cv=function(mat) {
	mean=apply(mat,2,mean)
	sd=apply(mat,2,sd)
	return(cv=sd/mean)
}
	
## Simple function for calculating root mean square distance
rmsDist=function(mat) {
	mean=colMeans(mat)
	rmsd=sqrt(sum(apply(mat,1,function(x) sum((x-mean)^2)))/nrow(mat))
	return(rmsd)
}

## Perform clustering
clust=function(QCInjs,QCFeats,modelNames=c('VVE'),G=seq(1,100,by=3),report=FALSE) {
	# modelNames='VVV'
	# modelNames=c('VEV','VVV')
	# modelNames=c('VEE','VEV','VVE','VVV')
	# modelNames=c('VEE','VVE')
	# modelNames=c('VII','VEI','VVI','VEE','VEV','VVE','VVV')
	# modelNames=c('EEE','EEV','EVE','EVV','VEE','VEV','VVE','VVV')
	startTime=proc.time()[3]
	mclBIC=mclustBIC(t(QCFeats),G=G,modelNames=modelNames)
	endTime=proc.time()[3]
	BICtime=endTime-startTime
	startTime=proc.time()[3]
	MC=summary(mclBIC,data=t(QCFeats))
	endTime=proc.time()[3]
	sumTime=endTime-startTime
	if (report==TRUE) {
		pdf(file=paste('cluster_BIC_',format(Sys.time(),format="%y%m%d_%H%M"),'.pdf',sep=''))
		plot(mclBIC)
		dev.off()
	}
	return(list(QCInjs=QCInjs,QCFeats=QCFeats,BIC=mclBIC,clust=MC,BICTime=BICtime,clustTime=sumTime))
}

## Calculate drift clusters
driftCalc=function(QCClust,smoothFunc=c('spline','loess'),spar=0.2,report=FALSE) {
	if (missing(smoothFunc)) smoothFunc='spline'
	MC=QCClust$clust
	QCInjs=QCClust$QCInjs
	QCFeats=QCClust$QCFeats
	# Extract classes
		nclass=MC$G # Take out total number of identified clusters/components/groups/classes/whatever you want to call them
		classes=MC$classification # Take out the classifications for the different variables
	# Allocate variables
		cvRaw=cvCorr=deltaDist=numeric(nclass) # allocate vector for effect size of drift correction per cluster
		injs=min(QCInjs):max(QCInjs) # Make injection list
		corMat=matrix(nrow=length(injs),ncol=nclass) # Allocate matrix with correction function (rows) per cluster (column)
		cvs=varClust=list()
		ratios=matrix(nrow=nclass,ncol=4)
		rownames(ratios)=paste('cluster',1:nclass,sep='')
		colnames(ratios)=c('raw.15','corr.15','raw.2','corr.2')
	if (report==TRUE) {
		pdf(file=paste('cluster_G',nclass,'_',format(Sys.time(),format="%y%m%d_%H%M"),'.pdf',sep=''))
		par(mfrow=c(2,1))
		par(mar=c(2,4,2,0))
	}
	### Calculate drift correction for each cluster
	# Calculate distance on scaled variables
	rmsdRaw=rmsDist(QCFeats)
		for (n in 1:nclass) { 
			QCFeatsCorr=QCFeats # Allocate matrix (for drift corrected variables) for later QC distance calculation
			vars=QCFeats[,classes==n] # Take out cluster variables
			varClust[[n]]=colnames(vars)
			V=as.data.frame(cbind(QCInjs,vars)) # Arrange and rearrange data
			V=melt(V,id.vars='QCInjs')
			V=V[order(V$QCInjs),]
			# Interpolations
				if (length(QCInjs)<=3) {  # 2nd degree polynomial if <4 data points
					Fit=lm(value ~ poly(QCInjs,2),data=V)
					Pred=predict(Fit,data.frame(QCInjs=injs))
					Pred=data.frame(x=injs,y=Pred)
				} else {
					if (smoothFunc=='spline') {
						splineFit=smooth.spline(V$QCInjs,V$value,spar=spar) # Cubic spline regression otherwise
						Pred=predict(splineFit,injs)  # Predict drift over all injections
					} else {
						loessFit=loess(value~QCInjs,data=V,span=spar)
						Pred=predict(loessFit,data.frame(QCInjs=injs))
						Pred=data.frame(x=injs,y=Pred)
					}
				}
			corFact=Pred$y[1]/Pred$y  # Calculate correction factors for all injections
			corMat[,n]=corFact # Store cluster correction factors in "master" matrix
			corQC=corFact[QCInjs-min(QCInjs)+1]  # Bring out correction factors for QC samples specifically
			QCFeatsCorr[,classes==n]=QCFeats[,classes==n]*corQC  # correct drift within cluster
			## Calculate rmsDist
			rmsdCorr=rmsDist(QCFeatsCorr)
			deltaDist[n]=rmsdCorr-rmsdRaw
			# deltaDist[n]=mean(dist(QCFeatsCorr))-meanDistQCFeats # Calculate change in average QC distance
			cvRaw[n]=mean(cv(QCFeats[,classes==n]))
			cvCorr[n]=mean(cv(QCFeatsCorr[,classes==n]))
			cvs[[n]]=data.frame(Raw=cv(QCFeats[,classes==n]),Corr=cv(QCFeatsCorr[,classes==n]))
			ratios[n,]=c(sum(cvs[[n]]$Raw<0.15)/nrow(cvs[[n]]),sum(cvs[[n]]$Corr<0.15)/nrow(cvs[[n]]),sum(cvs[[n]]$Raw<0.2)/nrow(cvs[[n]]),sum(cvs[[n]]$Corr<0.2)/nrow(cvs[[n]]))
			if (report==TRUE) {
				# Plot drift and drift function
				matplot(QCInjs,QCFeats[,classes==n],type='l',lty=1,col='grey',ylim=range(QCFeats[,classes==n]),main=paste('Cluster ',n,'; n=',sum(classes==n),'; Raw; Mean CV=',round(mean(cv(QCFeats[,classes==n])),3),sep=''),ylab='Scaled intensity',xlab='Injection number')
				lines(Pred,pch=2)
				matplot(QCInjs,QCFeatsCorr[,classes==n],type='l',lty=1,col='grey',ylim=range(QCFeats[,classes==n]),main=paste('Corrected; Mean CV=',round(mean(cv(QCFeatsCorr[,classes==n])),3),sep=''),ylab='Scaled intensity',xlab='Injection number')
			}
		}
	if (report==TRUE) dev.off() # Close pdf file
	clustComm=rep('None',nclass)
	actionInfo=data.frame(number=1:nclass,n=sapply(varClust,length),action=clustComm,CVRaw=cvRaw,CVCorr=cvCorr)
	QCClust$actionInfo=actionInfo
	QCClust$ratios=ratios
	QCClust$corMat=corMat
	QCClust$deltaDist=deltaDist
	QCClust$varClust=varClust
	QCDriftCalc=QCClust
	return(QCDriftCalc)
}

## Remove bad clusters
cleanClust=function(QCDriftCalc,report=FALSE,keepRatio=.2) {
	MC=QCDriftCalc$clust
	# Extract classes
		nclass=MC$G # Take out total number of identified clusters/components/groups/classes/whatever you want to call them
		classes=MC$classification # Take out the classifications for the different variables
	QCFeats=QCDriftCalc$QCFeats
	ratios=QCDriftCalc$ratios
	clustComm=as.character(QCDriftCalc$actionInfo$action)
	removeClust=which(ratios[,4]<keepRatio)
	clustComm[removeClust]='Removed'
	keepClust=which(ratios[,4]>=keepRatio)
	removeFeat=which(classes%in%removeClust)
	removeFeats=colnames(QCFeats)[removeFeat]
	QCFeatsClean=QCFeats[,-removeFeat]
	QCDriftCalc$keepClust=keepClust
	QCDriftCalc$removeFeat=removeFeat
	QCDriftCalc$removeFeats=removeFeats
	QCDriftCalc$QCFeatsClean=QCFeatsClean
	QCDriftCalc$actionInfo$action=clustComm
	if (report==TRUE) {
		pdf(file=paste('Hist_ClustClean_',format(Sys.time(),format="%y%m%d_%H%M"),'.pdf',sep=''))
		hist(cv(QCFeats),100,col=rgb(0,0,0,1),main='Cluster cleanup',xlab='CV (feature)')
		hist(cv(QCFeatsClean),col=rgb(1,1,1,.5),add=T)
		legend('topright',legend=c('Raw','Clean'),fill=c(rgb(0,0,0,1),rgb(1,1,1,0.5)))
		dev.off()
	}
	QCClean=QCDriftCalc
	return(QCClean)
}

## Perform drift correction for clusters IF rmsdRef is improved 
driftCorr2=function(QCClean,refList=NA,refType=c('none','one','many'),CorrObj=NA,report=FALSE) {
	if (missing(refType)) refType='none'
	if (refType=='many') {
		cat('\nMultiple reference samples not yet implemented\n')
		break
	}
	deltaDist=QCClean$deltaDist
	varClust=QCClean$varClust
	keepClust=QCClean$keepClust
	removeFeats=QCClean$removeFeats
	corrQCTemp=corrQC=QCClean$QCFeatsClean
	injQC=QCClean$QCInjs
	corMat=QCClean$corMat
	clustComm=as.character(QCClean$actionInfo$action)
	ordDist=order(deltaDist)
	ordDist=ordDist[ordDist%in%keepClust]
	if (refType=='one') {
  	refClean=refList$Feats[,!colnames(refList$Feats)%in%removeFeats]
	  injRef=refList$inj
	  corrRefTemp=corrRef=refClean
	}
	if (missing(CorrObj)) {
		injTest=injQC
		corrTest=QCClean$QCFeatsClean
		CorrObj=list(inj=injTest,Feats=corrTest)
	} else {
		injTest=CorrObj$inj
		corrTest=CorrObj$Feats
		corrTest=corrTest[,!colnames(corrTest)%in%removeFeats]
	}
	for (i in 1:length(keepClust)) {
		n=ordDist[i]
		corFact=corMat[,n] # take out cluster correction factors from "master" matrix
		corrFeats=varClust[[n]]
		if (refType=='none') { #Scheme for (suboptimal) situation without Ref samples
		  corQC=corFact[injQC-min(injQC)+1]  # Bring out correction factors for QC samples specifically
		  corrQCTemp[,colnames(corrQC)%in%corrFeats]=corrQC[,colnames(corrQC)%in%corrFeats]*corQC
		  if (rmsDist(corrQCTemp)<rmsDist(corrQC)) {
		    clustComm[n]='Corr_QC'
		    corrQC=corrQCTemp
		    corQC=corFact[injQC-min(injQC)+1]  # Bring out correction factors for QC samples specifically
		    corrQC[,colnames(corrQC)%in%corrFeats]=corrQC[,colnames(corrQC)%in%corrFeats]*corQC
		    corTest=corFact[injTest-min(injQC)+1]  # Bring out correction factors for Test samples specifically
		    corrTest[,colnames(corrTest)%in%corrFeats]=corrTest[,colnames(corrTest)%in%corrFeats]*corTest
		  } else {
		    corrQCTemp=corrQC
		  }
		}
		if (refType=='one') { # Scheme for situation with multiple injections of singular non-QC Ref samples
		  corRef=corFact[injRef-min(injQC)+1]
		  corrRefTemp[,colnames(corrRefTemp)%in%corrFeats]=corrRefTemp[,colnames(corrRefTemp)%in%corrFeats]*corRef
  		if (rmsDist(corrRefTemp)<rmsDist(corrRef)) {
  			clustComm[n]='Corr_1Ref'
  			corrRef=corrRefTemp
  			corQC=corFact[injQC-min(injQC)+1]  # Bring out correction factors for QC samples specifically
  			corrQC[,colnames(corrQC)%in%corrFeats]=corrQC[,colnames(corrQC)%in%corrFeats]*corQC
  			corTest=corFact[injTest-min(injQC)+1]  # Bring out correction factors for Test samples specifically
  			corrTest[,colnames(corrTest)%in%corrFeats]=corrTest[,colnames(corrTest)%in%corrFeats]*corTest
  		} else {
  			corrRefTemp=corrRef
  		}
		}
	}
	if (report==TRUE) {
		pdf(file=paste('Hist_Corrected_',format(Sys.time(),format="%y%m%d_%H%M"),'.pdf',sep=''))
		hist(cv(QCClean$QCFeatsClean),30,col=rgb(0,0,0,1),main='Cluster correction',xlab='CV (feature)')
		hist(cv(corrQC),20,col=rgb(1,1,1,.5),add=T)
		legend('topright',legend=c('Clean','Corrected'),fill=c(rgb(0,0,0,1),rgb(1,1,1,0.5)))
		dev.off()		
	}
	QCClean$actionInfo$action=clustComm
	QCClean$QCFeatsCorr=corrQC
	QCClean$RefType=refType
	if (refType=='one') {
  	QCClean$RefInjs=injRef
  	QCClean$RefFeats=refList$Feats
  	QCClean$RefFeatsClean=refClean
  	QCClean$RefFeatsCorr=corrRef
	}
	QCClean$TestInjs=injTest
	QCClean$TestFeats=CorrObj$Feats
	# QCClean$TestFeatsClean=
	QCClean$TestFeatsCorr=corrTest
	QCCorr=QCClean
	return(QCCorr)
}

## Remove individual variables with CV>0.2
cleanVar2=function(QCCorr,report=FALSE) {
	QCFeats=QCCorr$QCFeats
	QCFeatsClean=QCCorr$QCFeatsClean
	QCFeatsCorr=QCCorr$QCFeatsCorr
	cvIndex=which(cv(QCFeatsCorr)>.2)
	QCFeatsFinal=QCFeatsCorr[,-cvIndex]
	if (QCCorr$RefType=='one') {
	  RefFeatsFinal=QCCorr$RefFeatsCorr[,-cvIndex]
	}
	TestFeatsFinal=QCCorr$TestFeatsCorr[,-cvIndex]
	finalVars=colnames(QCFeatsFinal)
	cvFeats=mean(cv(QCFeats))      # 0.2276
	cvFeatsClean=mean(cv(QCFeatsClean))    # 0.1455
	cvFeatsCorr=mean(cv(QCFeatsCorr))   # 0.1144
	cvFeatsFinal=mean(cv(QCFeatsFinal)) # 0.1070
	QCcvs=data.frame(cvFeats=cvFeats,cvFeatsClean=cvFeatsClean,cvFeatsCorr=cvFeatsCorr,cvFeatsFinal=cvFeatsFinal)
	if (report==TRUE) {
		pdf(file=paste('Hist_Final_',format(Sys.time(),format="%y%m%d_%H%M"),'.pdf',sep=''))
		hist(cv(QCFeatsClean),30,col=rgb(0,0,0,1),main='Cluster cleanup',xlab='CV (feature)')
		hist(cv(QCFeatsFinal),20,col=rgb(1,1,1,.5),add=T)
		legend('topright',legend=c('Clean','Final'),fill=c(rgb(0,0,0,1),rgb(1,1,1,0.5)))
		dev.off()		
	}
	if (QCCorr$RefType=='one') {
	  rmsdRef=rmsDist(QCCorr$RefFeats)       # 5800887
	  rmsdRefClean=rmsDist(QCCorr$RefFeatsClean)     # 3670441
	  rmsdRefCorr=rmsDist(QCCorr$RefFeatsCorr)   # 2264153
	  rmsdRefFinal=rmsDist(RefFeatsFinal) # 2209906
	  RefRMSD=data.frame(rmsdRef=rmsdRef,rmsdRefClean=rmsdRefClean,rmsdRefCorr=rmsdRefCorr,rmsdRefFinal=rmsdRefFinal)
	  QCCorr$RefFeatsFinal=RefFeatsFinal
	  QCCorr$RefRMSD=RefRMSD
	}
	# Write algorithm to remove "removed" variables from clusters
	varClust=QCCorr$varClust
	N=length(varClust)
	CVAfter=nFeatAfter=numeric(N)
	for (n in 1:N) {
		vars=varClust[[n]]
		clustVars=QCFeatsFinal[,finalVars%in%vars]
		CVAfter[n]=mean(cv(clustVars))
		nFeatAfter[n]=ncol(clustVars)
	}
	actionInfo=QCCorr$actionInfo
	actionInfo=cbind(actionInfo[,1:2],nFeatAfter,actionInfo[,3:5],CVAfter)
	colnames(actionInfo)[2:3]=c('nBefore','nAfter')
	QCCorr$actionInfo=actionInfo
	QCCorr$QCFeatsFinal=QCFeatsFinal
	QCCorr$TestFeatsFinal=TestFeatsFinal
	QCCorr$finalVars=finalVars
	QCCorr$QCcvs=QCcvs
	QCFinal=QCCorr
	return(QCFinal)
}

## Perform between batch correction if average intensity > threshold value
twoBatch=function(B1,B2,thresh=1000) {
	b1vars=B1$finalVars
	b2vars=B2$finalVars
	bbvars=b1vars[b1vars%in%b2vars]
	b1Incl=b1vars%in%bbvars
	b2Incl=b2vars%in%bbvars
	b1Ref=B1$RefFeatsFinal[,b1Incl]
	b1RefMean=colMeans(b1Ref)
	b2Ref=B2$RefFeatsFinal[,b2Incl]
	b2RefMean=colMeans(b2Ref)
	corrFact=b1RefMean/b2RefMean
	subThresh=ifelse(b1RefMean<thresh,1,0)
	subThresh=subThresh+ifelse(b2RefMean<thresh,1,0)
	subThresh=ifelse(subThresh>=1,1,0)
	corrFact=ifelse(subThresh==1,subThresh,corrFact)
	b2RefCorr=t(t(b2Ref)*corrFact)
	bbRef=rbind(b1Ref,b2RefCorr)
	# matplot(bbRef[,1:10],type='l')
	b1Feats=B1$TestFeatsFinal[,b1Incl]
	b2Feats=B2$TestFeatsFinal[,b2Incl]
	b2FeatsCorr=t(t(b2Feats)*corrFact)
	bbFeats=rbind(b1Feats,b2FeatsCorr)
	return(list(finalVars=bbvars,RefFeatsFinal=bbRef,TestFeatsFinal=bbFeats))
}


## Wrapper function for grabbing QCs, reference and entire batch samples from XCMS-set
grabWrap=function(XS,batch,QC='QC',Ref='Ref') {
	QCObj=grabQC(XS,batch=batch,grp=QC)
	RefObj=grabRef(XS,QCObj,grp=Ref)
	BatchObj=grabBatch(XS,QCObj)
	return(list(QC=QCObj,Ref=RefObj,Batch=BatchObj))
}

## Wrapper function for all drift subfunctions (clust, driftCalc, cleanClust, driftCorr, cleanVar)
driftWrap=function(grabObj,report=FALSE) {
	QCObj=grabObj$QC
	RefObj=grabObj$Ref
	CorrObj=grabObj$Batch
	if (min(RefObj$inj)<min(QCObj$inj) | min(CorrObj$inj)<min(QCObj$inj)) {
		cat('\nReference or test injection outside drift calculation region: Before first QCs injection.')
		cat('\nCalculation aborted\n')
		break
	}
	if (max(RefObj$inj)>max(QCObj$inj) | max(CorrObj$inj)>max(QCObj$inj)) {
		cat('\nReference or test injection outside drift calculation region: After last QCs injection.')
		cat('\nCalculation aborted\n')
		break
	}
	A=driftList=clust(QCObj$inj,QCObj$Feats,report=report)
	B=driftList=driftCalc(driftList,report=report)
	C=driftList=cleanClust(driftList,report=report)
	D=driftList=driftCorr(driftList,refList=RefObj,CorrObj=CorrObj,report=report)
	E=driftList=cleanVar(driftList,report=report)
}

## Wrapper function for between batch correction
batchWrap=function(batchList) {
	n=length(batchList)
	B1=batchList[[1]]
	for (i in 2:n) {
		B2=batchList[[i]]
		B1=twoBatch(B1,B2)
	}
	return(B1)
}
