## remember to load libraries 'doParallel' and 'foreach' to use multiple cores
# library(doParallel)
# library(foreach)
# library(mixOmics)
# library(pROC)

vectSamp=function(vect,n=4,sampLen) { # Pick 'n' random samples within vector 'vect'
	# sampLen is a vector of number of observations within sample
	# If sampLen is not specified it is automatically calculated
	j=length(vect)
	if (missing(sampLen)) { # Calculation of sampLen vector
		sampLen=numeric(n)
		oddAdd=numeric(n)
		if (j%%n!=0) oddAdd[1:(j%%n)]=1
		sampLen=floor(j/n)+oddAdd
	}
	if (length(sampLen)!=n) 
		stop("Mismatch number of samples")
	if (sum(sampLen)>j) 
		stop("Exceeds maximum samples")
	if (sum(sampLen)<j) 
		cat('\n Warning: Undersampling \n')
	vectSamp=sample(vect) # Randomise vector vect
	samp=list()
	vectInd=1 # Index in randomised vector
	for (i in 1:n) { # Divide randomised vector into samples
		samp[[i]]=sort(vectSamp[vectInd:(vectInd+sampLen[i]-1)])  # Make sample
		vectInd=vectInd+sampLen[i]  # Count up index
	}
	return(samp)
}

mplsDAPar=function(X,nRep=25,nSeg=7,nInner=6,comps=5,featRatio=0.75,metric=c('miss','auc'),resampInner=c(TRUE,FALSE),perm=F) {  # Perform rdCV MPLS DA 
# X is only the "positive" matrix, the negative is calculated automatically.
	# Allocate list for function return 
	if (missing(metric)) metric='auc'
	if (missing(resampInner)) resampInner=FALSE
	cat('\nMultilevel PLS-DA')
	cat('\n ',nRep,' repetitions. ',nSeg,' outer segments. ',nInner,' inner segments. Max ',comps,' components.',sep='')
	cat('\n Feature ratio=',featRatio,'. Quality metric=',metric,'. Resample inner segments=',resampInner,'.\n',sep='')
	modelReturn=list()
	start.time=proc.time()[3]
	# Make binary vector for contrasts
	if (perm==F) {
		Y=rep(c(1,-1),each=nrow(X))
	} else {
		Y=sample(c(1,-1),size=nrow(X),replace=T)
		Y=c(Y,-Y)
	}
	# Find number of features in X matrix
	nFeat0=ncol(X)
	# count number of iterations for feature elimination
	cnt=0
	nFeat=nFeat0
	feat=numeric()
	while (nFeat>1) {  
		cnt=cnt+1
		feat=c(feat,nFeat)
		nFeat=floor(featRatio*nFeat)
	}
	# Make indices of observations for outer CV
	outInd=1:nrow(X) 
	# Allocate matrices for final models (min & max features) VIP ranks averaged over repetitions
	VIPMin=VIPMid=VIPMax=matrix(nrow=nFeat0,ncol=3) # Allocate matrix for variable importance averaged over repetition
	rownames(VIPMin)=rownames(VIPMid)=rownames(VIPMax)=colnames(X)
	colnames(VIPMin)=colnames(VIPMid)=colnames(VIPMax)=c('Mean','SD','CV')
	# Allocate matrix for prediction of outer segments Y per repetition
	yPredMin=yPredMid=yPredMax=matrix(nrow=length(Y),ncol=nRep)
	colnames(yPredMin)=colnames(yPredMid)=colnames(yPredMax)=paste('Rep',1:nRep,sep='')
	yPredMinR=yPredMidR=yPredMaxR=numeric(length(Y))
	# Allocate response vectors and matrices for feat's, nComp and VIP ranks over repetitions
	featRepMin=featRepMid=featRepMax=nCompRepMin=nCompRepMid=nCompRepMax=missRep=numeric(nRep)
	names(featRepMin)=names(featRepMid)=names(featRepMax)=names(nCompRepMin)=names(nCompRepMid)=names(nCompRepMax)=names(missRep)=paste(rep('rep',nRep),1:nRep,sep='')
	VIPRepMin=VIPRepMid=VIPRepMax=matrix(data=nFeat0,nrow=nFeat0,ncol=nRep)
	rownames(VIPRepMin)=rownames(VIPRepMid)=rownames(VIPRepMax)=colnames(X)
	colnames(VIPRepMin)=colnames(VIPRepMid)=colnames(VIPRepMax)=paste(rep('rep',nRep),1:nRep,sep='')
	# Perform rdCV Repetitions 
	reps=foreach(r=1:nRep, .packages=c('mixOmics','pROC'), .export='vectSamp') %dopar% {
	# for (r in 1:nRep) {  # Loop for all repetitions
		# r=1 # For testing
		cat('\nRepetition ',r,' of ',nRep,sep='')
		outSamp=vectSamp(outInd,n=nSeg)  # Draw random samples within total set
		# Allocate response vectors and matrices for feat's, nComp and VIP ranks over outer segments
		featOutMin=featOutMid=featOutMax=nCompOutMin=nCompOutMid=nCompOutMax=numeric(nSeg)
		names(featOutMin)=names(featOutMid)=names(featOutMax)=names(nCompOutMin)=names(nCompOutMax)=paste(rep('outSeg',nSeg),1:nSeg,sep='')
		VIPOutMin=VIPOutMid=VIPOutMax=matrix(data=nFeat0,nrow=nFeat0,ncol=nSeg)
		rownames(VIPOutMin)=rownames(VIPOutMid)=rownames(VIPOutMax)=colnames(X)
		colnames(VIPOutMin)=colnames(VIPOutMid)=colnames(VIPOutMax)=paste(rep('outSeg',nSeg),1:nSeg,sep='')
		# Perform outer CV over nSeg segments
		for (i in 1:nSeg) {   # Create 'nSeg' models within repetition
			# i=1 # for testing
			cat('\n Segment ',i,' (features): ',sep='') # Counter
			testInd=outSamp[[i]] # Draw out segment = test set
			inInd=outInd[-testInd]  # The rest are for training and validation
			# Allocate response vectors and matrices for feat's, nComp and VIP ranks over inner segments
			aucIn=missIn=nCompIn=matrix(nrow=nInner,ncol=cnt)
			rownames(aucIn)=rownames(missIn)=rownames(nCompIn)=paste(rep('inSeg',nInner),1:nInner,sep='')
			colnames(aucIn)=colnames(missIn)=colnames(nCompIn)=feat
			VIPInner=array(data=nFeat0,dim=c(nFeat0,cnt,nInner))
			rownames(VIPInner)=colnames(X)
			colnames(VIPInner)=feat
			dimnames(VIPInner)[[3]]=paste(rep('inSeg',nInner),1:nInner,sep='')
			nFeat=nFeat0
			inSamp=vectSamp(inInd,n=nInner)  # Draw random samples within current 2CV segment
			incFeat=colnames(X)
			# Perform (inner) 1CV over successively fewer features
			for (count in 1:cnt) {  # Build models with successively fewer feature. Quality metric = number of missclassifications for Validation set
				# count=1 # for testing
				# count=count+1
				cat(nFeat)
				if (resampInner==TRUE) inSamp=vectSamp(inInd,n=nInner)  # Resample inner segments for each feature elimination step
				comp=ifelse(nFeat<comps,nFeat,comps)
				for (j in 1:nInner) {
					# j=1 # for testing
					cat('.') # Counter
					valInd=inSamp[[j]] # Draw out segment = validation set
					xVal=rbind(X[valInd,],-X[valInd,]) 
					xVal=subset(xVal,select=incFeat)
					yVal=Y[c(valInd,valInd+nrow(X))]
					trainInd=inInd[-match(valInd,inInd)] # Define Training segment
					xTrain=rbind(X[trainInd,],-X[trainInd,])
					xTrain=subset(xTrain,select=incFeat)
					yTrain=Y[c(trainInd,trainInd+nrow(X))]
					plsInner=pls(xTrain,yTrain,ncomp=comp,mode="classic")
					yValInner=predict(plsInner,newdata=xVal)$predict[,,]  # Store class prediction probabilities per holdout (test) segment and repetition			
					if (metric=='miss') {
						# cat(' miss',count)
						yClassInner=ifelse(yValInner>0,1,-1)
						misClass=apply((yVal-yClassInner)*yVal/2,2,sum,na.rm=T)
						missIn[j,count]=min(misClass)
						nCompIn[j,count]=which.min(misClass)
					} else {
						# cat(' auc',count)
						auc=apply(yValInner,2,function(x) roc(yVal,x)$auc)
						aucIn[j,count]=max(auc)
						nCompIn[j,count]=which.max(auc)
					}
					VIPInner[match(names(vip(plsInner)[,nCompIn[j,count]]),rownames(VIPInner)),count,j]=rank(-vip(plsInner)[,nCompIn[j,count]])
					# VIPInner[5,,1:6]
				}
				# Eliminate features for subsequent model
				nFeat=floor(featRatio*nFeat)
				## NB!!! Average VIP ranks over inner segments before feature elimination!!!
				VIPInAve=apply(VIPInner[,count,],1,mean)
				incFeat=names(VIPInAve[order(VIPInAve)])[1:nFeat]
			}
			if (metric=='miss') {
				minIndex=max(which(apply(t(apply(missIn,1,rank)),2,mean)==min(apply(t(apply(missIn,1,rank)),2,mean))))
				maxIndex=min(which(apply(t(apply(missIn,1,rank)),2,mean)==min(apply(t(apply(missIn,1,rank)),2,mean))))
			} else {
				minIndex=max(which(apply(t(apply(aucIn,1,rank)),2,mean)==max(apply(t(apply(aucIn,1,rank)),2,mean))))
				maxIndex=min(which(apply(t(apply(aucIn,1,rank)),2,mean)==max(apply(t(apply(aucIn,1,rank)),2,mean))))
			}
			# Find a middle index | Either arithmetic or geometric mean
			# midIndex=which.min(abs(feat-mean(c(feat[minIndex],feat[maxIndex]))))  # Arithmetic
			midIndex=which.min(abs(feat-exp(mean(log(c(feat[minIndex],feat[maxIndex]))))))  # Geometric
			# Per outer segment: Average inner loop features, nComp and VIP ranks 
			featOutMin[i]=feat[minIndex]
			featOutMid[i]=feat[midIndex]
			featOutMax[i]=feat[maxIndex]
			nCompOutMin[i]=round(mean(nCompIn[,minIndex]))
			nCompOutMid[i]=round(mean(nCompIn[,midIndex]))
			nCompOutMax[i]=round(mean(nCompIn[,maxIndex]))
			VIPOutMin[,i]=apply(VIPInner[,minIndex,],1,mean)
			VIPOutMid[,i]=apply(VIPInner[,midIndex,],1,mean)
			VIPOutMax[,i]=apply(VIPInner[,maxIndex,],1,mean)
			# Build outer model for min and max nComp and predict YTEST
			xTest=rbind(X[testInd,],-X[testInd,])
			xNotTest=rbind(X[inInd,],-X[inInd,])
			yTrain=Y[c(inInd,inInd+nrow(X))]
			incFeatMin=rownames(VIPOutMin)[rank(VIPOutMin[,i])<=featOutMin[i]]
			xTrainMin=subset(xNotTest,select=incFeatMin)
			plsOutMin=pls(xTrainMin,yTrain,ncomp=nCompOutMin[i],mode="classic")
			xTestMin=subset(xTest,select=incFeatMin)
			yPredMinR[c(testInd,testInd+nrow(X))]=predict(plsOutMin,newdata=xTestMin)$predict[,,nCompOutMin[i]]  # 	
			incFeatMid=rownames(VIPOutMid)[rank(VIPOutMid[,i])<=featOutMid[i]]
			xTrainMid=subset(xNotTest,select=incFeatMid)
			plsOutMid=pls(xTrainMid,yTrain,ncomp=nCompOutMid[i],mode="classic")
			xTestMid=subset(xTest,select=incFeatMid)
			yPredMidR[c(testInd,testInd+nrow(X))]=predict(plsOutMid,newdata=xTestMid)$predict[,,nCompOutMid[i]]  # 	
			incFeatMax=rownames(VIPOutMax)[rank(VIPOutMax[,i])<=featOutMax[i]]
			xTrainMax=subset(xNotTest,select=incFeatMax)
			plsOutMax=pls(xTrainMax,yTrain,ncomp=nCompOutMax[i],mode="classic")
			xTestMax=subset(xTest,select=incFeatMax)
			yPredMaxR[c(testInd,testInd+nrow(X))]=predict(plsOutMax,newdata=xTestMax)$predict[,,nCompOutMax[i]]  # 	
		}
		# Per repetition: Average outer loop features, nComp and VIP ranks 
		featRepMinR=round(mean(featOutMin))
		nCompRepMinR=round(mean(nCompOutMin))
		VIPRepMinR=apply(VIPOutMin,1,mean)
		featRepMidR=round(mean(featOutMid))
		nCompRepMidR=round(mean(nCompOutMid))
		VIPRepMidR=apply(VIPOutMid,1,mean)
		featRepMaxR=round(mean(featOutMax))
		nCompRepMaxR=round(mean(nCompOutMax))
		VIPRepMaxR=apply(VIPOutMax,1,mean)
		parReturn=list(yPredMin=yPredMinR,featRepMin=featRepMinR,nCompRepMin=nCompRepMinR,VIPRepMin=VIPRepMinR,
			yPredMid=yPredMidR,featRepMid=featRepMidR,nCompRepMid=nCompRepMidR,VIPRepMid=VIPRepMidR,
			yPredMax=yPredMaxR,featRepMax=featRepMaxR,nCompRepMax=nCompRepMaxR,VIPRepMax=VIPRepMaxR)
		return(parReturn)
	}
	for (r in 1:nRep) {
		yPredMin[,r]=reps[[r]]$yPredMin
		featRepMin[r]=reps[[r]]$featRepMin
		nCompRepMin[r]=reps[[r]]$nCompRepMin
		VIPRepMin[,r]=reps[[r]]$VIPRepMin
		yPredMid[,r]=reps[[r]]$yPredMid
		featRepMid[r]=reps[[r]]$featRepMid
		nCompRepMid[r]=reps[[r]]$nCompRepMid
		VIPRepMid[,r]=reps[[r]]$VIPRepMid
		yPredMax[,r]=reps[[r]]$yPredMax
		featRepMax[r]=reps[[r]]$featRepMax
		nCompRepMax[r]=reps[[r]]$nCompRepMax
		VIPRepMax[,r]=reps[[r]]$VIPRepMax
	}
	yMin=yMid=yMax=matrix(nrow=nrow(X),ncol=3)
	colnames(yMin)=colnames(yMid)=colnames(yMax)=c('Mean','SD','CV')
	yMin[,1]=apply(yPredMin,1,mean)[1:nrow(X)]
	yMin[,2]=apply(yPredMin,1,sd)[1:nrow(X)]
	yMin[,3]=abs(100*yMin[,2]/yMin[,1])
	yMid[,1]=apply(yPredMid,1,mean)[1:nrow(X)]
	yMid[,2]=apply(yPredMid,1,sd)[1:nrow(X)]
	yMid[,3]=abs(100*yMid[,2]/yMid[,1])
	yMax[,1]=apply(yPredMax,1,mean)[1:nrow(X)]
	yMax[,2]=apply(yPredMax,1,sd)[1:nrow(X)]
	yMax[,3]=abs(100*yMax[,2]/yMax[,1])
	# Classify predictions
	yClassMin=ifelse(yMin[,1]>0,1,-1)
	yClassMid=ifelse(yMid[,1]>0,1,-1)
	yClassMax=ifelse(yMax[,1]>0,1,-1)
	# Bind together all Y data
	yMat=cbind(yMin[,1],yMid[,1],yMax[,1],yClassMin,yClassMax,Y)
	colnames(yMat)[1:3]=c('yMin','yMid','yMax')
	# Calculate misclassifications
	missMin=sum((1-yClassMin)/2)
	missMid=sum((1-yClassMid)/2)
	missMax=sum((1-yClassMax)/2)
	# Average VIP ranks over repetitions
	VIPMin[,1]=apply(VIPRepMin,1,mean)
	VIPMin[,2]=apply(VIPRepMin,1,sd)
	VIPMin[,3]=100*abs(VIPMin[,2]/VIPMin[,1])
	VIPMid[,1]=apply(VIPRepMid,1,mean)
	VIPMid[,2]=apply(VIPRepMid,1,sd)
	VIPMid[,3]=100*abs(VIPMid[,2]/VIPMid[,1])
	VIPMax[,1]=apply(VIPRepMax,1,mean)
	VIPMax[,2]=apply(VIPRepMax,1,sd)
	VIPMax[,3]=100*abs(VIPMax[,2]/VIPMax[,1])
	# Average nComp over repetitions
	nCompMin=c(mean(nCompRepMin),sd(nCompRepMin),abs(100*sd(nCompRepMin)/mean(nCompRepMin)))
	nCompMid=c(mean(nCompRepMid),sd(nCompRepMid),abs(100*sd(nCompRepMid)/mean(nCompRepMid)))
	nCompMax=c(mean(nCompRepMax),sd(nCompRepMax),abs(100*sd(nCompRepMax)/mean(nCompRepMax)))
	names(nCompMin)=names(nCompMid)=names(nCompMax)=c('Mean','SD','CV')
	# Average nFeat over repetitions
	nFeatMin=c(mean(featRepMin),sd(featRepMin),abs(100*sd(featRepMin)/mean(featRepMin)))
	nFeatMid=c(mean(featRepMid),sd(featRepMid),abs(100*sd(featRepMid)/mean(featRepMid)))
	nFeatMax=c(mean(featRepMax),sd(featRepMax),abs(100*sd(featRepMax)/mean(featRepMax)))
	names(nFeatMin)=names(nFeatMid)=names(nFeatMax)=c('Mean','SD','CV')
	# Stop timer
	end.time=proc.time()[3]
	# Generate report-list
	modelReturn$Y=Y[1:(length(Y)/2)]
	modelReturn$nFeat=list(minModel=nFeatMin,midModel=nFeatMid,maxModel=nFeatMax)
	modelReturn$nComp=list(minModel=nCompMin,midModel=nCompMid,maxModel=nCompMax)
	modelReturn$nFeatPerRep=list(minModel=featRepMin,midModel=featRepMid,maxModel=featRepMax)
	modelReturn$nCompPerRep=list(minModel=nCompRepMin,midModel=nCompRepMid,maxModel=nCompRepMax)
	modelReturn$yPred=list(minModel=yMin,midModel=yMid,maxModel=yMax)
	modelReturn$yClass=list(minModel=yClassMin[1:nrow(X)],midModel=yClassMid[1:nrow(X)],maxModel=yClassMax[1:nrow(X)])
	modelReturn$misClass=list(minModel=missMin,midModel=missMid,maxModel=missMax)
	modelReturn$VIPRank=list(minModel=VIPMin,midModel=VIPMid,maxModel=VIPMax)
	modelReturn$calcMins=(end.time-start.time)/60
	cat('\n Elapsed time',(end.time-start.time)/60,'mins \n')
	return(modelReturn)
}
