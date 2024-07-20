library(MASS)
#dyn.load("~/projects/def-chekouo/chekouo/BayesianFA/CallfromR/BayesianFA.so")

BIP <- function(dataList=dataList,IndicVar=IndicVar, groupList=NULL,Method=Method,nbrcomp=4, sample=5000, burnin=1000,nbrmaxmodels=50,
		  priorcompselv=c(1,1),priorcompselo=c(1,1),priorb0=c(2,2),priorb=c(1,1),priorgrpsel=c(1,1),probvarsel=0.05,chainNbr=1) {
       
if (sample<=burnin){
stop("Argument burnin must be smaller than sample: the number of MCMC iterations.")
}
if (sample<=20){
stop("Please specify a larger number of MCMC iterations")
}

if (is.null(nbrcomp)){
D=which(IndicVar==0);
mysvd=lapply(D, function(i)  svd(dataList[[i]]))
mysumsvd=lapply(1:length(D), function(i) cumsum(mysvd[[i]]$d)/max(cumsum(mysvd[[i]]$d))*100)
KMax=max(unlist(lapply(1:length(D), function(i) min(which(mysumsvd[[i]]>=80, arr.ind = TRUE))))) #chooses maximum from D cumulative proportions
nbrcomp=min(KMax+1,10);
}


#Method="GroupInfo"
if (Method=="BIP"){
meth=0;
} else if (Method=="BIPnet") {
meth=1;
} else {
stop("You must provide either BIP or BIPnet")
}

Np=length(dataList);P=NULL
n=nrow(dataList[[1]])
P=NULL;
MeanData=list();SD=list();
for (i in 1:Np){
	dataList[[i]]=as.matrix(dataList[[i]])
	        P[i]=ncol(dataList[[i]]);
	if ((IndicVar[i]==0)||(IndicVar[i]==2)){
	MeanData[[i]]=apply(dataList[[i]],2,mean);
	SD[[i]]=apply(dataList[[i]],2,sd);
         dataList[[i]]=scale(as.matrix(dataList[[i]]),T,T)
	}
	dataList[[i]]=t(dataList[[i]]);
}
datasets=unlist(dataList)

if (is.null(groupList) && (meth==1)){
stop("You must provide a list of grouping information")
} else if (is.null(groupList)){
for (ll in 1:length(IndicVar)){
groupList[[ll]]=matrix(1,P[ll],1)
}
}


ll=1;
K=NULL;
for (i in 1:Np){
if (IndicVar[i]!=0){
K[i]=1;
} else {
K[i]=ncol(groupList[[ll]]);
ll=ll+1;
}
}
groupList=lapply(groupList,t);Paths=unlist(groupList);

result <- .C("mainfunction",
               Method1=as.integer(meth),n1=as.integer(n),P=as.integer(P),r1=as.integer(nbrcomp),Np1=as.integer(Np), datasets=as.double(datasets),IndVar=as.integer(IndicVar), K=as.integer(K), Paths=as.integer(Paths),maxmodel1=as.integer(nbrmaxmodels),                nbrsample1=as.integer(sample), burninsample1=as.integer(burnin),
		CompoSelMean=as.double(rep(0,Np*nbrcomp)),VarSelMean=as.double(rep(0,nbrcomp*sum(P))),
		VarSelMeanGlobal=as.double(rep(0,sum(P))),GrpSelMean=as.double(rep(0,nbrcomp*sum(K))),
	       GrpEffectMean=as.double(rep(0,nbrcomp*sum(K))),IntGrpMean=as.double(rep(0,nbrcomp*Np)),
	       EstU=as.double(rep(0,n*nbrcomp)),EstSig2=as.double(rep(0,sum(P))),InterceptMean=as.double(rep(0,1)),
	       EstLoadMod=as.double(rep(0,nbrmaxmodels*nbrcomp*sum(P))),EstLoad=as.double(rep(0,nbrcomp*sum(P))),
	       nbrmodel1=as.integer(0),postgam=rep(0,nbrmaxmodels),priorcompsel=priorcompselv,
	       priorcompselo=priorcompselo,priorb0=priorb0,priorb=as.double(priorb),priorgrpsel=priorgrpsel,
		probvarsel=as.double(probvarsel),chainNbr=as.integer(chainNbr));



reseffect=result$EstLoadMod;
nbrmodel=result$nbrmodel1;
EstLoadModel=rep( list(list()), nbrmodel ) 
for (mo in 1:nbrmodel){
	for (m in 1:Np){
		x=sum(P[1:(m-1)]);y=sum(P[1:m]);
		if (m==1) {x=0}
		init=1+x*nbrcomp+(mo-1)*nbrcomp*sum(P);
		final=y*nbrcomp+(mo-1)*nbrcomp*sum(P);
	EstLoadModel[[mo]][[m]]=matrix(reseffect[init:final],nbrcomp,byrow=T);
	}
}

resvarsel1=result$VarSelMeanGlobal;
resvarsel2=result$VarSelMean;
resvarsel3=result$GrpSelMean;;
resvarsel4=result$GrpEffectMean
resvarsel5=result$EstLoad;
EstimateU=matrix(result$EstU,n,byrow=T);
CompoSelMean=matrix(result$CompoSelMean,Np,byrow=T);
IntGrpMean=matrix(result$IntGrpMean,Np,byrow=T);
VarSelMeanGlobal=list();
VarSelMean=list();
GrpSelMean=list();
GrpEffectMean=list();
EstLoad=list();
EstSig2=list();
m1=m2=m3=1;
for (m in 1:Np){
VarSelMeanGlobal[[m]]=resvarsel1[m1:(m1-1+P[m])]
VarSelMean[[m]]=matrix(resvarsel2[m2:(m2-1+P[m]*nbrcomp)],P[m],byrow=T)
GrpSelMean[[m]]=matrix(resvarsel3[m3:(m3-1+K[m]*nbrcomp)],nbrcomp,byrow=T)
GrpEffectMean[[m]]=matrix(resvarsel4[m3:(m3-1+K[m]*nbrcomp)],nbrcomp,byrow=T);
EstLoad[[m]]=matrix(resvarsel5[m2:(m2-1+P[m]*nbrcomp)],nbrcomp,byrow=T);
EstSig2[[m]]=result$EstSig2[m1:(m1-1+P[m])]
m1=m1+P[m];
m2=m2+P[m]*nbrcomp;
m3=m3+K[m]*nbrcomp;

}

return (list(EstU=EstimateU,VarSelMean=VarSelMean,VarSelMeanGlobal=VarSelMeanGlobal,CompoSelMean=CompoSelMean,GrpSelMean=GrpSelMean,GrpEffectMean=GrpEffectMean,IntGrpMean=IntGrpMean,EstLoad=EstLoad,EstLoadModel=EstLoadModel,nbrmodel=result$nbrmodel1,EstSig2=EstSig2,EstIntcp=result$InterceptMean,PostGam=result$postgam,IndicVar=IndicVar,nbrcomp=nbrcomp,MeanData=MeanData,SDData=SD))
}

