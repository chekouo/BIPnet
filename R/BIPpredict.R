BIPpredict=function(dataListNew=dataListNew,Result=Result,meth="BMA"){
meth=match.arg(meth,c("BMA","noBMA")) ## stops on anything else, e.g. "bma"
IndicVar=Result$IndicVar;nbrcomp=Result$nbrcomp;
X_new=list();
np=which(IndicVar==1)
Np=length(IndicVar);
## one matrix per data type used in BIP() except the response, in the same order
if (is.matrix(dataListNew) || is.data.frame(dataListNew) || !is.list(dataListNew)) dataListNew=list(dataListNew)
if (length(dataListNew)!=sum(IndicVar!=1)){
stop(sprintf("dataListNew must contain %d data sets: the data types used in BIP() except the response, in the same order.",sum(IndicVar!=1)))
}
j=1;
for (i in  1:Np){
X_new[[i]]=NULL;
if (IndicVar[i]!=1){ ## We do not center the response variable
X_new[[i]]=dataListNew[[j]];
if (is.null(dim(X_new[[i]]))) X_new[[i]]=matrix(X_new[[i]],nrow=1) ## a vector is one test sample
X_new[[i]]=as.matrix(X_new[[i]])
if (ncol(X_new[[i]])!=length(Result$MeanData[[i]])){
stop(sprintf("Data set %d of dataListNew has %d columns, but data type %d had %d columns in BIP().",j,ncol(X_new[[i]]),i,length(Result$MeanData[[i]])))
}
X_new[[i]]=t((t(X_new[[i]])-Result$MeanData[[i]])/Result$SDData[[i]])
j=j+1;
}
}

nbrmodel=Result$nbrmodel;
EstLoad=Result$EstLoadModel
ntest=nrow(X_new[[which(IndicVar!=1)[1]]]) ## from the first data type that is not the response, wherever the response is
if (any(sapply(X_new[IndicVar!=1],nrow)!=ntest)){
stop("All data sets in dataListNew must have the same number of rows (samples).")
}
Upredict=matrix(0,ntest,nbrcomp)
#Upredict=matrix(0,nrow(X_new[[1]]),ncol(Latent_new))
ypredict=rep(0,ntest);
nbrmodel1=nbrmodel;
if (meth!="BMA"){
nbrmodel1=1;
}

for (nb in 1:nbrmodel1){
SelLoadPacked=NULL;SelLoadPackedFull=NULL
SelVarXnew=NULL;SelVarXnewFull=NULL;
Sig2=NULL;Sig2Full=NULL;
if (meth=="BMA"){
EstL=EstLoad[[nb]];
pb=Result$PostGam[nb];
} else {
EstL=Result$EstLoad; pb=1;
}
for (m in 1:Np){

#if (m<=Np-1){
if (IndicVar[m]!=1){ # We exclude the response variable
nc=nrow(EstL[[m]]);
selvar=apply(EstL[[m]],2, function(x) length((which(x==0)))!=nc)
SelLoadPacked=cbind(SelLoadPacked,EstL[[m]][,selvar,drop=FALSE]); # stays nbrcomp x #selected, also when nbrcomp=1
Sig2=c(Sig2,Result$EstSig2[[m]][selvar])
SelVarXnew=cbind(SelVarXnew,X_new[[m]][,selvar,drop=FALSE]) # drop=FALSE keeps a single test sample as a 1-row matrix



}
}


SelLoadPackedFull=SelLoadPacked;## I added this
Sig2Full=Sig2;
SelVarXnewFull=SelVarXnew;
SelLoadPackedFull=as.matrix(SelLoadPackedFull);
## A %*% diag(1/Sig2) by scaling the columns of A: diag() of a single value builds an identity matrix
## of that size instead of a 1x1 matrix, and the p x p diagonal matrix is avoided
LoadInvSig2=SelLoadPackedFull*rep(1/Sig2Full,each=nrow(SelLoadPackedFull));
Sigma2Full=solve(LoadInvSig2%*%t(SelLoadPackedFull)+diag(nrow(SelLoadPackedFull)));
Upredict1=t(Sigma2Full%*%LoadInvSig2%*%t(SelVarXnewFull))
Upredict=Upredict+Upredict1*pb;
ypredict=ypredict+as.matrix(Upredict1%*%EstL[[np]])*pb;
}


#mse=mean((as.vector(ypredict)+Result$EstIntcp-as.vector(y_new))^2);
return (list(ypredict=as.vector(ypredict)+Result$EstIntcp,Upredtest=Upredict))
}

