library(MASS)
BIPpredict=function(dataListNew=dataListNew,Result=Result,meth="BMA"){
IndicVar=Result$IndicVar;nbrcomp=Result$nbrcomp;
X_new=list();
np=which(IndicVar==1)
Np=length(IndicVar);
j=1;
for (i in  1:Np){
X_new[[i]]=NULL;
if (IndicVar[i]!=1){ ## We do not center the response variable
X_new[[i]]=dataListNew[[j]];
X_new[[i]]=t((t(X_new[[i]])-Result$MeanData[[i]])/Result$SDData[[i]])
j=j+1;
}
}

nbrmodel=Result$nbrmodel;
EstLoad=Result$EstLoadModel
Upredict=matrix(0,nrow(X_new[[1]]),nbrcomp)
#Upredict=matrix(0,nrow(X_new[[1]]),ncol(Latent_new))
ypredict=rep(0,nrow(X_new[[1]]));
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
SelLoadPacked=cbind(SelLoadPacked,EstL[[m]][,selvar]);
Sig2=c(Sig2,Result$EstSig2[[m]][selvar])
SelVarXnew=cbind(SelVarXnew,X_new[[m]][,selvar])


#SelLoadPackedFull=cbind(SelLoadPackedFull,EstL[[m]][,selvar]);
#Sig2Full=c(Sig2Full,Result$EstSig2[[m]][selvar])
#SelVarXnewFull=cbind(SelVarXnewFull,X_new[[m]][,selvar])
}
}
#SelLoadPacked=as.matrix(SelLoadPacked);
#Sigma22=t(SelLoadPacked)%*%SelLoadPacked+diag(Sig2);
#AoA=t(Result$EstLoad[[Np]])%*%SelLoadPacked
#ypredict=(AoA%*%solve(Sigma22))%*%t(SelVarXnew)
#mse=mean((as.vector(ypredict)-as.vector(y_new))^2);

SelLoadPackedFull=SelLoadPacked;## I added this
Sig2Full=Sig2;
SelVarXnewFull=SelVarXnew;
SelLoadPackedFull=as.matrix(SelLoadPackedFull);
Sigma2Full=solve(SelLoadPackedFull%*%diag(1/Sig2Full)%*%t(SelLoadPackedFull)+diag(nrow(SelLoadPackedFull)));
Upredict1=t(Sigma2Full%*%SelLoadPackedFull%*%diag(1/Sig2Full)%*%t(SelVarXnewFull))
Upredict=Upredict+Upredict1*pb;
ypredict=ypredict+as.matrix(Upredict1%*%EstL[[np]])*pb;
}


#mse=mean((as.vector(ypredict)+Result$EstIntcp-as.vector(y_new))^2);
return (list(ypredict=as.vector(ypredict)+Result$EstIntcp,Upredtest=Upredict))
}

