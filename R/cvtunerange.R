cvtunerange=function(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,CovStructure="Iden",standardize=TRUE){

  #set defaults
  #check data sizes

  p=dim(Xdata1)[2]
  n=dim(Xdata2)[1]
  q=dim(Xdata2)[2]

  if(dim(Xdata1)[1]!=dim(Xdata2)[1]){
    stop('The datasets  have different number of observations')
  }


  if(is.null(ncancorr)){
    ncancorr=1
  }

  if(is.null(CovStructure)){
    CovStructure="Iden"
  }

  if(is.null(standardize)){
    standardize=TRUE
  }

  if(standardize==TRUE){
    Xdata1=scale(Xdata1)
    Xdata2=scale(Xdata2)
  }


  myfastcca=myfastnonsparsecca(Xdata1,Xdata2,CovStructure)
  tildeA=myfastcca$tildeA
  tildeB=myfastcca$tildeB
  Ux=myfastcca$Ux
  Sigma12r=myfastcca$Sigma12r
  Uy=myfastcca$Uy

  mymaxlambdaalpha=matrix(rep(NA,ncancorr),nrow=ncancorr,ncol=1)
  mymaxlambdabeta=matrix(rep(NA,ncancorr),nrow=ncancorr,ncol=1)

  for(j in 1:ncancorr){
    mymaxlambdaalpha[j,1]= norm(Ux%*%Sigma12r%*%t(Uy)%*%tildeB[,j],"I")
    mymaxlambdabeta[j,1] = norm(Uy%*%t(Sigma12r)%*%t(Ux)%*%tildeA[,j],"I")
  }

  if(CovStructure=="Iden"){
    lb=mymaxlambdaalpha
    TauX1range=cbind(sqrt(log(p)/n)*lb, mymaxlambdaalpha/1.2)
    TauX2range=cbind(sqrt(log(q)/n)*lb, mymaxlambdabeta/1.2)
  }else if(CovStructure=="Ridge"){
    lb= mymaxlambdaalpha/3
    TauX1range=cbind(sqrt(log(p)/n)*lb, mymaxlambdaalpha/3)
    TauX2range=cbind(sqrt(log(q)/n)*lb, mymaxlambdabeta/3)
  }
  return(list(TauX1range=TauX1range,TauX2range=TauX2range))
}
