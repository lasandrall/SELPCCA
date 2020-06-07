multiplescca=function(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,Tau=Tau,CovStructure="Iden",standardize=TRUE,maxiteration=20,thresh=1e-04){

  if(dim(Xdata1)[1]!=dim(Xdata2)[1]) stop("Xdata1 and Xdata2 have different number of observations")

  if(is.null(ncancorr)){
    ncancorr=1
  }

  if(is.null(CovStructure)){
    CovStructure="Iden"
  }

  if(is.null(standardize)){
    standardize=TRUE
  }

  if(is.null(Tau)){
    stop('Tuning parameter (s) missing. You must provide one' )
  }

  if(standardize==TRUE){
    Xdata1=scale(Xdata1)
    Xdata2=scale(Xdata2)
  }

  if(is.null(maxiteration)){
    maxiteration=20
  }

  if(is.null(thresh)){
    thresh=1e-04
  }

  myfastcca=myfastnonsparsecca(Xdata1,Xdata2,CovStructure)
  tildealpha=myfastcca$tildeA
  tildebeta=myfastcca$tildeB
  tilderho=myfastcca$tilderho
  Ux=myfastcca$Ux
  Sigma12r=myfastcca$Sigma12r
  Uy=myfastcca$Uy

  if(max(Tau[,1])==1){
    maxiteration=maxiteration
    Tau=do.call(rbind, replicate(maxiteration, Tau, simplify=FALSE))
    Tau[,1]=matrix(rep(1:maxiteration,each=ncancorr),nrow= 1)
  }else{
    maxiteration=max(Tau[,1])
  }


  mybeta=as.matrix(tildebeta[,1:ncancorr])
  myalpha=as.matrix(tildealpha[,1:ncancorr])
  myrho=tilderho[1:ncancorr]
  iter=0
  diffalpha=1
  diffbeta=1
  reldiff=1
  while((iter<maxiteration) && min(reldiff, max(diffalpha,diffbeta))>thresh){
    iter=iter+1
    print(paste("Current Iteration Is:", iter))
    mybetaold=mybeta
    myalphaold=myalpha
    tilderhoold=myrho

    myselps=selpscca(Xdata1,Xdata2,mybetaold,myalphaold,tilderhoold,ncancorr,Ux,Sigma12r,Uy,as.matrix(Tau[Tau[,1]==iter,2]),as.matrix(Tau[Tau[,1]==iter,3]),CovStructure)
    myalpha=myselps$myhatalpha
    mybeta=myselps$myhatbeta
    myrho=myselps$mycorrmat

    if( min(colSums(abs(myalpha)))==0 || min(colSums(abs(mybeta)))==0){
      myalpha=myalphaold
      mybeta=mybetaold
      break
    }

    diffalpha=norm(myalpha-myalphaold,'f')^2 / norm(myalphaold,'f')^2
    diffbeta=norm(mybeta-mybetaold,'f')^2 / norm(mybetaold,'f')^2
    sumnew=norm(myalpha-myalphaold,'f')^2 + norm(mybeta-mybetaold,'f')^2
    sumold=norm(myalphaold,'f')^2+norm(mybetaold,'f')^2
    reldiff=sumnew/sumold
  }

  maxcorr=round(diag(abs(cor(Xdata1%*%myalpha,Xdata2%*%mybeta))),3)
  hatalpha=myalpha
  hatbeta=mybeta

  NonZeroAlpha=colSums(hatalpha!=0)
  NonZeroBeta=colSums(hatbeta!=0)
  print(paste(c("Number of non-zero alphas:", NonZeroAlpha),collapse=" "))
  print(paste(c("Number of non-zero betas:", NonZeroBeta),collapse=" "))
  print(paste(c("Corr(Xdata1*alpha,Xdata2*beta):", t(maxcorr)),collapse=" "))
  print(paste("Sparse CCA CovStructure used is:", CovStructure))

  return(list(estCorr=maxcorr,myalpha=hatalpha,mybeta=hatbeta))
}
