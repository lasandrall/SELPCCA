cvselpscca=function(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,CovStructure="Iden",isParallel=TRUE,ncores=NULL,nfolds=5,ngrid=10,standardize=TRUE,thresh=1e-04,maxiteration=20){

  if(is.null(Xdata2)){
    stop('There must be two datasets')
  }

  nX=dim(Xdata1)[1]
  nY=dim(Xdata2)[1]
  n=nX
  if(nX!=nY){
    stop('Xdata1 and Xdata2 have different number of observations')
  }
  #set defaults

  if(is.null(ncancorr)){
    ncancorr=1
  }

  if(is.null(CovStructure)){
    CovStructure="Iden"
  }
  if(is.null(isParallel)){
    isParallel=TRUE
  }

  if(is.null(nfolds)){
    nfolds=5
  }

  if(is.null(ngrid)){
    ngrid=10
  }


  if(is.null(maxiteration)){
    maxiteration=20
  }

  if(is.null(thresh)){
    thresh=1e-04
  }

  if(is.null(standardize)){
    standardize=TRUE
  }

  if(standardize==TRUE){
    Xdata1=scale(Xdata1)
    Xdata2=scale(Xdata2)
  }

  Tauxvec=list()
  Tauyvec=list()
  tunerange=cvtunerange(Xdata1,Xdata2,ncancorr,CovStructure)
  Tauxrange=tunerange$TauX1range
  Tauyrange=tunerange$TauX2range

  for(j in 1:ncancorr){
    Tauxvec[[j]]=10^(seq(log(Tauxrange[j,1])/log(10),log(Tauxrange[j,2])/log(10), (log(Tauxrange[j,2])/log(10)-log(Tauxrange[j,1])/log(10))/ngrid))
    Tauyvec[[j]]=10^(seq(log(Tauyrange[j,1])/log(10),log(Tauyrange[j,2])/log(10), (log(Tauyrange[j,2])/log(10)-log(Tauyrange[j,1])/log(10))/ngrid))
  }
  Tauxvec=matrix(unlist(Tauxvec),nrow=length(Tauxvec),byrow=TRUE)
  Tauyvec=matrix(unlist(Tauyvec),nrow=length(Tauyvec),byrow=TRUE)

  myoptTau=list()

  set.seed(1)
  if((n%%nfolds)>=1){foldid=sample(c(rep(1:nfolds,floor(n/nfolds)),1:(n%%nfolds)),n,replace=FALSE)}else
  {foldid=sample(rep(1:nfolds,floor(n/nfolds)),n,replace=FALSE)}


  diff_alpha=1
  diff_beta=1
  iter=0
  reldiff=1

  #ncores=c()
  while((iter<maxiteration) && min(reldiff, max(diff_alpha,diff_beta))> thresh){
    iter=iter+1
    print(paste("Current iteration is", iter))
    mycorr2=list()
    if(dim(Tauyvec)[1]==1){
      Tauy=t(apply(t(Tauyvec[,1:ngrid]),1,median))
    }else{Tauy=as.matrix(apply(Tauyvec[,1:ngrid],1,median))}

    if(isParallel==TRUE){ #if using parallel

      registerDoParallel()
      if(is.null(ncores)){
        ncores=parallel::detectCores()
        ncores=ceiling(ncores/2)}
      cl=makeCluster(ncores)
      registerDoParallel(cl)

      mycorr2=foreach(jj=1:(dim(Tauxvec)[2]-1),.export=c('myfastnonsparsecca','mysqrtminv','selpscca','minv','fastsvd'),.packages=c('CVXR')) %dopar%{
        mycorr=list()
        if(dim(Tauxvec)[1]==1){
          Taux=t(Tauxvec[,jj])
        }else{Taux=as.matrix(Tauxvec[,jj])}

        for(ii in 1:nfolds){
          which.row=foldid==ii
          myfastcca=myfastnonsparsecca(Xdata1[-which(which.row),],Xdata2[-which(which.row),],CovStructure)
          tildeA=myfastcca$tildeA
          tildeB=myfastcca$tildeB
          tilderho=myfastcca$tilderho
          Ux=myfastcca$Ux
          Sigma12r=myfastcca$Sigma12r
          Uy=myfastcca$Uy

          if(iter==1){
            myalphaold=as.matrix(tildeA[,1:ncancorr])
            mybetaold=as.matrix(tildeB[,1:ncancorr])
            tilderhoold=tilderho[1:ncancorr]
          }

          myselps=selpscca(Xdata1[-which(which.row),],Xdata2[-which(which.row),],mybetaold,myalphaold, tilderhoold,ncancorr,Ux,Sigma12r,Uy,Taux,Tauy,CovStructure)
          myhatalpha=myselps$myhatalpha
          myhatbeta=myselps$myhatbeta

          myUtest=Xdata1[which(which.row),]%*%myhatalpha
          myVtest=Xdata2[which(which.row),]%*%myhatbeta

          myUtrain=Xdata1[-which(which.row),]%*%myhatalpha
          myVtrain=Xdata2[-which(which.row),]%*%myhatbeta

          mycorrTrain=prod(abs(diag(cor(myUtrain,myVtrain))))
          mycorrTest=prod(abs(diag(cor(myUtest,myVtest))))

          mycorr[[ii]]=c(t(Taux), ii ,mycorrTrain, mycorrTest)
        }

        mycorr=matrix(unlist(mycorr),nrow=length(mycorr),byrow=TRUE)
        isnanTrain=is.nan(mycorr[,ncancorr+2])
        isnanTest=is.nan(mycorr[,ncancorr+3])
        return(c(t(Taux),mean(mycorr[!isnanTrain,ncancorr+2]), mean(mycorr[!isnanTest,ncancorr+3])))
      }
      mycorr2=matrix(unlist(mycorr2),nrow=length(mycorr2),byrow=TRUE)
      mydiffmean=abs(mycorr2[,ncancorr+1]-mycorr2[,ncancorr+2])
      if(length(na.omit(mydiffmean[mydiffmean!=0]))==0){#If all are zeros
        optTaux=Tauxvec[,1]
      }else{
        row=max(which((mydiffmean==min(na.omit(mydiffmean[mydiffmean!=0]))),arr.ind=TRUE))
        optTaux=Tauxvec[,row]
      }

      mycorr2=list()
      Taux=as.matrix(optTaux)
      mycorr2=foreach(jj=1:(dim(Tauyvec)[2]-1),.export=c('myfastnonsparsecca','mysqrtminv','selpscca','minv','fastsvd'),.packages=c('CVXR')) %dopar%{
        mycorr=list()
        if(dim(Tauyvec)[1]==1){
          Tauy=t(Tauyvec[,jj])
        }else{Tauy=as.matrix(Tauyvec[,jj])}

        for(ii in 1:nfolds){
          which.row=foldid==ii
          myfastcca=myfastnonsparsecca(Xdata1[-which(which.row),],Xdata2[-which(which.row),],CovStructure)
          tildeA=myfastcca$tildeA
          tildeB=myfastcca$tildeB
          tilderho=myfastcca$tilderho
          Ux=myfastcca$Ux
          Sigma12r=myfastcca$Sigma12r
          Uy=myfastcca$Uy

          if(iter==1){
            myalphaold=as.matrix(tildeA[,1:ncancorr])
            mybetaold=as.matrix(tildeB[,1:ncancorr])
            tilderhoold=tilderho[1:ncancorr]
          }

          myselps=selpscca(Xdata1[-which(which.row),],Xdata2[-which(which.row),],mybetaold,myalphaold, tilderhoold,ncancorr,Ux,Sigma12r,Uy,Taux,Tauy,CovStructure)
          myhatalpha=myselps$myhatalpha
          myhatbeta=myselps$myhatbeta

          myUtest=Xdata1[which(which.row),]%*%myhatalpha
          myVtest=Xdata2[which(which.row),]%*%myhatbeta

          myUtrain=Xdata1[-which(which.row),]%*%myhatalpha
          myVtrain=Xdata2[-which(which.row),]%*%myhatbeta

          mycorrTrain=prod(abs(diag(cor(myUtrain,myVtrain))))
          mycorrTest=prod(abs(diag(cor(myUtest,myVtest))))

          mycorr[[ii]]=c(t(Tauy), ii ,mycorrTrain, mycorrTest)
        }
        mycorr=matrix(unlist(mycorr),nrow=length(mycorr),byrow=TRUE)
        isnanTrain=is.nan(mycorr[,ncancorr+2])
        isnanTest=is.nan(mycorr[,ncancorr+3])
        return(c(t(Tauy),mean(mycorr[!isnanTrain,ncancorr+2]), mean(mycorr[!isnanTest,ncancorr+3])))
      }
      mycorr2=matrix(unlist(mycorr2),nrow=length(mycorr2),byrow=TRUE)
      mydiffmean=abs(mycorr2[,ncancorr+1]-mycorr2[,ncancorr+2])
      if(length(na.omit(mydiffmean[mydiffmean!=0]))==0){#If all are zeros
        optTauy=Tauyvec[,1]
      }else{
        row=max(which((mydiffmean==min(na.omit(mydiffmean[mydiffmean!=0]))),arr.ind=TRUE))
        optTauy=Tauyvec[,row]}

      #stop cluster
      stopCluster(cl)
    }else{ #if not using parallel
      for(jj in 1:(dim(Tauxvec)[2]-1)){
        mycorr=list()
        if(dim(Tauxvec)[1]==1){
          Taux=t(Tauxvec[,jj])
        }else{Taux=as.matrix(Tauxvec[,jj])}

        for(ii in 1:nfolds){
          which.row=foldid==ii
          myfastcca=myfastnonsparsecca(Xdata1[-which(which.row),],Xdata2[-which(which.row),],CovStructure)
          tildeA=myfastcca$tildeA
          tildeB=myfastcca$tildeB
          tilderho=myfastcca$tilderho
          Ux=myfastcca$Ux
          Sigma12r=myfastcca$Sigma12r
          Uy=myfastcca$Uy

          if(iter==1){
            myalphaold=as.matrix(tildeA[,1:ncancorr])
            mybetaold=as.matrix(tildeB[,1:ncancorr])
            tilderhoold=tilderho[1:ncancorr]
          }

          myselps=selpscca(Xdata1[-which(which.row),],Xdata2[-which(which.row),],mybetaold,myalphaold, tilderhoold,ncancorr,Ux,Sigma12r,Uy,Taux,Tauy,CovStructure)
          myhatalpha=myselps$myhatalpha
          myhatbeta=myselps$myhatbeta

          myUtest=Xdata1[which(which.row),]%*%myhatalpha
          myVtest=Xdata2[which(which.row),]%*%myhatbeta

          myUtrain=Xdata1[-which(which.row),]%*%myhatalpha
          myVtrain=Xdata2[-which(which.row),]%*%myhatbeta

          mycorrTrain=prod(abs(diag(cor(myUtrain,myVtrain))))
          mycorrTest=prod(abs(diag(cor(myUtest,myVtest))))

          mycorr[[ii]]=c(t(Taux), ii ,mycorrTrain, mycorrTest)
        }

        mycorr=matrix(unlist(mycorr),nrow=length(mycorr),byrow=TRUE)
        isnanTrain=is.nan(mycorr[,ncancorr+2])
        isnanTest=is.nan(mycorr[,ncancorr+3])
        mycorr2[[jj]]=c(t(Taux),mean(mycorr[!isnanTrain,ncancorr+2]), mean(mycorr[!isnanTest,ncancorr+3]))
      }
      mycorr2=matrix(unlist(mycorr2),nrow=length(mycorr2),byrow=TRUE)
      mydiffmean=abs(mycorr2[,ncancorr+1]-mycorr2[,ncancorr+2])
      if(length(na.omit(mydiffmean[mydiffmean!=0]))==0){#If all are zeros
        optTaux=Tauxvec[,1]
      }else{
        row=max(which((mydiffmean==min(na.omit(mydiffmean[mydiffmean!=0]))),arr.ind=TRUE))
        optTaux=Tauxvec[,row]
      }

      mycorr2=list()
      Taux=as.matrix(optTaux)
      for(jj in 1:(dim(Tauyvec)[2]-1)){
        mycorr=list()
        if(dim(Tauyvec)[1]==1){
          Tauy=t(Tauyvec[,jj])
        }else{as.matrix(Tauyvec[,jj])}

        for(ii in 1:nfolds){
          which.row=foldid==ii
          myfastcca=myfastnonsparsecca(Xdata1[-which(which.row),],Xdata2[-which(which.row),],CovStructure)
          tildeA=myfastcca$tildeA
          tildeB=myfastcca$tildeB
          tilderho=myfastcca$tilderho
          Ux=myfastcca$Ux
          Sigma12r=myfastcca$Sigma12r
          Uy=myfastcca$Uy

          if(iter==1){
            myalphaold=as.matrix(tildeA[,1:ncancorr])
            mybetaold=as.matrix(tildeB[,1:ncancorr])
            tilderhoold=tilderho[1:ncancorr]
          }

          myselps=selpscca(Xdata1[-which(which.row),],Xdata2[-which(which.row),],mybetaold,myalphaold, tilderhoold,ncancorr,Ux,Sigma12r,Uy,Taux,Tauy,CovStructure)
          myhatalpha=myselps$myhatalpha
          myhatbeta=myselps$myhatbeta

          myUtest=Xdata1[which(which.row),]%*%myhatalpha
          myVtest=Xdata2[which(which.row),]%*%myhatbeta

          myUtrain=Xdata1[-which(which.row),]%*%myhatalpha
          myVtrain=Xdata2[-which(which.row),]%*%myhatbeta

          mycorrTrain=prod(abs(diag(cor(myUtrain,myVtrain))))
          mycorrTest=prod(abs(diag(cor(myUtest,myVtest))))

          mycorr[[ii]]=c(t(Tauy), ii ,mycorrTrain, mycorrTest)
        }
        mycorr=matrix(unlist(mycorr),nrow=length(mycorr),byrow=TRUE)
        isnanTrain=is.nan(mycorr[,ncancorr+2])
        isnanTest=is.nan(mycorr[,ncancorr+3])
        mycorr2[[jj]]=c(t(Tauy),mean(mycorr[!isnanTrain,ncancorr+2]), mean(mycorr[!isnanTest,ncancorr+3]))
      }
      mycorr2=matrix(unlist(mycorr2),nrow=length(mycorr2),byrow=TRUE)
      mydiffmean=abs(mycorr2[,ncancorr+1]-mycorr2[,ncancorr+2])
      if(length(na.omit(mydiffmean[mydiffmean!=0]))==0){#If all are zeros
        optTauy=Tauyvec[,1]
      }else{
        row=max(which((mydiffmean==min(na.omit(mydiffmean[mydiffmean!=0]))),arr.ind=TRUE))
        optTauy=Tauyvec[,row]}
    }

    Tauy=as.matrix(optTauy)
    myfastcca=myfastnonsparsecca(Xdata1 ,Xdata2,CovStructure)
    tildeA=myfastcca$tildeA
    tildeB=myfastcca$tildeB
    tilderho=myfastcca$tilderho

    myalphaold=as.matrix(tildeA[,1:ncancorr])
    mybetaold=as.matrix(tildeB[,1:ncancorr])
    tilderhoold=tilderho[1:ncancorr]

    if(iter>1){
      myalphaold=myalpha
      mybetaold=mybeta
      tilderhoold=mytilderho

      myalpha2old=myalpha2
      mybeta2old=mybeta2
    }

    myselpscca= selpscca(Xdata1,Xdata2,mybetaold, myalphaold, tilderhoold,ncancorr,myfastcca$Ux,myfastcca$Sigma12r,myfastcca$Uy,Taux,Tauy,CovStructure)
    myoptTau[[iter]]=cbind(iter*matrix(1,nrow=ncancorr,ncol=1),as.matrix(optTaux),as.matrix(optTauy))

    myalpha=myselpscca$myhatalpha
    mybeta=myselpscca$myhatbeta
    mytilderho=myselpscca$mycorrmat

    myalpha2=myselpscca$myalphaicon
    mybeta2=myselpscca$myalphaicon

    if( min(colSums(abs(myalpha)))==0 || min(colSums(abs(mybeta)))==0){
      myalpha=myalphaold
      mybeta=mybetaold
      break
    }

    myalpha2=myselpscca$myalphaicon
    mybeta2=myselpscca$myalphaicon

    diff_alpha=norm(myalpha-myalphaold,'f')^2 / norm(myalphaold,'f')^2 #
    diff_beta=norm(mybeta-mybetaold,'f')^2 / norm(mybetaold,'f')^2 #
    sumnew=norm(myalpha-myalphaold,'f')^2 + norm(mybeta-mybetaold,'f')^2 #
    sumold=norm(myalphaold,'f')^2+norm(mybetaold,'f')^2 #
    reldiff=sumnew/sumold #
  }



  print("Applying optimal tuning parameter on whole data")

  for(i in 1:length(myoptTau)){
    if(i==1){temp=matrix(unlist(myoptTau[[i]]),nrow=ncancorr,byrow=FALSE)}else{
      temp=rbind(temp,matrix(unlist(myoptTau[[i]]),nrow=ncancorr,byrow=FALSE))}}

  myoptTau=temp

  mymulticca=multiplescca(Xdata1,Xdata2,ncancorr,myoptTau,CovStructure,standardize,maxiteration,thresh)

  return(list(optTau=myoptTau,hatalpha=mymulticca$myalpha,
              hatbeta=mymulticca$mybeta,maxcorr=mymulticca$estCorr,CovStructure=CovStructure,
              tunerange=matrix(c(Tauxvec,Tauyvec),nrow = 2, byrow=TRUE)))
}
