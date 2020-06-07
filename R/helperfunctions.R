mysqrtminv=function(X){
  mysvd=svd(X)
  d=diag(mysvd$d^(-0.5))
  out=mysvd$u%*%d%*%t(mysvd$u)
  return(out)
}

minv=function(X){
  mysvd=svd(X)
  if(length(mysvd$d)==1){
    d=diag(as.matrix(mysvd$d^(-1)))
  }else{
    d=diag(mysvd$d^(-1))
  }
  out=mysvd$v%*%d%*%t(mysvd$u)
  return(out)
}

myfastnonsparsecca=function(Xdata1,Xdata2,CovStructure){

  n=dim(Xdata1)[1]
  p=dim(Xdata1)[1]
  q=dim(Xdata2)[2]

  if(CovStructure=="Iden"){
    SVD=svd(t(Xdata1))
    Ux=SVD$u
    Dx=diag(SVD$d)
    Vx=SVD$v
    Rx=Dx%*%t(Vx)
    cRx=Rx-rowMeans(Rx)

    SVD=svd(t(Xdata2))
    Uy=SVD$u
    Dy=diag(SVD$d)
    Vy=SVD$v
    Ry=Dy%*%t(Vy)
    cRy=Ry-rowMeans(Ry)

    Sigma12r=(cRx%*%t(cRy))/(n-1)
    tildeAr=svd(Sigma12r)$u
    tildeBr=svd(Sigma12r)$v

    tildeA=Ux%*%tildeAr
    tildeA=tildeA/sqrt(colSums(tildeA*tildeA))

    tildeB=Uy%*%tildeBr
    tildeB=tildeB/sqrt(colSums(tildeB*tildeB))
    tilderho=abs(diag(cor(Xdata1%*%tildeA, Xdata2%*%tildeB)))
  }else if(CovStructure=="Ridge"){
    SVD=svd(t(Xdata1))
    Ux=SVD$u
    Dx=diag(SVD$d)
    Vx=SVD$v
    Rx=Dx%*%t(Vx)
    cRx=Rx-rowMeans(Rx)

    SVD=svd(t(Xdata2))
    Uy=SVD$u
    Dy=diag(SVD$d)
    Vy=SVD$v
    Ry=Dy%*%t(Vy)
    cRy=Ry-rowMeans(Ry)

    Sigma12r=(cRx%*%t(cRy))/(n-1)
    lambda1=sqrt(log(p)/n)
    lambda2=sqrt(log(q)/n)

    Sigma11r=(cRx%*%t(cRx))/(n-1)
    Sigma22r=(cRy%*%t(cRy))/(n-1)

    AA1=Sigma11r +lambda1*diag(n)
    AA2=Sigma22r +lambda2*diag(n)

    out1=mysqrtminv(AA1);
    out2=mysqrtminv(AA2);

    tilde=svd(out1%*%Sigma12r%*%out2)
    tildeAr=tilde$u
    Dr=tilde$d
    tildeBr=tilde$v

    tilderho=Dr

    tildeA=Ux%*%out1%*%tildeAr
    tildeA=tildeA/sqrt(colSums(tildeA*tildeA))

    tildeB=Uy%*%out2%*%tildeBr
    tildeB=tildeB/sqrt(colSums(tildeB*tildeB))
  }
  return(list(tildeA=tildeA, tildeB=tildeB, tilderho=tilderho,Ux=Ux,Sigma12r=Sigma12r,Uy=Uy))
}




fastsvd=function(X){
  n=dim(X)[1]
  p=dim(X)[2]
  SVD=svd(t(X))
  U=SVD$u
  D=diag(SVD$d)
  V=SVD$v

  R=D%*%t(V)
  cR=R-rowMeans(R)
  Sigmar=(cR%*%t(cR))/(n-1)
  return(list(U=U,Sigmar=Sigmar))
}


selpscca=function(Xdata1=Xdata1,Xdata2=Xdata2,mybetaold, myalphaold, tilderhoold,ncancorr,Ux,Sigma12r,Uy,Taux,Tauy,CovStructure){
  myalphamat=list()
  mybetamat=list()
  mycorrmat=list()

  p=dim(Xdata1)[2]
  n=dim(Xdata2)[1]
  q=dim(Xdata2)[2]
  Xn=Xdata1
  Yn=Xdata2

  for(ii in 1:ncancorr){
    if(ii>1){
      alphamat=matrix(unlist(myalphamat), ncol=length(myalphamat),byrow =FALSE)
      betamat=matrix(unlist(mybetamat), ncol=length(mybetamat),byrow =FALSE)

      ProjmX=Xdata1%*%(as.matrix(alphamat[,1:(ii-1)])%*%minv(t(as.matrix(alphamat[,1:(ii-1)]))%*%as.matrix(alphamat[,1:(ii-1)])+0.001*diag(ii-1))%*%t(as.matrix(alphamat[,1:(ii-1)])))
      ProjmX[is.nan(ProjmX)]=0
      Xn=Xdata1-ProjmX

      ProjmY=Xdata2%*%(as.matrix(betamat[,1:(ii-1)])%*%minv(t(as.matrix(betamat[,1:(ii-1)]))%*%as.matrix(betamat[,1:(ii-1)])+0.001*diag(ii-1))%*%t(as.matrix(betamat[,1:(ii-1)])))
      ProjmY[is.nan(ProjmY)]=0
      Yn=Xdata2-ProjmY

      myfastcca=myfastnonsparsecca(Xn,Yn,CovStructure)
      Ux=myfastcca$Ux
      Sigma12r=myfastcca$Sigma12r
      Uy=myfastcca$Uy
    }

    Alphai=Variable(p)
    Objx=norm1(Alphai)
    if(CovStructure=="Iden") {constraint=list(norm_inf(Ux%*%Sigma12r%*%t(Uy)%*%mybetaold[,ii]- tilderhoold[ii]*Alphai) <= Taux[ii,1])
    }else if(CovStructure=="Ridge"){
      myfastsvd=fastsvd(Xn)
      UtrAlphai=t(myfastsvd$U)%*%Alphai
      constraint=list(norm_inf(Ux%*%Sigma12r%*%t(Uy)%*%mybetaold[,ii] - tilderhoold[ii]*(myfastsvd$U%*%myfastsvd$Sigmar%*%UtrAlphai + sqrt(log(p)/n)*Alphai)) <= Taux[ii,1])
    }
    prob=Problem(Minimize(Objx),constraint)
    result=solve(prob,solver="ECOS")
    alphai=result$getValue(Alphai)

    Betai=Variable(q)
    Objx=norm1(Betai)
    if(CovStructure=="Iden") {constraint=list(norm_inf(Uy%*%t(Sigma12r)%*%t(Ux)%*%myalphaold[,ii]- tilderhoold[ii]*Betai) <= Tauy[ii,1])
    }else if(CovStructure=="Ridge"){
      myfastsvd=fastsvd(Yn)
      UtrBetai=t(myfastsvd$U)%*%Betai
      constraint=list(norm_inf(Uy%*%t(Sigma12r)%*%t(Ux)%*%myalphaold[,ii] - tilderhoold[ii]*(myfastsvd$U%*%myfastsvd$Sigmar%*%UtrBetai + sqrt(log(q)/n)*Betai)) <= Tauy[ii,1])
    }
    prob=Problem(Minimize(Objx),constraint)
    result=solve(prob,solver="ECOS")
    betai=result$getValue(Betai)

    alphai[abs(alphai)<=10^(-5)]=0
    if(all(alphai==0) || sum(sum(is.nan(alphai)!=0))){
      myalpha=alphai
    }else{
      myalpha=alphai/norm(alphai,"2")
    }

    betai[abs(betai)<=10^(-5)]=0
    if(all(betai==0) || sum(sum(is.nan(betai)!=0))){
      mybeta=betai
    }else{
      mybeta=betai/norm(betai,"2")
    }

    if((all(myalpha==0) || sum(sum(is.nan(myalpha)!=0))) || (all(mybeta==0)|| sum(sum(is.nan(mybeta)!=0)))){
      myrho=0
    }else{
      myrho=abs(cor(Xdata1%*%myalpha, Xdata2%*%mybeta))
    }

    myalphamat[[ii]]=myalpha
    mybetamat[[ii]]=mybeta
    mycorrmat[[ii]]=myrho
  }

  myalphamat=matrix(unlist(myalphamat), ncol=length(myalphamat),byrow =FALSE)
  mybetamat=matrix(unlist(mybetamat), ncol=length(mybetamat),byrow=FALSE)
  mycorrmat=matrix(unlist(mycorrmat), ncol=1,byrow=FALSE)

  return(list(myhatalpha=myalphamat,myhatbeta=mybetamat,mycorrmat=mycorrmat))
}
