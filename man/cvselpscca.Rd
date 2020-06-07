\name{cvselpscca}
\alias{cvselpscca}
\title{
Cross validation for Sparse Canonical Correlation Analysis
}
\description{
Peforms nfolds cross validation  to select optimal tuning parameters for SELPCCA based on training data. If you want to apply optimal tuning parameters to testing data, you may also use multiplescca.
}
\usage{
cvselpscca(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,CovStructure="Iden",
          isParallel=TRUE,ncores=NULL,nfolds=5,ngrid=10,
          standardize=TRUE,thresh=0.0001,maxiteration=20)
}
\arguments{
  \item{Xdata1}{
A matrix of size \eqn{n \times p} for first dataset. Rows are samples and columns are variables.
}
 \item{Xdata2}{
A matrix of size \eqn{n \times q} for second dataset. Rows are samples and columns are variables.
}

 \item{ncancorr}{
Number of canonical correlation vectors. Default is 1.
}
\item{CovStructure}{
Covariance structure to use in estimating sparse canonical correlation vectors. Either "Iden" or "Ridge". Iden assumes the  covariance matrix for each dataset is identity. Ridge uses the sample covariance for each dataset. See reference article for more details.
}

 \item{isParallel}{
TRUE or FALSE for parallel computing. Default is TRUE.
}
\item{ncores}{
Number of cores to be used for parallel computing. Only used if isParallel=TRUE. If isParallel=TRUE and ncores=NULL, defaults to half the size of the number of system cores.
}

\item{nfolds}{
  Number of cross validation folds. Default is 5.
}
\item{ngrid}{
   Number of grid points for tuning parameters. Default is 10 for each dataset.
}
 \item{standardize}{
 TRUE or FALSE. If TRUE, data will be normalized to have mean zero  and variance one for each variable. Default is TRUE.
}
  \item{maxiteration}{
  Maximum iteration for the algorithm if not converged. Default is 20.
}

  \item{thresh }{
  Threshold for convergence. Default is 0.0001.
}
}
\details{
The function will return several R objects, which can be assigned to a variable. To see the results, use the ``$" operator.
}
\value{
	\item{hatalpha}{ Estimated sparse canonical correlation vectors for first dataset. }
	\item{hatbeta}{ Estimated sparse canonical correlation vectors for second dataset. }
\item{CovStructure}{Covariance structure used in estimating sparse canonical correlation vectors. Either "Iden" or "Ridge". }
	\item{optTau}{ Optimal tuning parameters for each dataset.}
		\item{maxcorr}{ Estimated canonical correlation coefficient.}
\item{tunerange}{ Grid values for each dataset used for searching optimal tuning paramters. }


}

\references{
Sandra E. Safo,  Jeongyoun Ahn, Yongho Jeon, and Sungkyu Jung (2018) , \emph{Sparse Generalized Eigenvalue Problem with Application
to Canonical Correlation Analysis for Integrative Analysis
of Methylation and Gene Expression Data}. \emph{Biometrics}%.\url{}
}

\seealso{
\code{\link{multiplescca} }
}
\examples{
library(SELPCCA)
##---- read in data
data(DataExample)

Xdata1=DataExample[[1]]
Xdata2=DataExample[[2]]


##---- call cross validation to estimate first canonical correlation vectors
ncancorr=1
mycv=cvselpscca(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,CovStructure="Iden",
                isParallel=FALSE,ncores=NULL,nfolds=5,ngrid=10,
                standardize=TRUE,thresh=0.0001,maxiteration=20)

#check output
train.correlation=mycv$maxcorr

optTau=mycv$optTau

hatalpha=mycv$hatalpha

hatbeta=mycv$hatbeta

#obtain correlation plot using training data
scoresX1=Xdata1\%*\% hatalpha
scoresX2=Xdata2\%*\% hatbeta
plot(scoresX1, scoresX2,lwd=3,
       ,xlab=paste(
         "First Canonical correlation variate for dataset", 1),
         ylab=paste("First Canonical correlation variate for dataset", 2),
         main=paste("Correlation plot for datasets",1, "and" ,2, ",", "\u03C1 =", mycv$maxcorr))


#obtain correlation plot using testing data

Xtestdata1=DataExample[[3]]
Xtestdata2=DataExample[[4]]
scoresX1=Xtestdata1\%*\%hatalpha
scoresX2=Xtestdata2\%*\%hatbeta
mytestcorr=round(abs(cor(Xtestdata1\%*\%hatalpha,Xtestdata2\%*\%hatbeta)),3)

plot(scoresX1, scoresX2,lwd=3,xlab=paste(
         "First Canonical correlation variate for dataset", 1),
         ylab=paste("First Canonical correlation variate for dataset", 2),
         main=paste("Correlation plot for datasets",1, "and" ,2, ",", "\u03C1 =", mytestcorr))

}
