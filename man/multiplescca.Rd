\name{multiplescca}
\alias{multiplescca}
\title{
Sparse canonical correlation vectors for fixed tuning paramters.
}
\description{
Obtain sparse canonical correlation vectors for fixed tuning parameters. It is recommended to use cvselpscca to choose optimal tuning paramters for each dataset, or use cvtunerange for range of tuning parameters.
}
\usage{
multiplescca(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,Tau=Tau,
            CovStructure="Iden",standardize=TRUE,maxiteration=20, thresh=0.0001)

}
\arguments{
  \item{Xdata1}{
A matrix of size \eqn{n \times p} for first dataset. Rows are samples and columns are variables.
}
 \item{Xdata2}{
A matrix of size \eqn{n \times q} for second dataset. Rows are samples and columns are variables.
}

 \item{ncancorr}{
Number of canonical correlation vectors. Default is one.
}
 \item{Tau}{
A vector or matrix of fixed tuning parameters for each dataset.
}

\item{CovStructure}{
Covariance structure to use in estimating sparse canonical correlation vectors. Either "Iden" or "Ridge". Iden assumes the  covariance matrix for each dataset is identity. Ridge uses the sample covariance for each dataset. See reference article for more details.
}

 \item{standardize}{
 TRUE or FALSE. If TRUE, data will be normalized to have mean zero  and variance one for each variable. Default is TRUE.
}
  \item{maxiteration}{
  Maximum iteration for the algorithm if not converged. Default is 20.
}

  \item{thresh}{
  Threshold for convergence. Default is 0.0001.
}
}
\details{
The function will return three R objects, which can be assigned to a variable. To see the results, use the ``$" operator.
}
\value{
	\item{hatalpha}{ Estimated sparse canonical correlation vectors for first dataset. }
	\item{hatbeta}{ Estimated sparse canonical correlation vectors for second dataset. }
		\item{maxcorr}{ Estimated correlation from canonical correlation vectors.}
}

\references{
Sandra E. Safo,  Jeongyoun Ahn, Yongho Jeon, and Sungkyu Jung (2018) , \emph{Sparse Generalized Eigenvalue Problem with Application
to Canonical Correlation Analysis for Integrative Analysis
of Methylation and Gene Expression Data}. \emph{Biometrics}%.\url{}
}

\seealso{
\code{\link{cvselpscca},\link{cvtunerange} }
}
\examples{
library(SELPCCA)
##---- read in data
data(DataExample)

Xdata1=DataExample[[1]]
Xdata2=DataExample[[2]]

##---- estimate first canonical correlation vectors
ncancorr=1

#use cvtunerange for range of tuning parameters
mytunerange=cvtunerange(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,
                        CovStructure="Iden",standardize=TRUE)
print(mytunerange)

#Fix Tau for first and second datasets as 1.1 and 1.0 respectively
Tau=matrix(c(1,1.2,1),nrow=1)
mysparsevectors=multiplescca(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,
                            Tau=Tau, CovStructure="Iden",standardize=TRUE,
                            maxiteration=20, thresh=0.0001)


#example with two canonical correlation vectors
#use cvselpscca to obtain optimal tuning parameters
mycv=cvselpscca(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,
                CovStructure="Iden",isParallel=FALSE,ncores=NULL,nfolds=5,
                ngrid=10, standardize=TRUE,thresh=0.0001,maxiteration=20)


Tau=mycv$optTau
mysparsevectors=multiplescca(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,
                Tau=Tau, CovStructure="Iden",standardize=TRUE,maxiteration=20,
                thresh=0.0001)


}
