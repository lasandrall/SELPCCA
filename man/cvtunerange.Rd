\name{cvtunerange}
\alias{cvtunerange}
\title{
Tuning parameter range
}
\description{
Obtain upper and lower bounds of tuning parameters for each canonical correlation vector. It is recommended to use cvselpscca to choose optimal tuning paramters for each dataset.
}
\usage{
cvtunerange(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,
            CovStructure="Iden",standardize=TRUE)

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

\item{CovStructure}{
Covariance structure to use in estimating sparse canonical correlation vectors. Either "Iden" or "Ridge". Iden assumes the  covariance matrix for each dataset is identity. Ridge uses the sample covariance for each dataset. See reference article for more details.
}

 \item{standardize}{
 TRUE or FALSE. If TRUE, data will be normalized to have mean zero  and variance one for each variable. Default is TRUE.
}

}
\details{
The function will return tuning ranges for sparse estimation of canonical correlation vectors. To see the results, use the ``$" operator.
}
\value{
	\item{TauX1range}{ A \eqn{ncancorr \times 2} matrix of upper and lower bounds of tuning parameters for each canonical correlation vector for first dataset. }
	\item{TauX2range}{A \eqn{ncancorr \times 2} matrix Upper and lower bounds of tuning parameters for each canonical correlation vector for second dataset.}
}

\references{
Sandra E. Safo,  Jeongyoun Ahn, Yongho Jeon, and Sungkyu Jung (2018) , \emph{Sparse Generalized Eigenvalue Problem with Application
to Canonical Correlation Analysis for Integrative Analysis
of Methylation and Gene Expression Data}. \emph{Biometrics}%.\url{}
}

\seealso{
\code{\link{cvselpscca}, \link{multiplescca} }
}
\examples{

#see example in multiplescca

}
