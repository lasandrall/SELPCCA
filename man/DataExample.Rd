\docType{data}
\name{DataExample}
\alias{DataExample}
\title{
Simulated data with one true canonical correlation vectors.
}

\format{A list with 7 elements
\describe{
	\item{Xdata1}{ A matrix of size \eqn{80 \times 200} for first dataset. Rows are samples and columns are variables. }
\item{Xdata2}{ A matrix of size \eqn{80 \times 150} for second dataset. Rows are samples and columns are variables. }
\item{Xtestdata1}{ A matrix of size \eqn{400 \times 200} for first dataset. Rows are samples and columns are variables. }
\item{Xtestdata2}{ A matrix of size \eqn{400 \times 150} for second dataset. Rows are samples and columns are variables. }
\item{TrueAlpha}{ The first  canonical correlation vector for Xdata1. }
\item{TrueBeta}{ The first  canonical correlation vector for Xdata2. }
\item{TrueCorr}{The first canonical correlation coefficient. }
}}

\usage{
data(DataExample)
}
\description{
Simulated data with one true canonical correlation vectors for first and second datasets. The first 20 and 15 variables are nonzero (i.e., signal variables) in the first canonical correlation vectors for the first and second datasets respectively.
}
\examples{

#see example in multiplescca or cvselpscca

}
\references{
Sandra E. Safo,  Jeongyoun Ahn, Yongho Jeon, and Sungkyu Jung (2018) , \emph{Sparse Generalized Eigenvalue Problem with Application
to Canonical Correlation Analysis for Integrative Analysis
of Methylation and Gene Expression Data}. \emph{Biometrics}
}




