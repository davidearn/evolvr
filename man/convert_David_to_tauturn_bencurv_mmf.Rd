% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_sigmoid.R
\name{convert_David_to_tauturn_bencurv_mmf}
\alias{convert_David_to_tauturn_bencurv_mmf}
\title{Converts David's sigmoid parametrization to Chai's.}
\usage{
convert_David_to_tauturn_bencurv_mmf(A, K, L)
}
\arguments{
\item{A, K, L}{David's parameters}
}
\description{
Converts David's sigmoid parametrization to the inflection point
\code{tauturn}, the curvature of the benefit function at
\code{taumax}, \code{benefitCurvatureAtTaumax}, and the maximal
marginal fitness \code{mmf}.
}
\note{
See easierMS appendix A.3 for derivation.
}
\examples{
pDavidV1 <- convert_Cornforth_to_David(a=100, b=0.2, beta=4.76, k=10, d=1)
pTauturnBencurvMMF <- convert_David_to_tauturn_bencurv_mmf(A=pDavidV1[["A"]], K=pDavidV1[["K"]], L = pDavidV1[["L"]])
pDavidV2 <- convert_tauturn_bencurv_mmf_to_David(tauturn=pTauturnBencurvMMF[["tauturn"]], benefitCurvatureAtTaumax=pTauturnBencurvMMF[["benefitCurvatureAtTaumax"]], mmf=pTauturnBencurvMMF[["mmf"]])
}
\references{
\insertRef{Corn+12}{evolvr}
}
