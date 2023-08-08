#' Convert Cornforth style sigmoid parameters to TauMMF
#'
#' Converts original (Cornforth et al. 2012) sigmoid parametrization
#' to the inflection point \code{tauturn}, the distance between
#' \code{taumin} and \code{taumax} and the maximal marginal fitness
#' \code{mmf}.
#' @param a,b,beta,k,d Cornforth et al parameters
#' @note See \code{SigmoidParametrization.pdf} for derivation.
#' @examples
#' pDavidV1 <- convert_Cornforth_to_David(a=100, b=0.2, beta=4.76, k=10, d=1)
#' pTauMMF <- convert_Cornforth_to_tau_mmf(a=100, b=0.2, beta=4.76, k=10, d=1)
#' pDavidV2 <- convert_tau_mmf_to_David(tauturn=pTauMMF[["tauturn"]], dtau=pTauMMF[["dtau"]], mmf=pTauMMF[["mmf"]])
#' @references
#' \insertRef{Corn+12}{evolvr}
#' @export
convert_Cornforth_to_tau_mmf <-function(a, b, beta, k, d)
{
  s <- log(d)+k
  tauturn <- (s-log(beta))/b
  mmb <- a*b/(4*beta) #maximal marginal benefit
  mmf <- mmb -1
  mmbInv <- 1/mmb
  dtau <- (-1/b)*log(1-2*sqrt(1-mmbInv)/((1-mmbInv/2)+ sqrt(1-mmbInv)))

  return(list(tauturn=tauturn,dtau=dtau, mmf=mmf))
}

#' Convert sigmoid parametrization based on the inflection point
#' tauturn, the distance between taumin and taumax and the maximal
#' marginal fitness mmf, to David's parametrization.
#' @param tauturn inflection point
#' @param dtau distance between \code{taumin} and \code{taumax}
#' @param mmf maximal marginal fitness
#' @note See easierMS appendix A.3 for derivation.
#' @examples
#' pDavidV1 <- convert_Cornforth_to_David(a=100, b=0.2, beta=4.76, k=10, d=1)
#' pTauMMF <- convert_Cornforth_to_tau_mmf(a=100, b=0.2, beta=4.76, k=10, d=1)
#' pDavidV2 <- convert_tau_mmf_to_David(tauturn=pTauMMF[["tauturn"]], dtau=pTauMMF[["dtau"]], mmf=pTauMMF[["mmf"]])
#' @references
#' \insertRef{Corn+12}{evolvr}
#' @export
convert_tau_mmf_to_David<-function(tauturn, dtau, mmf)
{

  A <- (1/dtau)*log((2*mmf+1+2*sqrt(mmf*(mmf+1)))/(2*mmf+1-2*sqrt(mmf*(mmf+1))))
  L <- exp(A*tauturn)
  K <- 4*(mmf+1)/(A*(L+1))

  return(list(A=A, K=K, L=L))
}

#' Convert original (Cornforth et al. 2012) sigmoid parametrization to
#' David's parametrization.
#'
#' @inheritParams convert_Cornforth_to_tau_mmf
#' @note I didn't write this conversion up anywhere.
#' @examples
#' pDavidV1 <- convert_Cornforth_to_David(a=100, b=0.2, beta=4.76, k=10, d=1)
#' pTauMMF <- convert_Cornforth_to_tau_mmf(a=100, b=0.2, beta=4.76, k=10, d=1)
#' pDavidV2 <- convert_tau_mmf_to_David(tauturn=pTauMMF[["tauturn"]], dtau=pTauMMF[["dtau"]], mmf=pTauMMF[["mmf"]])
#' @references
#' \insertRef{Corn+12}{evolvr}
#' @export
convert_Cornforth_to_David <-function(a, b, beta, k, d)
{

  A <- b
  K <- a/(beta+d*exp(k))
  L<- a/(beta*K) - 1

  return(list(A=A, K=K, L=L))
}

#' Convert Chai-style sigmoid parameters to David's.
#'
#' Converts the inflection point \code{tauturn}, the curvature of the
#' benefit function at \code{taumax}, \code{benefitCurvatureAtTaumax},
#' and the maximal marginal fitness \code{mmf} to David's sigmoid
#' parametrization.
#' @param tauturn inflection point
#' @param benefitCurvatureAtTaumax the curvature of the benefit
#'     function at \code{taumax}
#' @param mmf maximal marginal fitness
#' @note See easierMS appendix A.3 for derivation.
#' @examples
#' pDavidV1 <- convert_Cornforth_to_David(a=100, b=0.2, beta=4.76, k=10, d=1)
#' pTauturnBencurvMMF <- convert_David_to_tauturn_bencurv_mmf(A=pDavidV1[["A"]], K=pDavidV1[["K"]], L = pDavidV1[["L"]])
#' pDavidV2 <- convert_tauturn_bencurv_mmf_to_David(tauturn=pTauturnBencurvMMF[["tauturn"]], benefitCurvatureAtTaumax=pTauturnBencurvMMF[["benefitCurvatureAtTaumax"]], mmf=pTauturnBencurvMMF[["mmf"]])
#' @references
#' \insertRef{Corn+12}{evolvr}
#' @export
convert_tauturn_bencurv_mmf_to_David <-function(tauturn, benefitCurvatureAtTaumax, mmf)
{
  A <- benefitCurvatureAtTaumax*(1+ 1/(mmf-sqrt(mmf*(mmf+1))))
  L <- exp(A*tauturn)
  K <- 4*(mmf+1)/(A*(L+1))
  return(list(A=A, K=K, L=L))
}


#' Converts David's sigmoid parametrization to Chai's.
#'
#' Converts David's sigmoid parametrization to the inflection point
#' \code{tauturn}, the curvature of the benefit function at
#' \code{taumax}, \code{benefitCurvatureAtTaumax}, and the maximal
#' marginal fitness \code{mmf}.
#' @param A,K,L David's parameters
#' @note See easierMS appendix A.3 for derivation.
#' @examples
#' pDavidV1 <- convert_Cornforth_to_David(a=100, b=0.2, beta=4.76, k=10, d=1)
#' pTauturnBencurvMMF <- convert_David_to_tauturn_bencurv_mmf(A=pDavidV1[["A"]], K=pDavidV1[["K"]], L = pDavidV1[["L"]])
#' pDavidV2 <- convert_tauturn_bencurv_mmf_to_David(tauturn=pTauturnBencurvMMF[["tauturn"]], benefitCurvatureAtTaumax=pTauturnBencurvMMF[["benefitCurvatureAtTaumax"]], mmf=pTauturnBencurvMMF[["mmf"]])
#' @references
#' \insertRef{Corn+12}{evolvr}
#' @export
convert_David_to_tauturn_bencurv_mmf <-function(A, K, L)
{
  mmf <- -1 + K*A*(1+L)/4
  tauturn <- log(L)/A
  benefitCurvatureAtTaumax <- A *(1- 1/(1+ mmf-sqrt(mmf*(mmf+1))))
  return(list(tauturn=tauturn, benefitCurvatureAtTaumax=benefitCurvatureAtTaumax, mmf=mmf))
}


#' Convert taumax-style sigmoid parameters to David's.
#'
#' Converts the maximizing total good \code{taumax}, the curvature of the
#' benefit function at \code{taumax}, \code{benefitCurvatureAtTaumax},
#' and the maximal marginal fitness \code{mmf} to David's sigmoid
#' parametrization.
#' @param taumax maximizing total good
#' @param benefitCurvatureAtTaumax the curvature of the benefit
#'     function at \code{taumax}
#' @param mmf maximal marginal fitness
#' @note See easierMS appendix A.3 for derivation.
#' @examples
#' pDavidV1 <- convert_Cornforth_to_David(a=100, b=0.2, beta=4.76, k=10, d=1)
#' pTaumax <- convert_David_to_taumax(A=pDavidV1[["A"]], K=pDavidV1[["K"]], L = pDavidV1[["L"]])
#' pDavidV2 <- convert_taumax_to_David(taumax=pTaumax[["taumax"]], benefitCurvatureAtTaumax=pTaumax[["benefitCurvatureAtTaumax"]], mmf=pTaumax[["mmf"]])
#' @references
#' \insertRef{Corn+12}{evolvr}
#' @export
convert_taumax_to_David <-function(taumax, benefitCurvatureAtTaumax, mmf)
{
  A <- benefitCurvatureAtTaumax*(1+ 1/(mmf-sqrt(mmf*(mmf+1))))
  L <- exp(A*taumax)*(2*mmf +1 -2*sqrt(mmf*(mmf+1)))
  K <- 4*(mmf+1)/(A*(L+1))
  return(list(A=A, K=K, L=L))
}


#' Converts David's sigmoid parametrization to Taumax-style.
#'
#' Converts David's sigmoid parametrization to the maximizing total good
#' \code{taumax}, the curvature of the benefit function at
#' \code{taumax}, \code{benefitCurvatureAtTaumax}, and the maximal
#' marginal fitness \code{mmf}.
#' @param A,K,L David's parameters
#' @note See easierMS appendix A.3 for derivation.
#' @examples
#' pDavidV1 <- convert_Cornforth_to_David(a=100, b=0.2, beta=4.76, k=10, d=1)
#' pTaumax <- convert_David_to_taumax(A=pDavidV1[["A"]], K=pDavidV1[["K"]], L = pDavidV1[["L"]])
#' pDavidV2 <- convert_taumax_to_David(taumax=pTaumax[["taumax"]], benefitCurvatureAtTaumax=pTaumax[["benefitCurvatureAtTaumax"]], mmf=pTaumax[["mmf"]])
#' @references
#' \insertRef{Corn+12}{evolvr}
#' @export
convert_David_to_taumax <-function(A, K, L)
{
  mmf <- -1 + K*A*(1+L)/4
  taumax <- log(L/(2*mmf +1 -2*sqrt(mmf*(mmf+1))))/A
  benefitCurvatureAtTaumax <- A *(1- 1/(1+ mmf-sqrt(mmf*(mmf+1))))
  return(list(taumax=taumax, benefitCurvatureAtTaumax=benefitCurvatureAtTaumax, mmf=mmf))
}


#' Convert DC-style sigmoid parameters to David's.
#'
#' Converts the maximizing total good \code{taumax}, the payoff extrema difference (PED),
#' \code{PED},
#' and the maximal marginal fitness \code{mmf} to David's sigmoid
#' parametrization.
#' @param taumax maximizing total good
#' @param ped the payoff extrema difference
#' @param mmf maximal marginal fitness
#' @note See easierMS appendix A.3 for derivation.
#' @examples
#' pDavidV1 <- convert_Cornforth_to_David(a=100, b=0.2, beta=4.76, k=10, d=1)
#' pDC <- convert_David_to_DC(A=pDavidV1[["A"]], K=pDavidV1[["K"]], L = pDavidV1[["L"]])
#' pDavidV2 <- convert_DC_to_David(taumax=pDC[["taumax"]], ped=pDC[["ped"]], mmf=pDC[["mmf"]])
#' @references
#' \insertRef{Corn+12}{evolvr}
#' @export
convert_DC_to_David <-function(taumax, ped, mmf)
{
  A <- (4*sqrt(mmf*(mmf+1))- log((2*mmf+1+2*sqrt(mmf*(mmf+1)))/(2*mmf+1-2*sqrt(mmf*(mmf+1)))))/ped
  L <- exp(A*taumax)*(2*mmf +1 -2*sqrt(mmf*(mmf+1)))
  K <- 4*(mmf+1)/(A*(L+1))
  return(list(A=A, K=K, L=L))
}


#' Converts David's sigmoid parametrization to DC-style.
#'
#' Converts David's sigmoid parametrization to the maximizing total good
#' \code{taumax}, the payoff extrema difference, \code{PED}, and the maximal
#' marginal fitness \code{mmf}.
#' @param A,K,L David's parameters
#' @note See easierMS appendix A.3 for derivation.
#' @examples
#' pDavidV1 <- convert_Cornforth_to_David(a=100, b=0.2, beta=4.76, k=10, d=1)
#' pDC <- convert_David_to_DC(A=pDavidV1[["A"]], K=pDavidV1[["K"]], L = pDavidV1[["L"]])
#' pDavidV2 <- convert_DC_to_David(taumax=pDC[["taumax"]], ped=pDC[["ped"]], mmf=pDC[["mmf"]])
#' @references
#' \insertRef{Corn+12}{evolvr}
#' @export
convert_David_to_DC <-function(A, K, L)
{
  mmf <- -1 + K*A*(1+L)/4
  taumax <- log(L/(2*mmf +1 -2*sqrt(mmf*(mmf+1))))/A
  ped <- (4*sqrt(mmf*(mmf+1))- log((2*mmf+1+2*sqrt(mmf*(mmf+1)))/(2*mmf+1-2*sqrt(mmf*(mmf+1)))))/A
  return(list(taumax=taumax, ped=ped, mmf=mmf))
}


#' Convert original (Cornforth et al. 2012) sigmoid parametrization to
#' parametrization using the asymptotic benefit \code{L},
#' the maximal marginal fitness \code{mmf}, and the inflection point \code{tauturn}.
#'
#' @inheritParams convert_Cornforth_to_David
#' @note I didn't write this conversion up anywhere.
#' @examples
#' pDavidV1 <- convert_Cornforth_to_David(a=100, b=0.2, beta=4.76, k=10, d=1)
#' pLTauturnMMFV1 <- convert_Cornforth_to_L_tauturn_mmf(a=100, b=0.2, beta=4.76, k=10, d=1)
#' pDavidV2 <- convert_L_tauturn_mmf_to_David(L=pLTauturnMMFV1[["L"]], tauturn=pLTauturnMMFV1[["tauturn"]], mmf=pLTauturnMMFV1[["mmf"]])
#' pLTauturnMMFV2 <- convert_David_to_L_tauturn_mmf(A=pDavidV1[["A"]], K=pDavidV1[["K"]], L = pDavidV1[["L"]])
#' @references
#' \insertRef{Corn+12}{evolvr}
#' @export
convert_Cornforth_to_L_tauturn_mmf <-function(a, b, beta, k, d)
{
  K <- a/(beta+d*exp(k))
  L<- a/(beta*K) - 1
  s <- log(d)+k
  tauturn <- (s-log(beta))/b
  mmb <- a*b/(4*beta) #maximal marginal benefit
  mmf <- mmb -1


  return(list(L=L, tauturn=tauturn, mmf=mmf))
}

#' Convert David's sigmoid parametrization to a parametrization based on the
#' asymptotic benefit \code{L}
#' inflection point \code{tauturn},
#' and the maximal marginal fitness \code{mmf}.
#'
#' @param A,K,L David's parameters
#' @examples
#' pDavidV1 <- convert_Cornforth_to_David(a=100, b=0.2, beta=4.76, k=10, d=1)
#' pLTauturnMMFV1 <- convert_Cornforth_to_L_tauturn_mmf(a=100, b=0.2, beta=4.76, k=10, d=1)
#' pDavidV2 <- convert_L_tauturn_mmf_to_David(L=pLTauturnMMFV1[["L"]], tauturn=pLTauturnMMFV1[["tauturn"]], mmf=pLTauturnMMFV1[["mmf"]])
#' pLTauturnMMFV2 <- convert_David_to_L_tauturn_mmf(A=pDavidV1[["A"]], K=pDavidV1[["K"]], L = pDavidV1[["L"]])
#' @references
#' \insertRef{Corn+12}{evolvr}
#' @export
convert_David_to_L_tauturn_mmf<-function(A, L, K)
{

  mmf <- -1 + K*A*(1+L)/4
  tauturn <- log(L)/A

  return(list(L=L, tauturn=tauturn, mmf=mmf))
}

#' Convert parametrization based on the
#' asymptotic benefit \code{L}
#' inflection point \code{tauturn},
#' and the maximal marginal fitness \code{mmf} to David's sigmoid parametrization.
#' @param L asymptotic benefit
#' @param tauturn inflection point
#' @param mmf maximal marginal fitness
#' @examples
#' pDavidV1 <- convert_Cornforth_to_David(a=100, b=0.2, beta=4.76, k=10, d=1)
#' pLTauturnMMFV1 <- convert_Cornforth_to_L_tauturn_mmf(a=100, b=0.2, beta=4.76, k=10, d=1)
#' pDavidV2 <- convert_L_tauturn_mmf_to_David(L=pLTauturnMMFV1[["L"]], tauturn=pLTauturnMMFV1[["tauturn"]], mmf=pLTauturnMMFV1[["mmf"]])
#' pLTauturnMMFV2 <- convert_David_to_L_tauturn_mmf(A=pDavidV1[["A"]], K=pDavidV1[["K"]], L = pDavidV1[["L"]])
#' @references
#' \insertRef{Corn+12}{evolvr}
#' @export
convert_L_tauturn_mmf_to_David<-function(L, tauturn, mmf)
{
  A <- log(L)/tauturn
  K <- 4*(mmf+1)/(A*(L+1))

  return(list(A=A, K=K, L=L))
}


#' Mean fitness slope (between extrema) as a function of \code{mmf}, for tanh-style sigmoid
#'
#' \eqn{\alpha} is the ratio between the difference between the
#' heights of the two fitness extrema and the difference between the
#' maximizing and minimizing total goods, calculated for tanh-style sigmoid
#' @param mmf maximal marginal fitness
#' @export
mean_fitness_slope_tanh <- function( mmf ) {
    out <- 4*sqrt(mmf^2+mmf)/log((2*mmf+1 + 2*sqrt(mmf^2+mmf))/(2*mmf+1 - 2*sqrt(mmf^2+mmf))) - 1
    return(out)
}

#' Mean fitness slope (between extrema) as a function of \code{mmf}, for erf-style sigmoid
#'
#' \eqn{\alpha} is the ratio between the difference between the
#' heights of the two fitness extrema and the difference between the
#' maximizing and minimizing total goods, calculated for erf-style sigmoid
#' @param mmf maximal marginal fitness
#' @export
mean_fitness_slope_erf <- function( mmf ) {
  slog<-sqrt(log(mmf+1))
  out <- (sqrt(pi)/2)*(mmf+1)*erf(slog)/slog-1
  return(out)
}

