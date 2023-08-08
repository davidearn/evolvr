#' Set up snowdrift game with sigmoidal benefit function
#'
#' This function is called by \code{\link{create_snowdrift_game}} as
#' one of the basic setup options.  The benefit function is designed
#' using \code{\link[=gerf]{generalized error functions}}.
#' The cost function is quadratic in principle, \deqn{C(x) = c_1 x +
#' c_2 x^2 ,} but the default parameter values yield \eqn{C(x) = x}
#' and hence a \strong{\emph{natural snowdrift game}} (NSG) in the
#' sense of Molina and Earn (2018).
#'
#' @details The returned \code{benefitString} is a character string
#'     giving the full expression for the fitness \emph{benefit} as a
#'     function of the total group contribution \code{tau}.
#'     Similarly, the returned \code{costString} is a character string
#'     giving the full expression for the fitness \emph{cost} as a
#'     function of the focal individual's contribution \code{x}.
#'     These strings are used in \code{\link{create_snowdrift_game}}
#'     to construct functions and their derivatives.
#'
#' @param L benefit limit: the limit of \eqn{B(\tau)} as \eqn{\tau}
#'     approaches \eqn{\infty}
#' @param tauTurn inflection point
#' @param mmf maximal marginal fitness
#' @param nben the index of the generalized error function is \code{2*nben}
#' (FIX: think of a better name than \code{nben} for this...)
#' @param c1 coefficient of linear term of cost function
#' @param c2 coefficient of quadratic term of cost function
#' @return list containing benefit and cost functions (as strings) and
#'     their parameters as lists
#' @seealso \code{\link{gerf}}
#' @export
setup_gerf_snowdrift_game <-
    function(mmf=1, tauTurn=1, L=1, nben=1,
             c1=1, c2=0) {

        benefitString <- "L * gerf( (mmf+1)*gsl_sf_gamma(1/(2*nben))/(2*nben*L) * (tau - tauTurn), 2*nben )"
        ## Since gerf is odd (zero at tauTurn), and we need a
        ## non-negative benefit function, subtract value at 0:
        benefitString <- paste0(
            benefitString,
            " - ",
            gsub("tau ", "0 ", benefitString, fixed=TRUE))
        benefitParms <- list(mmf=mmf, L=L, tauTurn=tauTurn, nben=nben)
        costString <- "c1*x + c2*x^2"
        costParms <- list(c1=c1, c2=c2)

        ## FIX: compute derived parameters...  taumin, taumax, ...
        ## ... perhaps better to do this in a separate function that
        ## is called after setting up any snowdrift game, though for
        ## the natural snowdrift game we have convenient analytical
        ## formulae...
        benefitParmsDerived <- c()
        costParmsDerived <- c()

        benefitParms <- c(benefitParms, benefitParmsDerived)
        costParms <- c(costParms, costParmsDerived)

        return(list(benefitString=benefitString, costString=costString,
                    benefitParms=benefitParms, costParms=costParms))
}

#' Set up snowdrift game with hyperbolic tangent benefit function
#'
#' This function is called by \code{\link{create_snowdrift_game}} as
#' one of the basic setup options.  The benefit function is taken to
#' be that of Cornforth et al (2012), but parameterized following
#' Molina and Earn (2018), \deqn{B(\tau) = K [ (1+L)/(1+L
#' e^{-A\tau}) - 1 ] .}  The cost function is quadratic in
#' principle, \deqn{C(x) = c_1 x + c_2 x^2 ,} but the default
#' parameter values yield \eqn{C(x) = x} and hence a
#' \strong{\emph{natural snowdrift game}} (NSG) in the sense of Molina
#' and Earn (2018).
#'
#' @details The returned \code{benefitString} is a character string
#'     giving the full expression for the fitness \emph{benefit} as a
#'     function of the total group contribution \code{tau}.
#'     Similarly, the returned \code{costString} is a character string
#'     giving the full expression for the fitness \emph{cost} as a
#'     function of the focal individual's contribution \code{x}.
#'     These strings are used in \code{\link{create_snowdrift_game}}
#'     to construct functions and their derivatives.
#'
#' @param K benefit scale factor: converts benefit to units of
#'     fitness
#' @param L benefit limit: the limit of \eqn{B(\tau)/K} as \eqn{\tau}
#'     approaches \eqn{\infty}
#' @param A benefit exponential rate parameter
#' @param c1 coefficient of linear term of cost function
#' @param c2 coefficient of quadratic term of cost function
#' @return list containing benefit and cost functions (as strings) and
#'     their parameters as lists
#' @references
#' \insertRef{Corn+12}{evolvr}
#' @export
setup_tanh_snowdrift_game <-
    function(K=0.01, L=1000, A=1,
             c1=1, c2=0) {

        benefitString <- "K*( (1+L)/(1+L*exp(-A*tau)) - 1 )"
        benefitParms <- list(K=K, L=L, A=A)
        costString <- "c1*x + c2*x^2"
        costParms <- list(c1=c1, c2=c2)

        ## FIX: compute derived parameters...  taumin, taumax, ...
        ## ... perhaps better to do this in a separate function that
        ## is called after setting up any snowdrift game, though for
        ## the natural snowdrift game we have convenient analytical
        ## formulae...
        benefitParmsDerived <- c()
        costParmsDerived <- c()

        benefitParms <- c(benefitParms, benefitParmsDerived)
        costParms <- c(costParms, costParmsDerived)

        return(list(benefitString=benefitString, costString=costString,
                    benefitParms=benefitParms, costParms=costParms))
}

#' Quadratic cost coefficient
#'
#' @inheritParams setup_quadratic_snowdrift_game
#' @seealso \code{\link{setup_quadratic_snowdrift_game}}
#' @export
c2_fun <- function(Npop, nPlayers, nu, xi) {
  stopifnot((!is.null(xi))&(xi>nPlayers))
  stopifnot((nu >= 1)&(nu<2))
  ((xi-nPlayers)/(xi-1))*(1 + nu*(nPlayers-1)/xi)
}

#' Set up \code{snowdrift.game} fitness function based on quadratic
#' cost and benefit functions
#'
#' This is the class of snowdrift games used in the Invasion-Fixation
#' paper to illustrate the results.
#'
#' @param Npop population size
#' @param nPlayers number of players
#' @param nu non-negative real number that indexes the examples
#' @param b1,b2,c1,c2 coefficients of linear and quadratic terms in
#'     benefit and cost functions (defaults are chosen to yield the
#'     class of examples used in the Invasion-Fixation paper)
#' @return list containing benefit and cost functions (as strings) and
#'     their parameters as lists
#' @note The values of the arguments that are used to define `c2` by
#'     default (i.e., `N`, `nGroups`, `nu` and `xi`) are saved as attributes of
#'     the `costParms` list.
#' @seealso \code{\link{c2_fun}}
#' @export
setup_quadratic_snowdrift_game <-
    function(Npop=10, nPlayers=2, nu=3/2, xi=Npop,
             b1=0,
             b2=1,
             c1=2,
             c2=c2_fun(Npop, nPlayers, nu,xi)
             ) {
        stopifnot((nu >= 1)&(nu<2))
        stopifnot(xi>nPlayers)

        benefitString <- "b1*tau + b2*tau^2"
        benefitParms <- list(b1=b1, b2=b2)
        costString <- "c1*x + c2*x^2"
        costParms <- list(c1=c1, c2=c2)
        ## save the parameters that yielded the value of c2:
        attr(costParms[["c2"]],"Npop") <- Npop
        attr(costParms[["c2"]],"nPlayers") <- nPlayers
        attr(costParms[["c2"]],"nu") <- nu
        attr(costParms[["c2"]],"xi") <- xi

        return(list(benefitString=benefitString, costString=costString,
                    benefitParms=benefitParms, costParms=costParms))
    }




#' Create a \code{snowdrift.game} model
#'
#' A \strong{\emph{snowdrift game}} is a public goods game in which
#' the fitness of a focal individual can be written \deqn{W(x,\tau) =
#' B(\tau) - C(x) ,} where \eqn{x} is the focal individual's
#' contribution and \eqn{\tau} is the total group contribution.
#'
#' @aliases snowdrift.game
#'
#' @param type character string indicating type of fitness function
#'     (e.g., \code{"tanh"}, \code{"quadratic"})
#' @param nPlayers number of players in a single group
#' @param nGroups number of groups into which the population can be divided
#' @param nRepetitions number of repetitions of the game to be played between reproductive events (default=1)
#' @param simParms named list of simulation parameters
#' @param ... further arguments passed to
#'     \code{\link{setup_gerf_snowdrift_game}},
#'     \code{\link{setup_tanh_snowdrift_game}} or
#'     \code{\link{setup_quadratic_snowdrift_game}}
#'
#' @note The total population size is \code{Npop = nPlayers * nGroups}.
#'
#' @return a \code{snowdrift.game} object, which is a list containing
#'     many parameters, function definitions, and properties of the
#'     game (e.g., singular strategies).
#'
#' @importFrom Deriv Deriv
#'
#' @seealso \code{\link{setup_tanh_snowdrift_game}},
#'     \code{\link{setup_gerf_snowdrift_game}},
#'     \code{\link{setup_quadratic_snowdrift_game}},
#'     \code{\link{set_simParms}}
#'
#' @examples
#' sgtanh <- create_snowdrift_game()
#' summary(sgtanh)
#' plot(sgtanh)
#' sg <- create_snowdrift_game(type="gerf", nPlayers=4, nGroups=2)
#' summary(sg)
#' plot(sg)
#'
#' @export
create_snowdrift_game <- function(type="tanh", # type of fitness function
                                  nPlayers=2, nGroups=5, nRepetitions=1,
                                  ... ) {

  Npop <- nPlayers*nGroups
  message("type = ", type, ", nPlayers = ", nPlayers, ", nGroups = ", nGroups, ", Npop = ", Npop, ", nRepetitions = ", nRepetitions)

  ## basic setup
  basic <- switch(type,
                  gerf = setup_gerf_snowdrift_game(...),
                  tanh = setup_tanh_snowdrift_game(...),
                  quadratic = setup_quadratic_snowdrift_game(
                    Npop=Npop,
                    nPlayers=nPlayers,
                    ...),
                  NULL)
  if (is.null(basic)) stop("fitness function type '", type, "' unknown")

  ## EXTRACT BENEFIT AND COST FUNCTIONS AND PARAMETERS
  benefitString <- basic$benefitString
  costString <- basic$costString
  benefitParms <- basic$benefitParms
  costParms <- basic$costParms

  ## POPULATION PARAMETERS LIST
  popParms <- list(nPlayers = nPlayers, nGroups = nGroups, Npop = Npop, nRepetitions = nRepetitions)

  ## COMBINE PARMS LISTS
  parms <- c(benefitParms, costParms, popParms)

  ## SAVE PARAMETER TYPES AS ATTRIBUTES
  for (name in names(benefitParms)) attr(parms[[name]],"type") <- "benefit"
  for (name in names(costParms)) attr(parms[[name]],"type") <- "cost"
  for (name in names(popParms)) attr(parms[[name]],"type") <- "population"

  ## COST AND BENEFIT FUNCTIONS AND THEIR DERIVATIVES
  costStringList <- get_deriv_strings( costString, var="x" )
  costExprList <- string_to_expr( costStringList )
  benefitStringList <- get_deriv_strings( benefitString, var="tau" )
  benefitExprList <- string_to_expr( benefitStringList )

  ## STRING AND EXPR LISTS
  string <- list(cost=costStringList, benefit=benefitStringList)
  expr <- list(cost=costExprList, benefit=benefitExprList)

  ## CREATE MODEL
  model <- list(type = type,
                parms = parms,
                string = string,
                expr = expr)
  class(model) <- c("snowdrift.game", class(model))

  ## CALCULATE DERIVED PARAMETERS AND ADD TO MODEL
  derivedParms <- calculate_derived_parms( model )
  for (name in names(derivedParms)) attr(derivedParms[[name]],"type") <- "derived"
  parms <- c(parms, derivedParms)
  model$parms <- parms

  return(model)
}

##' Get symbolic derivatives as strings
##'
##' @param string function expression as character string (\code{f}
##'     argument of \code{\link[Deriv]{Deriv}})
##' @param var variable name (\code{x} argument of
##'     \code{\link[Deriv]{Deriv}})
##' @param nderiv vector of desired derivative orders (\code{nderiv}
##'     argument of \code{\link[Deriv]{Deriv}})
##' @return named list of character strings
##' @importFrom Deriv Deriv
##' @examples
##' get_deriv_strings("a*x^4")
##'
##' @export
get_deriv_strings <- function( string, var="x", nderiv=0:2 ) {
    derivList <- lapply(nderiv,
                        function(n) Deriv::Deriv(f=string, x=var,
                                                 cache.exp = FALSE,
                                                 nderiv=n)
                        )
    names(derivList) <- nderiv
    return(derivList)
}

#' Convert strings to expressions
#'
#' @param stringList a list of character strings
#' @export
string_to_expr <- function( stringList ) {
    if (is.null(names(stringList))) stop("argument must be a named list")
    exprList <- list()
    for (name in names(stringList)) {
        exprList[[name]] <- parse(text = stringList[[name]])
    }
    return(exprList)
}

## functions for gerf sigmoid:
tau_pm_gerf <- function(tauTurn, nben, L, mmf, pm) {
    tauTurn +
        pm*(2*nben*L/((mmf+1)*gsl::gsl_sf_gamma(1/(2*nben))) *
            (log(mmf+1))^(1/(2*nben)))
}
#' Finite population singular strategy for gerf sigmoid
#'
#' @inheritParams setup_gerf_snowdrift_game
#' @inheritParams XstarN_fun
#' @import gsl
#' @export
XstarN_gerf <- function(tauTurn, nben, L, mmf, pm, Npop, nPlayers) {
    (1/nPlayers)*(tauTurn +
        pm*(2*nben*L/((mmf+1)*gsl::gsl_sf_gamma(1/(2*nben))) *
            (log((mmf+1)*(Npop-nPlayers)/(Npop-1)))^(1/(2*nben))))
}

## functions for "tanh" type sigmoid:
## mmf: maximum marginal fitness; eq (68) in easierMS:
mmf_tanh <- function(K, L, A) {(K*A*(1+L)/4) - 1}
tau_pm <- function(K, L, A, pm) {
    m <- mmf_tanh(K, L, A)
    if (m < 0) warning("tau_pm: m < 0: A=",A, ", L=", L, ", K=", K, ", m=", m)
    out <- (1/A)*log(L/(2*m+1+pm*2*sqrt(m*(m+1))))
    return( out )
}
explore_tau_pm <- function(K, L, A, ...) {
    taumax <- tau_pm(K=K,L=L,A=A,pm=+1)
    taumin <- tau_pm(K=K,L=L,A=A,pm=-1)
    matplot(xx,cbind(taumin,taumax),type="l",las=1,
            ylab="taumin and taumax"
            )
}
K_finite <- function(model) {
    with(model$parms, {
        K*(Npop - nPlayers)/(Npop - 1)
    })
}

#' maximum marginal fitness for NSG
#'
#' @inheritParams dbenefit
#' @export
mmf_fun <- function(model) {
  if (model$type == "tanh") {
    mmf <- with(model$parms, mmf_tanh(K, L, A))
  }
  else {if (model$type == "gerf") {
    mmf <- model$parms$mmf
  }
  else
  {
    mmf<-NA
  }
  }
  return(mmf)
}

#' Infinite population singular strategy
#'
#' @inheritParams dbenefit
#' @export
XstarInf_fun <- function(model,
                         nPlayers=model$parms$nPlayers) {
    nPlayersArg <- nPlayers # so with() doesn't replace passed argument
    if (model$type == "quadratic") {
        ## eqns (37) and (38) in invasion-fixation ms:
        XstarInf <- with(model$parms, (c1-b1) / (2*(nPlayersArg*b2 - c2)))
    } else if (model$type == "tanh") {
        tauMax <- with(model$parms, tau_pm(K,L,A,-1))
        XstarInf <- tauMax / nPlayersArg
    }
    else if (model$type == "gerf") {
      tauMax <- with(model$parms, tau_pm_gerf(tauTurn, nben, L, mmf, +1))
      XstarInf <- tauMax / nPlayersArg
    }
    return(XstarInf)
}

#' Finite population singular strategy
#'
#' @inheritParams dbenefit
#' @param pm plus or minus: \code{+1/-1} for stable/unstable singular strategy
#' @param nPlayers number of players
#' @param Npop population size
#' @examples
#' sg <- create_snowdrift_game(type="gerf")
#' XstarN_fun(sg)
#' XstarN_fun(sg, -1)
#' @export
XstarN_fun <- function(model,
                       pm=+1,
                       nPlayers=model$parms$nPlayers,
                       Npop=model$parms$Npop) {
    ## NOTE: We must be careful to take the model parameters from the
    ## model object, but use the arguments for nPlayers and Npop.
    if (model$type == "quadratic") {
        fac <- (Npop - nPlayers)/(Npop - 1)
        b1 <- model$parms$b1
        b2 <- model$parms$b2
        c1 <- model$parms$c1
        c2 <- model$parms$c2
        XstarN <- (c1-fac*b1) / (2*(fac*nPlayers*b2 - c2))
    } else if (model$type == "tanh") {
        K <- model$parms$K
        Ktilde <- K*(Npop - nPlayers)/(Npop - 1)
        L <- model$parms$L
        A <- model$parms$A
        tauMaxN <- tau_pm(Ktilde,L,A,-1)
        XstarN <- tauMaxN / nPlayers
    } else if (model$type == "gerf") {
        tauTurn <- model$parms$tauTurn
        nben <- model$parms$nben
        L <- model$parms$L
        mmf <- model$parms$mmf
        XstarN <- XstarN_gerf(tauTurn, nben, L, mmf, pm, Npop, nPlayers)
    } else {
        warning("XstarN_fun: model type '", model$type, "' unknown")
        XstarN <- NA
    }
    ## as.numeric removes irrelevant attributes:
    return(as.numeric(XstarN))
}

#' Is this a Natural Snowdrift Game?
#'
#' @inheritParams dbenefit
#' @export
is_nsg <- function(model) {
    if ((model$type == "tanh")||(model$type == "gerf")) {
        return(TRUE)
    } else if (model$type =="quadratic") {
        return(FALSE)
    } else { # model unknown
        return(NA)
    }
}


#' Critical population size
#'
#' The minimum population size for which a cooperative ESSN exists for
#' the given number of players.
#' @inheritParams dbenefit
#' @param nPlayers number of players
#' @export
critical_Npop <- function(model, nPlayers=model$parms$nPlayers) {
    if (is_nsg(model)) {
        mmf <- mmf_fun(model)
        NpopCrit <- nPlayers + (nPlayers-1)/mmf # easierMS eq (4)
        return(NpopCrit)
    } else {
        return(NA)
    }
}

#' Critical population size
#'
#' The maximum population size for which a cooperative ESSN exists for
#' the given number of groups.
#' @inheritParams dbenefit
#' @param nGroups number of groups
#' @export
critical_Npop_from_nGroups <- function(model, nGroups=model$parms$nGroups) {
    if (is_nsg(model)) {
        mmf <- mmf_fun(model)
        NpopCrit <- nGroups / (1 - mmf*(nGroups-1)) # easierMS eq (7)
        return(NpopCrit)
    } else {
        return(NA)
    }
}

#' Critical number of players
#'
#' The maximum number of players for which a cooperative ESSN exists
#' for the given population size.
#' @inheritParams dbenefit
#' @param Npop population size
#' @export
critical_nPlayers <- function(model, Npop=model$parms$Npop) {
    if (is_nsg(model)) {
        mmf <- mmf_fun(model)
        nPlayersCrit <- (mmf*Npop + 1)/(mmf + 1) # easierMS eq (8)
        return(nPlayersCrit)
    } else {
        return(NA)
    }
}

#' Critical number of groups
#'
#' The minimum number of groups for which a cooperative ESSN exists
#' for the given population size.
#' @inheritParams dbenefit
#' @param Npop population size
#' @export
critical_nGroups <- function(model, Npop=model$parms$Npop) {
    if (is_nsg(model)) {
        mmf <- mmf_fun(model)
        nGroupsCrit <- (mmf + 1)/(mmf + 1/Npop) # easierMS eq (9)
        return(nGroupsCrit)
    } else {
        return(NA)
    }
}

#' Calculate derived parameters for a \code{snowdrift.game} model
#'
#' @param model a \code{snowdrift.game} object
#' @seealso \code{\link{create_snowdrift_game}}
#' @export
calculate_derived_parms <- function( model ) {

    stopifnot("snowdrift.game" %in% class(model))

    if (model$type == "quadratic") {
        derivedParms <- with(model$parms, {
            ## eqns (37) and (38) in invasion-fixation ms:
            ##XstarInf <- (c1-b1) / (2*(nPlayers*b2 - c2))
            XstarInf <- XstarInf_fun( model )
            ##fac <- (Npop - nPlayers)/(Npop - 1)
            ##XstarN <- (c1-fac*b1) / (2*(fac*nPlayers*b2 - c2))
            XstarN <- XstarN_fun( model )
            nPlayersCrit <- NA # FIX: sort this out
            nGroupsCrit <- NA # FIX: sort this out
            message("calculate_derived_parms: setting nPlayersCrit and nGroupsCrit to NA")
            return(list(XstarInf=XstarInf, XstarN=XstarN,
                        nPlayersCrit=nPlayersCrit))
        })
        return(derivedParms)
    } else if (model$type == "gerf") {
        derivedParms <- with(model$parms, {
            tauMin <- tau_pm_gerf(tauTurn, nben, L, mmf, -1)
            tauMax <- tau_pm_gerf(tauTurn, nben, L, mmf, +1)
            XstarInf <- tauMax / nPlayers
            XstarN <- XstarN_gerf(tauTurn, nben, L, mmf, +1, Npop, nPlayers)
            ## payoff extrema difference (ped):
            ped <- benefit(model, tauMax) - benefit(model, tauMin) - (tauMax-tauMin)
            XstarInf <- XstarInf_fun( model )
            NpopCrit <- critical_Npop( model )
            nPlayersCrit <- critical_nPlayers( model )
            if (nPlayersCrit < nPlayers) #too many players for cooperation to evolve
              warning("calculate_derived_parms:  nPlayersCrit = ",
                      nPlayersCrit, " < ", nPlayers, " = nPlayers")
            nGroupsCrit <- critical_nGroups( model )

            return(list(tauMin=tauMin, tauMax=tauMax,
                        XstarInf=XstarInf, XstarN=XstarN,
                        ped=ped, NpopCrit=NpopCrit, nPlayersCrit=nPlayersCrit,
                        nGroupsCrit=nGroupsCrit))
        })
        return(derivedParms)
    }

    if (model$type != "tanh")
        stop("model type '", model$type, "' unknown")

    derivedParms <- with(model$parms, {
        tauTurn <- log(L) / A
        mmf <- mmf_tanh(K, L, A)
        tauMin <- tau_pm(K,L,A,+1)
        tauMax <- tau_pm(K,L,A,-1)
        ##XstarInf <- tauMax / nPlayers
        XstarInf <- XstarInf_fun( model )
        NpopCrit <- critical_Npop( model )
        ##nPlayersCrit <- Npop - (4/(A*K*(1+L)))*(Npop-1)
        nPlayersCrit <- critical_nPlayers( model )
        if (nPlayersCrit < nPlayers)  #too many players for cooperation to evolve
            warning("calculate_derived_parms:  nPlayersCrit = ",
                    nPlayersCrit, " < ", nPlayers, " = nPlayers")
        nGroupsCrit <- critical_nGroups( model )
        Ktilde <- K_finite(model)
        tauMinN <- tau_pm(Ktilde,L,A,+1)
        tauMaxN <- tau_pm(Ktilde,L,A,-1)
        ##XstarN <- tauMaxN / nPlayers
        XstarN <- XstarN_fun(model, nPlayers=nPlayers, Npop=Npop)
        return(list(tauTurn=tauTurn,
                    mmf=mmf,
                    tauMin=tauMin,
                    tauMax=tauMax,
                    XstarInf=XstarInf,
                    NpopCrit=NpopCrit,
                    nPlayersCrit=nPlayersCrit,
                    nGroupsCrit=nGroupsCrit,
                    Ktilde=Ktilde,
                    tauMinN=tauMinN,
                    tauMaxN=tauMaxN,
                    XstarN=XstarN
                    ))
    })
    return(derivedParms)

    ## FIX: below is old notes from original version...

    ## FIX: following depends on NSG...
    ## Add etaMax and etaTurn to the model, since they're needed for simulation params
    etaDeriv <- list(etaMax = calc_eta_max(m = model),
                     etaTurn = calc_eta_turn(m = model))

    ## DERIVED SIMULATION PARAMS:
    HFinESS <- find_ESS(model) # finite population ESS
    if (is.null(hMax)) hMax <- floor(HFinESS + 2) # make sure we go significantly beyond the ESS
    if (is.null(hMin)) hMin <- hMax/simParms$hNum
    hVec <- seq(hMin, hMax, length = simParms$hNum) # sequence of mutant strategies
    HresVec <- seq(0, hMax*0.999, length = simParms$HresNum - 1) # 0.999 to avoid intersecting hVec
    HresVec <- sort(c(HresVec, HFinESS)) # sequence of resident strategies
    ##FIX: DE: making sure HFinESS is what is used:

    ##HresVec <- c(HFinESS)
    ## HresVec is used only when multiple processors in use ## FIX: true??

    return(c(etaDeriv, list(HFinESS)))
}

#' Set simulation parameters
#' @inheritParams dbenefit
#' @param maxGen maximum allowed generations
#' @param nReplicates number of replicate simulations
#' @param nMutantStrategies number of mutant strategies
#' @param nResidentStrategies number of resident strategies; enables
#'     multiple simultaneous jobs (should be \code{<=} number of
#'     processors on machine).  FIX: not currently set to run in
#'     parallel
#' @param selectionProcess character string indicating which selection
#'     process to use: Wright-Fisher (\code{"WF"}), Moran
#'     (\code{"M"}), or birth-death (\code{"BD"})
#' @return a \code{simParms} object, which is just a list of the
#'     parameters and their values
#' @export
set_simParms <- function(model,
                         maxGen = 100,
                         nReplicates = 100,   # was simNum
                         mutantMin = 1e-6,    # was hMin
                         mutantMax = 5,       # was hMax
                         nMutantStrategies = 100, # was hNum
                         residentMin = 0,
                         residentMax = model$parms$XstarInf,
                         nResidentStrategies = 2, # was HresNum
                         selectionProcess = "WF"
                         ) {

    stopifnot("snowdrift.game" %in% class(model))

    mutantStrategies <- seq( mutantMin, mutantMax, length = nMutantStrategies )
    residentStrategies <-
        if (nResidentStrategies > 1)
            seq( residentMin, residentMax, length = nResidentStrategies ) else residentMax

    simParms <- list(
        maxGen = maxGen,
        nReplicates = nReplicates,
        nMutantStrategies = nMutantStrategies,
        mutantStrategies = mutantStrategies,
        nResidentStrategies = nResidentStrategies,
        residentStrategies = residentStrategies,
        selectionProcess = selectionProcess
    )

    ## check the mutant and resident strategies are always different:
    nSame <- length(intersect(mutantStrategies,residentStrategies))
    if (nSame != 0) {
        print(mutantStrategies)
        print(residentStrategies)
        msg <- sprintf("ERROR: %d identical strategies", nSame)
        stop(msg)
    }

    class(simParms) <- c("simParms", class(simParms))
    return(simParms)
}

#' Summary method for \code{simParms} objects
#'
#' @param x a \code{\link[=set_simParms]{simParms}} object
#' @export
summary.simParms <- function( x ) {
    cat("maximum number of generations:  ", x$maxGen, "\n")
    cat("number of replicates:  ", x$nReplicates, "\n")
    cat("range of", x$nMutantStrategies, "mutant(s): ", range(x$mutantStrategies), "\n")
    cat("range of", x$nResidentStrategies, "resident(s): ", range(x$residentStrategies), "\n")
    cat("selection process:  ", x$selectionProcess, "\n")
}

#' Legend for \code{simParms} objects
#'
#' @param sp a \code{\link[=set_simParms]{simParms}} object
#' @param iRes index of desired resident strategy
#' @param ... further arguments passed to \code{\link{legend.default}}
#' @export
elegend.simParms <- function(sp, x="left", y=NULL, iRes=1,
                            ... ) {
    resident <- signif(sp$residentStrategies[iRes], 4)
    graphics::legend(x, y,
           bty="n",
           legend=c(paste0("maxGen: ",sp$maxGen),
                    paste0("nReplicates: ", sp$nReplicates),
                    paste0("resident: ", resident),
                    paste0("selection process: ", sp$selectionProcess)))
}

#' \emph{n}th derivative of benefit function for \code{snowdrift.game}
#'
#' @param model a \code{\link{snowdrift.game}} object
#' @param tau total group contribution
#' @details \code{nderiv=0} yields the benefit function itself
#' @seealso \code{\link{benefit}}
#' @export
dbenefit <- function(model, tau, nderiv=1) {
    stopifnot("snowdrift.game" %in% class(model))
    derivOrder <- as.character(nderiv)
    functionExpr <- model$expr$benefit[[derivOrder]]
    val <- eval(functionExpr, c(list(tau), model$parms))
    return(val)
}

#' Benefit function for \code{snowdrift.game}
#'
#' @inheritParams dbenefit
#' @seealso \code{\link{dbenefit}}
#' @export
benefit <- function(model, tau) {
    ## zeroth derivative of benefit function:
    val <- dbenefit(model, tau, nderiv=0)
    return(val)
}

## FIX: This is identical to dbenefit but with "benefit" replaced by
## "cost" everywhere.  Not ideal, but probably not worth getting
## fancy.
##
#' \emph{n}th derivative of cost function for \code{snowdrift.game}
#'
#' @inheritParams dbenefit
#' @param x focal individual's contribution
#' @export
dcost <- function(model, x, nderiv=1) {
    stopifnot("snowdrift.game" %in% class(model))
    derivOrder <- as.character(nderiv)
    functionExpr <- model$expr$cost[[derivOrder]]
    val <- eval(functionExpr, c(list(x), model$parms))
    return(val)
}

#' Cost function for \code{snowdrift.game}
#'
#' @inheritParams dcost
#' @export
cost <- function(model, x) {
    ## zeroth derivative of cost function:
    val <- dcost(model, x, nderiv=0)
    return(val)
}

#' Fitness function for \code{snowdrift.game}
#'
#' @inheritParams benefit
#' @inheritParams cost
#' @export
fitness <- function(model, x, tau) {
    stopifnot("snowdrift.game" %in% class(model))
    val <- benefit(model,tau) - cost(model,x)
    return(val)
}

#' Total group contribution in \code{snowdrift.game}
#'
#' @inheritParams benefit
#' @export
total_group_contribution <- function(model, resident, focal) {
    out <- (model$parms$nPlayers-1)*resident + focal
    return(out)
}

#' Compute reasonable upper limit for fitness function plots
#'
#' What is reasonable depends on whether an ESS exists.
#' If not, we used \code{xmax.default}.
#' If so, we use the ESS multiplied by \code{xmax.default}.
#' @inheritParams benefit
#' @export
reasonable_xmax <- function( model, xmax.default=2 ) {
    XstarInf <- model$parms$XstarInf
    if (is.null(XstarInf)) {
        xmax <- xmax.default
    } else {
        xmax <- xmax.default * XstarInf
    }
    if (xmax <= 0) stop("xmax = ", xmax)
    return(xmax)
}

#' Get singular strategy
#'
#' Either the infinite population singular strategy \code{XstarInf} if
#' it exists or the \code{defaultStrategy} otherwise.
#' @inheritParams benefit
#' @export
get_singular_strategy <- function (model, defaultStrategy=1) {
    if (is.null(model$parms$XstarInf)) {
        strat <- defaultStrategy
    } else {
        strat <- model$parms$XstarInf
    }
    return(strat)
}

#' Plot method for \code{snowdrift.game} objects
#'
#' @inheritParams dbenefit
#' @param what character vector: what to plot (\code{"cost"},
#'     \code{"benefit"}, \code{"fitness"}, \code{"dcost"},
#'     \code{"dbenefit"}, \code{"d2benefit"}, \code{"all"}, or any
#'     combination)
#' @param residentStrategy contribution level of residents
#' @param xmin,xmax domain limits for focal individual's strategy
#' @param xmax.group domain limit for total group contribution
#' @param ylim.fitness upper limit for fitness plot
#' @param xlab.focal,xlab.total \eqn{x} axis labels for individual
#'     vs. total group contributions
#' @param ... other \code{\link{graphical parameters}}
#'
#' @export
plot.snowdrift.game <- function(model,
                                add=FALSE,
                                residentStrategy=get_singular_strategy(model),
                                what="all",
                                xmin=0, xmax=reasonable_xmax(model),
                                xmax.group=NULL,
                                ylim.fitness=NULL,
                                show.main=TRUE,
                                show.leftlab=FALSE, cex.leftlab=2,
                                leftlab.fitness="Fitness",
                                show.annotation=TRUE,
                                draw.special.strategies=TRUE,
                                draw.XstarN=TRUE,
                                draw.XstarInf=TRUE,
                                col.resident="yellow",
                                lwd.resident=5,
                                lty.resident="dotted",
                                label.resident="resident strategy",
                                cex.label.resident=0.8,
                                lty.XstarInf="solid",
                                lty.XstarN="solid",
                                col.XstarInf="grey",
                                col.XstarN="black",
                                draw.grid=TRUE,
                                draw.benefit.tangents=FALSE,
                                xlab.focal="focal individual's contribution",
                                xlab.total="total group contribution",
                                ylab.cost="Fitness cost",
                                ylab.benefit="Fitness benefit",
                                ylab.fitness="Fitness",
                                col.cost="darkred",
                                col.benefit="darkblue",
                                col.fitness="darkgreen",
                                xpd.cost=FALSE,
                                xpd.benefit=NA,
                                xpd.fitness=NA,
                                lwd=5, las=1, xaxs="i", bty="L",
                                ... ) {
    if (add) plot <- lines

    if ("all" %in% what)
        what <- c("cost", "benefit", "fitness", "dcost", "dbenefit", "d2benefit")

    ## x axis values
    xx <- seq(xmin, xmax, length=1001)
    tau <- total_group_contribution(model, resident=residentStrategy, focal=xx)

    ## function values
    bb <- benefit(model, tau=tau)
    cc <- cost(model, x=xx)
    ff <- fitness(model, x=xx, tau=tau)

    ## COST
    if ("cost" %in% what) {
        plot(xx, cc, type="n",
             las=las, xaxs=xaxs, bty=bty,
             xlab=xlab.focal,
             ylab=ylab.cost, ...)
        if (draw.grid) my_grid()
        lines(xx, cc, lwd=lwd, col=col.cost, xpd=xpd.cost)
        if (show.annotation) {
            legend("topleft",model$string$cost$`0`, bty="n")
            annotate_parms(model, "cost")
        }
        if (show.main) title(main="Cost function")
        if (show.leftlab) legend("topleft", legend="Cost", bty="n", cex=cex.leftlab)
    }

    ## BENEFIT
    tauLowerLimit <- 0
    tauUpperLimit <-
        if (is.null(xmax.group)) {
            if (!is.null(model$parms$tauMin) && !is.null(model$parms$tauMax)) {
                model$parms$tauMax + model$parms$tauMin
            } else {
                max(tau, na.rm=TRUE)
            }
        } else {
            xmax.group
        }
    tauall <- seq(tauLowerLimit, tauUpperLimit, length=1001)
    bball <- benefit(model, tauall)
    if ("benefit" %in% what) {
        plot(tauall, bball, type="n", las=las, xaxs=xaxs, bty=bty,
             xlab=xlab.total,
             ylab=ylab.benefit, ...)
        if (draw.grid) my_grid()
        if (show.main) title(main="Benefit function")
        if (model$type %in% c("tanh","gerf")) {
            draw_tauMinMax(model, "benefit", lwd=lwd, draw.benefit.tangents)
        }
        lines(tauall, bball, lwd=lwd, col=col.benefit, xpd=xpd.benefit)
        if (show.annotation) {
            legend("topleft",model$string$benefit$`0`, bty="n", xpd=NA)
            annotate_parms(model, "benefit")
        }
        if (show.leftlab) legend("topleft", legend="Benefit", bty="n", cex=cex.leftlab)
    }

    ## FITNESS
    if ("fitness" %in% what) {
        plot(xx, ff, type="n", las=las, xaxs=xaxs, bty=bty, ylim=ylim.fitness,
             xlab=xlab.focal,
             ylab=ylab.fitness, ...)
        if (draw.grid) my_grid()
        if (model$type %in% c("tanh","gerf") && residentStrategy==0) {
            draw_tauMinMax(model, "fitness", lwd=lwd, resident=residentStrategy)
        }
        if (show.main) title(main="Fitness function")
        if (draw.special.strategies) {
            tauRes <- total_group_contribution(model, resident=residentStrategy, focal=residentStrategy)
            fitnessAtRes <- fitness(model,x=residentStrategy,tau=tauRes)
            ##abline(v=residentStrategy, lwd=lwd.resident, col=col.resident)
            ymin <- par("usr")[3]
            lines(x=c(residentStrategy,residentStrategy),
                  y=c(ymin,fitnessAtRes), lwd=lwd.resident, col=col.resident, lty=lty.resident,
                  xpd=NA)
            yposition <- ymin + 0.6*(fitnessAtRes - ymin)
            text(label=label.resident, x=residentStrategy, y=yposition, srt=90, cex=cex.label.resident)
            with(model$parms, {
                if (draw.XstarInf) abline(v=XstarInf, lty=lty.XstarInf, col=col.XstarInf)
                if (draw.XstarN) abline(v=XstarN, lty=lty.XstarN, col=col.XstarN)
            })
        }
        lines(xx, ff, type="l", lwd=lwd, col=col.fitness, xpd=xpd.fitness)
        if (show.annotation) {
            ## annotate_parms(model, "population", "bottomleft")
            legend("topleft",
                   legend=paste0("residentStrategy = ",
                                 signif(residentStrategy,4)), bty="n")
            legend("bottomleft",
                   legend=paste0("nPlayers = ", model$parms$nPlayers), bty="n")
        }
        if (show.leftlab) legend("topleft", legend=leftlab.fitness, bty="n", cex=cex.leftlab)
    }

    ## C'(x)
    dc <- dcost(model, x=xx)
    if ("dcost" %in% what) {
        plot(xx, dc, type="l", xaxs=xaxs, bty=bty,
             lwd=lwd, las=las, col=col.cost,
             xlab=xlab.focal,
             ylab="", ...)
        if (show.main) title(main="Derivative 1 of cost function")
    }

    ## B'(tau)
    db <- dbenefit(model, tau=tauall)
    if ("dbenefit" %in% what) {
        plot(tauall, db, type="l", xaxs=xaxs, bty=bty,
             lwd=lwd, las=las, col=col.benefit,
             xlab=xlab.total,
             ylab="", ...)
        if (show.main) title(main="Derivative 1 of benefit function")
    }

    ## B''(tau)
    db2 <- dbenefit(model, tau=tauall, nderiv=2)
    if (length(db2) == 1) {
        ##FIX: this happens if the derivative is zero... not sure why...
        message("plot.snowdrift.game: hacking db2 to have length(tau)...")
        db2 <- rep(db2, length(tau))
    }
    if ("d2benefit" %in% what) {
        plot(tauall, db2, type="l", col=col.benefit,
             lwd=lwd, las=las, xaxs=xaxs, bty=bty,
             xlab=xlab.total,
             ylab="", ...)
        if (show.main) title(main="Derivative 2 of benefit function")
    }
}

draw_tauMinMax <- function(model, funName, lwd,
                           draw.benefit.tangents,
                           resident=model$parms$XstarInf) {
    with(model$parms, {
        if (funName=="benefit") {
            ftauMin <- benefit(model, tau=tauMin)
            ftauMax <- benefit(model, tau=tauMax)
            if (draw.benefit.tangents) {
                draw_tangent(tauMin, ftauMin, slope=1, width=2)
                draw_tangent(tauMax, ftauMax, slope=1, width=2)
            }
        } else if (funName=="fitness") {
            tau <- total_group_contribution(model, resident=resident, focal=tauMin)
            ftauMin <- fitness(model, x=tauMin, tau=tau)
            tau <- total_group_contribution(model, resident=resident, focal=tauMax)
            ftauMax <- fitness(model, x=tauMax, tau=tau)
        } else {
            stop("function f unavailable")
        }
        ymin <- par("usr")[3]
        lines(c(tauMin,tauMin), c(ymin,ftauMin),
              lty="solid", col="grey", lwd=lwd/2)
        lines(c(tauMax,tauMax), c(ymin,ftauMax),
              lty="solid", col="grey", lwd=lwd/2)
        tauMinLabel <- if (dev_is_tikz()) "$\\tau_{\\rm min}$" else expression(tau[{min}])
        tauMaxLabel <- if (dev_is_tikz()) "$\\tau_{\\rm max}$" else expression(tau[{max}])
        ## FIX: mtext() works nicely when the positions of tauMin and tauMax
        ##      don't overlap with axis labels, but when they do, the text()
        ##      approach is better.  In principle, some fancy calculations
        ##      could work out what's best automatically.
        ##mtext(text=tauMinLabel, side=1, line=0.25, at=tauMin)
        ##mtext(text=tauMaxLabel, side=1, line=0.25, at=tauMax)
        ymeanMin <- ymin + 0.5*(ftauMin - ymin)
        text(label=tauMinLabel, x=tauMin, y=ymeanMin, xpd=NA, adj=0.2)
        ymeanMax <- ymeanMin # looks best to have them at same height
        text(label=tauMaxLabel, x=tauMax, y=ymeanMax, xpd=NA, adj=0.2)
    })
}

annotate_parms <- function(model, parmType, location="bottomright", signif=4) {
    out <- c()
    for (parmName in names(model$parms)) {
        parmValue <- model$parms[[parmName]]
        if (attr(parmValue,"type") == parmType) {
            out <-c(out, paste0(parmName, " = ", signif(parmValue,signif)))
        }
    }
    legend(location, out, bty="n")
    return(invisible(out))
}

#' Summary method for \code{snowdrift.game} objects
#'
#' @inheritParams benefit
#' @export
summary.snowdrift.game <- function( model ) {
    cat("type: ",  model$type, "\n")
    p <- model$parms
    cat("SPECIAL CONTRIBUTION LEVELS:\n")
    cat("  tauMin:  ", p$tauMin, "\n")
    cat("  tauTurn: ", p$tauTurn, "\n")
    cat("  tauMax:  ", p$tauMax, "\n")
    Deltatau <- p$tauMax - p$tauMin
    cat("  Deltatau:", Deltatau, "\n")
    cat("  ped:     ", p$ped, "\n")
    cat("  ped/Deltatau: ", p$ped/Deltatau, "\n")
    cat("  mmf = B'(tauTurn)-1:  ", p$mmf, "\n")
    cat("SINGULAR STRATEGIES:\n")
    cat("  XstarInf:   ", p$XstarInf, "\n")
    cat("  XstarN:     ", p$XstarN, "\n")
    cat("  Difference: ", p$XstarInf - p$XstarN, "\n")
    cat("POPULATION:\n")
    cat("  Npop: ", p$Npop, "\n")
    cat("  nPlayers: ", p$nPlayers, "\n")
    cat("  nGroups: ", p$nGroups, "\n")
    cat("  Minimum population size for which a cooperative ESSN exists:\n")
    cat("    NpopCrit: ", p$NpopCrit, "\n")
    cat("  Maximum number of players for which a cooperative ESSN exists:\n")
    cat("    nPlayersCrit: ", p$nPlayersCrit, "\n")
    cat("  Minimum number of groups for which a cooperative ESSN exists:\n")
    cat("    nGroupsCrit: ", p$nGroupsCrit, "\n")
}

#' \code{elegend} S3 method for \code{snowdrift.game} objects
#'
#' @inheritParams dbenefit
#' @param ... further arguments passed to \code{\link{legend}}
#' @export
elegend.snowdrift.game <- function(model, x="topleft", y=NULL,
                                  digits=4,
                                  ... ) {
    p <- model$parms
    if (model$type == "tanh"){
      parmString <- sprintf("K = %g, L = %g, A = %g", p$K, p$L, p$A)
      } else if (model$type == "gerf"){
        parmString <- sprintf("mmf = %g, tauturn = %g, L = %g, nben=%d", p$mmf, p$tauTurn, p$L, p$nben)
      } else parmString <- ""

    graphics::legend(x=x, y=y,
           bty="n",
           legend=c(paste0("type: ", model$type),
                    paste0("XstarInf: ", signif(p$XstarInf,digits)),
                    paste0("XstarN:   ", signif(p$XstarN,digits)),
                    paste0("Npop:     ", p$Npop),
                    paste0("nPlayers: ", p$nPlayers),
                    paste0("NpopCrit: ", p$NpopCrit),
                    paste0("nPlayersCrit: ", signif(p$nPlayersCrit,digits)),
                    paste0("nRepetitions: ", p$nRepetitions),
                    parmString
                    ))
}

#' Run info legend
#'
#' @param ... further arguments passed to \code{\link{legend}}
#' @export
legend_run_info <- function(x, xwhere="topright", ywhere=NULL,
                            ##xjust=NULL, yjust=NULL,
                            digits=4,
                            ... ) {
    ##loc <- legend_location( x=xwhere, y=ywhere )
    ##xwhere <- loc$x
    ##ywhere <- loc$y
    ##if (is.null(xjust)) xjust <- loc$xjust
    ##if (is.null(yjust)) yjust <- loc$yjust
    pt <- attr(x,"procTime")
    ptSum <- pt[["user"]] + pt[["system"]]
    nodename <- attr(x,"sysInfo")["nodename"]
    ## extract run number if the working directory was saved:
    wd <- attr(x,"wd")
    runNum <- if (is.null(wd)) NULL else dplyr::last(unlist(strsplit(wd,"/")))
    legend(xwhere, ywhere,
           ##xjust=xjust, yjust=yjust,
           bty="n",
           legend=c(paste0("CPU time: ", time_string(ptSum)),
                    paste0("Node: ", nodename),
                    paste0("Run:  ", runNum)),
           ...)
}

#' Calculate fitness for each member of a population playing a
#' \code{\link{snowdrift.game}}
#'
#' @inheritParams dbenefit
#' @param pop population matrix; FIX: document this in one place
#' @return population matrix
#' @export
compute_fitnesses <- function(model, pop) {
  nPlayers <- model$parms$nPlayers
  nGroups <- model$parms$nGroups

  #reinitialize fitnesses to 0
  pop[,"fitness"] <- rep(0,dim(pop)[1])
  for (repetition in 1:(model$parms$nRepetitions))
  {
    # Randomizes group assignment by permuting strategies (along with cumulative fitnesses).
    # Weirdly, sticky causes column names to be lost if the code
    # pop <- pop[sample(1:model$parms$Npop),]
    # is used, but adding [,] somehow manages to avoid this.
    pop[,]<- pop[sample(1:model$parms$Npop),]

    ## calculate group contribution, assuming that the population is
    ## ordered by group.  FIX: is this assumption always valid?
    groupCont <-
      colSums(matrix(pop[,"strategy"], nrow = nPlayers, ncol = nGroups),
              dims = 1)

    ## calculate fitnesses
    ## FIX: why not use fitness(model, x, tau) ???
    pop[,"fitness"] <- pop[,"fitness"]+
      benefit(model, tau = rep(groupCont, each = nPlayers)) -
      cost(model, x = pop[,"strategy"])
  }
  return(pop)
}


#' Merge several \code{fixationData} objects
#' @param x a list of \code{fixationData} objects
#' @import dplyr
#' @export
merge_fixationData <- function( x ) {
    ##FIX: check x is a list of compatible fixationData objects
    fdd <- list()  # fixationData list of data frames
    for (i in seq_along(x)) {
        fdd[[i]] <- as.data.frame(x[[i]]) # comes as matrix; easier to work with data frames
        colnames(fdd[[i]]) <- colnames(x[[i]])
    }
    names(fdd) <- "fdd"
    
    fixationData <- bind_rows(fdd) %>% arrange(mutant)
    attr(fixationData, "model") <- attr(x[[1]],"model")
    attr(fixationData, "simParms") <- attr(x[[1]],"simParms")
    class(fixationData) <- c("fixationData", class(fixationData))
    return(fixationData)
}


#' Plot fixation data as a function of mutant strategy
#'
#' @param x a \code{fixationData} object
#' @param what what to plot, i.e., which column of \code{x}
#' @param show.resident logical; if \code{TRUE}, show vertical line at
#'     resident strategy
#' @param show.XstarInf logical; if \code{TRUE}, show vertical line at
#'     infinite population ESS
#' @param show.XstarN logical; if \code{TRUE}, show vertical line at
#'     ESSN
#' @param show.model.info logical; if \code{TRUE}, show information
#'     about model (game parameters)
#' @param show.sim.info logical; if \code{TRUE}, show information
#'     about simulations that produced the data
#' @param show.run.info logical; if \code{TRUE}, show information
#'     about the run that created the data (e.g., nodename, cpu time)
#' @param ... \code{\link{graphical parameters}}
#' @export
plot.fixationData <- function(x,
                              what="mutantFixationProbability",
                              show.resident=TRUE,
                              lwd.resident=5,
                              show.XstarInf=TRUE,
                              col.XstarInf="black",
                              lwd.XstarInf=2,
                              show.XstarN=TRUE,
                              col.XstarN="black",
                              lwd.XstarN=2,
                              show.model.info=TRUE,
                              show.sim.info=TRUE,
                              show.run.info=TRUE,
                              ## standard parameters:
                              xlim=NULL, ylim=NULL,
                              bg="darkred", cex=1, las=1,
                              xlab="Mutant strategy",
                              ylab=what,
                              bty="L",
                              ...) {

    stopifnot("fixationData" %in% class(x))

    xall <- x[,"mutant"]
    yall <- x[,what]

    if (is.null(xlim)) xlim <- range(xall, na.rm=TRUE)
    if (is.null(ylim)) ylim <- range(yall[is.finite(yall)], na.rm=TRUE)
    iResAll <- unique(x[,"iRes"])

    for (iRes in iResAll) {
        resident <- signif(attr(x,"simParms")$residentStrategies[iRes], 4)
        xsub <- x[x[,"iRes"]==iRes,]
        if (is.null(dim(xsub))) {
            ## there is only one mutant so we need to explicitly
            ## remake it as a matrix:
            xsub <- matrix(data=xsub, nrow=1)
            colnames(xsub) <- colnames(x)
        }
        xx <- xsub[,"mutant"]
        yy <- xsub[,what]
        if (!any(is.finite(yy))) {
            message("For iRes=", iRes, " (resident=", resident,
                    ") all elements of ", what,
                    " are NaN or Inf; not plotting...")
        } else {
            plot(xx, yy, xlim=xlim, ylim=ylim, bty=bty, las=1,
                 xlab=xlab, ylab=ylab,
                 type="n",
                 ...)
            Npop <- attr(x,"model")$parms$Npop
            if (show.resident) abline(v=resident, col="yellow", lwd=lwd.resident)
            abline(h=1/Npop, col="grey")
            if (show.XstarInf) abline(v=attr(x,"model")$parms$XstarInf,
                                      lty="dotted",
                                      lwd=lwd.XstarInf, col=col.XstarInf)
            if (show.XstarN) abline(v=attr(x,"model")$parms$XstarN,
                                    lty="solid", lwd=lwd.XstarN, col=col.XstarN)
            points(xx, yy, type="o", pch=21, bg=bg, cex=cex)
            ##legend("topleft", legend=paste0("resident = ", resident), bty="n")
            if (show.sim.info) elegend(sp=attr(x,"simParms"), iRes=iRes)
            if (show.model.info) elegend(attr(x,"model"))
            if (show.run.info) legend_run_info(x)
        }
    }
}

#' Plot fixation time as a function of mutant strategy
#'
#' @param x a \code{fixationData} object
#' @param ... \code{\link{graphical parameters}}
#' @export
plot_fixation_time <- function(x,
                              xlim=NULL, ylim=NULL,
                              bg=c(mutant="darkred",resident="darkblue"),
                              cex=1, las=1,
                              xlab="Mutant strategy",
                              ylab=" fixation time",
                              bty="L",
                              ...) {

    stopifnot("fixationData" %in% class(x))

    xall <- x[,"mutant"]
    yall <- x[,-(1:5)] ## remove iRes,iMut,resident,mutant,probability

    if (is.null(xlim)) xlim <- range(xall, na.rm=TRUE)
    if (is.null(ylim)) ylim <- range(yall[is.finite(yall)], na.rm=TRUE)
    iResAll <- unique(x[,"iRes"])

    for (iRes in iResAll) {
        resident <- signif(attr(x,"simParms")$residentStrategies[iRes], 4)
        xsub <- x[x[,"iRes"]==iRes,]
        if (is.null(dim(xsub))) {
            ## there is only one mutant so we need to explicitly
            ## remake it as a matrix:
            xsub <- matrix(data=xsub, nrow=1)
            colnames(xsub) <- colnames(x)
        }
        for (strategy in c("mutant", "resident")) {
            xx <- xsub[,"mutant"]
            ymean <- xsub[,paste0(strategy,"FixationTimeMean")]
            ysd <- xsub[,paste0(strategy,"FixationTimeSD")]
            ymin <- xsub[,paste0(strategy,"FixationTimeMin")]
            ymax <- xsub[,paste0(strategy,"FixationTimeMax")]
            if (!any(is.finite(ymean))) {
                message("For iRes=", iRes, " (resident=", resident,
                        ") all elements of mean are NaN or Inf; not plotting...")
            } else {
                plot(xx, ymean, xlim=xlim, ylim=ylim, bty=bty, las=las,
                     xlab=xlab, ylab=paste0(strategy, " ", ylab),
                     type="n",
                     ...)
                ## XstarInf, XstarN, and resident strategy:
                abline(v=resident, col="yellow", lwd=5)
                abline(v=attr(x,"model")$parms$XstarInf, lty="dotted")
                abline(v=attr(x,"model")$parms$XstarN, lty="solid")
                ## error bars:
                arrows(xx, ymean-ysd, xx, ymean+ysd, length=0.05, angle=90, code=3)
                ## means:
                points(xx, ymean, pch=21, bg=bg[strategy], cex=cex)
                ## extremes:
                points(xx, ymin, pch=21, bg=bg[strategy], cex=cex/2)
                points(xx, ymax, pch=21, bg=bg[strategy], cex=cex/2)
                ## annotation:
                elegend(sp=attr(x,"simParms"), iRes=iRes)
                elegend(attr(x,"model"))
                legend_run_info(x)
            }
        }
    }
}
