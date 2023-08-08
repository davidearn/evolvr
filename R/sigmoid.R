#' Generalized error function via incomplete gamma function
#'
#' \deqn{E_n(x) = (1/\sqrt\pi) * \Gamma(n) * ( \Gamma(1/n) - \Gamma(1/n, x^n) ) }
#' See \url{https://en.wikipedia.org/wiki/Error_function#Generalized_error_functions}.
#' Unlike \eqn{E_n}, \code{gerf} is normalized so that its maximum
#' is \code{1} regardless of the value of \code{n}.
#' @param x argument
#' @param n order of generalized erf; \eqn{n=2} yields the normal erf
#' @seealso \code{\link[gsl]{gamma_inc}}, \code{\link{gerf_norm}}
#' @import gsl
#' @export
#' @examples
#' x <- seq(-3,3,length=1000)
#' g2 <- gerf(x,2)
#' g4 <- gerf(x,4)
#' g8 <- gerf(x,8)
#' g16 <- gerf(x,16)
#' g32 <- gerf(x,32)
#' y <- cbind(g2,g4,g8,g16,g32)
#' plot(NA, NA, xlim=range(x), ylim=c(-1,1), type="n", bty="L", las=1,
#'      xlab=expression(x), ylab="gerf(x,n)",
#'      main="Generalized Error Function")
#' grid()
#' matlines(x,y, lty=1, lwd=2)
#' plot(NA, NA, xlim=c(0.75,1.25), ylim=c(0.8,1), type="n", bty="L", las=1,
#'      xlab=expression(x), ylab="gerf(x,n)",
#'      main="Generalized Error Function (zoomed in)")
#' grid()
#' matlines(x,y, lty=1, lwd=2)
#'
gerf <- function(x, n) {
    ## FIX: the warning should handle vector n
    if (any((n %% 2) == 1)) warning("gerf: n = ", n, " is odd")
    g1n <- gsl::gsl_sf_gamma(1/n) # Gamma(1/n)
    g <- 1 - (gamma_inc(1/n, x^n)/g1n)
    out <- ifelse( x >= 0, g, -g)
    return( out )
}

#' Derivative of gerf
#' @import gsl
#' @seealso \code{\link{gerf}}
#' @importFrom Deriv drule
#' @export
dgerf <- function(x, n) {
    g1n <- gsl::gsl_sf_gamma(1/n)
    dg <- n * exp(-x^n) / g1n
    return( dg )
}

#' Derivative of dgerf
#' @import gsl
#' @importFrom Deriv drule
#' @seealso \code{\link{gerf}}
#' @export
d2gerf <- function(x, n) {
    g1n <- gsl::gsl_sf_gamma(1/n)
    d2g <- -n^2 * x^(n-1) * exp(-x^n) / g1n
    return( d2g )
}

#' Mean fitness slope with sigmoid via gerf
#'
#' Mean fitness slope with sigmoid via generalized error functions
#' @inheritParams setup_gerf_snowdrift_game
#' @seealso \code{\link{gerf}}, \code{\link{benefit}},
#'     \code{\link{fitness}},
#'     \code{\link{setup_gerf_snowdrift_game}},
#'     \code{\link{mean_fitness_slope_erf}},
#'     \code{\link{mean_fitness_slope_tanh}}
#' @examples
#' nvec <- 1:50
#' mfsg <- mean_fitness_slope_gerf(mmf=1,nben=nvec)
#' plot(NA, NA, xlim=range(nvec), ylim=c(0,1), bty="L", las=1,
#'      xlab="nben", ylab="Mean fitness slope")
#' abline(h=1, col="grey")
#' lines(nvec, mfsg, type="o", pch=21, bg="darkred", cex=0.5)
#' nn <- length(nvec)
#' message("For n=", nvec[nn], ", mfsg=", mfsg[nn])
#' @export
mean_fitness_slope_gerf <- function(mmf,nben) {
    n <- nben
    mmb <- mmf +1 # maximal marginal benefit
    tmp <- (log(mmb))^(1/(2*n))
    g1o2n <- gsl::gsl_sf_gamma(1/(2*n)) # Gamma(1/(2n))
    meanFitnessSlope <- -1 + 
        mmb*(g1o2n/(2*n))*gerf(tmp,2*n)/tmp
    return(meanFitnessSlope)
}
