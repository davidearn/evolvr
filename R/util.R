#' Human-friendly time string
#' @param seconds a real number
#' @param sigdig number of significant digits to retain
#' @importFrom lubridate seconds_to_period
#' @export
time_string <- function(seconds,sigdig=4) {
  return(as.character(signif(lubridate::seconds_to_period(seconds),sigdig)))
}

#' Draw the tangent to a curve at a given point
#'
#' @param x,y coordinates
#' @param slope slope of tangent line
#' @param width length of segment centered at the point \eqn{(x,y)}
#'     projected onto \eqn{x} axis
#' @param ... further \link{graphical parameters}
#' @export
draw_tangent <- function(x,y,slope,width=1, ...) {
    x1 <- x - width/2
    x2 <- x + width/2
    y1 <- y + slope*(x1 - x)
    y2 <- y + slope*(x2 - x)
    lines(c(x1,x2), c(y1,y2), ...)
}

#' More flexible grid function
#'
#' @inheritParams graphics::grid
#' @inheritParams base::pretty
#' @param ... additional parameters passed to \code{\link{abline}}
#' @seealso \code{\link{grid}}
#' @note See the examples in the help for \code{\link{plot.window}}
#' @export
my_grid <- function(nx = NULL, ny = nx, col = "lightgray", lty = "dotted", 
                    lwd = par("lwd"),
                    eps.correct = 1,
                    avoid.axes = TRUE,
                    ...) {
    if (is.null(nx)) {
        ny <- nx <- 5
    }
    rx <- par("usr")[1:2]
    ry <- par("usr")[3:4]
    message("my_grid: rx = ", rx, ",  ry = ", ry)
    vvals <- pretty(rx, nx, eps.correct=eps.correct)
    hvals <- pretty(ry, ny, eps.correct=eps.correct)
    if (avoid.axes) {
        ## it is a bit awkward to be sure to avoid the
        ## axes, but this seems to work OK:
        hvals <- hvals[hvals > pretty(ry[1])[2]]
        vvals <- vvals[vvals > pretty(rx[1])[2]]
    }
    abline(h = hvals,
           v = vvals,
           col = col,
           lty = lty,
           lwd = lwd,
           ...)
    return(invisible(list(h=hvals, v=vvals)))
}

#' Add system attributes to an R object
#'
#' For long jobs, it is helpful to retain various system information
#' in simulation output.  We retain the output of
#' \code{\link{Sys.info}}, \code{\link{Sys.getpid}},
#' \code{\link{Sys.time}}, \code{\link{getwd}},
#' \code{\link{sessionInfo}}, and the value of
#' \code{\link{R.version}}.
#' @param x an R object
#' @return an object of the same time with system attributes attached
#' @export
add_system_attributes <- function( x ) {
    attr(x,"sysInfo") <- Sys.info()
    attr(x,"pid") <- Sys.getpid()
    attr(x,"time") <- Sys.time()
    wd <- getwd()
    attr(x,"wd") <- wd
    bn <- basename(wd)
    attr(x,"basename") <- bn # this is often a run number
    attr(x,"sessionInfo") <- sessionInfo()
    attr(x,"R.version") <- R.version
    return(x)
}

#' Is \code{\link{tikzDevice}} currently in use?
#'
#' @description Convenient for dealing with text in graphics,
#' which can be rendered much more professionally with \code{\link{tikz}}.
#' @return logical
#' @export
dev_is_tikz <- function() {
  return(names(dev.cur()) == "tikz output")
}

#' Legend location
#'
#' This function adds additional named locations to those allowed by
#'     \code{\link{legend}}.  In addition to what \code{\link{legend}}
#'     allows, you can set \code{x} to any of \code{"middleleft"},
#'     \code{"middleright"}, \code{"topmiddle"},
#'     \code{"bottommiddle"}.
#' @note It turns out \code{\link{legend}} already has these options
#'     built in, but called \code{"top"}, \code{"bottom"} \emph{etc.}
#'     See \code{example(legend) (one of the last examples)}.
#' @inheritParams graphics::legend
#' @return list of revised versions of arguments suitable for passing
#'     to \code{\link{legend}}
#' @seealso \code{\link{legend_simParms}}, \code{\link{legend_snowdrift_game}}
#' @examples
#' x <- 1:100
#' y <- rnorm(100)
#' plot(x,y)
#' legloc <- legend_location("middleleft")
#' legend(legloc$x, legloc$y, xjust=legloc$xjust, yjust=legloc$yjust, legend="middle left legend")
#' legloc <- legend_location("middleright")
#' legend(legloc$x, legloc$y, xjust=legloc$xjust, yjust=legloc$yjust, legend="middle right legend")
#' legloc <- legend_location("topmiddle")
#' legend(legloc$x, legloc$y, xjust=legloc$xjust, yjust=legloc$yjust, legend="top middle legend")
#' legloc <- legend_location("bottommiddle")
#' legend(legloc$x, legloc$y, xjust=legloc$xjust, yjust=legloc$yjust, legend="bottom middle legend")
#' @export
legend_location <- function(x, y = NULL, xjust = 0, yjust = 1) {
    if (!is.numeric(x)) {
        if (!is.character(x)) stop("class(x) = ", class(x))
        standard.options <- c("topleft", "bottomleft", "topright", "bottomright")
        if (!(x %in% standard.options)) {
            if (x=="middleleft") {
                x <- par("usr")[1]
                y <- mean(par("usr")[3:4])
                yjust <- 0.5
            } else if (x=="middleright") {
                x <- par("usr")[2]
                y <- mean(par("usr")[3:4])
                xjust <- 1
                yjust <- 0.5
            } else if (x=="topmiddle") {
                x <- mean(par("usr")[1:2])
                y <- par("usr")[4]
                xjust <- 0.5
            } else if (x=="bottommiddle") {
                x <- mean(par("usr")[1:2])
                y <- par("usr")[3]
                xjust <- 0.5
                yjust <- 0
            }
        }
    }
    return(list(x=x, y=y, xjust=xjust, yjust=yjust))
}

#' Default legend method for \pkg{evolvr} objects
#'
#' @inheritParams graphics::legend
#' @seealso \code{\link[graphics]{legend}}
#' @export
elegend <- function (object, x, y=NULL, ...) {
    UseMethod("elegend")
}

#' Transparent colour
#' @description
#' Given a colour, change it to a transparent version with a given \code{alpha} level
#' @param col as in \code{\link{col2rgb}}: vector of any of the three kinds of R color specifications, i.e., either a colour name (as listed by \code{colors()}), a hexadecimal string of the form \code{"#rrggbb"} or \code{"#rrggbbaa"} (see \code{\link{rgb}}), or a positive integer \code{i} meaning \code{\link{palette}}\code{()[i]}.
#' @param alpha opacity level on scale from \code{0} to \code{1} (i.e., \code{1} means completely opaque)
#' @seealso \code{\link{col2rgb}}, \code{\link{rgb}}
#' @export
transparent_colour <- function(col,alpha=150/255) {
    alpha <- round(alpha * 255)
    v <- col2rgb(col)[,1] # color as rgb vector
    tcol <- rgb(v["red"],v["green"],v["blue"],alpha=alpha,maxColorValue = 255)
    return(tcol)
}

