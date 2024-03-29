% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snowdrift.R
\name{plot.fixationData}
\alias{plot.fixationData}
\title{Plot fixation data as a function of mutant strategy}
\usage{
\method{plot}{fixationData}(
  x,
  what = "mutantFixationProbability",
  show.resident = TRUE,
  lwd.resident = 5,
  show.XstarInf = TRUE,
  col.XstarInf = "black",
  lwd.XstarInf = 2,
  show.XstarN = TRUE,
  col.XstarN = "black",
  lwd.XstarN = 2,
  show.model.info = TRUE,
  show.sim.info = TRUE,
  show.run.info = TRUE,
  xlim = NULL,
  ylim = NULL,
  bg = "darkred",
  cex = 1,
  las = 1,
  xlab = "Mutant strategy",
  ylab = what,
  bty = "L",
  ...
)
}
\arguments{
\item{x}{a \code{fixationData} object}

\item{what}{what to plot, i.e., which column of \code{x}}

\item{show.resident}{logical; if \code{TRUE}, show vertical line at
resident strategy}

\item{show.XstarInf}{logical; if \code{TRUE}, show vertical line at
infinite population ESS}

\item{show.XstarN}{logical; if \code{TRUE}, show vertical line at
ESSN}

\item{show.model.info}{logical; if \code{TRUE}, show information
about model (game parameters)}

\item{show.sim.info}{logical; if \code{TRUE}, show information
about simulations that produced the data}

\item{show.run.info}{logical; if \code{TRUE}, show information
about the run that created the data (e.g., nodename, cpu time)}

\item{...}{\code{\link{graphical parameters}}}
}
\description{
Plot fixation data as a function of mutant strategy
}
