% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/susie_utils.R
\name{susie_plot}
\alias{susie_plot}
\title{Plot per variable summary in SuSiE CSs}
\usage{
susie_plot(
  model,
  y,
  add_bar = FALSE,
  pos = NULL,
  b = NULL,
  max_cs = 400,
  add_legend = FALSE,
  ...
)
}
\arguments{
\item{model}{a susie fit, the output of `susieR::susie()`.
It has to contain `z`, `PIP` and optionally `sets`.
It is also possible to take in a vector of z-score or PIP,
in order to plot data from other software program.}

\item{y}{a string indicating what to plot: z (for z-score), PIP, log10PIP
or a random label to plot input data as is.}

\item{add_bar}{add horizontal bar to signals in credible interval.}

\item{pos}{index of variables to plot, default to all variables}

\item{b}{for simulated data, specify b = true effects (highlights in red).}

\item{max_cs}{the biggest CS to display, based on purity (set max_cs in between 0 and 1) or size (>1).}

\item{add_legend}{if TRUE, add a legend to annotate the size and purity of each CS discovered.}
}
\description{
Plot per variable summary in SuSiE CSs
}
