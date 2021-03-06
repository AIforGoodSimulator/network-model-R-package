% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_times}
\alias{plot_times}
\title{Plotting time-metric distributions}
\usage{
plot_times(net_build_object, interactive = TRUE)
}
\arguments{
\item{net_build_object}{object returned by \code{net_simulate}}

\item{interactive}{TRUE/FALSE whether output plot should be interactive (plotly) or
static (ggplot)}
}
\description{
Plotting function to visually evaluate the overall goodness of fit
of the network epidemic simulation with respect to well known
external metrics (incubation period, time-to-recovery)
}
\examples{
plot_times(net_example)
}
