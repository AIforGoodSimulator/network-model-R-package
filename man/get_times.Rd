% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{get_times}
\alias{get_times}
\title{Get distribution of time-metrics}
\usage{
get_times(net_build_object)
}
\arguments{
\item{net_build_object}{object returned by \code{net_simulate}}
}
\description{
For diagnostic purposes, the user may wish to explore metrics such as the
distribution of the incubation period (transition between exposed and infected), or recovery time,
or time-to-failure/fatality and compare this to the real world estimates. For example,
it is known that the COVID19 incubation period is anywhere between 5 and 14 days.
}
\details{
The network_simulation object contain several time features for each individual,
namely, the time to exposure, the time to infection, etc. This function is used
to calculate common metric that may be chosen to be evaluated. It is mostly
for internal usage.

Copied from Tim Churches' blogpost.

For more info on tidy evaluation check https://dplyr.tidyverse.org/articles/programming.html#eliminating-r-cmd-check-notes
(i.e. why there is so much .data in this function source)
}
\examples{
times_df = get_times(net_example)
}
