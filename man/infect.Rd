% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SEIQHRFNetModules.R
\name{infect}
\alias{infect}
\title{Transition between exposed to infected state (one-to-one transition)}
\usage{
infect(dat, at)
}
\arguments{
\item{dat}{native Epimodel object}

\item{at}{simulation timestep}
}
\description{
Whether the transition occurs is governed by a Bernouilli trial with probability equal to ei.rate
}
