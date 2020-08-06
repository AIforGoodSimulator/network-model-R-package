#' Age distribution of individuals in Moria Camp (in bins)
#'
#' Age distribution in 10 year bins from individuals in Moria camp. 
#' For each interval the lower range is inclusive and the upper range
#' is exclusive.
#'
#' @format A vector with 12883 entries containing the age bin 
#' of each individual
#' @source internal
"agedistribution"

#' Hospitalisation and Fatality rate for each age bin
#' 
#' @format For each 10-year age bin a hospitalisation/fatality 
#' rate is known (units: number of individuals/100/day)
#' \describe{
#'   \item{AgeGroup}{bins of age with open lower end and closed upper end}
#'   \item{hosp.rate}{Hospitalisation rate by age group (units: ?number of fatalities per case/1000)}
#'   \item{fat.rate}{Fatality rate by age group (units: ?number of fatalities per case/1000)}
#' }
"ratesbyage"

#' Example network object
#' 
#' Object that is returned from running `net_simulate()`, available
#' for plot examples 
#' 
"net_example"
