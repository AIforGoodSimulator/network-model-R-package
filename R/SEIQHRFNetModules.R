#' Network Modules
#' Based majorly in https://timchurches.github.io/blog/posts/2020-03-10-modelling-the-effects-of-public-health-interventions-on-covid-19-transmission-part-1/
#' S, Susceptible
#' E, Exposed,
#' I, Infected
#' Q, Quarantined
#' H, Need hospitalization
#' R, Recovered
#' F, Fatality
#'
#' In this work we are going to redefine practically all base modules from EpiModel 


#' Transition from susceptible to exposed upon contact with transmitting individuals
#'
#' Transmitting individuals are exposed, infected and quarantined. Whether individuals that come upon contact
#' with transmitting individuals is drawn from a Bernoulli distribution with P(transmitted) = __finalProb__. To
#' see how finalProb is calculate please see the documentation for the __exposure__ function
#' This model assumes an exposure period before an infection period for any invidividual exposed to the
#' virus. 
#' 
#' This script replaces infection.net
#' check ?infection.net for specification of what the base module does
#' 
#' @details This function is called in the __exposure__ function
#' 
#' @param dat Native EpiModel object
#' @param at simulation timestep
#' @param finalProb Probability of transition between the susceptible and the infected state upon contact with individuals
#' from the state __contact_with__, can be a vector of lenght equal to the number of simulation timesteps.
#' @param contact_with State of the individuals that the susceptible come into contact with
#' @export
FromXToExposed = function(dat, at, finalProb, contact_with = "i"){
  # utility function for exposure module
  
  # get active individuals  
  active <- dat$attr$active
  # get their status (s, i, r...)
  status <- dat$attr$status
  
  # Get infected nodes and count them ---------------------------------------
  # get the ids of the infected
  idsFrom <- which(active == 1 & status == contact_with)
  # get the number of active nodes
  nActive <- sum(active == 1)
  
  # get the number of infected
  nElig <- length(idsFrom)
  nExpFromX <- 0
  
  # Get discordant_edgelist and draw those that get exposed FROM INFECTED --------------
  if (nElig > 0 && nElig < nActive) {
    #discordant_edgelist – a matrix of ID numbers of active dyads 
    # in the network in which one member of the dyad is susceptible
    # and the other is infected -- 
    del <- custom_discord_edgelist(dat, at, from = contact_with)
    
    if (!(is.null(del))) { # if del exists
      # draw samples from binomial (Bernouilli) distribution with 
      # final prob
      transmit <- stats::rbinom(nrow(del), 1, finalProb)
      # select those that have transmitted according to the trial
      del <- del[which(transmit == 1), ]
      # Get new infected
      idsNewExp <- unique(del$sus)
      # update number of infected
      nExpFromX <- length(idsNewExp)
      
      if (nExpFromX > 0) {
        # set new infected as "e", for exposed
        dat$attr$status[idsNewExp] <- "e"
        # store infection time
        dat$attr$expTime[idsNewExp] <- at
      }
    }
  }
  return(list(dat, nExpFromX))
}


#' Transition to exposed (from susceptible) module
#' 
#' Models transition to the exposed state. Makes use of the FromXtoExposed function.
#' 
#' @details The final probability is calculated based on the contact rate (act.rate) and infection probability (inf.prob)
#' using the following equation: 
#' \deqn{finalProb = 1 - (1 - inf.prob)^act.rate}
#' Read more at the pages 65-66 of [this book](http://courses.washington.edu/b578a/readings/bookchap4.pdf)
#' This corresponding to a Binomial trial with a p of success (transmission equal) to __inf.prob__ 
#' and __act.rate__ trials
#' 
#' @details Currently, the transition are processed sequentially. The transition in question are one-to-multiple state transitions]
#' which should be modelled concurrently using transition matrices. The current implementation introduces bias as
#' the probability of the transmission not being succesful is greater than it would be if this was to be processed
#' concurrently. Read more [here](https://github.com/statnet/EpiModel/issues/417) for further discussion or 
#' get in touch if you implement this
#' 
#' @param dat native Epimodel object
#' @param at simulation timestep
#' @export
exposure <- function(dat, at) {
  
  transProb.se <- ifelse(!is.null(dat$param$inf.prob.se), 
                         dat$param$inf.prob.se, 
                         ifelse(!is.null(dat$param$inf.prob.si),
                                dat$param$inf.prob.si,
                                dat$param$inf.prob))
  
  transProb.sq <- ifelse(!is.null(dat$param$inf.prob.sq), 
                         dat$param$inf.prob.sq, 
                         ifelse(!is.null(dat$param$inf.prob.si),
                                dat$param$inf.prob.si,
                                dat$param$inf.prob))
  
  transProb.si = ifelse(!is.null(dat$param$inf.prob.si), 
                        dat$param$inf.prob.si,
                        dat$param$inf.prob)
  
  # get act rate (~contact rate) from relevant sub-state
  actRate.se <- ifelse(!is.null(dat$param$act.rate.se), 
                       dat$param$act.rate.se, 
                       ifelse(!is.null(dat$param$act.rate.si),
                              dat$param$act.rate.si,
                              dat$param$act.rate))
  
  actRate.sq <- ifelse(!is.null(dat$param$act.rate.sq), 
                       dat$param$act.rate.sq, 
                       ifelse(!is.null(dat$param$act.rate.si),
                              dat$param$act.rate.si,
                              dat$param$act.rate))
  
  actRate.si <- ifelse(!is.null(dat$param$act.rate.si), 
                       dat$param$act.rate.si, 
                       dat$param$act.rate)
  
  # Binomial trials
  finalProb.se <- 1 - (1 - transProb.se)^actRate.se
  finalProb.sq <- 1 - (1 - transProb.sq)^actRate.sq
  finalProb.si <- 1 - (1 - transProb.si)^actRate.si
  
  result = FromXToExposed(dat, at, finalProb = finalProb.se, contact_with = "e")
  dat = result[[1]]
  nExpFromExp = result[[2]]
  result = FromXToExposed(dat, at, finalProb = finalProb.sq, contact_with ="q")
  dat = result[[1]]
  nExpFromQuar = result[[2]]
  result = FromXToExposed(dat, at, finalProb = finalProb.si, contact_with = "i")
  dat = result[[1]]
  nExpFromInf = result[[2]]
  
  status = dat$attr$status
  active = dat$attr$active
  
  nExp = nExpFromQuar+nExpFromExp+nExpFromInf
  
  
  # Store results from this iteration and carry on --------------------------
  if (at == 2) {
    # set flow of susceptible to infected for at 1 and at 2
    # remember at 1 the model is initialized nothing actually happens
    dat$epi$se.flow <- c(0, nExp)
    dat$epi$e.num <- c(0, sum(active == 1 & status == "e"))
  }
  else {
    dat$epi$se.flow[at] <- nExp
    dat$epi$e.num[at] <- sum(active == 1 & status == "e")
  }
  return(dat)
}

#' Transition between exposed to infected state (one-to-one transition)
#' 
#' Whether the transition occurs is governed by a Bernouilli trial with probability equal to ei.rate
#' 
#' @param dat native Epimodel object
#' @param at simulation timestep
#' @export
infect = function(dat, at){
  # Exposed to infected transition
  
  # Get active nodes and their status ---------------------------------------
  
  # active nodes
  active <- dat$attr$active
  # get status too (s, i, r, e)
  status <- dat$attr$status
  
  # get the ei rate and ir rate (Exposed to infected and infected to recovered)
  ei.rate <- dat$param$ei.rate
  
  ## E to I progression
  nInf <- 0 # initialise number of infected variable
  # get IDs of those that are eligible for infection (aka active nodes that are exposed)
  # in this case exposure is a prerequisite to be infected
  idsEligInf <- which(active == 1 & status == "e")
  # get number of people that are eligible for infection
  nEligInf <- length(idsEligInf)
  
  if (nEligInf > 0) { # if anyone is eligible
    # pick those that will be infected accordging to a bernouilli trial
    vecInf <- which(stats::rbinom(nEligInf, 1, ei.rate) == 1) # like the toss of a coin
    if (length(vecInf) > 0) { # if by the trial results somepeople were picked
      idsInf <- idsEligInf[vecInf] 
      # get the id for those and transition them from "e" to "i"
      nInf <- length(idsInf)
      status[idsInf] <- "i"
      dat$attr$infTime[idsInf] = at
    }
  }
  
  dat$attr$status <- status
  
  if (at == 2) {
    dat$epi$ei.flow <- c(0, nInf)
    dat$epi$i.num = c(0, sum(active == 1 & status == "i"))
  }
  else {
    dat$epi$ei.flow[at] <- nInf
    dat$epi$i.num[at] = sum(active == 1 & status == "i")
  }
  
  return(dat)
}

#' Transition bettween infected to quarantined state (one-to-one transition)
#' 
#' Whether the transition occurs is governed by a Bernouilli trial with probability equal to ei.rate
#' 
#' This transition is useful to model population compliance with Public Health campaigns and their consequences.
#' 
#' @param dat native Epimodel object
#' @param at simulation timestep
#' @export
quarantining = function(dat,at){
  # Exposed to infected transition
  
  # Get active nodes and their status ---------------------------------------
  
  # active nodess
  active <- dat$attr$active
  # get status too (s, i, r, e)
  status <- dat$attr$status
  
  # get the ei rate and ir rate (Exposed to infected and infected to recovered)
  iq.rate <- dat$param$iq.rate
  
  ## E to I progression
  nQuar <- 0 # initialise number of infected variable
  # get IDs of those that are eligible for infection (aka active nodes that are exposed)
  # in this case exposure is a prerequisite to be infected
  idsEligQuar <- which(active == 1 & status == "i")
  # get number of people that are eligible for infection
  nEligQuar <- length(idsEligQuar)
  
  if (nEligQuar > 0) { # if anyone is eligible
    # pick those that will be infected accordging to a bernouilli trial
    vecQuar <- which(stats::rbinom(nEligQuar, 1, iq.rate) == 1) # like the toss of a coin
    if (length(vecQuar) > 0) { # if by the trial results somepeople were picked
      idsQuar <- idsEligQuar[vecQuar] 
      # get the id for those and transition them from "e" to "i"
      nQuar <- length(idsQuar)
      status[idsQuar] <- "q"
      dat$attr$quarTime[idsQuar] = at
    }
  }
  
  dat$attr$status <- status
  if (at == 2) {
    dat$epi$iq.flow <- c(0, nQuar)
    dat$epi$q.num <- c(0, sum(active == 1 & status == "q"))
  }
  else {
    dat$epi$iq.flow[at] <- nQuar
    dat$epi$q.num[at] <- sum(active == 1 & status == "q")
  }
  
  return(dat)
}

#' Transition between the infected or quarantined state to requiring hospitalisation (multi-to-one transition)
#' 
#' This can be seen as analog to some sort of complication rate.
#' Whether the transition occurs is governed by a Bernouilli trial with probability equal to ih.rate or qh.rate
#' 
#' @param dat native Epimodel object
#' @param at simulation timestep
#' @export
RequireHospitalization = function(dat,at){
  # Quarantined to Hosp and Infected to Hosp
  
  # active nodess
  active <- dat$attr$active
  # get status too (s, i, r, e)
  status <- dat$attr$status
  
  # get the ei rate and ir rate (Exposed to infected and infected to recovered)
  ih.rate <- dat$param$ih.rate
  qh.rate <- dat$param$qh.rate
  
  ## I to H progression
  nHosp <- 0 
  # get IDs of those that are eligible for infection (aka active nodes that are exposed)
  # in this case exposure is a prerequisite to be infected
  idsElig <- which(active == 1 & status == "i")
  # get number of people that are eligible for infection
  nElig <- length(idsElig)
  
  if (!is.null(dat$attr$age)){
    lookupHosp = stats::setNames(dat$param$ratesbyAge$hosp.rate/100, # divided by 100 because the things
                                 # is given as proportion
                                 dat$param$ratesbyAge$AgeGroup)
    ages = dat$attr$age
    ih.rate = lookupHosp[ages] # overwrite single value for hf.rate by age specific one
    ih.rate = ih.rate[idsElig]
  }
  
  if (nElig > 0) { # if anyone is eligible
    # pick those that will be infected accordging to a bernouilli trial
    vec <- which(stats::rbinom(nElig, 1, ih.rate) == 1) # like the toss of a coin
    if (length(vec) > 0) { # if by the trial results somepeople were picked
      ids <- idsElig[vec] 
      # get the id for those and transition them from "e" to "i"
      nHosp<- nHosp + length(ids)
      status[ids] <- "h"
      # store hosp time
      dat$attr$hospTime[ids] <- at
    }
  }
  
  ## Q to H progression
  # get IDs of those that are eligible for infection (aka active nodes that are exposed)
  # in this case exposure is a prerequisite to be infected
  idsElig <- which(active == 1 & status == "q")
  # get number of people that are eligible for infection
  nElig <- length(idsElig)
  
  if (nElig > 0) { # if anyone is eligible
    # pick those that will be infected accordging to a bernouilli trial
    vec <- which(stats::rbinom(nElig, 1, qh.rate) == 1) # like the toss of a coin
    if (length(vec) > 0) { # if by the trial results somepeople were picked
      ids <- idsElig[vec] 
      # get the id for those and transition them from "e" to "i"
      nHosp<- nHosp + length(ids)
      status[ids] <- "h"
      # store hosp time
      dat$attr$hospTime[ids] <- at
    }
  }
  
  dat$attr$status <- status
  
  if (at == 2) {
    dat$epi$h.flow <- c(0, nHosp)
    dat$epi$h.num <- c(0, sum(active == 1 & status == "h"))
  }
  else {
    dat$epi$h.flow[at] <- nHosp
    dat$epi$h.num[at] <- sum(active == 1 & status == "h")
  }
  
  return(dat)
}

#' Transition from quarantined, infected, requiring hospitalisation to recovered (multi-to-one transition)
#' 
#' Whether the transition occurs is governed by a Bernouilli trial with probability equal to qr.rate or ir.rate or
#' hr.rate. Again this happens sequentially which results in same bias explained in the __exposure__ model 
#' (see `help(exposure)`)
#' @param dat native Epimodel object
#' @param at simulation timestep
#' @export
recover = function(dat,at){
  # Quarantined to Hosp and Infected to Hosp
  
  # active nodess
  active <- dat$attr$active
  # get status too (s, i, r, e)
  status <- dat$attr$status
  
  # get the rates to recovery
  qr.rate <- dat$param$qr.rate
  hr.rate <- dat$param$hr.rate
  ir.rate <- dat$param$ir.rate
  
  
  # H -> R progression ------------------------------------------------------
  
  
  nRec <- 0 
  # get IDs of those that are eligible for infection (aka active nodes that are exposed)
  # in this case exposure is a prerequisite to be infected
  idsElig <- which(active == 1 & status == "h")
  # get number of people that are eligible for infection
  nElig <- length(idsElig)
  
  if (nElig > 0) { # if anyone is eligible
    # pick those that will be infected accordging to a bernouilli trial
    vec <- which(stats::rbinom(nElig, 1, hr.rate) == 1) # like the toss of a coin
    if (length(vec) > 0) { # if by the trial results somepeople were picked
      ids <- idsElig[vec] 
      # get the id for those and transition them from "e" to "i"
      nRec<- nRec + length(ids)
      status[ids] <- "r"
      # store hosp time
      dat$attr$dischTime[ids] <- at
    }
  }
  
  
  # Q -> R transition -------------------------------------------------------
  # get IDs of those that are eligible for infection (aka active nodes that are exposed)
  # in this case exposure is a prerequisite to be infected
  idsElig <- which(active == 1 & status == "q")
  # get number of people that are eligible for infection
  nElig <- length(idsElig)
  
  if (nElig > 0) { # if anyone is eligible
    # pick those that will be infected accordging to a bernouilli trial
    vec <- which(stats::rbinom(nElig, 1, qr.rate) == 1) # like the toss of a coin
    if (length(vec) > 0) { # if by the trial results somepeople were picked
      ids <- idsElig[vec] 
      # get the id for those and transition them from "e" to "i"
      nRec<- nRec + length(ids)
      status[ids] <- "r"
      # store hosp time
      dat$attr$recTime[ids] <- at
    }
  }
  
  # I -> R transition -------------------------------------------------------
  # get IDs of those that are eligible for infection (aka active nodes that are exposed)
  # in this case exposure is a prerequisite to be infected
  idsElig <- which(active == 1 & status == "i")
  # get number of people that are eligible for infection
  nElig <- length(idsElig)
  if (nElig > 0) { # if anyone is eligible
    # pick those that will be infected accordging to a bernouilli trial
    vec <- which(stats::rbinom(nElig, 1, ir.rate) == 1) # like the toss of a coin
    if (length(vec) > 0) { # if by the trial results somepeople were picked
      ids <- idsElig[vec] 
      # get the id for those and transition them from "e" to "i"
      nRec<- nRec + length(ids)
      status[ids] <- "r"
      # store hosp time
      dat$attr$recTime[ids] <- at
    }
  }
  
  # end of transitions ------------------------------------------------------
  
  
  dat$attr$status <- status
  
  if (at == 2) {
    dat$epi$r.flow <- c(0, nRec)
    dat$epi$r.num <- c(0, sum(active == 1 & status == "r"))
  }
  else {
    dat$epi$r.flow[at] <- nRec
    dat$epi$r.num[at] <- sum(active == 1 & status == "r")
  }
  
  return(dat)
}

#' Transition from requiring hospitalisation to fatality (death from infection)
#' 
#' Whether the transition occurs is governed by a Bernouilli trial with probability equal to hf.rate. 
#' Further more if someone has spent several days in the H compartment that results in a linear increase
#' in their probability of dying (with a slope of __hosp.tcoeff__). Additionally, if the hospital is at 
#' full capacity or beyond the chance of becoming a fatality increases for everyone (baseline increased. 
#' Finally, the hf.rate is dependent on age as per the data in `ratesbyage`.
#' to __hf.rate.overcap__)
#' @param dat native Epimodel object
#' @param at simulation timestep
#' @export
fatality = function(dat,at){
  # active nodess
  active <- dat$attr$active
  # get status too (s, i, r, e)
  status <- dat$attr$status
  
  # get the rates
  hf.rate <- dat$param$hf.rate
  hf.rate.overcap = dat$param$hf.rate.overcap
  
  # hospital parameters
  hosp.cap = dat$param$hosp.cap
  hosp.tcoeff = dat$param$hosp.tcoeff
  
  ## H to R progression
  nFat <- 0 
  # get IDs of those that are eligible for infection (aka active nodes that are exposed)
  # in this case exposure is a prerequisite to be infected
  idsElig <- which(active == 1 & status == "h")
  # get number of people that are eligible for infection
  nElig <- length(idsElig)
  
  ## by age rates
  if (!is.null(dat$attr$age)){
    lookupFat = setNames(dat$param$ratesbyAge$fat.rate/1000,
                         dat$param$ratesbyAge$AgeGroup)
    ages = dat$attr$age
    hf.rate = lookupFat[ages] # overwrite single value for hf.rate by age specific one
    hf.rate = hf.rate[idsElig]
  }
  
  if (nElig > 0) { # if anyone is eligible
    
    timeInHospElig = at - dat$attr$hospTime[idsElig]
    h.num.yesterday = dat$epi$h.num[at-1]
    if (h.num.yesterday>hosp.cap){
      # weighted average
      hf.rate = ((hosp.cap * hf.rate) + 
                   ((h.num.yesterday - hosp.cap) * hf.rate.overcap)) / 
        h.num.yesterday
    }
    hf.rate <- hf.rate + timeInHospElig*hosp.tcoeff*hf.rate
    
    # pick those that will be infected according to a Bernoulli trial
    vec <- which(stats::rbinom(nElig, 1, hf.rate) == 1) # like the toss of a coin
    if (length(vec) > 0) { # if by the trial results some people were picked
      ids <- idsElig[vec] 
      # get the id for those and transition them from "e" to "i"
      nFat<- nFat + length(ids)
      status[ids] <- "f"
      # store hosp time
      dat$attr$fatTime[ids] <- at
      # deactivate these nodes so that edges for this node are not 
      # counted anymore
      dat$attr$active[ids] <- 0
      dat$attr$exitTime[ids] <- at
      dat$nw[[1]] <- networkDynamic::deactivate.vertices(dat$nw[[1]], onset = at, terminus = Inf,
                                                    v = ids, deactivate.edges = TRUE)
    }
  }
  
  dat$attr$status <- status
  
  if (at == 2) {
    dat$epi$f.flow <- c(0, nFat)
    dat$epi$f.num <- c(0, sum(status == "f"))
  }
  else {
    dat$epi$f.flow[at] <- nFat
    dat$epi$f.num[at] <- sum(status == "f")
  }
  
  return(dat)
}

#' @title Departures: netsim Module
#'
#' @description This function simulates departure for use in \link{netsim} simulations.
#'
#' @param dat Master list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @seealso \code{\link{netsim}}
#'
#' @export
#' @keywords netMod internal
#' @source EpiModel 2.0 source code
custom.departures.net <- function(dat, at) {
  
  # Conditions --------------------------------------------------------------
  vital <- EpiModel::get_param(dat, "vital")
  if (vital == FALSE){
    return(dat)
  }
  
  type <- EpiModel::get_control(dat, "type")
  active <- EpiModel::get_attr(dat, "active")
  status <- EpiModel::get_attr(dat, "status")
  exitTime <- EpiModel::get_attr(dat, "exitTime")
  rates.sus <- EpiModel::get_param(dat, "ds.rate")
  rates.inf <- EpiModel::get_param(dat, "di.rate")
  
  
  rates.rec <- EpiModel::get_param(dat, "dr.rate")
  
  
  # Susceptible departures ------------------------------------------------------
  nDepartures.sus <- 0
  idsElig.sus <- which(active == 1 & status == "s")
  nElig.sus <- length(idsElig.sus)
  if (nElig.sus > 0) {
    vecDepartures.sus <- which(stats::rbinom(nElig.sus, 1, rates.sus) == 1)
    if (length(vecDepartures.sus) > 0) {
      idsDpt.sus <- idsElig.sus[vecDepartures.sus]
      nDepartures.sus <- length(idsDpt.sus)
      active[idsDpt.sus] <- 0
      exitTime[idsDpt.sus] <- at
    }
  }
  
  # Infected departures ---------------------------------------------------------
  nDepartures.inf <- 0
  idsElig.inf <- which(active == 1 & status == "i")
  nElig.inf <- length(idsElig.inf)
  if (nElig.inf > 0) {
    vecDepartures.inf <- which(stats::rbinom(nElig.inf, 1, rates.inf) == 1)
    if (length(vecDepartures.inf) > 0) {
      idsDpt.inf <- idsElig.inf[vecDepartures.inf]
      nDepartures.inf <- length(idsDpt.inf)
      active[idsDpt.inf] <- 0
      exitTime[idsDpt.inf] <- at
    }
  }
  
  # Recovered departures --------------------------------------------------------
  nDepartures.rec <- 0
  idsElig.rec <- which(active == 1 & status == "r")
  nElig.rec <- length(idsElig.rec)
  if (nElig.rec > 0) {
    vecDepartures.rec <- which(stats::rbinom(nElig.rec, 1, rates.rec) == 1)
    if (length(vecDepartures.rec) > 0) {
      idsDpt.rec <- idsElig.rec[vecDepartures.rec]
      nDepartures.rec <- length(idsDpt.rec)
      active[idsDpt.rec] <- 0
      exitTime[idsDpt.rec] <- at
    }
  }
  
  
  # Output ------------------------------------------------------------------
  
  dat <- EpiModel::set_attr(dat, "active", active)
  dat <- EpiModel::set_attr(dat, "exitTime", exitTime)
  
  dat <- EpiModel::set_epi(dat, "ds.flow", at, nDepartures.sus)
  dat <- EpiModel::set_epi(dat, "di.flow", at, nDepartures.inf)
  dat <- EpiModel::set_epi(dat, "dr.flow", at, nDepartures.rec)
  return(dat)
}


