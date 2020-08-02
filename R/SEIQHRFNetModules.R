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

#' Find edges that are in touch and are in different status
#'
#' The EpiModel discordedgelist function finds susceptible individuals that are in touch with exposed individuals.
#' In our case, because susceptible individuals can be in touch with exposed, quarantined or infected individuals
#' we have extended the function to be able to match susceptible individuals to individuals in either of those states.
#'
#' @param dat native simulation object from EpiModel
#' @param at timestep (also native to EpimModel)
#' @param from that susceptible individuals will be matched to (options: i, q, e)
#' @export
custom_discord_edgelist <- function(dat, at, from = "i") {
  
  status <- dat$attr$status
  el <- networkDynamic::get.dyads.active(dat$nw, at = at)
  
  del <- NULL
  if (nrow(el) > 0) {
    el <- el[sample(1:nrow(el)), , drop = FALSE]
    stat <- matrix(status[el], ncol = 2)
    isInf <- matrix(stat %in% from, ncol = 2)
    isSus <- matrix(stat %in% "s", ncol = 2)
    SIpairs <- el[isSus[, 1] * isInf[, 2] == 1, , drop = FALSE]
    ISpairs <- el[isSus[, 2] * isInf[, 1] == 1, , drop = FALSE]
    pairs <- rbind(SIpairs, ISpairs[, 2:1])
    if (nrow(pairs) > 0) {
      sus <- pairs[, 1]
      inf <- pairs[, 2]
      del <- data.frame(at, sus, inf)
    }
  }
  
  return(del)
}

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
    hf.rate=  lookupFat[ages] # overwrite single value for hf.rate by age specific one
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
      dat$nw <- networkDynamic::deactivate.vertices(dat$nw, onset = at, terminus = Inf,
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

#' Modification of the EpiModel::netsim function
#' 
#' The performed modifications are clearly highlighted.
#' 
#' INTERNAL PURPOSE MOSTLY
#' @keywords internal
#' @importFrom foreach %dopar%
custom.netsim <- function(x, param, init, control) {
  
  EpiModel::crosscheck.net(x, param, init, control)
  if (!is.null(control[["verbose.FUN"]])) {
    do.call(control[["verbose.FUN"]], list(control, type = "startup"))
  }
  
  nsims <- control$nsims
  ncores <- ifelse(nsims == 1, 1, min(parallel::detectCores(), control$ncores))
  control$ncores <- ncores
  
  if (is.null(control$depend)) {
    control$depend <- FALSE
  }
  
  if (ncores == 1) {
    for (s in 1:control$nsims) {
      
      ## Initialization Module
      if (!is.null(control[["initialize.FUN"]])) {
        dat <- do.call(control[["initialize.FUN"]], list(x, param, init, control, s))
      }
      
      
      ### TIME LOOP
      if (control$nsteps > 1) {
        for (at in max(2, control$start):control$nsteps) {
          
          ## Module order
          morder <- control$module.order
          if (is.null(morder)) {
            lim.bi.mods <- control$bi.mods[-which(control$bi.mods %in%
                                                    c("initialize.FUN", "verbose.FUN"))]
            morder <- c(control$user.mods, lim.bi.mods)
          }
          
          ## Evaluate modules
          for (i in seq_along(morder)) {
            dat <- do.call(control[[morder[i]]], list(dat, at))
          }
          
          ## Verbose module
          if (!is.null(control[["verbose.FUN"]])) {
            do.call(control[["verbose.FUN"]], list(dat, type = "progress", s, at))
          }
          
        }
      }
      
      # Set output
      if (s == 1) {
        ##################################################################
        ##                         Changed here                         ##
        ##################################################################
        ## out <- saveout.net(dat, s)
        out = custom.saveout.net(dat, s)
        ##################################################################
        ##                         Changed here                         ##
        ##################################################################
        
      } else {
        ##################################################################
        ##                         Changed here                         ##
        ##################################################################
        # out <- saveout.net(dat, s, out)
        out = custom.saveout.net(dat, s, out)
        ##################################################################
        ##                         Changed here                         ##
        ##################################################################
      }
      class(out) <- "netsim"
    }
  }
  
  if (ncores > 1) {
    doParallel::registerDoParallel(ncores)
    
    sout <- foreach::foreach(s = 1:nsims) %dopar% {
      
      control$nsims <- 1
      control$currsim <- s
      
      ## Initialization Module
      if (!is.null(control[["initialize.FUN"]])) {
        dat <- do.call(control[["initialize.FUN"]], list(x, param, init, control, s))
      }
      
      
      ### TIME LOOP
      if (control$nsteps > 1) {
        for (at in max(2, control$start):control$nsteps) {
          
          ## Module order
          morder <- control$module.order
          if (is.null(morder)) {
            lim.bi.mods <- control$bi.mods[-which(control$bi.mods %in%
                                                    c("initialize.FUN", "verbose.FUN"))]
            morder <- c(control$user.mods, lim.bi.mods)
          }
          
          ## Evaluate modules
          for (i in seq_along(morder)) {
            dat <- do.call(control[[morder[i]]], list(dat, at))
          }
          
          ## Verbose module
          if (!is.null(control[["verbose.FUN"]])) {
            do.call(control[["verbose.FUN"]], list(dat, type = "progress", s, at))
          }
          
        }
      }
      
      # Set output
      ## out <- saveout.net(dat, s = 1)
      ##################################################################
      ##                         Changed here                         ##
      ##################################################################
      out = custom.saveout.net(dat, s = 1)
      ##################################################################
      ##                         Changed here                         ##
      ##################################################################
      
      
      class(out) <- "netsim"
      return(out)
    }
    
    #################################################################
    ##                           my code                           ##
    #################################################################
    
    collected_times = list()
    # collect the times from sout then delete them
    for (i in 1:length(sout)) {
      collected_times[[paste0("sim", i)]] <- sout[[i]]$times$sim1 
      sout[[i]]$times <- NULL
    }
    
    #################################################################
    ##                           my code                           ##
    #################################################################
    
    merged.out <- sout[[1]]
    for (i in 2:length(sout)) {
      merged.out <- merge(merged.out, sout[[i]], param.error = FALSE)
    }
    out <- merged.out
    
    #################################################################
    ##                           my code                           ##
    #################################################################
    out$times = collected_times
    #################################################################
    ##                           my code                           ##
    #################################################################
    
    class(out) <- "netsim"
  } # end of parallel execution
  
  return(out)
}

#' Modification of the EpiModel::saveout.net function
#' 
#' The performed modifications are clearly highlighted.
#' 
#' This function has been modified with the purpose of modifying how the net object is saved
#' so that all of our custom states and their statistics can be saved.
#' 
#' INTERNAL PURPOSE MOSTLY
#' @keywords internal
custom.saveout.net <- function(dat, s, out = NULL) {
  
  # Counts number of simulated networks
  num.nw <- ifelse(any(class(dat$nw) == "network"), 1, length(dat$nw))
  
  if (s == 1) {
    out <- list()
    out$param <- dat$param
    out$control <- dat$control
    out$nwparam <- dat$nwparam
    out$control$num.nw <- num.nw
    
    out$epi <- list()
    for (j in 1:length(dat$epi)) {
      out$epi[[names(dat$epi)[j]]] <- data.frame(dat$epi[j])
    }
    
    ##########################################################
    # get transitions times
    ##########################################################
    out$times[[paste0("sim",s)]] = data.frame(expTime = dat$attr$expTime,
                                              infTime = dat$attr$infTime,
                                              quarTime = dat$attr$quarTime,
                                              recTime = dat$attr$recTime,
                                              hospTime = dat$attr$hospTime,
                                              dischTime = dat$attr$dischTime,
                                              fatTime = dat$attr$fatTime,
                                              exitTime = dat$attr$exitTime)
    
    ##########################################################
    #       #       #       #       #       #       #       #
    ##########################################################
    
    out$stats <- list()
    if (dat$control$save.nwstats == TRUE) {
      out$stats$nwstats <- list(dat$stats$nwstats)
    }
    if (dat$control$save.transmat == TRUE) {
      if (!is.null(dat$stats$transmat)) {
        row.names(dat$stats$transmat) <- 1:nrow(dat$stats$transmat)
        out$stats$transmat <- list(dat$stats$transmat)
      } else {
        out$stats$transmat <- list(data.frame())
      }
      class(out$stats$transmat) <- c("transmat", class(out$stats$transmat))
    }
    
    if (dat$control$save.network == TRUE) {
      out$network <- list(dat$nw)
    }
    
    if (!is.null(dat$control$save.other)) {
      for (i in 1:length(dat$control$save.other)) {
        el.name <- dat$control$save.other[i]
        out[[el.name]] <- list(dat[[el.name]])
      }
    }
  }
  
  if (s > 1) {
    for (j in 1:length(dat$epi)) {
      out$epi[[names(dat$epi)[j]]][, s] <- data.frame(dat$epi[j])
    }
    
    ##########################################################
    # get transitions times
    ##########################################################
    out$times[[paste0("sim",s)]] = data.frame(expTime = dat$attr$expTime,
                                              infTime = dat$attr$infTime,
                                              quarTime = dat$attr$quarTime,
                                              recTime = dat$attr$recTime,
                                              hospTime = dat$attr$hospTime,
                                              dischTime = dat$attr$dischTime,
                                              fatTime = dat$attr$fatTime,
                                              exitTime = dat$attr$exitTime)
    
    ##########################################################
    #       #       #       #       #       #       #       #
    ##########################################################
    
    if (dat$control$save.nwstats == TRUE) {
      out$stats$nwstats[[s]] <- dat$stats$nwstats
    }
    if (dat$control$save.transmat == TRUE) {
      if (!is.null(dat$stats$transmat)) {
        row.names(dat$stats$transmat) <- 1:nrow(dat$stats$transmat)
        out$stats$transmat[[s]] <- dat$stats$transmat
      } else {
        out$stats$transmat[[s]] <- data.frame()
      }
    }
    if (dat$control$save.network == TRUE) {
      out$network[[s]] <- dat$nw
    }
    
    if (!is.null(dat$control$save.other)) {
      for (i in 1:length(dat$control$save.other)) {
        el.name <- dat$control$save.other[i]
        out[[el.name]][[s]] <- dat[[el.name]]
      }
    }
  }

  
  ## Final processing
  if (s == dat$control$nsims) {
    
    # Set names for out
    simnames <- paste0("sim", 1:dat$control$nsims)
    for (i in as.vector(which(lapply(out$epi, class) == "data.frame"))) {
      colnames(out$epi[[i]]) <- simnames
    }
    if (dat$control$save.nwstats == TRUE) {
      names(out$stats$nwstats) <- simnames
    }
    if (dat$control$save.transmat == TRUE) {
      names(out$stats$transmat) <- simnames[1:length(out$stats$transmat)]
    }
    if (dat$control$save.network == TRUE) {
      names(out$network) <- simnames
    }
    
    # Remove functions from control list
    ftodel <- grep(".FUN", names(out$control), value = TRUE)
    out$control[ftodel] <- NULL
    out$control$currsim <- NULL
    environment(out$control$nwstats.formula) <- NULL
    
    if (!("temp" %in% dat$control$save.other)) {
      out$temp <- NULL
    }
    
  }
  
  return(out)
}

#' Modification of the EpiModel::initialize.net function
#' 
#' The performed modifications are clearly highlighted.
#' 
#' This function has been modified with the purpose of modifying how the net object is saved
#' so that all of our custom states and their statistics can be saved.
#' 
#' INTERNAL PURPOSE MOSTLY
#' @keywords internal
custom.initialize.net <- function(x, param, init, control, s) {
  
  if (control$start == 1) {
    # Master Data List --------------------------------------------------------
    dat <- list()
    dat$param <- param
    dat$init <- init
    dat$control <- control
    
    dat$attr <- list()
    dat$stats <- list()
    dat$temp <- list()
    
    
    # Network Simulation ------------------------------------------------------
    nw <- simulate(x$fit, basis = x$fit$newnetwork,
                   control = control$set.control.ergm)
    modes <- ifelse(nw %n% "bipartite", 2, 1)
    if (control$depend == TRUE) {
      if (class(x$fit) == "stergm") {
        nw <- networkDynamic::network.collapse(nw, at = 1)
      }
      nw <- EpiModel::sim_nets(x, nw, nsteps = 1, control)
    }
    if (control$depend == FALSE) {
      nw <- EpiModel::sim_nets(x, nw, nsteps = control$nsteps, control)
    }
    nw <- networkDynamic::activate.vertices(nw, onset = 1, terminus = Inf)
    
    
    # Network Parameters ------------------------------------------------------
    dat$nw <- nw
    dat$nwparam <- list(x[-which(names(x) == "fit")])
    dat$param$modes <- modes
    
    
    # Initialization ----------------------------------------------------------
    
    ## Initialize persistent IDs
    if (control$use.pids == TRUE) {
      dat$nw <- EpiModel::init_pids(dat$nw, dat$control$pid.prefix)
    }
    
    ## Pull network val to attr
    form <- EpiModel::get_nwparam(dat)$formation
    
    #fterms <- get_formula_term_attr(form, nw)
    fterms <- setdiff(network::list.vertex.attributes(nw), c("active", "vertex.names", "na"))
    if(length(fterms) == 0) fterms <- NULL
    
    dat <- EpiModel::copy_toall_attr(dat, at = 1, fterms)
    
    ## Store current proportions of attr
    dat$temp$fterms <- fterms
    dat$temp$t1.tab <- EpiModel::get_attr_prop(dat$nw, fterms)
    
    
    ## Infection Status and Time Modules
    ##################################################################
    ##                       This is modified                       ##
    ##################################################################
    # dat <- init_status.net(dat)
    dat = custom.init_status.net(dat)
    ##################################################################
    ##                       This is modified                       ##
    ##################################################################
    
    ## Get initial prevalence
    dat <- custom.get_prev.net(dat, at = 1)
  } else {
    dat <- list()
    
    dat$nw <- x$network[[s]]
    dat$param <- x$param
    dat$control <- control
    dat$nwparam <- x$nwparam
    dat$epi <- sapply(x$epi, function(var) var[s])
    names(dat$epi) <- names(x$epi)
    dat$attr <- x$attr[[s]]
    dat$stats <- sapply(x$stats, function(var) var[[s]])
    dat$temp <- list()
  }
  
  return(dat)
}

#' Modification of the EpiModel::init_status.net function
#' 
#' The performed modifications are clearly highlighted.
#' 
#' This function has been modified with the purpose of modifying how the times spent in each 
#' state
#' 
#' INTERNAL PURPOSE MOSTLY
#' @keywords internal
custom.init_status.net <- function(dat) {

# Variables ---------------------------------------------------------------
tea.status <- dat$control$tea.status
i.num <- dat$init$i.num
i.num.m2 <- dat$init$i.num.m2
r.num <- dat$init$r.num
r.num.m2 <- dat$init$r.num.m2

status.vector <- dat$init$status.vector
num <- network::network.size(dat$nw)
statOnNw <- "status" %in% dat$temp$fterms

modes <- dat$param$modes
if (modes == 1) {
  mode <- rep(1, num)
} else {
  mode <- EpiModel::idmode(dat$nw)
}

type <- dat$control$type


# Status ------------------------------------------------------------------

## Status passed on input network
if (statOnNw == TRUE) {
  status <- network::get.vertex.attribute(dat$nw, "status")
} else {
  if (!is.null(status.vector)) {
    status <- status.vector
  } else {
    status <- rep("s", num)
    status[sample(which(mode == 1), size = i.num)] <- "i"
    if (modes == 2) {
      status[sample(which(mode == 2), size = i.num.m2)] <- "i"
    }
    if (type == "SIR") {
      status[sample(which(mode == 1 & status == "s"), size = r.num)] <- "r"
      if (modes == 2) {
        status[sample(which(mode == 2 & status == "s"), size = r.num.m2)] <- "r"
      }
    }
  }
}
dat$attr$status <- status


## Save out other attr
dat$attr$active <- rep(1, length(status))
dat$attr$entrTime <- rep(1, length(status))
dat$attr$exitTime <- rep(NA, length(status))
if (tea.status == TRUE) {
  dat$nw <- networkDynamic::activate.vertex.attribute(dat$nw,
                                      prefix = "testatus",
                                      value = status,
                                      onset = 1,
                                      terminus = Inf)
}


# Infection Time ----------------------------------------------------------
## Set up inf.time vector
idsInf <- which(status == "i")
infTime <- rep(NA, length(status))

if (!is.null(dat$init$infTime.vector)) {
  infTime <- dat$init$infTime.vector
} else {
  # If vital dynamics, infTime is a geometric draw over the duration of infection
  if (dat$param$vital == TRUE && dat$param$di.rate > 0) {
    if (dat$control$type == "SI") {
      infTime[idsInf] <- -rgeom(n = length(idsInf), prob = dat$param$di.rate) + 2
    } else {
      infTime[idsInf] <- -rgeom(n = length(idsInf),
                                prob = dat$param$di.rate +
                                  (1 - dat$param$di.rate)*mean(dat$param$rec.rate)) + 2
    }
  } else {
    if (dat$control$type == "SI" || mean(dat$param$rec.rate) == 0) {
      # if no recovery, infTime a uniform draw over the number of sim time steps
      infTime[idsInf] <- EpiModel::ssample(1:(-dat$control$nsteps + 2),
                                 length(idsInf), replace = TRUE)
    } else {
      infTime[idsInf] <- -rgeom(n = length(idsInf), prob = mean(dat$param$rec.rate)) + 2
    }
  }
}

#################################################################
##              This is my bit the rest is source              ##
#################################################################
# Infection Time ----------------------------------------------------------
## Set up all time vectors
dat$attr$quarTime = rep(NA, length(status))
dat$attr$hospTime = rep(NA, length(status))
dat$attr$recTime = rep(NA, length(status))
dat$attr$dischTime = rep(NA, length(status))
dat$attr$fatTime = rep(NA, length(status))
dat$attr$expTime = rep(NA, length(status))
#################################################################
##              This is my bit the rest is source              ##
#################################################################

dat$attr$infTime <- infTime

return(dat)
}

#' Modification of the EpiModel::get_prev.net function
#' 
#' The performed modifications are clearly highlighted.
#' 
#' This function has been modified with to be able to store the different number 
#' by vertex attribute (as per epi.by)
#' 
#' INTERNAL PURPOSE MOSTLY
#' @keywords internal
custom.get_prev.net <- function(dat, at) {
  
  active <- dat$attr$active
  modes <- dat$param$modes
  
  #################################################################
  ##                        my commenting                        ##
  #################################################################
  # Subset attr to active == 1
  l <- lapply(1:length(dat$attr), function(x) dat$attr[[x]][active == 1 | active == 0])
  names(l) <- names(dat$attr)
  l$active <- l$infTime <- NULL
  #################################################################
  ##                        my commenting                        ##
  #################################################################
  
  status <- l$status
  
  if (modes == 2) {
    mode <- EpiModel::idmode(dat$nw)[active == 1]
  }
  
  ## Subsetting for epi.by control
  eb <- !is.null(dat$control$epi.by)
  if (eb == TRUE) {
    ebn <- dat$control$epi.by
    ebv <- dat$temp$epi.by.vals
    ebun <- paste0(".", ebn, ebv)
    assign(ebn, l[[ebn]])
  }
  
  ## One mode networks
  if (modes == 1) {
    if (at == 1) {
      dat$epi <- list()
      dat$epi$s.num <- sum(status == "s")
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("s.num", ebun[i])]] <- sum(status == "s" &
                                                       get(ebn) == ebv[i])
        }
      }
      dat$epi$i.num <- sum(status == "i")
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("i.num", ebun[i])]] <- sum(status == "i" &
                                                       get(ebn) == ebv[i])
        }
      }
      #################################################################
      ##                           my code                           ##
      #################################################################
      
      dat$epi$h.num <- sum(status == "h")
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("h.num", ebun[i])]] <- sum(status == "h" &
                                                       get(ebn) == ebv[i])
        }
      }
      dat$epi$q.num <- sum(status == "q")
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("q.num", ebun[i])]] <- sum(status == "q" &
                                                       get(ebn) == ebv[i])
        }
      }
      dat$epi$f.num <- sum(status == "f")
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("f.num", ebun[i])]] <- sum(status == "f" &
                                                       get(ebn) == ebv[i])
        }
      }
      dat$epi$r.num <- sum(status == "r")
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("r.num", ebun[i])]] <- sum(status == "r" &
                                                       get(ebn) == ebv[i])
        }
      }
      
      dat$epi$e.num <- sum(status == "e")
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("e.num", ebun[i])]] <- sum(status == "e" &
                                                       get(ebn) == ebv[i])
        }
      }
      
      
      #################################################################
      ##                           my code                           ##
      #################################################################
      
      if (dat$control$type == "SIR") {
        dat$epi$r.num <- sum(status == "r")
        if (eb == TRUE) {
          for (i in 1:length(ebun)) {
            dat$epi[[paste0("r.num", ebun[i])]] <- sum(status == "r" &
                                                         get(ebn) == ebv[i])
          }
        }
      }
      ## CHANGE 
      # dat$epi$num[at] <- length(status)
      dat$epi$num <- sum(active == 1)
      ## CHANGE
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("num", ebun[i])]] <- sum(get(ebn) == ebv[i] &
                                                     active == 1) # I added the active == 1
        }
      }
    } else {
      # at > 1
      dat$epi$s.num[at] <- sum(status == "s")
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("s.num", ebun[i])]][at] <- sum(status == "s" &
                                                           get(ebn) == ebv[i])
        }
      }
      dat$epi$i.num[at] <- sum(status == "i")
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("i.num", ebun[i])]][at] <- sum(status == "i" &
                                                           get(ebn) == ebv[i])
        }
      }
      
      #################################################################
      ##                           my code                           ##
      #################################################################
      
      dat$epi$h.num[at] <- sum(status == "h")
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("h.num", ebun[i])]][at] <- sum(status == "h" &
                                                       get(ebn) == ebv[i])
        }
      }
      dat$epi$q.num[at] <- sum(status == "q")
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("q.num", ebun[i])]][at] <- sum(status == "q" &
                                                       get(ebn) == ebv[i])
        }
      }
      # there is a bug here, this line below is overwriting the actual number of fatalities
      # as the deactivated vertices are being removed hence their status is deleted
      ## - Fixed by not filtering out inactive nodes at the top oof this function
      dat$epi$f.num[at] <- sum(status == "f")
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("f.num", ebun[i])]][at] <- sum(status == "f" &
                                                       get(ebn) == ebv[i])
        }
      }
      dat$epi$r.num[at] <- sum(status == "r")
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("r.num", ebun[i])]][at] <- sum(status == "r" &
                                                       get(ebn) == ebv[i])
        }
      }
      
      dat$epi$e.num[at] <- sum(status == "e")
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("e.num", ebun[i])]][at] <- sum(status == "e" &
                                                       get(ebn) == ebv[i])
        }
      }
      
      
      #################################################################
      ##                           my code                           ##
      #################################################################
      
      if (dat$control$type == "SIR") {
        dat$epi$r.num[at] <- sum(status == "r")
        if (eb == TRUE) {
          for (i in 1:length(ebun)) {
            dat$epi[[paste0("r.num", ebun[i])]][at] <- sum(status == "r" &
                                                             get(ebn) == ebv[i])
          }
        }
      }
      
      ## CHANGE 
      # dat$epi$num[at] <- length(status)
      dat$epi$num[at] <- sum(active == 1)
      ## CHANGE
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("num", ebun[i])]][at] <- sum(get(ebn) == ebv[i] &
                                                         active == 1) # I added the active = 1
        }
      }
    }
    
  } else {
    stop("Bipartite networks not currently supported in this custom function")
  }
  
  return(dat)
}


