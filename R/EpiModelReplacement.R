#' Find edges that are in touch and are in different status
#'
#' The EpiModel discordedgelist function finds susceptible individuals that are in touch with exposed individuals.
#' In our case, because susceptible individuals can be in touch with exposed, quarantined or infected individuals
#' we have extended the function to be able to match susceptible individuals to individuals in either of those states.
#'
#' @param dat native simulation object from EpiModel
#' @param at timestep (also native to EpimModel)
#' @param network (not applicable) In case of models with multiple networks, the network to pull the current edgelist from. Default of network = 1.
#' @param from that susceptible individuals will be matched to (options: i, q, e)
#' @source EpiModel 2.0 source code
#' @export
custom_discord_edgelist <- function(dat, at, network = 1, from = "i") {
  
  status <- EpiModel::get_attr(dat, "status")
  active <- EpiModel::get_attr(dat, "active")
  tergmLite <- EpiModel::get_control(dat, "tergmLite")
  
  if (tergmLite == TRUE) {
    el <- dat$el[[network]]
  } else {
    el <- networkDynamic::get.dyads.active(dat$nw[[network]], at = at)
  }
  
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
      
      # Check for active status
      keep <- rowSums(matrix(c(active[del$sus], active[del$inf]), ncol = 2)) == 2
      del <- del[keep, ]
      if (nrow(del) < 1) {
        del <- NULL
      }
    }
  }
  
  return(del)
}

#' Modification of the EpiModel::netsim function
#' 
#' DEPRECATED (use custom.netsim instead)
#' The performed modifications are clearly highlighted.
#' 
#' INTERNAL PURPOSE MOSTLY
#' @keywords internal
#' @source EpiModel2.0
#' @importFrom foreach %dopar%
custom.netsim <- function(x, param, init, control) {
  
  EpiModel::crosscheck.net(x, param, init, control)
  if (!is.null(control[["verbose.FUN"]])) {
    do.call(control[["verbose.FUN"]], list(control, type = "startup"))
  }
  
  nsims <- control$nsims
  ncores <- ifelse(nsims == 1, 1, min(parallel::detectCores(), control$ncores))
  control$ncores <- ncores
  
  if (is.null(control$resimulate.network)) {
    control$resimulate.network <- FALSE
  }
  
  s <- NULL
  if (ncores == 1) {
    sout <- lapply(seq_len(control$nsims), function(s) {
      # Run the simulation
      netsim_loop(x, param, init, control, s)
    })
  }
  
  if (ncores > 1) {
    doParallel::registerDoParallel(ncores)
    sout <- foreach::foreach(s = 1:nsims) %dopar% {
      # Run the simulation
      netsim_loop(x, param, init, control, s)
    }
  }
  
  # Process the outputs unless `control$raw.output` is `TRUE`
  if (!is.null(control$raw.output) && control$raw.output == TRUE) {
    out <- sout
  } else {
    out <- custom.process_out.net(sout)
  }
  
  return(out)
}

#' Processes how output is saved 
#' 
#' This function was introduced in EpiModel 2 to simplify the code in netsim
#' @param dat_list (From `EpiModel` docs) A list of Master objects in netsim simulations.
custom.process_out.net <- function(dat_list) {
  for (s in seq_along(dat_list)) {
    # Set output
    if (s == 1) {
      out <- custom.saveout.net(dat_list[[s]], s)
    } else {
      out <- custom.saveout.net(dat_list[[s]], s, out)
    }
  }
  class(out) <- "netsim"
  
  return(out)
}

#' @title Internal function running the network simulation loop
#'
#' @description This function run the initialization and simulation loop for one
#'              simulation. CNS note: this is not an exported object from the EpiModel package hence
#'              I am pasting it here
#' @source EpiModel2.0
#' @inheritParams custom.initialize.net
#' @keywords internal
netsim_loop <- function(x, param, init, control, s) {
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
  
  return(dat)
}

#' Modification of the EpiModel::saveout.net function
#' 
#' The performed modifications are clearly highlighted.
#' 
#' This function has been modified with the purpose of modifying how the net object is saved
#' so that all of our custom states and their statistics can be saved.
#' 
#' @source EpiModel 2.0 source code
#' INTERNAL PURPOSE MOSTLY
#' @keywords internal
#' 
custom.saveout.net <- function(dat, s, out = NULL) {
  
  # Counts number of simulated networks
  num.nw <- length(dat$nw)
  
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
    
    if (dat$control$tergmLite == FALSE) {
      if (!is.null(dat$stats$transmat)) {
        row.names(dat$stats$transmat) <- 1:nrow(dat$stats$transmat)
        out$stats$transmat <- list(dat$stats$transmat)
      } else {
        out$stats$transmat <- list(data.frame())
      }
      class(out$stats$transmat) <- c("transmat", class(out$stats$transmat))
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
    
    if (dat$control$tergmLite == FALSE) {
      if (!is.null(dat$stats$transmat)) {
        row.names(dat$stats$transmat) <- 1:nrow(dat$stats$transmat)
        out$stats$transmat[[s]] <- dat$stats$transmat
      } else {
        out$stats$transmat[[s]] <- data.frame()
      }
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
    
    if (dat$control$tergmLite == FALSE) {
      names(out$stats$transmat) <- simnames[1:length(out$stats$transmat)]
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
#' @source EpiModel 2.0
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
    
    # Initial Network Simulation ----------------------------------------------
    if (x$edapprox == TRUE) {
      nw <- stats::simulate(x$fit, basis = x$fit$newnetwork,
                            control = control$set.control.ergm)
    } else {
      nw <- x$fit$network
    }
    if (control$resimulate.network == TRUE) {
      nw <- EpiModel::sim_nets(x, nw, nsteps = 1, control)
    } else {
      nw <- EpiModel::sim_nets(x, nw, nsteps = control$nsteps, control)
    }
    nw <- networkDynamic::activate.vertices(nw, onset = 1, terminus = Inf)
    dat$nw[[1]] <- nw
    
    # Network Parameters ------------------------------------------------------
    dat$nwparam <- list(x[-which(names(x) == "fit")])
    groups <- length(unique(tergmLite::get_vertex_attribute(nw, "group")))
    dat <- EpiModel::set_param(dat, "groups", groups)
    
    # Nodal Attributes --------------------------------------------------------
    
    # Standard attributes
    num <- network::network.size(nw)
    dat <- EpiModel::set_attr(dat, "active", rep(1, num), override.length.check = TRUE)
    dat <- EpiModel::set_attr(dat, "entrTime", rep(1, num))
    dat <- EpiModel::set_attr(dat, "exitTime", rep(NA, num))
    
    ## Pull attr on nw to dat$attr
    dat <- EpiModel::copy_nwattr_to_datattr(dat)
    
    ## Store current proportions of attr
    nwterms <- EpiModel::get_network_term_attr(nw)
    if (!is.null(nwterms)) {
      dat$temp$nwterms <- nwterms
      dat$temp$t1.tab <- EpiModel::get_attr_prop(dat, nwterms)
    }
    
    ## Infection Status and Time
    dat <- custom.init_status.net(dat)
    
    # Conversions for tergmLite
    if (control$tergmLite == TRUE) {
      dat <- tergmLite::init_tergmLite(dat)
    }
    
    ## Infection Status and Time Modules
    ##################################################################
    ##                       This is modified                       ##
    ##################################################################
    # dat <- init_status.net(dat)
    dat = custom.init_status.net(dat)
    ##################################################################
    ##                       This is modified                       ##
    ##################################################################
    
    # Summary Stats -----------------------------------------------------------
    dat <- do.call(control[["prevalence.FUN"]],list(dat, at = 1))
    
    
    # Restart/Reinit Simulations ----------------------------------------------
  } else if (control$start > 1) {
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
#' 
#' @source EpiModel 2.0
custom.init_status.net <- function(dat) {
  
  type <- EpiModel::get_control(dat, "type", override.null.error = TRUE)
  nsteps <- EpiModel::get_control(dat, "nsteps")
  tergmLite <- EpiModel::get_control(dat, "tergmLite")
  vital <- EpiModel::get_param(dat, "vital")
  groups <- EpiModel::get_param(dat, "groups")
  status.vector <- EpiModel::get_init(dat, "status.vector", override.null.error = TRUE)
  if (type %in% c("SIS", "SIR") && !is.null(type)) {
    rec.rate <- EpiModel::get_param(dat, "rec.rate")
  }
  if (vital == TRUE) {
    di.rate <- EpiModel::get_param(dat, "di.rate")
  }
  
  # Variables ---------------------------------------------------------------
  i.num <- EpiModel::get_init(dat, "i.num", override.null.error = TRUE)
  if (type  == "SIR" && is.null(status.vector) && !is.null(type)) {
    r.num <- EpiModel::get_init(dat, "r.num")
  }
  
  num <- sum(EpiModel::get_attr(dat, "active") == 1)
  
  if (groups == 2) {
    group <- EpiModel::get_attr(dat, "group")
    i.num.g2 <- EpiModel::get_init(dat, "i.num.g2")
    if (type  == "SIR" && is.null(status.vector) && !is.null(type)) {
      r.num.g2 <- EpiModel::get_init(dat, "r.num.g2", override.null.error = TRUE)
    }
  } else {
    group <- rep(1, num)
  }
  
  statOnNw <- "status" %in% dat$temp$nwterms
  
  # Status ------------------------------------------------------------------
  
  ## Status passed on input network
  if (statOnNw == FALSE) {
    if (!is.null(status.vector)) {
      status <- status.vector
    } else {
      status <- rep("s", num)
      status[sample(which(group == 1), size = i.num)] <- "i"
      if (groups == 2) {
        status[sample(which(group == 2), size = i.num.g2)] <- "i"
      }
      if (type == "SIR"  && !is.null(type)) {
        status[sample(which(group == 1 & status == "s"), size = r.num)] <- "r"
        if (groups == 2) {
          status[sample(which(group == 2 & status == "s"), size = r.num.g2)] <- "r"
        }
      }
    }
    dat <- EpiModel::set_attr(dat, "status", status)
  } else {
    status <- tergmLite::get_vertex_attribute(dat$nw[[1]], "status")
    dat <- EpiModel::set_attr(dat, "status", status)
  }
  
  
  ## Set up TEA status
  if (tergmLite == FALSE) {
    if (statOnNw == FALSE) {
      dat$nw[[1]] <- tergmLite::set_vertex_attribute(dat$nw[[1]], "status", status)
    }
    dat$nw[[1]] <- networkDynamic::activate.vertex.attribute(dat$nw[[1]],
                                                             prefix = "testatus",
                                                             value = status,
                                                             onset = 1,
                                                             terminus = Inf)
  }
  
  
  # Infection Time ----------------------------------------------------------
  ## Set up inf.time vector
  if (!is.null(type)) {
    idsInf <- which(status == "i")
    infTime <- rep(NA, length(status))
    infTime.vector <- EpiModel::get_init(dat, "infTime.vector", override.null.error = TRUE)
    
    if (!is.null(infTime.vector)) {
      infTime <- infTime.vector
    } else {
      # If vital dynamics, infTime is a geometric draw over the duration of infection
      if (vital == TRUE && di.rate > 0) {
        if (type == "SI") {
          infTime[idsInf] <- -stats::rgeom(n = length(idsInf), prob = di.rate) + 2
        } else {
          infTime[idsInf] <- -stats::rgeom(n = length(idsInf),
                                           prob = di.rate +
                                             (1 - di.rate)*mean(rec.rate)) + 2
        }
      } else {
        if (type == "SI" || mean(rec.rate) == 0) {
          # if no recovery, infTime a uniform draw over the number of sim time steps
          infTime[idsInf] <- EpiModel::ssample(1:(-nsteps + 2),
                                               length(idsInf), replace = TRUE)
        } else {
          infTime[idsInf] <- -stats::rgeom(n = length(idsInf), prob = mean(rec.rate)) + 2
        }
      }
    }
    
    dat <- EpiModel::set_attr(dat, "infTime", infTime)
  } else {
    infTime <- rep(NA, num)
    idsInf <- idsInf <- which(status == "i")
    infTime[idsInf] <- 1
    dat <- EpiModel::set_attr(dat, "infTime", infTime)
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
#' 
custom.get_prev.net <- function(dat, at) {
  
  active <- EpiModel::get_attr(dat, "active")
  type <- EpiModel::get_control(dat, "type")
  groups <- EpiModel::get_param(dat, "groups")
  
  # Subset attr to active == 1
  l <- lapply(1:length(dat$attr), function(x) dat$attr[[x]][active == 1])
  names(l) <- names(dat$attr)
  l$active <- l$infTime <- NULL
  
  status <- l$status
  
  ## Subsetting for epi.by control
  eb <- !is.null(dat$control$epi.by)
  if (eb == TRUE) {
    ebn <- EpiModel::get_control(dat, "epi.by")
    ebv <- dat$temp$epi.by.vals
    ebun <- paste0(".", ebn, ebv)
    assign(ebn, l[[ebn]])
  }
  
  if (at == 1) {
    dat$epi <- list()
  }
  
  if (groups == 1) {
    
    dat <- EpiModel::set_epi(dat, "s.num", at, sum(status == "s"))
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        ebn.temp <- paste0("s.num", ebun[i])
        dat <- EpiModel::set_epi(dat, ebn.temp, at, sum(status == "s" &
                                                          get(ebn) == ebv[i]))
      }
    }
    
    dat <- EpiModel::set_epi(dat, "i.num", at, sum(status == "i"))
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        ebn.temp <- paste0("i.num", ebun[i])
        dat <- EpiModel::set_epi(dat, ebn.temp, at, sum(status == "i" &
                                                          get(ebn) == ebv[i]))
      }
    }
    
    # REMOVED AS TYPE WILL ALWAYS BE NULL
    #if (!is.null(type)) {
    #   if (type == "SIR") {
    #     dat <- EpiModel::set_epi(dat, "r.num", at, sum(status == "r"))
    #     if (eb == TRUE) {
    #       for (i in 1:length(ebun)) {
    #         ebn.temp <- paste0("r.num", ebun[i])
    #         dat <- EpiModel::set_epi(dat, ebn.temp, at, sum(status == "r" &
    #                                                           get(ebn) == ebv[i]))
    #       }
    #     }
    #   }
    # }
    dat <- EpiModel::set_epi(dat, "num", at, length(status))
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        ebn.temp <- paste0("num", ebun[i])
        dat <- EpiModel::set_epi(dat, ebn.temp, at, sum(get(ebn) == ebv[i]))
      }
    }
    
    #################################################################
    ##                           my code                           ##
    #################################################################
    
    dat <- EpiModel::set_epi(dat, "h.num", at, sum(status == "h"))
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        ebn.temp <- paste0("h.num", ebun[i])
        dat <- EpiModel::set_epi(dat, ebn.temp, at, sum(status == "h" &
                                                          get(ebn) == ebv[i]))
      }
    }
    
    dat <- EpiModel::set_epi(dat, "q.num", at, sum(status == "q"))
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        ebn.temp <- paste0("q.num", ebun[i])
        dat <- EpiModel::set_epi(dat, ebn.temp, at, sum(status == "q" &
                                                          get(ebn) == ebv[i]))
      }
    }
    
    dat <- EpiModel::set_epi(dat, "f.num", at, sum(status == "f"))
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        ebn.temp <- paste0("f.num", ebun[i])
        dat <- EpiModel::set_epi(dat, ebn.temp, at, sum(status == "f" &
                                                          get(ebn) == ebv[i]))
      }
    }
    
    dat <- EpiModel::set_epi(dat, "r.num", at, sum(status == "r"))
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        ebn.temp <- paste0("r.num", ebun[i])
        dat <- EpiModel::set_epi(dat, ebn.temp, at, sum(status == "r" &
                                                          get(ebn) == ebv[i]))
      }
    }
    
    dat <- EpiModel::set_epi(dat, "e.num", at, sum(status == "e"))
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        ebn.temp <- paste0("e.num", ebun[i])
        dat <- EpiModel::set_epi(dat, ebn.temp, at, sum(status == "e" &
                                                          get(ebn) == ebv[i]))
      }
    }
    
    
    #################################################################
    ##                           my code                           ##
    #################################################################
    
    
  } else {
    stop("Bipartite networks not currently supported in this custom function")
  }
  return(dat)
}




