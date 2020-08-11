#' Find edges that are in touch and are in different status
#'
#'
#' DEPRECATED
#' The EpiModel discordedgelist function finds susceptible individuals that are in touch with exposed individuals.
#' In our case, because susceptible individuals can be in touch with exposed, quarantined or infected individuals
#' we have extended the function to be able to match susceptible individuals to individuals in either of those states.
#'
#' @param dat native simulation object from EpiModel
#' @param at timestep (also native to EpimModel)
#' @param from that susceptible individuals will be matched to (options: i, q, e)
#' @export
custom_discord_edgelist.old <- function(dat, at, from = "i") {
  
  .Deprecated('custom.netsim',
              msg = 'This function should only be used if you are 
              using EpiModel 1.x for EpiModel 2.x use custom.netsim')
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

#' Modification of the EpiModel::netsim function
#' 
#' DEPRECATED (use custom.netsim instead)
#' The performed modifications are clearly highlighted.
#' 
#' INTERNAL PURPOSE MOSTLY
#' @keywords internal
#' @importFrom foreach %dopar%
custom.netsim.old <- function(x, param, init, control) {
  
  .Deprecated('custom.netsim',
              msg = 'This function should only be used if you are 
              using EpiModel 1.x for EpiModel 2.x use custom.netsim')
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
custom.saveout.net.old <- function(dat, s, out = NULL) {
  
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
#' DEPRECATED!!
#' The performed modifications are clearly highlighted.
#' 
#' This function has been modified with the purpose of modifying how the net object is saved
#' so that all of our custom states and their statistics can be saved.
#' 
#' INTERNAL PURPOSE MOSTLY
#' @keywords internal
custom.initialize.net.old <- function(x, param, init, control, s) {
  
  
  .Deprecated('custom.initialize.net',
              msg = 'This function should only be used if you are 
              using EpiModel 1.x for EpiModel 2.x use custom.initialize.net')
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

custom.init_status.net.old <- function(dat) {
  
  .Deprecated('custom.init_status.net',
              msg = 'This function should only be used if you are 
              using EpiModel 1.x for EpiModel 2.x use custom.init_status.net')
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
#' 
custom.get_prev.net.old <- function(dat, at) {
  
  .Deprecated('custom.get_prev.net',
              msg = 'This function should only be used if you are 
              using EpiModel 1.x for EpiModel 2.x use custom.get_prev.net')
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


