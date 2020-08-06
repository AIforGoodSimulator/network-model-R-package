#' Quick fn to calculate max.possible amount of edges given some number
#' of nodes n
#' @param n number of nodes
calculate.max.edges = function(n){
  return (n*(n-1)/2)
}

#' Network simulation wrapper
#' 
#' This function simulates a SEIQHRF epidemic on top of a 
#' network of n individuals with a specific list of parameters as provided.
#' S - Susceptible
#' E - Exposed (also Infected)
#' I - Infected 
#' Q - Quarantined (also Infected)
#' H - Requiring hospitalisation (also Infected)
#' R - Recovered
#' F - Fatality (due to infection)
#' 
#' @param n numbers of individuals, a random sample of this size will be taken from the Camp's population
#' @param nsims number of network simulations
#' @param nsteps number of timesteps (days, hours, weeks, arbitrary choice), the rates need to be defined in units
#' relative to this metric (e.g number of cases/day or /hour or /week)
#' @param prop.isobox Proportion of individuals in iso boxes, in this first implementation only two kinds of housing units 
#' are considered (isoboxes and tent) -- can (and should) be changed in the future
#' @param iso.capacity Capacity of iso boxes
#' @param tent.capacity Capacity of tents 
#' @param verbose whether function should be verbose in its command line output
#' @param cores number of cores to be used (default: number of available cores - 1)
#' @param d.rate departure rate (includes external deaths and actual departures) -- if set to 0 no invidual is consider to leave
#' or die from causes different than the epidemic agent 
#' @param external.contacts Average number of external contacts that individuals make every day (external meaning from outside
#'  their tent)
#' @param act.rate.se Act rate (contact rate) between susceptible and exposed individuals 
#' -- default: 10 interactions of S-to-E/timestep/person
#' @param inf.prob.se infection probability of susceptible individuals upon contact with exposed individuals, these value
#' is often considered to be lower than that upon contact with infected individuals (though literature needs to be 
#' checked for this) -- default: 0.02 (Out of every 100 contacts between S and E we expect 2 infections)
#' @param act.rate.si Act rate (contact rate) between susceptible and infected individuals 
#' -- default: 10 interactions/timestep /person
#' @param inf.prob.si infection probability of susceptible individuals upon contact with susceptible individuals -- default: 0.05
#' @param act.rate.sq Act rate (contact rate) between susceptible and quarantined individuals 
#' -- default: 2.5 interactions/timestep/person
#' @param inf.prob.sq infection probability of susceptible individuals upon contact with quarantined individuals -- default: 0.05
#' @param ei.rate Rate (can be seen as flow rate or probability) between exposed and infected individuals,
#' can take the value of a vector of length nsteps to simulate timevarying scenarios
#' @param iq.rate Rate between infected and quarantined states (~chances of quarantining when infected)
#' @param ih.rate Rate between infected and hospitalised states (~chances of requiring 
#' hospitalisation when infected)
#' @param qh.rate Rate between quarantined and hospitalised states (~ chances of requiring hospitalisation 
#' when quarantining)
#' @param hr.rate Rate between requiring hospitalisation and recovered states (~chances of recovering
#' when requiring hospitalisation)
#' @param qr.rate Rate between quarantined and recovered states (~chances of recovering when
#' in quarantined state)
#' @param ir.rate Rate between infected and recovered states (~chances of recovering when
#' in infected state)
#' @param hf.rate Rate between hospitalised and fatality state (~chances of dying from the infectious agent
#' given requiring hospitalisation)
#' @param hosp.cap Hospital capacity 
#' @param hf.rate.overcap Same as hf.rate but when the hospital is at full capacity (As per hosp.cap)
#' @param hosp.tcoeff Increase in the hr.rate for every timestep spent in the hospital
#' @param a.rate arrival rate -> if set to 0, inflow of people is not considered (i.e. no births, no 
#' immigration)
#' @param di.rate departure rate for infected individuals (defaults to d.rate)
#' @param ds.rate departure rate for susceptible individuals (defaults to d.rate)
#' @param dr.rate departure rate for susceptible individuals (defaults to d.rate)
#' @param i.num Initial number of infected individuals (defanetults to 1)
#' @param r.num Initial number of recovered individuals (defaults to 0)
#' @param e.num Initial number of exposed individuals (default to 0)
#' @param s.num Initial number of susceptible individuals (defaults to n-1)
#' @param f.num Initial number of fatalities (defaults to 0)
#' @param h.num Initial number of individuals requiring hospitalization (defaults to 0)
#' @param q.num Initial number of quarantined individuals (defaults to 0)
#' @param agedistribution vector containing the age distribution of 
#' individuals in the camp. By default this will be loaded internally,
#' see `CampNetworkSimulator::agedistribution`
#' @param ratesbyage dataframe containing the hospitalisation rate and
#' fatality rate by age, by default `CampNetworkSimulator::ratesbyage` is used
#'
#' @examples 
#' 
#' net_simulate(n= 100, nsims = 1, nsteps = 10)
#' @importFrom stats rgeom setNames simulate time
#' @importFrom network %n%
#' @export
net_simulate = function(
  n = 100,
  nsims = 3,
  nsteps = 90,
  prop.isobox = 0.43,
  iso.capacity = 10,
  tent.capacity = 4,
  verbose = F,
  cores = 1,
  d.rate = 0,
  external.contacts = 4,
  act.rate.se = 10,
  inf.prob.se = 0.02,
  act.rate.si = 10,
  inf.prob.si = 0.05,
  act.rate.sq = 2.5,
  inf.prob.sq = 0.05,
  ei.rate = 1/10,
  iq.rate = 1/30, #c(rep(1/30, 60), rep(15/30, 120)), # time varying works
  ih.rate = 1/100,
  qh.rate = 1/100,
  hr.rate = 1/15,
  qr.rate = 1/20,
  ir.rate = 1/20,
  hf.rate = 1/50,
  hf.rate.overcap = 1/25,
  hosp.cap = 5,
  hosp.tcoeff = 0.5,
  a.rate = 0,
  di.rate = d.rate,
  ds.rate = d.rate,
  dr.rate = d.rate,
  i.num = 1,
  r.num = 0,
  e.num = 0,
  s.num = n - 1,
  f.num = 0,
  h.num = 0,
  q.num = 0,
  agedistribution = CampNetworkSimulator::agedistribution,
  ratesbyage = CampNetworkSimulator::ratesbyage
) {
  
  # initialise log_build, in which all the build information will be stored
  log_build = list()
  
  # checking if params inputs make sense
  if(i.num + r.num + e.num + s.num + f.num + h.num + q.num > n){
    stop('The sum of individuals in each initial state is greater than the number of individuals
         in the network (n). Please check the numbers and re-run the function')
  }
  if (any(c(i.num, r.num, e.num, s.num, f.num, h.num, q.num)> n)){
    stop('At least one of the initial number in one of the states is greater than the number
         of individuals in the network.')
  }
  
  # Population parameters ####
  
  # Housing parameters ####
  prop.tent = 1 - prop.isobox
  
  num_of_iso = round(n*prop.isobox/iso.capacity) # number of isoboxes
  num_of_tents = round(n*prop.tent/tent.capacity) # number of tents
  
  if (n > num_of_iso*iso.capacity+num_of_tents*tent.capacity){
    num_of_tents = num_of_tents+1
  }
  
  num_in_iso = num_of_iso*iso.capacity
  num_in_tents = n - num_in_iso # num_of_tents*tent.capacity 
  
  msg1 = paste0(
    "Number of tents: ", num_of_tents, "\n",
    "Number of isoboxes: ", num_of_iso, "\n",
    "Number of people in tents: ", num_in_tents, "\n",
    "Number of people in isoboxes: ", num_in_iso, "\n\n",
    "Number of housing units in total: ", num_of_iso+num_of_tents)
  if (verbose){
    cat(msg1)
  }
  log_build$messages = msg1
  
  # Housing allocation #### 
  iso_ids = 1:num_of_iso
  tent_ids = 1:num_of_tents 
  
  housing_iso = EpiModel::apportion_lr(
    vector.length = num_in_iso,
    values = iso_ids,
    proportions = rep(1/(num_of_iso), num_of_iso)
  )
  housing_iso = paste0("iso", housing_iso)
  
  
  housing_tents = EpiModel::apportion_lr(
    vector.length = num_in_tents,
    values = tent_ids,
    proportions = rep(1/(num_of_tents), num_of_tents)
  )
  housing_tents = paste0("tent", housing_tents)
  
  
  housing = c(housing_iso,housing_tents)
  
  # residence vector to keep track of those in tents and isos
  residence = c(rep("iso", num_in_iso), rep("tent", num_in_tents))
  
  # initialise network and assign attributes
  nw = network::network.initialize(n = n, directed = FALSE)
  
  # housing
  nw = network::set.vertex.attribute(nw, "housing", housing)
  # age
  nw = network::set.vertex.attribute(nw, "age",
                                     sample(as.vector(agedistribution),n))
  
  # Setting network formation dynamics ####
  formation = ~edges + 
    offset(nodematch("housing", diff = F))
  
  ## max in-house ties
  max.inhouse.edges = 0
  for (num_in_house in table(housing)){
    max.inhouse.edges = max.inhouse.edges +
      calculate.max.edges(num_in_house)
  }
  
  # default degrees in housing units as per current occupancy
  ## The degree of a node in a network is the number of connections it has to other nodes
  # the -1 is tou account fr the lack of connection to one-self
  iso.default.degree = mean(table(housing)[iso_ids]-1)
  tent.default.degree = mean(table(housing)[tent_ids]-1)
  
  # number of external contacts per person in average
  # external.contacts = 4
  
  mean_degree.iso =  iso.default.degree + external.contacts
  mean_degree.tent =  tent.default.degree + external.contacts
  
  # calculate mean degree
  mean_degree = (num_in_iso*mean_degree.iso+
                   num_in_tents*mean_degree.tent)/n
  
  expected.edges = n*mean_degree/2
  
  target.stats = c(expected.edges)
  
  coef.diss = EpiModel::dissolution_coefs(
    dissolution = ~offset(edges)+
      offset(nodematch("housing", diff = F)),
    duration = c(2, 1e9),
    d.rate = d.rate)
  
  # Estimate network dynamics #### 
  cat("Network estimation step")
  est <- EpiModel::netest(nw,
                          formation,
                          target.stats,
                          coef.diss,
                          coef.form = Inf,
                          set.control.ergm = ergm::control.ergm(MCMLE.maxit = 500)
  )
  
  log_build$network_estimation = est
  
  # Network diagnostics
  cat('Getting network diagnostics\n')
  dx <- EpiModel::netdx(est,
                        nsims = 1e3,
                        nsteps = 90,
                        ncores = cores,
                        dynamic = FALSE,
                        nwstats.formula = ~edges + nodematch("housing", diff = FALSE),
                        set.control.ergm = ergm::control.simulate.ergm(MCMC.burnin = 1e6),
                        keep.tnetwork = T)
  log_build$network_diagnostics = dx
  
  # Setting up epidemic parameters and controls ####
  param = EpiModel::param.net(act.rate.se = act.rate.se,
                              inf.prob.se = inf.prob.se,
                              act.rate.si = act.rate.si,
                              inf.prob.si = inf.prob.si,
                              act.rate.sq = act.rate.sq,
                              inf.prob.sq = inf.prob.sq,
                              ei.rate = ei.rate,
                              iq.rate = iq.rate, #c(rep(1/30, 60), rep(15/30, 120)), # time varying works
                              ih.rate = ih.rate,
                              qh.rate = qh.rate,
                              hr.rate = hr.rate,
                              qr.rate = qr.rate,
                              ir.rate = ir.rate,
                              hf.rate = hf.rate,
                              hf.rate.overcap = hf.rate.overcap,
                              hosp.cap = hosp.cap,
                              hosp.tcoeff = hosp.tcoeff,
                              a.rate = a.rate,
                              di.rate = d.rate,
                              ds.rate = d.rate,
                              dr.rate = d.rate,
                              ratesbyAge = ratesbyage
  ) 
  
  init = EpiModel::init.net(
    i.num = i.num,
    r.num = r.num,
    e.num = e.num,
    s.num = s.num,
    f.num = f.num,
    h.num = h.num,
    q.num = q.num
  )
  
  control = EpiModel::control.net(
    nsims = nsims, 
    nsteps = nsteps,
    # delete.nodes = T,  this does not work for now
    ncores = cores,
    initialize.FUN = custom.initialize.net, # this bit is just so that I can extract time
    exposure.FUN = exposure,
    infect.FUN = infect,
    epi.by = "age",
    quarantine.FUN = quarantining,
    hospitalize.FUN = RequireHospitalization,
    recover.FUN = recover,
    fatality.FUN = fatality,
    recovery.FUN = NULL,
    infection.FUN = NULL,
    departures.FUN = EpiModel::departures.net,
    get_prev.FUN = custom.get_prev.net,
    skip_check = FALSE,
    depend = T
  )
  
  cat('Running simulation, this may take a while...\n')
  t0 = Sys.time()
  sim = custom.netsim(est, param, init, control)
  
  cat(Sys.time()-t0) # roughly 30 minutes for n = 1000 and nsteps = 60
  log_build$network_simulation_object = sim
  # res = as.data.frame(sim)
  
  return(log_build)
  
}


