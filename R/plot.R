#' Plot epidemic curves for each state
#' 
#' @param net_build_object object returned by `net_simulate`
#' 
#' @examples
#' data(net_example)
#' plot_states(net_example)
#' @importFrom rlang .data
#' @export
plot_states = function(net_build_object){
  res = as.data.frame(net_build_object$network_simulation)
  plotly::ggplotly(
    res[, c('s.num', 'e.num', 'i.num', 'q.num',
            'h.num', 'r.num', 'f.num', 'num', 'time')] %>% 
      dplyr::group_by(.data$time) %>%
      dplyr::summarise_all(~mean(.)) %>% 
      tidyr::pivot_longer(-c(.data$time)) %>%
      ggplot2::ggplot(
        ggplot2::aes(x = .data$time,
                     y = .data$value,
                     color = .data$name))+
      ggplot2::geom_line(size = 1)+
      ggplot2::scale_color_brewer(palette = "Set1")+
      ggplot2::xlab('Simulation time-steps')+
      ggplot2::ylab('Number of individuals')+
      ggplot2::labs(color = 'State')
  )
}

#' Get distribution of time-metrics
#' 
#' For diagnostic purposes, the user may wish to explore metrics such as the
#' distribution of the incubation period (transition between exposed and infected), or recovery time,
#' or time-to-failure/fatality and compare this to the real world estimates. For example, 
#' it is known that the COVID19 incubation period is anywhere between 5 and 14 days.
#' 
#' The network_simulation object contain several time features for each individual, 
#' namely, the time to exposure, the time to infection, etc. This function is used
#' to calculate common metric that may be chosen to be evaluated. It is mostly 
#' for internal usage.
#' 
#' Copied from Tim Churches' blogpost.
#' 
#' For more info on tidy evaluation check https://dplyr.tidyverse.org/articles/programming.html#eliminating-r-cmd-check-notes
#' (i.e. why there is so much .data in this function source)
#' 
#' @param net_build_object object returned by `net_simulate`
#' @examples 
#' times_df = get_times(net_example)
#' @importFrom rlang .data
#' @export
get_times <- function(net_build_object) {
  
  sim <- net_build_object$network_simulation
  
  for (s in 1:sim$control$nsims) {
    if (s == 1) {
      times <- sim$times[[paste0("sim", s)]]
      times <- times %>% dplyr::mutate(s = s)
    } else {
      times <- times %>% 
        dplyr::bind_rows(
          sim$times[[paste("sim", s, sep = "")]] %>%
            dplyr::mutate(s = s))
    }
  }
  
  times <- times %>%
    dplyr::mutate(infTime = ifelse(.data$infTime < 0, -5, .data$infTime),
                  expTime = ifelse(.data$expTime < 0, -5, .data$expTime)) %>% 
    dplyr::mutate(incubation_period = .data$infTime - .data$expTime,
                  illness_duration  = .data$recTime - .data$expTime,
                  illness_duration_hosp = .data$dischTime - .data$expTime, 
                  hosp_los              = .data$dischTime - .data$hospTime,
                  quarantine_delay      = .data$quarTime  - .data$infTime,
                  survival_time         = .data$fatTime   - .data$infTime) %>% 
    dplyr::select(.data$s,
                  .data$incubation_period,
                  .data$quarantine_delay,
                  .data$illness_duration, 
                  .data$illness_duration_hosp,
                  .data$hosp_los,
                  .data$survival_time) %>% 
    tidyr::pivot_longer(-.data$s,
                        names_to = "period_type",
                        values_to = "duration") %>% 
    dplyr::mutate(period_type = factor(
        .data$period_type,
        levels = c("incubation_period", 
                   "quarantine_delay",
                   "illness_duration",
                   "illness_duration_hosp", 
                   "hosp_los",
                   "survival_time"),
        labels = c("Incubation period", 
                   "Delay entering isolation",
                   "Illness duration",
                   "Illness duration (hosp)", 
                   "Hospital care required duration",
                   "Survival time of case fatalities"), 
        ordered = TRUE))
  return(times)
}


#' Plotting time-metric distributions 
#' 
#' Plotting function to visually evaluate the overall goodness of fit 
#' of the network epidemic simulation with respect to well known 
#' external metrics (incubation period, time-to-recovery)
#' @param net_build_object object returned by `net_simulate`
#' @param interactive TRUE/FALSE whether output plot should be interactive (plotly) or
#' static (ggplot)
#' @examples 
#' plot_times(net_example)
#' @export
plot_times = function(net_build_object, interactive = TRUE){
  
  times = CampNetworkSimulator::get_times(net_build_object)
  
  plt = times %>% 
    dplyr::filter(.data$duration <= 30) %>% 
    ggplot2::ggplot(ggplot2::aes(x = .data$duration)) + 
    ggplot2::geom_density() + 
    ggplot2::facet_wrap(period_type ~ ., scales = "free_y") + 
    ggplot2::labs(title = "Duration frequency distributions",
                  subtitle = "Baseline simulation")
  
  return(if (interactive) plotly::ggplotly(plt) 
         else plt)
}

#' Get network heatmap at a given time step
#' 
#' @param net_build_object objected returned by `net_simulate`
#' @param at timestep
#' @examples
#' data(net_example)
#' network.heatmap(net_example, at = 2)
#' @export
network.heatmap = function(net_build_object, at = 1){
  
  nw_object = EpiModel::get_network(net_build_object$network_simulation)
  net_at_1 = networkDynamic::network.collapse(nw_object, at = 2)
  
  graph = intergraph::asIgraph(net_at_1)
  
  adj = igraph::as_adjacency_matrix(graph, sparse = F)
  colnames(adj) = igraph::V(graph)$housing
  rownames(adj) = igraph::V(graph)$housing
  
  adj = adj[order(rownames(adj)), order(colnames(adj))]
  n = nrow(adj)
  pheatmap::pheatmap(adj,
                     color = c("grey50","black"),
                     border_color = "white",
                     angle_col = 45,
                     angle_row = 45,
                     fontsize = 6,
                     legend_breaks = c(0,1),
                     legend = F,
                     cluster_rows = F,
                     cluster_cols = F,
                     show_rownames = ifelse(n<100, T, F),
                     show_colnames = ifelse(n<100, T, F))
  
}


#' Get network heatmap at all timesteps as a GIF
#' 
#' @param net_build_object objected returned by `net_simulate`
#' @param interval interval of timesteps to produce GIF in
#' @param save TRUE/FALSE whether or not to save the GIF, if TRUE the GIF will be saved as
#' network.gif in the current working directory.
#' @examples
#' \dontrun{
#' data(net_example)
#' network.heatmap(net_example, at = 2)
#' }
#' @export
network.heatmap.GIF = function(net_build_object, interval = NULL, save = FALSE){
  
  nw_object = EpiModel::get_network(net_build_object$network_simulation)
  animation::ani.record(reset = TRUE)  # clear history before recording
  
  if (is.null(interval)){ # as many frames as simulation timesteps
    interval = nrow(as.data.frame(net_build_object$network_simulation))
  }
  
  for (at in interval){
    net_at = networkDynamic::network.collapse(nw_object, at = at)
    graph = intergraph::asIgraph(net_at)
    adj = igraph::as_adjacency_matrix(graph, sparse = F)
    
    colnames(adj) = igraph::V(graph)$housing
    rownames(adj) = igraph::V(graph)$housing
    
    adj = adj[order(rownames(adj)), order(colnames(adj))]
    n = nrow(adj)
    pheatmap::pheatmap(
      adj,
      color = c("grey50","black"),
      border_color = "white",
      angle_col = 45,
      angle_row = 45,
      fontsize = 6,
      legend_breaks = c(0,1),
      legend = F,
      cluster_rows = F,
      cluster_cols = F,
      show_rownames = ifelse(n<100, T, F),
      show_colnames = ifelse(n<100, T, F)
    )
    animation::ani.record()
  }
  
  oopts = animation::ani.options(interval = 0.5)
  
  if (save){
    animation::saveGIF(animation::ani.replay(),
                       movie.name = 'network.gif')
    print(paste0('GIF created at ', getwd(), 'network.gif'))
  } else {
    animation::ani.replay()
  }
  
}


