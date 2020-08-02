#' Plot epidemic curves for each state
#' 
#' @param net_build_object objected returned by `net_simulate`
#' @export
plot_states = function(net_build_object){
  res = as.data.frame(net_build_object$network_simulation)
  plotly::ggplotly(
    res[, c('s.num', 'e.num', 'i.num', 'q.num',
            'h.num', 'r.num', 'f.num', 'num', 'time')] %>% 
      dplyr::group_by('time') %>%
      dplyr::summarise_all(~mean(.)) %>% 
      tidyr::pivot_longer(-c('time')) %>%
      ggplot2::ggplot(ggplot2::aes(x = 'time', y = 'value', color = 'name'))+
      ggplot2::geom_line(size = 1)+
      ggplot2::scale_color_brewer(palette = "Set1")
  )
}

#' Get network heatmap at a given time step
#' 
#' @param net_build_object objected returned by `net_simulate`
#' @param at timestep
#' @export
network.heatmap = function(net_build_object, at = 1){
  
  nw_object = EpiModel::get_network(net_build_object$network_simulation_object)
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
#' @export
network.heatmap.GIF = function(net_build_object, interval = NULL){
  
  nw_object = EpiModel::get_network(net_build_object$network_simulation)
  animation::ani.record(reset = TRUE)  # clear history before recording
  
  if (is.null(interval)){
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
  animation::saveGIF(animation::ani.replay(),
                     movie.name = 'network.gif')
  print('GIF created.')
}


