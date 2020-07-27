#' Utility plotting functions


#' Plot epidemic curves for each state
#' 

plot.states = function(net_build_object){
  res = as.data.frame(net_build_object$network_simulation)
  ggplotly(
    res %>% select(s.num, e.num, i.num, q.num, h.num, r.num, f.num, num, time) %>%
      group_by(time) %>% summarise_all(~mean(.)) %>% 
      pivot_longer(-time) %>% ggplot(aes(x = time, y = value, color = name))+
      geom_line(size = 1)+scale_color_brewer(palette = "Set1")
  )
}

#' Get network heatmap at a given time step
#' 
#' 

network.heatmap = function(network_build_object, at = 1){
  
  nw_object = get_network(network_build_object$network_simulation_object)
  net_at_1 = network.collapse(nw_object, at = 2)
  
  graph = intergraph::asIgraph(net_at_1)
  
  adj = igraph::as_adjacency_matrix(graph, sparse = F)
  colnames(adj) = igraph::V(graph)$housing
  rownames(adj) = igraph::V(graph)$housing
  
  adj = adj[order(rownames(adj)), order(colnames(adj))]
  
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

network.heatmap.GIF = function(){
  library(animation)
  
  nw_object = get_network(sim)
  ani.record(reset = TRUE)  # clear history before recording
  for (at in c(1:30)){
    
    net_at = network.collapse(nw_object, at = at)
    graph = intergraph::asIgraph(net_at)
    adj = igraph::as_adjacency_matrix(graph, sparse = F)
    
    colnames(adj) = igraph::V(graph)$housing
    rownames(adj) = igraph::V(graph)$housing
    
    adj = adj[order(rownames(adj)), order(colnames(adj))]
    
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
    ani.record()
  }
  
  oopts = ani.options(interval = 0.5)
  saveGIF(ani.replay(), img.name = 'networkGIF.gif')## this line of code is wrong 
}


