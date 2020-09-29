## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=F, warning=F----------------------------------------------
library(CampNetworkSimulator)

## -----------------------------------------------------------------------------
plot_states(net_example)

## ---- eval = F----------------------------------------------------------------
#  net_obj = net_simulate()
#  plot_states(net_obj)

## -----------------------------------------------------------------------------
network.heatmap(net_example, at = 5)

## -----------------------------------------------------------------------------
network.heatmap.GIF(net_example) # from first time-step to end

## -----------------------------------------------------------------------------
network.heatmap.GIF(net_example, interval = 5:10) # from timestep 5 to timestep 10

