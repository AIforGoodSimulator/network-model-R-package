## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, warning=F, message=F----------------------------------------------
library(CampNetworkSimulator)

## -----------------------------------------------------------------------------
results_df = as.data.frame(net_example$network_simulation)
head(results_df)

## ----eval = F-----------------------------------------------------------------
#  write.csv(results_df, 'network_simulation_results.csv')

