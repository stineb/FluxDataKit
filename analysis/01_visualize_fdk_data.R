library(tidyverse)
library(FluxDataKit)

source("R/fdk_plot.R")


files <- list.files("/data/scratch/PLUMBER_X/","*Flux.nc", full.names = TRUE)

lapply(files, function(file){
  try(fdk_plot(
    file = file
    )
  )
  invisible()
})