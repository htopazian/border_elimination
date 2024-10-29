# packages ----

# spatial
library(sf)
library(geodata)           # gadm
library(malariaAtlas)      # malaria atlas project
library(units)
library(rmapshaper)
library(gdistance)
library(raster)
library(exactextractr)

# data manipulation
library(data.table)
library(readxl)
library(splines)
library(tidyverse)
library(patchwork)
library(MetBrewer)         # devtools::install_github("BlakeRMills/MetBrewer")
library(netz)              # https://github.com/mrc-ide/netz
library(igraph)
library(survey)
library(MASS)
library(quantreg)

# DHS
library(rdhs)              # https://github.com/ropensci/rdhs

# models
library(malariasimulation) # update documentation - devtools::document()
                           # then rebuild package - devtools::build()
# or download branch: (commit a591eb0)
# devtools::install_github("mrc-ide/malariasimulation@feat/metapop_asym", force = TRUE)

# site files
# library(foresite)          # https://github.com/mrc-ide/foresite
# library(site)              # https://github.com/mrc-ide/site

# or download branch:
# devtools::install_github("mrc-ide/site@feat/update_vax_function", force = TRUE)

# colors ----
# met.brewer("Hiroshige")

# paths ----
# high performance computing cluster
HPCpath <- "M:/Hillary/border_elimination/"

# local / project location
path <- "C:/Users/htopazia/OneDrive - Imperial College London/Github/border_elimination_private/"

# HPC packages
library(hipercow) # https://mrc-ide.github.io/hipercow/articles/hipercow.html
# hipercow_init(driver = "windows")
# windows_check()

# print session info
devtools::session_info()
writeLines(capture.output(devtools::session_info()), "sessionInfo.txt")
