# Run sites through malariasimulation.

# libraries
source("./02_code/packages_paths.R")

# pull in shapefiles from A_shapefiles_centroids.R
centroids_border <- readRDS("./03_output/centroids_border.rds")

# HPC set-up -------------------------------------------------------------------
setwd(HPCpath)

# set source function
hipercow_environment_create(sources = "function_run_metapopulation_sims.R")

# install packages from pkgdepends.txt
hipercow_provision()

# Set up your job --------------------------------------------------------------
cluster <- c(1:nrow(centroids_border))
mix_type <- 2 # 1 = distance, 2 = travel time
draw <- readRDS(paste0(HPCpath, "draws.rds")) # read in parameter draws
strategy <- c("baseline", "test&treat")

scenarios <- crossing(cluster, mix_type, draw, strategy)

# Run tasks --------------------------------------------------------------------
index <- scenarios

saveRDS(index, paste0(HPCpath, "index_main_runs.rds"))

# remove ones that have already been run
index <- index |> mutate(f = paste0("./results/", cluster, "_", mix_type, "_", draw, "_", strategy, ".rds")) |>
  mutate(exist = case_when(file.exists(f) ~ 1, !file.exists(f) ~ 0)) |>
  filter(exist == 0) |>
  dplyr::select(-f, -exist)

# submit jobs
id <- task_create_expr(run_metapopulation_sim(1, 2, 4, "baseline"))
task_status(id)
task_log_show(id)

# in bulk
bundle <- task_create_bulk_expr(run_metapopulation_sim(cluster, mix_type, draw, strategy), index[1:100,])
hipercow_bundle_status(bundle)

