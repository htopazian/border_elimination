# Calibrating PfPRs and EIRs for basic units with no interventions. 

# libraries
source("./02_code/packages_paths.R")

# HPC set-up -------------------------------------------------------------------
setwd(HPCpath)

# to edit HPC username and password below
# usethis::edit_r_environ()

src <- conan::conan_sources(c(
  "github::mrc-ide/malariasimulation@feat/metapopulation_model",
  "github::mrc-ide/cali"))

ctx <- context::context_save(path = paste0(HPCpath, "contexts"),
                             sources = c(paste0(HPCpath, "function_run_metapopulation_casestudies.R")),
                             packages = c("dplyr", "malariasimulation", "cali"),
                             package_sources = src)

share <- didehpc::path_mapping("malaria", "M:", "//fi--didenas1/malaria", "M:")

config <- didehpc::didehpc_config(credentials = list(
                          username = Sys.getenv("DIDE_USERNAME"),
                          password = Sys.getenv("DIDE_PASSWORD")),
                          shares = share,
                          use_rrq = FALSE,
                          cores = 1,
                          cluster = "fi--didemrchnb", # "fi--dideclusthn", "fi--didemrchnb"
                          # "GeneralNodes", "12Core", "16Core", 
                          # "12and16Core", "20Core", "24Core", "32Core"
                          template = "GeneralNodes", 
                          parallel = FALSE)

obj <- didehpc::queue_didehpc(ctx, config = config)

# Set up your job --------------------------------------------------------------

# choose pfpr values of interest
pfpr <- tibble(pfpr = seq(0.05, 0.80, 0.05))

# choose 50 unique parameter draws
set.seed(1234)
rsample <- sample(x = c(1:1000), size = 50, replace = FALSE, prob = NULL)
rsample <- tibble(draw = rsample)

# save random draws
saveRDS(rsample, paste0(path, "03_output/draws.rds"))
saveRDS(rsample, paste0(HPCpath, "draws.rds"))


# Run tasks --------------------------------------------------------------------
# create combination of pfpr and draws
index <- crossing(pfpr, rsample)

# remove ones that have already been run
index <- index |> mutate(f = paste0("./calibrate/PRmatch_", pfpr, "_", draw, ".rds")) |>
  mutate(exist = case_when(file.exists(f) ~ 1, !file.exists(f) ~ 0)) |>
  filter(exist == 0) |>
  select(-f, -exist)

# submit all remaining tasks
t <- obj$enqueue_bulk(index, cali_pfpr)
t$status()


# Processing -------------------------------------------------------------------
# read in results
files <- list.files(path = paste0(HPCpath, "calibrate"), 
                    pattern = "*.rds", full.names = TRUE)
dat_list <- lapply(files, function (x) readRDS(x))

# concatenate
match <-  do.call("bind_rows", dat_list) |> as_tibble()

summary(match$starting_EIR)

# save EIR estimates, both on local machine and on shared drive
saveRDS(match, paste0(path, "03_output/EIRestimates.rds"))
saveRDS(match, paste0(HPCpath, "EIRestimates.rds"))

