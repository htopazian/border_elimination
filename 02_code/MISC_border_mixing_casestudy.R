# Running a two unit example, varying mixing values. 

# Before running this case study, MISC_calibrate.R must be run to generate 
# a set of EIR values for each PfPR and parameter draw set

# libraries
source("./02_code/packages_paths.R")

# HPC set-up -------------------------------------------------------------------
setwd(HPCpath)

# set source function
hipercow_environment_create(sources = "function_run_metapopulation_casestudies.R")

# install packages from pkgdepends.txt
hipercow_provision()


# Set up your job --------------------------------------------------------------
# choose PfPR combinations

# international mixing summaries from Results_section.R
# > summary(international_mixing * 100)
# travel         
# Min.   : 0.000002  
# 1st Qu.: 0.656090  
# Median : 1.274075  
# Mean   : 1.671202  
# 3rd Qu.: 2.231578  
# Max.   :18.572697  

d1 <- as_tibble(t(combn(x = seq(0.2, 0.80, 0.2), m = 2)))
colnames(d1) <- c("pfpr1", "pfpr2")

d2 <- tibble(pfpr1 = seq(0.2, 0.80, 0.2),
             pfpr2 = seq(0.2, 0.80, 0.2))

# read in parameter draws
rsample <- readRDS(paste0(HPCpath, "draws.rds"))

# combine
d3 <- crossing(bind_rows(d1, d2), rsample)

strategy <- c("baseline", "test&treat")

mixing_i <- c(0.0066, 0.0127, 0.0223)

scenarios <- crossing(d3, mixing_i, strategy)


# Run tasks --------------------------------------------------------------------
index <- scenarios

# remove ones that have already been run
index <- index |> 
  mutate(f = paste0("./mixing_casestudy/", pfpr1, "_", pfpr2, "_", draw, "_", strategy, "_", mixing_i, ".rds")) |>
  mutate(exist = case_when(file.exists(f) ~ 1, !file.exists(f) ~ 0)) |>
  filter(exist == 0) |>
  select(-f, -exist)

# submit jobs
id <- task_create_expr(run_metapopulation_mixing_casestudy(0.2, 0.2, 4, 0.0066, "baseline"))
task_status(id)
task_log_show(id)

# in bulk
bundle <- task_create_bulk_expr(run_metapopulation_mixing_casestudy(pfpr1, pfpr2, draw,  mixing_i, strategy), index[2001:3000,])
hipercow_bundle_status(bundle)



# Processing -------------------------------------------------------------------
# read in results
files <- list.files(path = paste0(HPCpath, "mixing_casestudy"), 
                      pattern = "*.rds", full.names = TRUE)
dat_list <- lapply(files, function (x) readRDS(x))

# concatenate
results <- do.call("bind_rows", dat_list) |> as_tibble()

results2 <- results |>
  filter(year >= 20 & year <= 29) |> # intervention starts at 20 years
  # summarize over years
  group_by(ID, strategy, mixing_i, pfpr, time, draw) |>
  summarize(A_count = mean(A_count, na.rm = T),
            D_count = mean(D_count, na.rm = T),
            U_count = mean(U_count, na.rm = T),
            
            n_0_73000 = mean(n_0_73000, na.rm = T),
            n_730_3650 = mean(n_730_3650, na.rm = T),
            
            p_inc_clinical_0_73000 = sum(p_inc_clinical_0_73000, na.rm = T),
            p_detect_730_3650 = mean(p_detect_730_3650, na.rm = T)) |>
  rowwise() |>
  mutate(pfpr_2_10 = p_detect_730_3650 / n_730_3650,
         case_total = p_inc_clinical_0_73000,
         infectious_prev = sum(D_count, A_count, U_count) / n_0_73000,
         logit_prev = log(infectious_prev / (1 - infectious_prev)),
         logit_rdt = -0.968 + 1.186 * logit_prev,
         detect = exp(logit_rdt) / (1 + exp(logit_rdt))) |>
  # remove unnecessary variables
  dplyr::select(ID, strategy, mixing_i, pfpr, time, draw, pfpr_2_10, case_total, infectious_prev, detect)

# separate out baseline and intervention data
baseline <- results2 |> filter(strategy == "baseline") |>
  ungroup() |>
  rename(case_total_baseline = case_total,
         pfpr_2_10_baseline = pfpr_2_10) |>
  dplyr::select(-strategy, -time, -infectious_prev, -detect)

intervention <- results2 |> 
  ungroup() |> filter(strategy == "test&treat") |>
  dplyr::select(-strategy, -time)

# merge and calculate cases averted
avert <- intervention |> left_join(baseline, by = c("ID", "pfpr", "draw", "mixing_i")) |>
  rowwise() |>
  mutate(case_avert = case_total_baseline - case_total,
         case_avert_p = case_avert / case_total_baseline * 100,
         detect = median(detect))

avert_median <- avert |> ungroup() |> group_by(ID, pfpr, mixing_i) |>
  summarize(case_avert_p = median(case_avert_p),
            detect = median(detect),
            infectious_prev = median(infectious_prev))

summary(avert_median$case_avert_p); summary(avert_median$detect); summary(avert_median$infectious_prev)


# save output
saveRDS(avert, paste0(path, "03_output/mixing_casestudy.rds"))
saveRDS(avert, paste0(HPCpath, "mixing_casestudy.rds"))

