# Running a two unit example, varying PfPRs. 

# Before running this case study, MISC_calibrate.R must be run to generate 
# a set of EIR values for each PfPR and parameter draw set

# libraries
source("./02_code/packages_paths.R")

# HPC set-up -------------------------------------------------------------------
setwd(HPCpath)

# to edit HPC username and password below
# usethis::edit_r_environ()

src <- conan::conan_sources(c("github::mrc-ide/malariasimulation@feat/metapopulation_model"))

ctx <- context::context_save(path = paste0(HPCpath, "contexts"),
                             sources = c(paste0(HPCpath, "function_run_metapopulation_casestudies.R")),
                             packages = c("dplyr", "malariasimulation"),
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
# choose PfPR combinations
d1 <- as_tibble(t(combn(x = seq(0.05, 0.80, 0.05), m = 2)))
colnames(d1) <- c("pfpr1", "pfpr2")

d2 <- tibble(pfpr1 = seq(0.05, 0.80, 0.05),
             pfpr2 = seq(0.05, 0.80, 0.05))

# read in parameter draws
rsample <- readRDS(paste0(HPCpath, "draws.rds"))

# combine
d3 <- crossing(bind_rows(d1, d2), rsample)

strategy <- c("baseline", "test&treat")

diagnostic <- c("RDT", "PCR")

scenarios <- crossing(d3, strategy, diagnostic)


# Run tasks --------------------------------------------------------------------
index <- scenarios

# remove ones that have already been run
index <- index |> filter(diagnostic == "RDT") |> 
  mutate(f = paste0("./PfPR_casestudy/", pfpr1, "_", pfpr2, "_", draw, "_", strategy, ".rds")) |>
  mutate(exist = case_when(file.exists(f) ~ 1, !file.exists(f) ~ 0)) |>
  filter(exist == 0) |>
  select(-f, -exist)

index <- index |> filter(diagnostic == "PCR") |> 
  mutate(f = paste0("./PfPR_casestudy2/", pfpr1, "_", pfpr2, "_", draw, "_", strategy, ".rds")) |>
  mutate(exist = case_when(file.exists(f) ~ 1, !file.exists(f) ~ 0)) |>
  filter(exist == 0) |>
  select(-f, -exist)

# submit jobs, 100 as a time
sjob <- function(x, y){
  
  t <- obj$enqueue_bulk(index[x:y,], run_metapopulation_pfpr_casestudy)
  print(paste0(x, " to ", y))
  
}

map2_dfr(seq(0, nrow(index) - 100, 100),
         seq(99, nrow(index), 100),
         sjob)

# submit all remaining tasks
t <- obj$enqueue_bulk(index, run_metapopulation_pfpr_casestudy)
t$status()


# Processing -------------------------------------------------------------------
# read in results
files_a <- list.files(path = paste0(HPCpath, "PfPR_casestudy"), 
                    pattern = "*.rds", full.names = TRUE)
dat_list_a <- lapply(files_a, function (x) readRDS(x))

files_b <- list.files(path = paste0(HPCpath, "PfPR_casestudy2"), 
                      pattern = "*.rds", full.names = TRUE)
dat_list_b <- lapply(files_b, function (x) readRDS(x))

# concatenate
results_a <- do.call("bind_rows", dat_list_a) |> as_tibble()
results_b <- do.call("bind_rows", dat_list_b) |> as_tibble()

# results_a <- results_a |> mutate(diagnostic = "RDT")

results <- results_a |> full_join(results_b)

results2 <- results |>
  filter(year >= 20 & year <= 29) |> # intervention starts at 20 years
  # summarize over years
  group_by(ID, strategy, diagnostic, pfpr, time, draw) |>
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
  select(ID, strategy, diagnostic, pfpr, time, draw, pfpr_2_10, case_total, infectious_prev, detect)

# separate out baseline and intervention data
baseline <- results2 |> filter(strategy == "baseline") |>
  ungroup() |>
  rename(case_total_baseline = case_total,
         pfpr_2_10_baseline = pfpr_2_10) |>
  select(-strategy, -time, -infectious_prev, -detect)

intervention <- results2 |> 
  ungroup() |> filter(strategy == "test&treat") |>
  select(-strategy, -time)

# merge and calculate cases averted
avert <- intervention |> left_join(baseline, by = c("ID", "pfpr", "draw", "diagnostic")) |>
  rowwise() |>
  mutate(case_avert = case_total_baseline - case_total,
         case_avert_p = case_avert / case_total_baseline * 100,
         detect = median(detect))

avert_median <- avert |> ungroup() |> group_by(ID, pfpr, diagnostic) |>
  summarize(case_avert_p = median(case_avert_p),
            detect = median(detect),
            infectious_prev = median(infectious_prev))

summary(avert_median$case_avert_p); summary(avert_median$detect); summary(avert_median$infectious_prev)


# save output
saveRDS(avert, paste0(path, "03_output/PfPR_casestudy.rds"))
saveRDS(avert, paste0(HPCpath, "PfPR_casestudy.rds"))

