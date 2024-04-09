# Process malariasimulation runs.

# libraries
source("./02_code/packages_paths.R")

# read in all parameter draw runs and process
files <- list.files(path = paste0(HPCpath, "results/"), full.names = TRUE)

# turn off summarize warning in function below
options(dplyr.summarize.inform = FALSE)

# cycle through all clusters
process_results <- function(x){ # index of results file

  # read in file
  results <- readRDS(files[x])
  
  # calculate outcomes averted
  results_agg <- results |> ungroup() |>
    # summarize over years
    group_by(ID_0, NAME_1, pop, seed_point, mix_type, strategy, draw) |>
    filter(year >= 23 & year <= 32) |> # just capture intervention period
    summarize(time = mean(time), 
              S_count = mean(S_count, na.rm = T),
              A_count = mean(A_count, na.rm = T),
              D_count = mean(D_count, na.rm = T),
              U_count = mean(U_count, na.rm = T),
              Tr_count = mean(Tr_count, na.rm = T),
              
              n_730_3649 = mean(n_730_3649, na.rm = T),
              p_detect_730_3649 = mean(p_detect_730_3649, na.rm = T),
              
              n_total = mean(n_total, na.rm = T),
              
              p_inc_clinical_total = sum(p_inc_clinical_total, na.rm = T),
              p_inc_severe_total = sum(p_inc_severe_total, na.rm = T)) |>
    # weight cases by total population in each admin
    rowwise() |>
    mutate(pfpr_2_10 = p_detect_730_3649 / n_730_3649,
           case_total = (p_inc_clinical_total / n_total) * pop,
           severe_total = (p_inc_severe_total / n_total) * pop,
           infectious_prev = sum(D_count, A_count, U_count) / n_total,
           logit_prev = log(infectious_prev / (1 - infectious_prev)),
           logit_rdt = -0.968 + 1.186 * logit_prev,
           detect = exp(logit_rdt) / (1 + exp(logit_rdt))) |>
    
    # remove unnecessary variables
    dplyr::select(-n_730_3649, -p_detect_730_3649, -p_inc_clinical_total, 
           -p_inc_severe_total)
  
  print(x)
  return(results_agg)
  
}

results_all <- map_dfr(.x = c(1:length(files)), .f = process_results)

saveRDS(results_all, paste0(HPCpath, "output_aggregated.rds"))

# separate out baseline and intervention data
baseline <- results_all |> filter(strategy == "baseline") |>
  rename(case_total_baseline = case_total,
         severe_total_baseline = severe_total,
         pfpr_2_10_baseline = pfpr_2_10,
         infectious_prev_baseline = infectious_prev,
         detect_baseline = detect) |>
  dplyr::select(-strategy, -logit_prev, -logit_rdt)

intervention <- results_all |> filter(strategy == "test&treat") |>
  dplyr::select(-strategy, -logit_prev, -logit_rdt)

# merge and calculate cases averted
avert <- intervention |> left_join(baseline, by = c("ID_0", "NAME_1", "pop", "seed_point", "mix_type", "draw")) |>
  rowwise() |>
  mutate(case_avert = case_total_baseline - case_total,
         case_avert_p = case_avert / case_total_baseline * 100,
         severe_avert = severe_total_baseline - severe_total)

avert_median <- avert |> ungroup() |> group_by(ID_0, NAME_1, seed_point) |>
  summarize(case_avert = median(case_avert),
            case_avert_p = median(case_avert_p),
            detect = median(detect),
            infectious_prev = median(infectious_prev))

summary(avert_median$case_avert); summary(avert_median$case_avert_p); summary(avert_median$detect); summary(avert_median$infectious_prev)


# save output
saveRDS(avert, paste0(HPCpath, "avert_output.rds"))
saveRDS(avert, paste0(path, "03_output/avert_output.rds"))


# look at distribution of run times
run_times <- avert |> ungroup() |> 
  dplyr::select(ID_0, NAME_1, seed_point, draw, time.x, time.y) |>
  pivot_longer(cols = c(time.x, time.y), names_to = "xy", values_to = "time") |>
  group_by(seed_point, draw, xy) |>
  summarize(time = mean(time)) |> # remove duplicates - each point in cluster has the same run time
  rowwise() |>
  mutate(time = as.numeric(mean(time)))
  
summary(as.numeric(run_times$time))

ggplot(data = run_times) + 
  geom_histogram(aes(x = time), binwidth = .3, color = "#376795", 
                 fill = "#528fad") +
  geom_density(aes(x = time, y = .3 * after_stat(count)), linewidth = 1.3, color = "#e76254") + 
  labs(x = "time (hours)", y = "number of runs") + 
  theme_classic()

ggsave(filename = "./03_output/run_times.pdf", width = 5, height = 3)


