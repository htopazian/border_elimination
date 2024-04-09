# This script contains the function necessary to run the metapopulation model case studies


cali_pfpr <- function(pfpr,    # target pfpr value
                      draw) {  # model parameter draw
  
  # calibration ref: https://mrc-ide.github.io/cali/articles/Basic_calibration.html
  # define target: PfPR2-10 value
  target <- pfpr
  
  # write function for time points to average PfPR - years 4-6
  # 6 years is enough to match by PfPR - 21 years is needed to match by c inc
  year <- 365
  
  # define timesteps and parameter sets
  p <- get_parameters(overrides = list(
    human_population = 100000,
    prevalence_rendering_min_ages = 2 * year,
    prevalence_rendering_max_ages = 10 * year))
  
  p <- set_parameter_draw(
    parameters = p,
    draw = draw)
  
  p$timesteps <- 9 * year # simulation run time = 9 years
  
  summary_mean_pfpr_2_10_6y9y <- function(x){
    
    x$year <- ceiling(x$timestep / year)
    x <- x |> filter(year >= 7)
    prev_2_10 <- mean(x$n_detect_730_3650 / x$n_730_3650)
    return(prev_2_10)
    
  }
  
  # run calibration model
  set.seed(123)
  out <- cali::calibrate(parameters = p,
                         target = target,
                         summary_function = summary_mean_pfpr_2_10_6y9y,
                         tolerance = 0.02,
                         low = 0.001,
                         high = 2000)
  
  # store init_EIR results as an .rds file to be read in later
  PR <- tibble(starting_EIR = out, 
               pfpr = pfpr,
               draw = draw)

  saveRDS(PR, paste0("./calibrate/PRmatch_", pfpr, "_", draw, ".rds"))
  
}

# ------------------------------------------------------------------------------


run_metapopulation_pfpr_casestudy <- function(pfpr1,       # PfPR for unit 1
                                              pfpr2,       # PfPR for unit 2
                                              draw,        # parameter draw
                                              strategy,    # baseline or test&treat
                                              diagnostic){ # RDT or PCR
  
  start.time <- Sys.time()
  
  # define EIRs from pfprs
  match <- readRDS("./EIRestimates.rds")
  
  EIR1 <- as.numeric(match[round(match$pfpr, 2) == round(pfpr1, 2) & match$draw == draw, 1])
  EIR2 <- as.numeric(match[round(match$pfpr, 2) == round(pfpr2, 2) & match$draw == draw, 1])
  
  
  # define mixing matrix
  mixing_matrix <- matrix(data = c(0.95, 0.05, 
                                   0.05, 0.95),
                          nrow = 2,
                          ncol = 2,
                          byrow = T)
  
  # generic values
  year <- 365
  month <- year / 12
  run_time <-  30 * year
  intervention_start <- 20 * year
  coef <- 0.80 # test-and-treat coverage
  count <- nrow(mixing_matrix)
  
  # define timesteps and parameter sets
  params <- get_parameters(overrides = list(
    human_population = 100000,
    prevalence_rendering_min_ages = 2 * year,
    prevalence_rendering_max_ages = 10 * year,
    clinical_incidence_rendering_min_ages = c(0, 0) * year, 
    clinical_incidence_rendering_max_ages = c(5, 200) * year,
    severe_incidence_rendering_min_ages = c(0, 0) * year,
    severe_incidence_rendering_max_ages = c(5, 200) * year))
  
  # modify per parameter draw
  params <- set_parameter_draw(
    parameters = params,
    draw = draw)
  
  # set RDT or PCR diagnostics. RDT = default
  if(diagnostic == "PCR"){
    params$rdt_coeff = 1
    params$rdt_intercept = 0
  }
  
  params1 <- set_equilibrium(params, init_EIR = EIR1)
  params2 <- set_equilibrium(params, init_EIR = EIR2)
  
  paramslist <- list(params1, params2)
  
  # baseline, national cooperation, regional cooperation, synchronization
  if(strategy == "baseline"){
    mixing <- list(mixing_matrix)
    mixing_tt <- c(1)
    
    p_captured <- list(
      matrix(data = rep(0, count * count),
             nrow = count,
             ncol = count,
             byrow = T))
    
    p_captured_tt <- c(1)
  }
  
  # test-and-treat
  if(strategy == "test&treat"){
    mixing <- list(mixing_matrix, mixing_matrix)
    mixing_tt <- c(1, intervention_start)
    
    p_captured <- list(
      matrix(data = rep(0, count * count),
             nrow = count,
             ncol = count,
             byrow = T),
      matrix(data = c(rep(
        c(rep(0, count / 2), 
          rep(coef, count / 2)), 
        count / 2),
        
        rep(
          c(rep(coef, count / 2),
            rep(0, count / 2)), 
          count / 2)),
        
        nrow = count,
        ncol = count,
        byrow = T) 
    )
    p_captured_tt <- c(1, intervention_start)
  }
  
  # run model
  set.seed(123)
  output <- run_metapop_simulation(timesteps = run_time,
                                   parameters = paramslist,
                                   correlations = NULL,
                                   mixing_tt = mixing_tt, 
                                   export_mixing = mixing, # rows = destinations, cols = origins
                                   import_mixing = mixing, # rows = origins, cols = destinations
                                   p_captured = p_captured,
                                   p_captured_tt = p_captured_tt,
                                   p_success = 0.95) # AL efficacy: https://github.com/mrc-ide/malariasimulation/blob/master/man/AL_params.Rd#L19
  
  results <- do.call("bind_rows", output)
  results$ID <- paste0(pfpr1, "_", pfpr2)
  results$pfpr <- unlist(lapply(c(pfpr1, pfpr2), function(x){rep(x, run_time)}))
  results$draw <- draw
  results$strategy <- strategy
  results$diagnostic <- diagnostic
  
  results <- results |>
    mutate(year = ceiling(timestep/year)) |>
    group_by(ID, draw, strategy, diagnostic, pfpr, year) |>
    summarize(S_count = mean(S_count, na.rm = T),
              A_count = mean(A_count, na.rm = T),
              D_count = mean(D_count, na.rm = T),
              U_count = mean(U_count, na.rm = T),
              Tr_count = mean(Tr_count, na.rm = T),
              
              n_0_73000 = mean(n_0_73000, na.rm = T),
              p_inc_clinical_0_73000 = sum(p_inc_clinical_0_73000, na.rm = T),
              
              n_0_1825 = mean(n_0_1825, na.rm = T),
              p_inc_clinical_0_1825 = sum(p_inc_clinical_0_1825, na.rm = T),
              
              n_730_3650 = mean(n_730_3650, na.rm = T),
              p_detect_730_3650 = mean(p_detect_730_3650, na.rm = T),
              
              EIR_All = mean(EIR_All, na.rm = T),
              ica_mean = mean(ica_mean, na.rm = T),
              ib_mean = mean(ib_mean, na.rm = T))
  
  # add time elapsed
  end.time <- Sys.time()
  results$time <- end.time - start.time
  
  # save
  if(diagnostic == "RDT"){
    saveRDS(results, paste0("./PfPR_casestudy/", pfpr1, "_", pfpr2, "_", draw, "_", strategy, ".rds"))
  }
  
  if(diagnostic == "PCR"){
    saveRDS(results, paste0("./PfPR_casestudy2/", pfpr1, "_", pfpr2, "_", draw, "_", strategy, ".rds"))
  }
  
}


# ------------------------------------------------------------------------------


run_metapopulation_coverage_casestudy <- function(pfpr1,     # PfPR for unit 1
                                                  pfpr2,     # PfPR for unit 2
                                                  coverage,  # border post coverage
                                                  draw,      # parameter draw
                                                  strategy){ # baseline or test&treat
  
  start.time <- Sys.time()
  
  # define EIRs from pfprs
  match <- readRDS("./EIRestimates.rds")
  
  EIR1 <- as.numeric(match[round(match$pfpr, 2) == round(pfpr1, 2) & match$draw == draw, 1])
  EIR2 <- as.numeric(match[round(match$pfpr, 2) == round(pfpr2, 2) & match$draw == draw, 1])
  
  
  # define mixing matrix
  mixing_matrix <- matrix(data = c(0.95, 0.05, 
                                   0.05, 0.95),
                          nrow = 2,
                          ncol = 2,
                          byrow = T)
  
  # generic values
  year <- 365
  month <- year / 12
  run_time <-  30 * year
  intervention_start <- 20 * year
  coef <- coverage # test-and-treat coverage
  count <- nrow(mixing_matrix)
  
  # define timesteps and parameter sets
  params <- get_parameters(overrides = list(
    human_population = 100000,
    prevalence_rendering_min_ages = 2 * year,
    prevalence_rendering_max_ages = 10 * year,
    clinical_incidence_rendering_min_ages = c(0, 0) * year, 
    clinical_incidence_rendering_max_ages = c(5, 200) * year,
    severe_incidence_rendering_min_ages = c(0, 0) * year,
    severe_incidence_rendering_max_ages = c(5, 200) * year))
  
  # modify per parameter draw
  params <- set_parameter_draw(
    parameters = params,
    draw = draw)
  
  params1 <- set_equilibrium(params, init_EIR = EIR1)
  params2 <- set_equilibrium(params, init_EIR = EIR2)
  
  paramslist <- list(params1, params2)
  
  # baseline, national cooperation, regional cooperation, synchronization
  if(strategy == "baseline"){
    mixing <- list(mixing_matrix)
    mixing_tt <- c(1)
    
    p_captured <- list(
      matrix(data = rep(0, count * count),
             nrow = count,
             ncol = count,
             byrow = T))
    
    p_captured_tt <- c(1)
  }
  
  # test-and-treat
  if(strategy == "test&treat"){
    mixing <- list(mixing_matrix, mixing_matrix)
    mixing_tt <- c(1, intervention_start)
    
    p_captured <- list(
      matrix(data = rep(0, count * count),
             nrow = count,
             ncol = count,
             byrow = T),
      matrix(data = c(rep(
        c(rep(0, count / 2), 
          rep(coef, count / 2)), 
        count / 2),
        
        rep(
          c(rep(coef, count / 2),
            rep(0, count / 2)), 
          count / 2)),
        
        nrow = count,
        ncol = count,
        byrow = T) 
    )
    p_captured_tt <- c(1, intervention_start)
  }
  
  # run model
  set.seed(123)
  output <- run_metapop_simulation(timesteps = run_time,
                                   parameters = paramslist,
                                   correlations = NULL,
                                   mixing_tt = mixing_tt, 
                                   export_mixing = mixing, # rows = destinations, cols = origins
                                   import_mixing = mixing, # rows = origins, cols = destinations
                                   p_captured = p_captured,
                                   p_captured_tt = p_captured_tt,
                                   p_success = 0.95) # AL efficacy: https://github.com/mrc-ide/malariasimulation/blob/master/man/AL_params.Rd#L19
  
  results <- do.call("bind_rows", output)
  results$ID <- paste0(pfpr1, "_", pfpr2)
  results$pfpr <- unlist(lapply(c(pfpr1, pfpr2), function(x){rep(x, run_time)}))
  results$draw <- draw
  results$strategy <- strategy
  results$coverage <- coverage
  
  results <- results |>
    mutate(year = ceiling(timestep/year)) |>
    group_by(ID, draw, strategy, coverage, pfpr, year) |>
    summarize(S_count = mean(S_count, na.rm = T),
              A_count = mean(A_count, na.rm = T),
              D_count = mean(D_count, na.rm = T),
              U_count = mean(U_count, na.rm = T),
              Tr_count = mean(Tr_count, na.rm = T),
              
              n_0_73000 = mean(n_0_73000, na.rm = T),
              p_inc_clinical_0_73000 = sum(p_inc_clinical_0_73000, na.rm = T),
              
              n_0_1825 = mean(n_0_1825, na.rm = T),
              p_inc_clinical_0_1825 = sum(p_inc_clinical_0_1825, na.rm = T),
              
              n_730_3650 = mean(n_730_3650, na.rm = T),
              p_detect_730_3650 = mean(p_detect_730_3650, na.rm = T),
              
              EIR_All = mean(EIR_All, na.rm = T),
              ica_mean = mean(ica_mean, na.rm = T),
              ib_mean = mean(ib_mean, na.rm = T))
  
  # add time elapsed
  end.time <- Sys.time()
  results$time <- end.time - start.time
  
  # save
  saveRDS(results, paste0("./coverage_casestudy/", pfpr1, "_", pfpr2, "_", draw, "_", strategy, "_", coverage, ".rds"))
  
}

# ------------------------------------------------------------------------------

mixing_example <- function(x){
  
  print(x)
  
  # set variables
  year <- 365
  sim_length <- 20 * year
  EIR_vector <- c(0.001, 1, 20)
  
  # get parameters
  ms_parameterize <- function(x){ # index of EIR and ITN vector
    
    params <- get_parameters(list(human_population = 200000,
                                  model_seasonality = TRUE,
                                  g0 = c(0.285505),
                                  g = c(-0.325352, -0.0109352, 0.0779865),
                                  h = c(-0.132815, 0.104675, -0.013919),
                                  individual_mosquitoes = FALSE))
    
    params <- set_equilibrium(params, init_EIR = EIR_vector[x])
    
    print(paste0("Set params: ", x))
    
    return(params)
    
  }
  
  paramslist <- lapply(seq(1, length(EIR_vector), 1), ms_parameterize)
  
  # run models
  set.seed(123)
  
  # isolated
  output <- run_metapop_simulation(timesteps = sim_length,
                                   parameters = paramslist,
                                   correlations = NULL,
                                   mixing_tt = 1,
                                   export_mixing = list(diag(length(EIR_vector))), # FOIM, transposed EIR
                                   import_mixing = list(diag(length(EIR_vector))), # EIR
                                   p_captured = list(matrix(0,3,3)),
                                   p_captured_tt = 1,
                                   p_success = 0)
  
  isolate <- data.table::rbindlist(output, idcol = "EIR") |>
    mutate(EIR = c(sort(rep(EIR_vector, sim_length))))
  
  # perfectly mixed
  output <- run_metapop_simulation(timesteps = sim_length,
                                   parameters = paramslist,
                                   correlations = NULL,
                                   mixing_tt = c(1, 15 * year),
                                   export_mixing = list(diag(length(EIR_vector)),
                                                        matrix(rep(.33, 9),
                                                               nrow = 3, ncol = 3)), # FOIM, transposed EIR
                                   import_mixing = list(diag(length(EIR_vector)),
                                                        matrix(rep(.33, 9),
                                                               nrow = 3, ncol = 3)), # EIR
                                   p_captured = list(matrix(0,3,3)),
                                   p_captured_tt = 1,
                                   p_success = 0)
  
  mixed <- data.table::rbindlist(output, idcol = "EIR") |>
    mutate(EIR = c(sort(rep(EIR_vector, sim_length))))
  
  # slightly mixed
  output <- run_metapop_simulation(timesteps = sim_length,
                                   parameters = paramslist,
                                   correlations = NULL,
                                   mixing_tt = c(1, 15 * year),
                                   # export_mixing = list(diag(length(EIR_vector)),
                                   #                      t(matrix(c(0.9, 0.05, 0.05,
                                   #                               0.05, 0.9, 0.05,
                                   #                               0.05, 0.05, 0.9),
                                   #                             nrow = 3, ncol = 3))), # FOIM, transposed EIR
                                   # import_mixing = list(diag(length(EIR_vector)),
                                   #                      matrix(c(0.9, 0.05, 0.05,
                                   #                               0.05, 0.9, 0.05,
                                   #                               0.05, 0.05, 0.9),
                                   #                             nrow = 3, ncol = 3)), # EIR
                                   export_mixing = list(diag(length(EIR_vector)),
                                                        t(matrix(c(0.8, 0.1, 0.1,
                                                                   0.05, 0.85, 0.1,
                                                                   0.05, 0.2, 0.75),
                                                                 nrow = 3, ncol = 3))), # FOIM, transposed EIR
                                   import_mixing = list(diag(length(EIR_vector)),
                                                        matrix(c(0.8, 0.1, 0.1,
                                                                 0.05, 0.85, 0.1,
                                                                 0.05, 0.2, 0.75),
                                                               nrow = 3, ncol = 3)), # EIR
                                   p_captured = list(matrix(0,3,3)),
                                   p_captured_tt = 1,
                                   p_success = 0)
  
  middle <- data.table::rbindlist(output, idcol = "EIR") |>
    mutate(EIR = c(sort(rep(EIR_vector, sim_length))))
  
  
  # merge dataframes
  mdat <- isolate |> mutate(model = "isolated") |>
    full_join(mixed |> mutate(model = "perfectly mixed")) |>
    full_join(middle |> mutate(model = "semi-mixed")) |>
    mutate(model = factor(model, levels = c("isolated", "semi-mixed", "perfectly mixed"))) |>
    mutate(prev2to10 = p_detect_730_3650 / n_730_3650) |>
    mutate(year = ceiling(timestep / 365),
           month = ceiling(timestep / (365 / 12))) 
  
  # save
  saveRDS(mdat, paste0("./mixing_example.rds"))
  
  }

