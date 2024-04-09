# This script contains the function necessary to run the metapopulation model mixing example

# ------------------------------------------------------------------------------


mixing_example <- function(x){ # x = model name
  
  # set variables
  year <- 365
  sim_length <- 40 * year
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
  if (x == "isolated"){
  output <- run_metapop_simulation(timesteps = sim_length,
                                   parameters = paramslist,
                                   correlations = NULL,
                                   mixing_tt = 1,
                                   export_mixing = list(diag(length(EIR_vector))), # FOIM, transposed EIR
                                   import_mixing = list(diag(length(EIR_vector))), # EIR
                                   p_captured = list(matrix(0,3,3)),
                                   p_captured_tt = 1,
                                   p_success = 0)
  
  output <- data.table::rbindlist(output, idcol = "EIR") |>
    mutate(EIR = c(sort(rep(EIR_vector, sim_length))))
  }
  
  # perfectly mixed
  else if (x == "perfectly mixed"){
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
  
  output <- data.table::rbindlist(output, idcol = "EIR") |>
    mutate(EIR = c(sort(rep(EIR_vector, sim_length))))
  }
  

  # semi-mixed
  else if (x == "semi-mixed"){
  output <- run_metapop_simulation(timesteps = sim_length,
                                   parameters = paramslist,
                                   correlations = NULL,
                                   mixing_tt = c(1, 15 * year),
                                   export_mixing = list(diag(length(EIR_vector)),
                                                        t(matrix(c(0.8, 0.1, 0.1,
                                                                   0.05, 0.85, 0.1,
                                                                   0.05, 0.2, 0.75),
                                                                 nrow = 3, ncol = 3, 
                                                                 byrow = TRUE))), # FOIM, transposed EIR
                                   import_mixing = list(diag(length(EIR_vector)),
                                                        matrix(c(0.8, 0.1, 0.1,
                                                                 0.05, 0.85, 0.1,
                                                                 0.05, 0.2, 0.75),
                                                               nrow = 3, ncol = 3, 
                                                               byrow = TRUE)), # EIR
                                   p_captured = list(matrix(0,3,3)),
                                   p_captured_tt = 1,
                                   p_success = 0)
  
  output <- data.table::rbindlist(output, idcol = "EIR") |>
    mutate(EIR = c(sort(rep(EIR_vector, sim_length))))
  }
  
  # merge dataframes
  mdat <- output |>
    mutate(model = x) |>
    mutate(prev2to10 = p_detect_730_3650 / n_730_3650) |>
    mutate(year = ceiling(timestep / 365),
           month = ceiling(timestep / (365 / 12))) 
  
  # save
  saveRDS(mdat, paste0("./3unit_mixing_example_", x, ".rds"))
  
}

