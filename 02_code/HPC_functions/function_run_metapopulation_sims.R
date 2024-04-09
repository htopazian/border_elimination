# This script contains the function necessary to run the metapopulation model

run_metapopulation_sim <- function(cluster,   # cluster number
                                   mix_type,  # mixing matrix type
                                   draw,      # parameter draw
                                   strategy){ # intervention
  
  library(malariasimulation)
  library(dplyr)
  library(tidyr)
  library(sf)
  
  start.time <- Sys.time()
  
  # read in .rds files - cluster set and mixing matrix of interest
  cluster_set <- readRDS(paste0("./clusters/cluster_", cluster, ".rds")) |> st_drop_geometry()
  mixing_matrix <- readRDS(paste0("./mixing_matrices_grid/mixing_matrix_", cluster, "_", mix_type, ".rds"))
  
  # read in site parameters
  site_params <- readRDS("./site_params.rds")
  
  # read in population
  # pull in 2023 population estimates from B_population.R
  population <- readRDS("./population.rds") 
  
  # extract appropriate parameter sets based on cluster set values
  cluster_set <- cluster_set |> left_join(site_params) |>
    left_join(population, by = c("ID_0" = "iso3c", "NAME_1" = "name_1"))
  
  # generic values
  year <- 365
  month <- year / 12
  intervention_start <- 23 * 365 # intervention start date = year 2023
  coef <- 0.80 # test-and-treat coverage
  count <- nrow(cluster_set)
  
  # define timesteps and parameter sets
  timesteps <- cluster_set$params[[1]]$timesteps
  paramslist <- lapply(1:count, 
                       function(x){cluster_set$params[[x]]})
  
  # modify parameters per parameter draw
  add_draw <- function(x){
    p <- paramslist[[x]]

    p <- set_parameter_draw(
      parameters = p,
      draw = draw)
    
    paramslist[[x]] <- p
  }
  
  paramslist <- lapply(c(1:length(paramslist)), add_draw)
  
  # baseline, national cooperation, regional cooperation, synchronization
  if(strategy == "baseline"){
    import_mixing <- list(mixing_matrix)
    export_mixing <- list(t(mixing_matrix))
    mixing_tt <- c(1)
    
    p_captured <- list(
      matrix(data = rep(0, count * count),
             nrow = count,
             ncol = count,
             byrow = T))
    
    p_captured_tt <- c(1)
  }
  
  # border closure
  if(strategy == "closure"){
    mixing_matrix2 <- mixing_matrix
    
    # keep people traveling outside of country A at origin
    for(col in 1:4){
      mixing_matrix2[col,col] <- sum(mixing_matrix2[col, c(col,5:8)])
      mixing_matrix2[col,5:8] <- 0
    }
    
    # keep people traveling outside of country B at origin
    for(col in 5:8){
      mixing_matrix2[col,col] <- sum(mixing_matrix2[col, c(col,1:4)])
      mixing_matrix2[col,1:4] <- 0
    }
    
    import_mixing <- list(mixing_matrix, mixing_matrix2)
    export_mixing <- list(t(mixing_matrix), t(mixing_matrix2))
    mixing_tt <- c(1, intervention_start)
    
    p_captured <- list(
      matrix(data = rep(0, count * count),
             nrow = count,
             ncol = count,
             byrow = T),
      matrix(data = rep(0, count * count),
             nrow = count,
             ncol = count,
             byrow = T))
    p_captured_tt <- c(1, intervention_start)
  }
  
  
  # test-and-treat
  if(strategy == "test&treat"){
    import_mixing <- list(mixing_matrix, mixing_matrix)
    export_mixing <- list(t(mixing_matrix), t(mixing_matrix))
    
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
  output <- run_metapop_simulation(timesteps = timesteps,
                                   parameters = paramslist,
                                   correlations = NULL,
                                   mixing_tt = mixing_tt, 
                                   export_mixing = export_mixing, # rows = destinations, cols = origins
                                   import_mixing = import_mixing, # rows = origins, cols = destinations
                                   p_captured = p_captured,
                                   p_captured_tt = p_captured_tt,
                                   p_success = 0.95) # AL efficacy: https://github.com/mrc-ide/malariasimulation/blob/master/man/AL_params.Rd#L19
  
  results <- do.call("bind_rows", output)
  results$ID_0 <- unlist(lapply(cluster_set$ID_0, function(x){rep(x, timesteps)}))
  results$NAME_1 <- unlist(lapply(cluster_set$NAME_1, function(x){rep(x, timesteps)}))
  results$pop <- unlist(lapply(cluster_set$pop, function(x){rep(x, timesteps)}))
  results$seed_point <- cluster
  results$mix_type <- mix_type
  results$draw <- draw
  results$strategy <- strategy

  results <- results |>
    mutate(year = ceiling(timestep/year),
           month = ceiling(timestep/month)) |>
    group_by(ID_0, NAME_1, year, month, pop, seed_point, mix_type, draw, strategy) |>
    summarize(S_count = mean(S_count, na.rm = T),
              A_count = mean(A_count, na.rm = T),
              D_count = mean(D_count, na.rm = T),
              U_count = mean(U_count, na.rm = T),
              Tr_count = mean(Tr_count, na.rm = T),

              n_0_1824 = mean(n_0_1824, na.rm = T),
              n_1825_5474 = mean(n_1825_5474, na.rm = T),
              n_5475_36499 = mean(n_5475_36499, na.rm = T),
              
              n_730_3649 = mean(n_730_3649, na.rm = T),
              p_detect_730_3649 = mean(p_detect_730_3649, na.rm = T),
              
              p_inc_clinical_0_1824 = sum(p_inc_clinical_0_1824, na.rm = T),
              p_inc_clinical_1825_5474 = sum(p_inc_clinical_1825_5474, na.rm = T),
              p_inc_clinical_5475_36499 = sum(p_inc_clinical_5475_36499, na.rm = T),
              
              p_inc_severe_0_1824 = sum(p_inc_severe_0_1824, na.rm = T),
              p_inc_severe_1825_5474 = sum(p_inc_severe_1825_5474, na.rm = T),
              p_inc_severe_5475_36499 = sum(p_inc_severe_5475_36499, na.rm = T),
              
              n_treated = sum(n_treated, na.rm = T),
              EIR_arabiensis = mean(EIR_arabiensis, na.rm = T),
              EIR_funestus = mean(EIR_funestus, na.rm = T),
              EIR_gambiae = mean(EIR_gambiae, na.rm = T),
              ica_mean = mean(ica_mean, na.rm = T),
              ib_mean = mean(ib_mean, na.rm = T)) |>
    
    rowwise() |>
    mutate(n_total = sum(n_0_1824, n_1825_5474, n_5475_36499, na.rm = T),
           p_inc_clinical_total = sum(p_inc_clinical_0_1824, p_inc_clinical_1825_5474, p_inc_clinical_5475_36499, na.rm = T),
           p_inc_severe_total = sum(p_inc_severe_0_1824, p_inc_severe_1825_5474, p_inc_severe_5475_36499, na.rm = T),
           detect = sum(A_count, D_count, U_count))
  
  # add time elapsed
  end.time <- Sys.time()
  results$time <- end.time - start.time
  
  # save
  saveRDS(results, paste0("./results/", cluster, "_", mix_type, "_", draw, "_", strategy, ".rds"))

}

