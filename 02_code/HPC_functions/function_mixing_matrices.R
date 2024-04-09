# This script contains the function necessary to generate mixing matrices
# here we are using the gravity model from Marshall et al. 2018
# https://www.nature.com/articles/s41598-018-26023-1
# with trip duration data from Marshall et al. 2016
# https://malariajournal.biomedcentral.com/articles/10.1186/s12936-016-1252-3#MOESM2
# 12936_2016_1252_MOESM2_ESM.xlsx Additional file 2
# MAP friction surface from (with refs): https://malariaatlas.org/project-resources/accessibility-to-healthcare/

# ------------------------------------------------------------------------------


create_mix_matrix <- function(x, y){ # x = index of seed point; y = distance (1) or travel time (2) used to calculate mixing
  
  library(malariasimulation)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(sf)
  library(gdistance)
  library(raster)
  library(reshape2)
  library(exactextractr)
  library(units)
  library(quantreg)
  
  # read in cluster
  random_sample <- readRDS(paste0("clusters/cluster_", x, ".rds")) |>
    mutate(match = paste0(ID_0, "_", NAME_1))
  
  # pull in probability of staying values from MISC_travel_prob.R
  admin1s_rematch <- readRDS("admin1s_travel.rds")
  
  # pull in gravity model fits from MISC_trip_duration_travel_time_fits.R
  model_duration <- readRDS("dist_duration_quant_model.rds")
  model_tt <- readRDS("traveltime_duration_quant_model.rds")
  
  # load friction surface from MAP
  # tutorial from Amelia B-V: https://medium.com/@abertozz/mapping-travel-times-with-malariaatlas-and-friction-surfaces-f4960f584f08
  # friction <- malariaAtlas::getRaster(
  #   surface = "A global friction surface enumerating land-based travel speed for a nominal year 2015",
  #   file_path = "./01_data/MAP_friction_surface/",
  #   shp = sf::as_Spatial(countries))
  
  # malariaAtlas::autoplot_MAPraster(friction) # visualize
  
  # load raster from process ^
  raster_layer_tt <- raster::raster("01_data/MAP_friction_surface/getRaster/2015_friction_surface_v1_Decompressed_latest_.17.5_.34.8_51.41_27.29_2023_03_07.tiff")
  
  # import worldpop pop count data: unconstrained global mosaic 2020 (1km resolution)
  # https://hub.worldpop.org/project/categories?id=3
  raster_layer_pop <- raster::raster("01_data/ppp_2020_1km_Aggregated.tif")

  # crop the two rasterlayers to the cluster polygons
  # the buffer is needed because the shortest travel time might involve a route not included in the bbox of the admin units
  admin_sample <- admin1s_rematch |> 
    mutate(match = paste0(ID_0, "_", NAME_1)) |>
    filter(match %in% random_sample$match)

  bbox <- bbox(sf::as_Spatial(admin_sample))
  bbox_sp <- as(extent(bbox), "SpatialPolygons")
  buffer_size <- 1 # set the buffer size in arc_degree
  bbox_buffered <- st_buffer(st_as_sfc(bbox_sp), dist = buffer_size)

  raster_layer_pop_cropped <- crop(raster_layer_pop, st_as_sf(bbox_buffered))
  
  if (y == 2) {
    raster_layer_tt_cropped <- crop(raster_layer_tt, st_as_sf(bbox_buffered))
  }
  
  # make grids: cell size = .1 degree ~ 11 km 
  # st_crs(admin_sample)$proj4string
  
  # if x is in this set of numbers, the function results in an error because of degenerate edges
  if (x %in% c(13,38,39,70,76,109,131,142,151,178,186,187,189,190,195,211,240,244,274,278,282,314,357,374)) {
  grid <- st_make_grid(x = bbox_sp, cellsize = .1, what = "polygons", square = TRUE) |>
    st_as_sf(crs = st_crs(4326)) |> 
    st_join(admin_sample[,1:4], largest = FALSE) # assign to admin1s
  }
  
  # otherwise default to assigning admin1s to the largest cell
  if (!(x %in% c(13,38,39,70,76,109,131,142,151,178,186,187,189,190,195,211,240,244,274,278,282,314,357,374))) {
    grid <- st_make_grid(x = bbox_sp, cellsize = .1, what = "polygons", square = TRUE) |>
      st_as_sf(crs = st_crs(4326)) |> 
      st_join(admin_sample[,1:4], largest = TRUE) # assign to admin1s
  }
  
  # manually set some grid cells for admins which are too small to contain a grid centroid
  if (sum(str_detect(random_sample$ID_0, "BWA")) > 0 & sum(str_detect(random_sample$NAME_1, "Sowa")) > 0) {
    grid <- grid |>
      st_join(admin_sample[,1:4] |> filter(ID_0 == "BWA" & NAME_1 == "Sowa")) |> # assign to admin1s
      mutate(ID_0 = case_when(!is.na(ID_0.y) ~ ID_0.y,
                              is.na(ID_0.y) ~ ID_0.x),
             COUNTRY = case_when(!is.na(COUNTRY.y) ~ COUNTRY.y,
                                 is.na(COUNTRY.y) ~ COUNTRY.x),
             ID_1 = case_when(!is.na(ID_1.y) ~ ID_1.y,
                              is.na(ID_1.y) ~ ID_1.x),
             NAME_1 = case_when(!is.na(NAME_1.y) ~ NAME_1.y,
                                is.na(NAME_1.y) ~ NAME_1.x)) |>
      dplyr::select(ID_0, COUNTRY, ID_1, NAME_1, geometry)
  }
  
  if (sum(str_detect(random_sample$ID_0, "COG")) > 0 & sum(str_detect(random_sample$NAME_1, "Brazzaville")) > 0) {
    grid <- grid |>
      st_join(admin_sample[,1:4] |> filter(ID_0 == "COG" & NAME_1 == "Brazzaville")) |> # assign to admin1s
      mutate(ID_0 = case_when(!is.na(ID_0.y) ~ ID_0.y,
                              is.na(ID_0.y) ~ ID_0.x),
             COUNTRY = case_when(!is.na(COUNTRY.y) ~ COUNTRY.y,
                                 is.na(COUNTRY.y) ~ COUNTRY.x),
             ID_1 = case_when(!is.na(ID_1.y) ~ ID_1.y,
                              is.na(ID_1.y) ~ ID_1.x),
             NAME_1 = case_when(!is.na(NAME_1.y) ~ NAME_1.y,
                                is.na(NAME_1.y) ~ NAME_1.x)) |>
      dplyr::select(ID_0, COUNTRY, ID_1, NAME_1, geometry)
  }
  
  if (sum(str_detect(random_sample$ID_0, "MWI")) > 0 & sum(str_detect(random_sample$NAME_1, "Likoma")) > 0) {
    grid <- grid |>
      st_join(admin_sample[,1:4] |> filter(ID_0 == "MWI" & NAME_1 == "Likoma")) |> # assign to admin1s
      mutate(ID_0 = case_when(!is.na(ID_0.y) ~ ID_0.y,
                              is.na(ID_0.y) ~ ID_0.x),
             COUNTRY = case_when(!is.na(COUNTRY.y) ~ COUNTRY.y,
                                 is.na(COUNTRY.y) ~ COUNTRY.x),
             ID_1 = case_when(!is.na(ID_1.y) ~ ID_1.y,
                              is.na(ID_1.y) ~ ID_1.x),
             NAME_1 = case_when(!is.na(NAME_1.y) ~ NAME_1.y,
                                is.na(NAME_1.y) ~ NAME_1.x)) |>
      dplyr::select(ID_0, COUNTRY, ID_1, NAME_1, geometry)
  }
  
  if (sum(str_detect(random_sample$ID_0, "BWA")) > 0 & sum(str_detect(random_sample$NAME_1, "Francistown")) > 0) {
    grid <- grid |>
      st_join(admin_sample[,1:4] |> filter(ID_0 == "BWA" & NAME_1 == "Francistown")) |> # assign to admin1s
      mutate(ID_0 = case_when(!is.na(ID_0.y) ~ ID_0.y,
                              is.na(ID_0.y) ~ ID_0.x),
             COUNTRY = case_when(!is.na(COUNTRY.y) ~ COUNTRY.y,
                                 is.na(COUNTRY.y) ~ COUNTRY.x),
             ID_1 = case_when(!is.na(ID_1.y) ~ ID_1.y,
                              is.na(ID_1.y) ~ ID_1.x),
             NAME_1 = case_when(!is.na(NAME_1.y) ~ NAME_1.y,
                                is.na(NAME_1.y) ~ NAME_1.x)) |>
      dplyr::select(ID_0, COUNTRY, ID_1, NAME_1, geometry)
  }
  
  if (sum(str_detect(random_sample$ID_0, "BWA")) > 0 & sum(str_detect(random_sample$NAME_1, "Selibe Phikwe")) > 0) {
    grid <- grid |>
      st_join(admin_sample[,1:4] |> filter(ID_0 == "BWA" & NAME_1 == "Selibe Phikwe")) |> # assign to admin1s
      mutate(ID_0 = case_when(!is.na(ID_0.y) ~ ID_0.y,
                              is.na(ID_0.y) ~ ID_0.x),
             COUNTRY = case_when(!is.na(COUNTRY.y) ~ COUNTRY.y,
                                 is.na(COUNTRY.y) ~ COUNTRY.x),
             ID_1 = case_when(!is.na(ID_1.y) ~ ID_1.y,
                              is.na(ID_1.y) ~ ID_1.x),
             NAME_1 = case_when(!is.na(NAME_1.y) ~ NAME_1.y,
                                is.na(NAME_1.y) ~ NAME_1.x)) |>
      dplyr::select(ID_0, COUNTRY, ID_1, NAME_1, geometry)
  }
  
  if (sum(str_detect(random_sample$ID_0, "BWA")) > 0 & sum(str_detect(random_sample$NAME_1, "Jwaneng")) > 0) {
    grid <- grid |>
      st_join(admin_sample[,1:4] |> filter(ID_0 == "BWA" & NAME_1 == "Jwaneng")) |> # assign to admin1s
      mutate(ID_0 = case_when(!is.na(ID_0.y) ~ ID_0.y,
                              is.na(ID_0.y) ~ ID_0.x),
             COUNTRY = case_when(!is.na(COUNTRY.y) ~ COUNTRY.y,
                                 is.na(COUNTRY.y) ~ COUNTRY.x),
             ID_1 = case_when(!is.na(ID_1.y) ~ ID_1.y,
                              is.na(ID_1.y) ~ ID_1.x),
             NAME_1 = case_when(!is.na(NAME_1.y) ~ NAME_1.y,
                                is.na(NAME_1.y) ~ NAME_1.x)) |>
      dplyr::select(ID_0, COUNTRY, ID_1, NAME_1, geometry)
  }  
  
  if (sum(str_detect(random_sample$ID_0, "BWA")) > 0 & sum(str_detect(random_sample$NAME_1, "Lobatse")) > 0) {
    grid <- grid |>
      st_join(admin_sample[,1:4] |> filter(ID_0 == "BWA" & NAME_1 == "Lobatse")) |> # assign to admin1s
      mutate(ID_0 = case_when(!is.na(ID_0.y) ~ ID_0.y,
                              is.na(ID_0.y) ~ ID_0.x),
             COUNTRY = case_when(!is.na(COUNTRY.y) ~ COUNTRY.y,
                                 is.na(COUNTRY.y) ~ COUNTRY.x),
             ID_1 = case_when(!is.na(ID_1.y) ~ ID_1.y,
                              is.na(ID_1.y) ~ ID_1.x),
             NAME_1 = case_when(!is.na(NAME_1.y) ~ NAME_1.y,
                                is.na(NAME_1.y) ~ NAME_1.x)) |>
      dplyr::select(ID_0, COUNTRY, ID_1, NAME_1, geometry)
  }
  
  # identify grid cell centroids and assign to admin1s
  grid_centroids <- st_centroid(grid) |> 
    st_join(admin_sample[,1:4])
  
  # extract population size
  pop <- exactextractr::exact_extract(raster_layer_pop_cropped, grid, fun = "sum", progress = FALSE)
  
  # create distance matrix for either euclidean distance or travel time
  if (y == 1) {
    
    # distance matrix - in km
    d <- set_units(st_distance(x = grid_centroids,
                               y = grid_centroids), km) # in km
    
    # remove units
    clean_units <- function(x){
      attr(x, "units") <- NULL
      class(x) <- setdiff(class(x), "units")
      x
    }
    
    d <- clean_units(d)
    
  }
  
  if (y == 2) {
    
    # calculate travel time from cell to cell
    # make and geocorrect the transition matrix 
    T <- gdistance::transition(raster_layer_tt_cropped, function(x) 1/mean(x), 8) # RAM intensive, can be very slow for large areas
    T.GC <- gdistance::geoCorrection(T)                    
    
    # run accumulated cost algorithm to calculate travel times from points to all other points
    d <- gdistance::costDistance(T.GC, as_Spatial(grid_centroids), as_Spatial(grid_centroids))
    
  }
  
  # Gravity model parameters
  if (y == 1) {
    # define parameters [Marshall et al. 2018]
    a <- 1.91
    log_p <- 4.29
    t <- 1.22
    
    # select duration model
    model <- model_duration
  }
  
  if (y == 2) {
    # define parameters [code from Marshall et al. 2018 run for travel times]
    a <- 2.67
    log_p <- 5.11
    t <- 1.13
    
    # select duration model
    model <- model_tt
  }
  
  
  # destination population matrix
  Nj <- t(matrix(data = c(rep(pop, nrow(grid_centroids))),
                 nrow = nrow(grid_centroids),
                 ncol = nrow(grid_centroids)))
  
  # calculate the probability of trips to each destination
  totals <- (Nj ^ t) * ((1 + (d / exp(log_p)))^(-a))
  # set diagonals to 0 (staying within own admin)
  
  # label row and columns with admin1 IDs
  rownames(totals) <- grid$ID_1
  colnames(totals) <- grid$ID_1
  
  
  # collapse grid to summarize by admin1 location and destination
  # empty matrix
  m <- matrix(data = NA, 
              nrow = nrow(admin_sample),
              ncol = nrow(admin_sample))
  
  # melt matrix into dataframe of all combinations
  df <- reshape2::melt(totals) |> filter(!is.na(Var1) & !is.na(Var2))
  df <- df |> group_by(Var1, Var2) |> summarize(value = sum(value, na.rm = TRUE)) # summarize values
  
  # find the row number of each admin in the 8 group sample
  # then assign row numbers on to the gravity model dataframe
  admin_row <- admin_sample |> 
    st_drop_geometry() |> 
    dplyr::select(ID_1) |> mutate(rn = row_number()) 
  
  df <- df |> left_join(admin_row, by = c("Var1" = "ID_1")) |> 
    left_join(admin_row, by = c("Var2" = "ID_1")) |>
    rename(row = rn.x, col = rn.y) |>
    arrange(row, col) # arrange by row / column 
  
  # fill empty matrix m with corresponding values from gravity model
  for(z in 1:nrow(df)){
    m[df[z,]$row, df[z,]$col] <- df[z,]$value
  }
  
  # set diagonals to zero
  diag(m) <- 0
  
  
  # calculation duration of travel (days)
  
  if (y == 1) {
    # add in duration of travel (length of stay)
    f <- as_tibble(d) |> # convert distance matrix to long tibble
      pivot_longer(cols = everything(), values_to = "Distance", names_to = NULL)
    
    pred <- predict(model, newdata = f) # predict duration of travel with Marshall et al. model fits
    
    pred_dat <- f |> cbind(pred) |> # bind predictions to data
      mutate(pred = ifelse(Distance == 0, NA, pred), # remove pred for admin1 with itself
             pred = ifelse(pred < 0, 0, pred))        # turn negative durations to 0
    
    pred_mat <- matrix(pred_dat$pred, ncol = nrow(totals), nrow = nrow(totals)) # transform predictions back into matrix
  }
  
  if (y == 2) {
    # add in duration of travel (length of stay)
    f <- as_tibble(d) |> # convert distance matrix to long tibble
      pivot_longer(cols = everything(), values_to = "travel_time", names_to = NULL)
    
    pred <- predict(model, newdata = f) # predict duration of travel with Marshall et al. model fits
    
    pred_dat <- f |> cbind(pred) |> # bind predictions to data
      mutate(pred = ifelse(travel_time == 0, NA, pred), # remove pred for admin1 with itself
             pred = ifelse(pred < 0, 0, pred))        # turn negative durations to 0
    
    pred_mat <- matrix(pred_dat$pred, ncol = nrow(totals), nrow = nrow(totals)) # transform predictions back into matrix
  }
  
  # label row and columns with admin1 IDs
  rownames(pred_mat) <- grid$ID_1
  colnames(pred_mat) <- grid$ID_1
  
  # collapse grid to summarize by admin1 location and destination
  # empty matrix
  m2 <- matrix(data = NA, 
               nrow = nrow(admin_sample),
               ncol = nrow(admin_sample))
  
  # melt matrix into dataframe of all combinations
  df <- reshape2::melt(pred_mat) |> filter(!is.na(Var1) & !is.na(Var2))
  df <- df |> group_by(Var1, Var2) |> summarize(value = mean(value, na.rm = TRUE)) # summarize values
  
  # find the row number of each admin in the 8 group sample
  # then assign row numbers on to the gravity model dataframe
  admin_row <- admin_sample |> 
    st_drop_geometry() |> 
    dplyr::select(ID_1) |> mutate(rn = row_number()) 
  
  df <- df |> left_join(admin_row, by = c("Var1" = "ID_1")) |> 
    left_join(admin_row, by = c("Var2" = "ID_1")) |>
    rename(row = rn.x, col = rn.y) |>
    arrange(row, col) # arrange by row / column 
  
  # fill empty matrix m with corresponding values from gravity model
  for(z in 1:nrow(df)){
    m2[df[z,]$row, df[z,]$col] <- df[z,]$value
  }
  
  # set diagonal to zero
  diag(m2) <- 0
  
  # add p_stay travel info onto random_sample
  p_stay_sample <- left_join(random_sample, admin1s_rematch |> st_drop_geometry(), by = c("ID_0", "NAME_1"))
  
  p_travelled <- (1 - p_stay_sample$p_stay) 
  
  p_travelled <- matrix(p_travelled , ncol = nrow(random_sample), nrow = nrow(random_sample)) # transform predictions back into matrix
  
  # get number of trips
  med_trips <- matrix(p_stay_sample$med_trips , ncol = nrow(random_sample), nrow = nrow(random_sample)) # transform predictions back into matrix
  med_trips
 
  # CALCULATE overall P(travel) from i to j
  # align with notation in supplement
  qi <- p_travelled 
  mi <- med_trips
  ui <- m2
  ai <- qi * mi * ui

  pji <- m / rowSums(m, na.rm = TRUE)
  
  Mij <- ai * pji
  diag(Mij) <- 365 - rowSums(ai* pji) 

  # rowSums(Mij) = 1
  Mij <- Mij / 365 # normalize
 
  
  # save
  saveRDS(Mij, paste0("mixing_matrices_grid/mixing_matrix_", x, "_", y, ".rds"))
  
}
