# Results section

# libraries
source("./02_code/packages_paths.R")

# site descriptions ------------------------------------------------------------
# pull in shapefiles from A_shapefiles_centroids.R
centroids_all <- readRDS("./03_output/centroids_all.rds")
centroids_border <- readRDS("./03_output/centroids_border.rds")

# pull in population estimates from B_population.R
population <- readRDS("./03_output/population.rds") 

# load raster from process described in D_mixing matrices
raster_layer <- raster("./01_data/MAP_friction_surface/getRaster/2015_friction_surface_v1_Decompressed_latest_.17.5_.34.8_51.41_27.29_2023_03_07.tiff")

# sites
nrow(centroids_all)

# border sites
nrow(centroids_border)

# countries
unique(centroids_all$ID_0) |> length()
unique(centroids_border$ID_0) |> length()

# travel time
calc_tt <- function(x){ # x = index of seed point
  
  # read in cluster
  random_sample <- readRDS(paste0(HPCpath, "clusters/cluster_", x, ".rds"))
  
  # crop the friction surface rasterlayer to the random_points
  bbox <- bbox(sf::as_Spatial(random_sample))
  bbox_sp <- as(extent(bbox), "SpatialPolygons")

  # draw a buffer zone around the bounding box to expand it
  buffer_size <- 1 # set the buffer size in arc_degree
  bbox_buffered <- st_buffer(st_as_sfc(bbox_sp), dist = buffer_size)

  # assign crs to bounding box
  st_crs(bbox_buffered) <- st_crs(raster_layer)

  # crop rasterlayer
  raster_layer_cropped <- crop(raster_layer, st_as_sf(bbox_buffered))

  # make and geocorrect the transition matrix 
  T <- gdistance::transition(raster_layer_cropped, function(x) 1/mean(x), 8) # RAM intensive, can be very slow for large areas
  T.GC <- gdistance::geoCorrection(T)                    

  # extract lat / longs of random sample points
  points <- as.matrix(st_coordinates(random_sample)) 
  
  # run accumulated cost algorithm to calculate travel times from points to all other points
  d <- gdistance::costDistance(T.GC, points, points)
  
  colnames(d) = random_sample$NAME_1
  row.names(d) = random_sample$NAME_1
 
  d[lower.tri(d > 0, diag = TRUE)] <- 0
  
  d <- tibble(col = rownames(d)[row(d)],
              row = colnames(d)[col(d)], 
              value = c(d)) |> 
    filter(value != 0)
  
  print(paste0("finished cluster ", x))
  
  return(d)

}
  
# cycle through all clusters
travel_time <- map_dfr(.x = c(1:nrow(centroids_border)), .f = calc_tt)

# save
saveRDS(travel_time, "./03_output/traveltimes.rds")

travel_time <- readRDS("./03_output/traveltimes.rds")

travel_time |> ungroup() |>
  summarize(min = min(value), 
            median = median(value),
            mean = mean(value),
            max = max(value))


# population size
population |> ungroup() |>
  summarize(min = min(pop), 
            median = median(pop),
            mean = mean(pop),
            max = max(pop))

# cluster representation
# read in all cluster files and merge
files <- list.files(path = paste0(HPCpath, "clusters/"), full.names = TRUE)
dat_list <- lapply(files, function (x) readRDS(x))
output <- data.table::rbindlist(dat_list, fill = TRUE, idcol = "identifier")

output |> group_by(ID_0, NAME_1) |>
  summarize(n = n()) |>
  ungroup() |>
  summarize(min = min(n), 
            median = median(n),
            mean = mean(n),
            max = max(n))


# SSA cluster analysis ---------------------------------------------------------
# read in cases averted dataset from G_processing.R
avert <- readRDS("./03_output/avert_output.rds")

# attach seed point indices to each admin unit
seed_point <- centroids_border |> st_drop_geometry() |>
  mutate(seed = row_number())

output <- avert |> ungroup() |>
  left_join(seed_point) |>
  filter(seed_point == seed) |>
  group_by(ID_0, NAME_1, pop, seed_point, draw, mix_type) |>
  summarize(pfpr_2_10 = mean(pfpr_2_10, na.rm = T),
            case_total = sum(case_total, na.rm = T),
            severe_total = sum(severe_total, na.rm = T),
            pfpr_2_10_baseline = mean(pfpr_2_10_baseline, na.rm = T),
            case_total_baseline = sum(case_total_baseline, na.rm = T),
            severe_total_baseline = sum(severe_total_baseline, na.rm = T)) |>
  rowwise() |>
  mutate(case_avert = case_total_baseline - case_total,
         case_avert_p = case_avert / case_total_baseline * 100,
         cinc_baseline = (case_total_baseline / pop),
         cinc_diff = (case_total_baseline / pop) - (case_total / pop),
         pfpr_diff = pfpr_2_10_baseline - pfpr_2_10,
         pfpr_diff_p = (pfpr_2_10_baseline - pfpr_2_10) / pfpr_2_10_baseline * 100) |>
  # find median over draws
  ungroup() |> group_by(ID_0, NAME_1, pop, seed_point, mix_type) |>
  summarize(n = n(),
            case_avert = median(case_avert),
            case_avert_p = median(case_avert_p),
            cinc_baseline = median(cinc_baseline),
            cinc_diff = median(cinc_diff),
            pfpr_diff = median(pfpr_diff),
            pfpr_diff_p = median(pfpr_diff_p))

summary(output$case_avert); summary(output$case_avert_p); summary(output$cinc_diff); summary(output$pfpr_diff_p)

# calculate nnt
# pulling out seed points
seed_points <- centroids_border |> st_drop_geometry() |> 
  rename(ID_0_seed = ID_0,
         NAME_1_seed = NAME_1) |>
  mutate(seed_point = row_number())

# linking with results and calculating nnt
output <- avert |> ungroup() |>
  mutate(n_inf_detect = infectious_prev * detect, 
         # the probability of being infectious * the probability of detection through RDT
         # = the overall probability of testing + through RDT when screened
         # the number of people needing to be tested to result in 1 positive RDT
         nnt_detect = 1 / n_inf_detect,
         # the number of people needing to be tested to successfully treat 1 positive RDT
         nnt = 1 / (n_inf_detect * .95)) |>
  dplyr::select(ID_0, NAME_1, seed_point, draw, nnt_detect, n_inf_detect, nnt, pfpr_2_10, pfpr_2_10_baseline)

summary(output$nnt)

# calculate pfpr of domestic vs international neighbors in clusters
output <- output |> left_join(seed_points, by = c("seed_point")) |>
  mutate(domestic = case_when(ID_0 == ID_0_seed ~ "domestic",
                              ID_0 != ID_0_seed ~ "international"),
         pfpr_seed = case_when(ID_0 == ID_0_seed & NAME_1 == NAME_1_seed ~ pfpr_2_10,
                               ID_0 != ID_0_seed ~ pfpr_2_10, 
                               TRUE ~ NA_real_),
         nnt_seed = case_when(ID_0 == ID_0_seed & NAME_1 == NAME_1_seed ~ nnt,
                              TRUE ~ NA_real_)) |>
  # summarize across draws
  group_by(seed_point, domestic) |>
  summarize(pfpr_seed_2 = mean(pfpr_seed, na.rm = TRUE),
            nnt_seed_2 = mean(nnt_seed, na.rm = TRUE)) |>
  pivot_wider(names_from = "domestic", values_from = pfpr_seed_2) |>
  group_by(seed_point) |>
  summarize(across(c(nnt_seed_2, domestic, international), \(x) mean(x, na.rm = TRUE)))

summary(output$nnt_seed_2)  
summary(output[output$nnt_seed_2 < Inf,]$nnt_seed_2)


# PfPR case study --------------------------------------------------------------
output <- readRDS("./03_output/PfPR_casestudy.rds") |> 
  ungroup() |>
  #filter(diagnostic == "RDT") |>
  rowwise() |>
  mutate(case_avert = case_total_baseline - case_total,
         case_avert_p = case_avert / case_total_baseline * 100,
         n_inf_detect = infectious_prev * detect, 
         nnt = 1/ n_inf_detect,
         sens = detect / infectious_prev) |>
  group_by(ID, pfpr, diagnostic) |>
  # summarize over draws
  summarize(ca_25 = quantile(case_avert_p, p = 0.25),
            ca_median = median(case_avert_p),
            ca_75 = quantile(case_avert_p, p = 0.75),
            # sensitivity measures only interpretable for RDT
            sens_25 = quantile(sens, p = 0.25),
            sens_median = median(sens),
            sens_75 = quantile(sens, p = 0.75),
            nnt = median(nnt))

summary(output$ca_median); summary(output$sens_median); summary(output$nnt)

# cases averted outcomes
summary(output$ca_median);

# for pairings of 5 to 15% with 50 to 80%
a <- seq(0.05, 0.10, 0.05)
b <- seq(0.5, 0.8, 0.05)
  
ID_pair <- crossing(a, b) |>
  mutate(ID = paste0(a, "_", b))

output |> 
  filter(pfpr %in% c(0.05, 0.10) & ID %in% ID_pair$ID) |>
  group_by(diagnostic) |>
  summarize(l = min(ca_median),
            m = median(ca_median),
            u = max(ca_median))
  
# sensitivity measures only interpretable for RDT
output |> filter(ID == "0.05_0.05")
output |> filter(ID == "0.8_0.8")


# coverage case study ----------------------------------------------------------
output <- readRDS("./03_output/coverage_casestudy.rds") |> ungroup() |>
  mutate(sens = detect / infectious_prev)

tilep1 <- output |>
  separate(ID, c("pfpr1", "pfpr2"), "_") |>
  filter(pfpr1 == pfpr) |>
  select(pfpr1, pfpr2, coverage, draw, case_avert_p, detect, infectious_prev, sens)

timepall <- tilep1 |> distinct() |> 
  group_by(pfpr1, pfpr2, coverage) |>
  summarize(case_avert_p = median(case_avert_p),
            detect = median(detect),
            infectious_prev = median(infectious_prev),
            sens = median(sens))

test <- timepall |> filter(case_avert_p < 1)
table(test$coverage, test$pfpr2)

