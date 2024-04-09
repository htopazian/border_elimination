# Re-fitting gravity models from Marshall et al. 2016, 2018 with travel time data from MAP.

# libraries
source("./02_code/packages_paths.R")

# pull in supplementary data files from Marshall et al. 2018
# https://doi.org/10.1038/s41598-018-26023-1
datapath <- "./01_data/Marshall et al. 2018/Dataset 1.xlsx"

# origin and destination point coordinates
d1 <- readxl::read_excel(
  path = datapath, sheet = 1)[, c(1, 5:7)] |> 
  mutate(COUNTRY = "Mali")

d2 <- readxl::read_excel(
  path = datapath, sheet = 2)[, c(1, 5:7)] |> 
  mutate(COUNTRY = "Burkina Faso")

d3 <- readxl::read_excel(
  path = datapath, sheet = 3)[, c(1, 6:8)] |> 
  mutate(COUNTRY = "Zambia")

d4 <- readxl::read_excel(
  path = datapath, sheet = 4)[, c(1, 6:8)] |> 
  mutate(COUNTRY = "Tanzania")

coords <- bind_rows(d1, d2, d3, d4)

saveRDS(coords, "./01_data/Marshall et al. 2018/coords.rds")

# trip data
import_trips <- function(x, country){ # x = sheet number, country = country name
  data <- readxl::read_excel(
    path = datapath, sheet = x)[, c(1:3)] |>
    mutate(COUNTRY = country)
  
  names(data) <- c("ORIG", "DEST", "CLUSTER", "COUNTRY")
  
  return(data)
}

trips <- map2_dfr(.x = c(5:8), .y = c("Mali", "Burkina Faso", "Zambia", "Tanzania"), .f = import_trips)


# Distance and travel times ----------------------------------------------------
# merge Marshall et al. 2018 datasets

# coordinates
coords <- coords |>
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = st_crs(4326))

# origins
origins <- left_join(trips, coords, by = c("COUNTRY", "ORIG" = "INDEX")) |>
  dplyr::select(-DEST, -POPULATION) |>
  st_as_sf(crs = st_crs(4326))
  
  
# destinations
destinations <- left_join(trips, coords, by = c("COUNTRY", "DEST" = "INDEX")) |>
  dplyr::select(-ORIG) |>
  st_as_sf(crs = st_crs(4326))

# find Euclidean distances between origins and destinations
distances <- st_distance(x = origins, y = destinations, by_element = TRUE) # meters
distances <- tibble(distances = set_units(distances, km)); summary(distances)

# calculate distances via Marshall et a. 2016 process
orig_coords <- as_tibble(st_coordinates(origins)) # X Y | long lat
dest_coords <- as_tibble(st_coordinates(destinations))
LongRad1 <- orig_coords$X*pi/180
LatRad1 <- orig_coords$Y*pi/180
LongRad2 <- dest_coords$X*pi/180
LatRad2 <- dest_coords$Y*pi/180
LongDiff <- abs(LongRad2-LongRad1)
DistRad <- acos(sin(LatRad2)*sin(LatRad1)+cos(LatRad2)*cos(LatRad1)*cos(LongDiff))
dresults <- DistRad*3437.74677*1.852
head(distances); head(dresults) # same

# find travel times between origins and destinations
# load raster from D_mixing_matrices.R
raster_layer <- raster("./01_data/MAP_friction_surface/getRaster/2015_friction_surface_v1_Decompressed_latest_.17.5_.34.8_51.41_27.29_2023_03_07.tiff")

# function to calculate travel times by country instead of all at once
find_times <- function(country){
  
  print(paste0("Running country ", country))
  
  # crop the friction surface rasterlayer to the random_points
  sample <- coords |> filter(COUNTRY == country)
  bbox <- bbox(sf::as_Spatial(sample))
  bbox_sp <- as(extent(bbox), "SpatialPolygons")
  
  # draw a buffer zone around the bounding box to expand it
  buffer_size <- 1 # set the buffer size in arc_degree
  bbox_buffered <- st_buffer(st_as_sfc(bbox_sp), dist = buffer_size)
  
  # assign crs to bounding box
  st_crs(bbox_buffered) <- st_crs(raster_layer)
  
  # crop rasterlayer
  raster_layer_cropped <- crop(raster_layer, st_as_sf(bbox_buffered))
  
  # make and geocorrect the transition matrix 
  T <- gdistance::transition(raster_layer_cropped, function(x) 1/mean(x), 8) 
  T.GC <- gdistance::geoCorrection(T)                    
  
  # extract lat / longs of random sample points
  points_origins <- as.matrix(st_coordinates(origins |> filter(COUNTRY == country)))
  points_destination <- as.matrix(st_coordinates(destinations |> filter(COUNTRY == country)))
  
  # run accumulated cost algorithm to calculate travel times from points to all other points
  access.matrix <- gdistance::costDistance(x = T.GC, 
                                           fromCoords = points_origins, 
                                           toCoords = points_destination)
  
  # pull out only origin - destination points of interest
  travel_times <- tibble(travel_times = diag(access.matrix)); summary(travel_times)
  
  return(travel_times)
  
}

# run function
travel_times <- map_dfr(.x = unique(coords$COUNTRY), .f = find_times)

# combine datasets & add destination population size for each trip
trip_data <- bind_cols(trips, distances, travel_times) |> 
  left_join(st_drop_geometry(coords), by = c("COUNTRY", "DEST" = "INDEX"))

saveRDS(trip_data, "./01_data/Marshall et al. 2018/trip_data.rds")

# visualize correlation between distances and travel times
ggplot(data = trip_data) + 
  geom_point(aes(x = units::drop_units(distances), y = travel_times, color = COUNTRY), alpha = 0.5, size = 0.8) + 
  geom_abline(slope = 1, intercept = 0, lty = 2) + 
  labs(x = "distance (km)", y = "travel time (min)", color = "country") + 
  scale_color_met_d(name = "Hiroshige") + 
  theme_classic() 


# Travel times matrices to export to Matlab ------------------------------------
# function to calculate travel times by country instead of all at once
find_times2 <- function(country){
  
  print(paste0("Running country ", country))
  
  # crop the friction surface rasterlayer to the random_points
  sample <- coords |> filter(COUNTRY == country)
  bbox <- bbox(sf::as_Spatial(sample))
  bbox_sp <- as(extent(bbox), "SpatialPolygons")
  
  # draw a buffer zone around the bounding box to expand it
  buffer_size <- 1 # set the buffer size in arc_degree
  bbox_buffered <- st_buffer(st_as_sfc(bbox_sp), dist = buffer_size)
  
  # assign crs to bounding box
  st_crs(bbox_buffered) <- st_crs(raster_layer)
  
  # crop rasterlayer
  raster_layer_cropped <- crop(raster_layer, st_as_sf(bbox_buffered))
  
  # make and geocorrect the transition matrix 
  T <- gdistance::transition(raster_layer_cropped, function(x) 1/mean(x), 8) 
  T.GC <- gdistance::geoCorrection(T)                    
  
  # extract lat / longs of random sample points
  points_origins <- as.matrix(st_coordinates(coords |> filter(COUNTRY == country)))
  points_destination <- as.matrix(st_coordinates(coords |> filter(COUNTRY == country)))
  
  # run accumulated cost algorithm to calculate travel times from points to all other points
  access.matrix <- gdistance::costDistance(x = T.GC, 
                                           fromCoords = points_origins, 
                                           toCoords = points_destination)
  
  # save matrix as .csv in Matlab folder
  write_csv(as_tibble(access.matrix), paste0(HPCpath, "./Marshall et al. 2018 code/traveltime_", country, ".csv"), col_names = F)
  
  print(paste0("finished calculating travel times for ", country))
  
}

# run function
map(.x = unique(coords$COUNTRY), .f = find_times2)


# continue model fits in Matlab



