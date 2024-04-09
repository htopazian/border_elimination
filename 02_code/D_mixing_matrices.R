# Creating mixing matrices between clusters.

# libraries
source("./02_code/packages_paths.R")

# pull in shapefiles from A_shapefiles_centroids.R
countries <- readRDS("./03_output/countries.rds")
admin1s <- readRDS("./03_output/admin1s.rds")
centroids_border <- readRDS("./03_output/centroids_border.rds")
borders_inner <- readRDS("03_output/borders_inner.rds") 
border_admins <- readRDS("./03_output/border_admins.rds")

# pull in shapefiles from B_population.R
population <- readRDS("./03_output/population.rds")

# pull in gravity model fits from MISC_trip_duration_travel_time_fits.R
model_duration <- readRDS("./03_output/dist_duration_quant_model.rds")
model_tt <- readRDS("./03_output/traveltime_duration_quant_model.rds")

# pull in probability of staying values from MISC_travel_prob.R
admin1s_rematch <- readRDS("03_output/admin1s_travel.rds")

# load friction surface from MAP
# tutorial from Amelia B-V: https://medium.com/@abertozz/mapping-travel-times-with-malariaatlas-and-friction-surfaces-f4960f584f08
# friction <- malariaAtlas::getRaster(
#   surface = "A global friction surface enumerating land-based travel speed for a nominal year 2015",
#   file_path = "./01_data/MAP_friction_surface/",
#   shp = sf::as_Spatial(countries))

# malariaAtlas::autoplot_MAPraster(friction) # visualize

# load raster from process ^
raster_layer_tt <- raster("./01_data/MAP_friction_surface/getRaster/2015_friction_surface_v1_Decompressed_latest_.17.5_.34.8_51.41_27.29_2023_03_07.tiff")

# import worldpop pop count data: unconstrained global mosaic 2020 (1km resolution)
# https://hub.worldpop.org/project/categories?id=3
raster_layer_pop <- raster("./01_data/ppp_2020_1km_Aggregated.tif")


# HPC set-up -------------------------------------------------------------------
setwd(HPCpath)

# set source function
hipercow_environment_create(sources = "function_mixing_matrices.R")

# install packages from pkgdepends.txt
hipercow_provision()

# Set up your job --------------------------------------------------------------
# x = index of seed point; y = distance (1) or travel time (2) used to calculate mixing
x <- rep(c(1:nrow(centroids_border)), 2)
y <- c(rep(1, nrow(centroids_border)), rep(2, nrow(centroids_border)))

index <- tidyr::crossing(x, y)

# Run tasks --------------------------------------------------------------------
# submit test
id <- task_create_expr(create_mix_matrix(13, 1))
task_status(id)
task_log_show(id)

# in bulk
bundle <- task_create_bulk_expr(create_mix_matrix(x, y), index)
hipercow_bundle_status(bundle)

# remove ones that have already been run
index <- index |> mutate(f = paste0(HPCpath, "mixing_matrices_grid/mixing_matrix_", x, "_", y, ".rds")) |>
  mutate(exist = case_when(file.exists(f) ~ 1, !file.exists(f) ~ 0)) |>
  filter(exist == 0) |>
  dplyr::select(-f, -exist)


# Process ----------------------------------------------------------------------

# read in diagonals 
read_matrix <- function(x, y){
  m <- readRDS(paste0(HPCpath, "mixing_matrices_grid/mixing_matrix_", x, "_", y, ".rds"))
  as.data.frame(diag(m))
}

# cycle through all clusters
diagonals <- map2_dfr(index$x, 2, .f = read_matrix)

summary(diagonals)


# Plot mixing matrix for select cluster ----------------------------------------

# read in cluster
random_sample <- readRDS(paste0(HPCpath, "clusters/cluster_180.rds"))
mix <- readRDS(paste0(HPCpath, "mixing_matrices/mixing_matrix_180_2.rds"))

# remove mixing within admin unit
I <- diag(1, ncol = nrow(mix), nrow = nrow(mix))
mix[I == 1] <- 0

mix <- mix |>
  as.data.frame() |>
  pivot_longer(cols = everything(), names_to = "destination", values_to = "flow")

mix$country_origin <- unlist(lapply(random_sample$ID_0, 
                                    function(x){rep(x, nrow(random_sample))}))
mix$admin1_origin <- unlist(lapply(random_sample$NAME_1, 
                                   function(x){rep(x, nrow(random_sample))}))
mix$country_destination <- rep(random_sample$ID_0, nrow(random_sample))
mix$admin1_destination <- rep(random_sample$NAME_1, nrow(random_sample))

mix <- mix |>
  # connect geometry of origin points
  left_join(random_sample, by = c("country_origin" = "ID_0", "admin1_origin" = "NAME_1")) |>
  rename(centroid_origin = geometry) |>
  # connect geometry of destination points
  left_join(random_sample, by = c("country_destination" = "ID_0", "admin1_destination" = "NAME_1")) |>
  rename(centroid_destination = geometry) |>
  
  # add population 
  left_join(population, by = c("country_destination" = "iso3c", "admin1_destination" = "name_1")) |>
  
  # set lat and long coordinates for origin and destination and calculate change
  mutate(long_origin = unlist(map(centroid_origin, 1)),
         lat_origin = unlist(map(centroid_origin, 2)),
         long_destination = unlist(map(centroid_destination, 1)),
         lat_destination = unlist(map(centroid_destination, 2)),
         d_lat = (lat_destination - lat_origin),
         d_long = (long_destination - long_origin)) |>
  
  # remove points matched with themselves
  filter(flow > 0) |>
  # categorize flow
  mutate(flow_s = flow / max(mix$flow) * 100, 
         flow_f = case_when(flow > 0 & flow <= 0.001 ~ "less movement",
                            flow > 0.001 & flow <= 0.1 ~ "more movement",
                            flow > 0.1 ~ "max movement"),
         flow_f = factor(flow_f)) 

summary(mix$flow); summary(mix$flow_s)
table(mix$flow_f)

# plot
ggplot() + 
  geom_sf(data = admin1s, fill = "cornsilk2", color = "cornsilk3") + 
  geom_sf(data = border_admins) + 
  geom_sf(data = borders_inner, size = 1.5, color = "navy") + 
  geom_sf(data = random_sample |> left_join(population, 
                                         by = c("ID_0" = "iso3c", "NAME_1" = "name_1")), aes(size = pop)) + 
  geom_segment(data = mix,
               aes(x = long_origin, y = lat_origin,
                   xend = long_destination, yend = lat_destination,
                   color = flow_s),
               linewidth = 1,
               alpha = 0.2,
               arrow = grid::arrow(length = unit(0.05, "cm"))) +
  scale_color_met_c(palette_name = "Hiroshige") + 
  labs(color = "normalized mixing \nprobability", size = "population", x = NULL, y = NULL) +
  coord_sf(xlim = st_bbox(random_sample)[c(1,3)], 
           ylim = st_bbox(random_sample)[c(2,4)]) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank())

ggsave(filename = "./03_output/seed_sample_flow_diagram.pdf", 
       width = 5, height = 5); beepr::beep(1)


# plot raster example
bbox <- bbox(sf::as_Spatial(random_sample))
bbox_sp <- as(extent(bbox), "SpatialPolygons")
buffer_size <- 1 
bbox_buffered <- st_buffer(st_as_sfc(bbox_sp), dist = buffer_size)
st_crs(bbox_buffered) <- st_crs(raster_layer_tt)
raster_layer_cropped <- crop(raster_layer_tt, st_as_sf(bbox_buffered))

# change into dataframe to use ggplot
as(raster_layer_cropped, "SpatialPixelsDataFrame")
raster_layer_cropped <- as.data.frame(raster_layer_cropped, xy = TRUE, na.rm = TRUE)

# plot
ggplot() + 
  geom_raster(data = raster_layer_cropped, aes(x = x, y = y, fill = raster_layer_cropped[,3])) + 
  geom_sf(data = random_sample) + 
  scale_fill_met_c(palette_name = "Hiroshige") + 
  labs(fill = "travel times \n(pixel, minutes)", x = NULL, y = NULL) +
  coord_sf(xlim = st_bbox(random_sample)[c(1,3)], 
           ylim = st_bbox(random_sample)[c(2,4)]) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank())

ggsave(filename = "./03_output/seed_sample_travel_times.pdf", 
       width = 5, height = 5); beepr::beep(1)
