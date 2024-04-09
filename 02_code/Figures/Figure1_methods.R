# Figure 1. Methods diagrams

# libraries
source("./02_code/packages_paths.R")

# pull in shapefiles from A_shapefiles_centroids.R
countries_africa <- readRDS("./03_output/countries_africa.rds")
countries <- readRDS("./03_output/countries.rds")
admin1s <- readRDS("./03_output/admin1s.rds")
border_admins <- readRDS("./03_output/border_admins.rds")
centroids_all <- readRDS("./03_output/centroids_all.rds")
centroids_border <- readRDS("./03_output/centroids_border.rds")
borders <- readRDS("./03_output/borders.rds")
borders_inner <- readRDS("03_output/borders_inner.rds") 

# pull in population estimates from B_population.R
population <- readRDS("./03_output/population.rds") 


# border polygons
A <- ggplot() +
  geom_sf(data = countries_africa) + 
  geom_sf(data = countries, fill = "cornsilk2", color = "cornsilk3") + 
  geom_sf(data = border_admins, fill = "tomato") + 
  geom_sf(data = borders_inner, color = "navy", size = 1) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank())

# centroids SSA
B <- ggplot() + 
  geom_sf(data = countries_africa) + 
  geom_sf(data = countries, fill = "cornsilk2", color = "cornsilk3") + 
  geom_sf(data = borders_inner, color = "navy", size = 1) + 
  geom_sf(data = centroids_all, color = "goldenrod1", size = 0.8) + 
  geom_sf(data = centroids_border, color = "tomato", size = 0.8) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank())

# seed point cluster
random_sample <- readRDS(paste0(HPCpath, "clusters/cluster_180.rds"))
x = 180; seed_point <- centroids_border[x, ]

C <- ggplot() + 
  geom_sf(data = admin1s, fill = "cornsilk2", color = "cornsilk3") + 
  geom_sf(data = border_admins) + 
  geom_sf(data = borders_inner, size = 1.5, color = "navy") + 
  geom_sf(data = random_sample[1:4,], color = "orange") +
  geom_sf(data = seed_point, color = "red") + 
  geom_sf(data = random_sample[5:8,], color = "blue") + 
  coord_sf(xlim = st_bbox(random_sample)[c(1,3)], 
           ylim = st_bbox(random_sample)[c(2,4)]) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank())

# plot travel time raster example
raster_layer <- raster("./01_data/MAP_friction_surface/getRaster/2015_friction_surface_v1_Decompressed_latest_.17.5_.34.8_51.41_27.29_2023_03_07.tiff")

bbox <- bbox(sf::as_Spatial(random_sample))
bbox_sp <- as(extent(bbox), "SpatialPolygons")
buffer_size <- 1 
bbox_buffered <- st_buffer(st_as_sfc(bbox_sp), dist = buffer_size)
st_crs(bbox_buffered) <- st_crs(raster_layer)
raster_layer_cropped <- crop(raster_layer, st_as_sf(bbox_buffered))

# change into dataframe to use ggplot
as(raster_layer_cropped, "SpatialPixelsDataFrame")
raster_layer_cropped <- as.data.frame(raster_layer_cropped, xy = TRUE, na.rm = TRUE)

grid <- st_make_grid(x = bbox_buffered, cellsize = .1, what = "polygons", square = TRUE) |>
  st_as_sf(crs = st_crs(4326))

# identify grid cell centroids
grid_centroids <- st_centroid(grid) |> 
  # assign centroids their admin1 location
  st_join(admin1s[,1:4])

# plot
D <- ggplot() + 
  geom_raster(data = raster_layer_cropped, aes(x = x, y = y, fill = raster_layer_cropped[,3])) + 
  # geom_sf(data = grid, fill = NA, color = "navy") + 
  # geom_sf(data = grid_centroids, color = "navy", size = .5) + 
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


# calculate population at grid centroids
raster_layer_pop <- raster("./01_data/ppp_2020_1km_Aggregated.tif")
raster_layer_pop_cropped <- crop(raster_layer_pop, st_as_sf(bbox_buffered))
pop <- exactextractr::exact_extract(raster_layer_pop_cropped, grid, fun = "sum", progress = FALSE)

# identify grid cell centroids
grid_centroids2 <- cbind(grid_centroids, pop) |>
  mutate(pop_f = case_when(pop < 5000 ~ "< 5,000",
                           pop >= 5000 & pop < 50000 ~ "5,000 to 50,000",
                           pop > 50000 ~ "> 50,000"),
         pop_f = factor(pop_f, levels = c("< 5,000", "5,000 to 50,000", "> 50,000")))

E <- ggplot() + 
  geom_sf(data = admin1s, fill = "cornsilk2", color = "cornsilk3") + 
  geom_sf(data = border_admins) + 
  geom_sf(data = borders_inner, size = 1.5, color = "navy") + 
  geom_sf(data = grid, fill = NA, color = "navy") + 
  geom_sf(data = grid_centroids2, color = "#376795", aes(size = pop_f)) + 
  # scale_size_discrete(name = "population",
  #                      breaks = seq(0, 30000, 1190335),
  #                     labels = scales::comma,
  #                     range = c(0.2, 1.6)) +
  scale_size_discrete(name = "population",
                      range = c(0.2, 1)) +
  labs(size = "population", x = NULL, y = NULL) +
  coord_sf(xlim = st_bbox(random_sample)[c(1,3)], 
           ylim = st_bbox(random_sample)[c(2,4)]) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank())


# plot mixing flow and population size 
Pij_f <- readRDS(paste0(HPCpath, "mixing_matrices_grid/mixing_matrix_180_2.rds"))

# visualize flow between selected centroids
mix <- Pij_f 

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


F <- ggplot() + 
  geom_sf(data = admin1s, fill = "cornsilk2", color = "cornsilk3") + 
  geom_sf(data = border_admins) + 
  geom_sf(data = borders_inner, size = 1.5, color = "navy") + 
  geom_sf(data = random_sample, color = "#528fad") + 
  geom_segment(data = mix,
               aes(x = long_origin, y = lat_origin,
                   xend = long_destination, yend = lat_destination,
                   color = flow_s),
               linewidth = 1,
               alpha = 0.2,
               arrow = grid::arrow(length = unit(0.05, "cm"))) +
  scale_color_met_c(palette_name = "Hiroshige") +
  labs(color = "normalized mixing \nprobability", x = NULL, y = NULL) +
  coord_sf(xlim = st_bbox(random_sample)[c(1,3)], 
           ylim = st_bbox(random_sample)[c(2,4)]) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank())


# combine plots
plot_total <- A + C + E + B + D + F + 
  plot_layout(widths = c(1.5, 1, 1, 1.5, 1, 1), 
              nrow = 2, ncol = 3, guides = "keep") + 
  plot_annotation(tag_levels = list(c("A", "C", "E", "B", "D", "F")))

# save
ggsave(plot = plot_total, filename = "./03_output/methods.pdf", 
       width = 11, height = 6); beepr::beep(1)


