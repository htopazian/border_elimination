# Generating shapefiles for admin1 units, identifying border areas and centroids.

# libraries
source("./02_code/packages_paths.R")


# Shapefiles SSA ---------------------------------------------------------------

# define countries by ISO
# also downloadable at: https://geodata.ucdavis.edu/gadm/gadm4.0/pck/
Africa <- c("DZA","AGO","BEN","BWA","BFA","BDI","CMR","CPV","CAF","TCD","COM",
            "COD","DJI","EGY","GNQ","ERI","ETH","GAB","GMB","GHA","GIN","GNB",
            "CIV","KEN","LSO","LBR","LBY","MDG","MWI","MLI","MRT","MUS","MYT",
            "MAR","MOZ","NAM","NER","NGA","COG","REU","RWA","SHN","STP","SEN",
            "SYC","SLE","SOM","ZAF","SSD","SDN","SWZ","TZA","TGO","TUN","UGA",
            "ESH","ZMB","ZWE")

# save copies of gadm SpatVector as .rds files in data folder
import_gadm <- function(ISO, level){
  geodata::gadm(country = ISO, level = level, path = "./01_data/gadm", version = "4.0")
}

map2(Africa, 0, import_gadm)    # admin0
map2(Africa, 1, import_gadm)    # admin1


# create a list of all countries
countries_africa <- lapply(c(Africa), 
                           function(x){
                             list.files(path = "./01_data/gadm",
                                        pattern = paste0("*", x , "_0_pk.rds"), 
                                        full.names = TRUE)
                           }) |> 
  unlist()

# create a list of all admin1s
admin1s_africa <- lapply(c(Africa), 
                         function(x){
                           list.files(path = "./01_data/gadm",
                                      pattern = paste0("*", x , "_1_pk.rds"), 
                                      full.names = TRUE)
                         }) |> 
  unlist()

# unpack gadms
unpack_gadm <- function(file){
  
  object <- readRDS(file) # read in object
  object <- terra::vect(object) # unpack SpatVector
  st_as_sf(object) # transform to sf object
  
}

countries_africa <- map_dfr(countries_africa, unpack_gadm) # loop over each country
admin1s_africa <- map_dfr(admin1s_africa, unpack_gadm) # loop over each admin1
st_crs(countries_africa); st_crs(countries_africa) # view CRS


# select out African countries with malaria (i.e. countries in foresite package)
countries <- countries_africa[!(countries_africa$ID_0 %in% c("DZA", "CPV", "EGY", "LSO", "LBY", "MUS", "MYT", "MAR", "REU", "SHN", "STP", "SYC", "TUN", "ESH")), ]

admin1s <- admin1s_africa[!(admin1s_africa$ID_0 %in% c("DZA", "CPV", "EGY", "LSO", "LBY", "MUS", "MYT", "MAR", "REU", "SHN", "STP", "SYC", "TUN", "ESH")), ]

# select only mainland admin units
# border_graph <- st_intersects(admin1s, admin1s)
# islands <- components(graph_from_adj_list(border_graph))
# mainland <- which.max(islands$csize)
# admin1s <- admin1s[islands$membership == mainland,]

saveRDS(countries, "./03_output/countries.rds") # save
# countries <- readRDS("./03_output/countries.rds")
saveRDS(countries_africa, "./03_output/countries_africa.rds") # save
saveRDS(admin1s, "./03_output/admin1s.rds") # save
saveRDS(admin1s, paste0(HPCpath, "admin1s.rds")) # save for HPC
# admin1s <- readRDS("./03_output/admin1s.rds")


# Country borders --------------------------------------------------------------

# define borders between country polygons
# borders <- st_boundary(countries) # will find all outside and inside borders
# border_admins <- admin1s[unique(unlist(st_touches(borders, admin1s))),]

borders_inner <- countries |> # will find inner borders
  rmapshaper::ms_innerlines() |>
  # as_tibble() |>
  st_as_sf()

saveRDS(borders_inner, "03_output/borders_inner.rds") # save

# find intersections between admin1 units and inner borders
admin1s_intersect <- st_intersects(admin1s, borders_inner)

index <- tibble(index = as.numeric(admin1s_intersect %>% lengths > 0))

# filter to only those admin1s which are on the borders
border_admins <- bind_cols(admin1s, index) |> 
  mutate(nrow = row_number()) |> 
  # remove admins not touching the inner borders
  filter(index == 1) |> 
  # remove admin1s missed by the function (Madagascar, South Africa, Nigeria)
  filter(!(nrow %in% c(482, 490, 305, 307, 424))) 

saveRDS(border_admins, "03_output/border_admins.rds") # save


# visualize
p <- ggplot() +
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
 
ggsave(plot = p, filename = "./03_output/borders_Africa.pdf", 
       width = 5, height = 5); beepr::beep(1)



# Centroids --------------------------------------------------------------------
# find centroids for each admin1
centroids_all <- sf::st_centroid(admin1s) |>
  dplyr::select(ID_0, NAME_1, geometry) 

# find centroids for each border admin1
centroids_border <- sf::st_centroid(border_admins) |>
  dplyr::select(ID_0, NAME_1, geometry) 

saveRDS(centroids_all, "03_output/centroids_all.rds") # save
saveRDS(centroids_border, "03_output/centroids_border.rds") # save

# visualize
p <- ggplot() + 
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

ggsave(plot = p, filename = "./03_output/centroids_Africa.pdf", 
       width = 5, height = 5); beepr::beep(1)



# Distance between centroids ---------------------------------------------------

# examine distribution of distances between admin1 units within countries
find_distances <- function(ISO){
  
  centroids <- centroids_all |> filter(ID_0 == {{ISO}}) 
  
  # find distance between centroids
  m_distance <- st_distance(centroids)
  km_distance <- set_units(m_distance, km)
  
  # replace the distance of 0 m with NA
  km_distance[upper.tri(km_distance, diag = T)] <- NA
  
  # mean(km_distance, na.rm = T)
  # hist(km_distance)
  
  distance <- na.omit(as.numeric(km_distance))
  distance_table <- tibble(ISO = {{ISO}}, distance = distance)
  
  return(distance_table)
  
}

distance <- map_dfr(unique(centroids_all$ID_0), find_distances)

# plot histogram of distances
ggplot(data = distance) + 
  geom_histogram(aes(x = distance, fill = ISO), alpha = 0.4, position = "stack") + 
  theme_classic()
