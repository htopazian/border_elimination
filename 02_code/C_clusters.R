# Creating random clusters from border area seed points.

# libraries
source("./02_code/packages_paths.R")

# pull in shapefiles from A_shapefiles_centroids.R
centroids_all <- readRDS("./03_output/centroids_all.rds")
centroids_border <- readRDS("./03_output/centroids_border.rds")
borders <- readRDS("./03_output/borders.rds")
borders_inner <- readRDS("03_output/borders_inner.rds") 
border_admins <- readRDS("./03_output/border_admins.rds")

# Select clusters --------------------------------------------------------------

select_cluster <- function(x){ # x = index of seed point
  
  # randomly select a seed point from the list of all admin1 border centroids
  seed_point <- centroids_border[x, ]
  
  # select out all admin1 centroids in the same country as the seed_point (A) and 
  # all admin1 centroids in countries other than the seed_point (B)
  centroids_seed_A <- centroids_all[centroids_all$ID_0 == seed_point$ID_0,]
  centroids_seed_B <- centroids_all[centroids_all$ID_0 != seed_point$ID_0,]
  
  # calculate the distance from seed point to all other centroids in set A and set B
  dist_to_seed_A <-  set_units(st_distance(seed_point, centroids_seed_A), km)
  dist_to_seed_B <-  set_units(st_distance(seed_point, centroids_seed_B), km)
  
  # select the 7 points nearest to the seed point
  # three from set A (same country as seed_point) and four from set B
  nearest_points_A <- centroids_seed_A[dist_to_seed_A <= sort(dist_to_seed_A)[4] & dist_to_seed_A  != set_units(0, km), ] # select all points with a distance equal to or less than the 4th nearest value and above 0 (since the first value is always 0 (point with itself)). 
  nearest_points_B <- centroids_seed_B[dist_to_seed_B <= sort(dist_to_seed_B)[4] & dist_to_seed_B != set_units(0, km), ]
  
  # combine seed point with nearest neighbors for a group of 8
  random_sample <- st_as_sf(rbind(seed_point, nearest_points_A, nearest_points_B))
  
  # save
  saveRDS(random_sample, paste0(HPCpath, "clusters/cluster_", x, ".rds"))
  
  print(x)

}

# loop through, selecting each border centroid as a seed point
lapply(c(1:nrow(centroids_border)), select_cluster)

# visualize sample points - run the same steps as above for one point
random_sample <- readRDS(paste0(HPCpath, "clusters/cluster_180.rds"))
x = 180; seed_point <- centroids_border[x, ]

# quick look
ggplot() +
  geom_sf(data = random_sample, color = "blue") +
  geom_sf(data = seed_point, color = "red") +
  coord_sf(xlim = st_bbox(random_sample)[c(1,3)],
           ylim = st_bbox(random_sample)[c(2,4)]) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank())

ggsave(filename = "./03_output/seed_sample_points.pdf", 
       width = 5, height = 5); beepr::beep(1)


# visualize points with country borders
ggplot() + 
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

ggsave(filename = "./03_output/seed_sample_admins.pdf", 
       width = 5, height = 5); beepr::beep(1)

