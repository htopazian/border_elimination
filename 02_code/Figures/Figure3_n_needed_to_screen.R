# Figure 3. Number needed to screen and PAR density

# libraries
source("./02_code/packages_paths.R")

# pull in shapefiles from A_shapefiles_centroids.R
countries <- readRDS("./03_output/countries.rds")
centroids_border <- readRDS("./03_output/centroids_border.rds")
centroids_all <- readRDS("./03_output/centroids_all.rds")

# pull in population estimates from B_population.R
population <- readRDS("./03_output/population.rds") 

# read in cases averted dataset from G_processing.R
avert <- readRDS("./03_output/avert_output.rds")


# PLOT: nnt --------------------------------------------------------------------
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
  dplyr::select(ID_0, NAME_1, seed_point, draw, n_inf_detect, nnt_detect, nnt, pfpr_2_10, pfpr_2_10_baseline)

summary(output$nnt_detect); summary(output$nnt)

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
max <- max(output[output$nnt_seed_2 < Inf,]$nnt_seed_2)

output <- output |>
  mutate(nnt_bin = case_when(nnt_seed_2 > 100000000 ~ 100000001,
                             TRUE ~ nnt_seed_2))

summary(output$nnt_bin)  

# plot
p1 <- ggplot(output, aes(x = domestic, y = international)) +
  geom_point(aes(color = nnt_bin)) +
  scale_color_stepsn(colors = met.brewer("Hiroshige", n = 12)[seq(6, 11, 1)],
                    breaks = c(10, 1000, 10000, 1000000, 100000000),
                    values = scales::rescale(c(10, 1000, 10000, 1000000, 100000000)),
                    labels = c("10", "1,000", "10,000", "1,000,000", "100,000,000"),
                    guide = "colorsteps",
                    limits = c(0, 100000001)) +
  labs(x = expression(italic(Pf)~PR[2-10]~" seed point"), 
       y = expression(italic(Pf)~PR[2-10]~" international neighbors"), 
       color = "number screened \nto prevent one case") + 
  theme_classic() + 
  theme(legend.text = element_text(size = 7))

p1

ggsave(plot = p1, filename = "./03_output/nnt.pdf", width = 6, height = 4)



# PLOT: density weighted by par by pfpr differences ----------------------------

# select population for each admin1 unit from site files
pop_pfpr_lookup <- function(ISO){
  
  # pull country sites
  site <- eval(parse(text = paste0("foresite::", ISO)))
  
  # extract population from year 2023
  site_pop <- site$population[site$population$year == 2023, ]
  
  # extract pfpr from year 2023
  site_pfpr <- site$prevalence[site$prevalence$year == 2019, c(2, 3, 4, 6)]
  
  # sum by admin1 across urban / rural strata
  pop <- site_pop |> group_by(iso3c, name_1) |> 
    summarize(pop = sum(pop, na.rm = T),
              par = sum(par, na.rm = T))
  pfpr <- site_pfpr |> group_by(iso3c, name_1) |> slice(1) |> dplyr::select(-urban_rural)
  
  # merge
  dat <- left_join(pfpr, pop, by = c("iso3c", "name_1"))
  
  return(dat)
  
}

# loop through function for all malaria endemic countries
pop_pfpr_dat <- map_dfr(unique(countries$ID_0), pop_pfpr_lookup)

# extract each pair-wise combination of nearest intl neighbouring admin with each border site
select_cluster <- function(x){ # x = index of seed point
  
  # randomly select a seed point from the list of all admin1 border centroids
  seed_point <- centroids_border[x, ]
  
  # select out all admin1 centroids in countries other than the seed_point (B)
  centroids_seed_B <- centroids_all[centroids_all$ID_0 != seed_point$ID_0,]
  
  # calculate the distance from seed point to all other centroids in set B
  dist_to_seed_B <-  set_units(st_distance(seed_point, centroids_seed_B), km)
  
  # select the point nearest to the seed point but across the border
  nearest_points_B <- centroids_seed_B[dist_to_seed_B <= sort(dist_to_seed_B)[1] & dist_to_seed_B != set_units(0, km), ]
  
  # combine seed point with nearest neighbor for a group of 2
  random_sample <- rbind(seed_point, nearest_points_B) |> 
    st_drop_geometry() |> 
    mutate(seed = x)
  
  print(x)
  
  return(random_sample)

}

# loop through, selecting each border centroid as a seed point
clusters <- map_dfr(c(1:nrow(centroids_border)), select_cluster)

# estimate PfPR Δ and PAR
dat <- clusters |> 
  left_join(pop_pfpr_dat, by = c("ID_0" = "iso3c", "NAME_1" = "name_1")) |>
  group_by(seed) |>
  mutate(ID_0_seed = ID_0[1],
         NAME_1_seed = NAME_1[1],
         pop_total = sum(pop), 
         par_total = sum(par),
         pop_seed = case_when(row_number() == 1 ~ pop, TRUE ~ NA_real_),
         par_seed = case_when(row_number() == 1 ~ par, TRUE ~ NA_real_),
         pfpr_abs = max(pfpr) - min(pfpr),
         pfpr_diff = pfpr - lead(pfpr)) |>
  group_by(seed, ID_0_seed, NAME_1_seed) |>
  summarize(across(pop_total:pfpr_diff, \(x) mean(x, na.rm = TRUE)))

p2 <- ggplot(data = dat) + 
  labs(x = expression("difference in "~italic(Pf)~PR[2-10]~" between seed and the nearest international neighbor"),
       y = "density, weighted by population at risk") + 
  geom_density(data = dat, aes(x = pfpr_diff, weight = par_seed), color = "#e76254", fill = "#ef8a47", alpha = 0.8) + 
  theme_classic()

ggsave(plot = p2, filename = "./03_output/par_pfpr.pdf", width = 6, height = 4)


# save plots together ----------------------------------------------------------
p3 <- p1 + 
  
  (p2 + labs(
  x = expression(Delta~italic(Pf)~PR[2-10]~" between seed and nearest border neighbor"),
  y = "density, weighted by population at risk")) + 
  
  plot_layout(guides = "keep") + 
  plot_annotation(tag_levels = "A")

ggsave(plot = p3, filename = "./03_output/nnt_par.pdf", 
       width = 11, height = 4)

