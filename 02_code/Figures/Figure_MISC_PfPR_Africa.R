# Figure MISC. PfPR SSA 2022

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

# make a list of all foresite package site ISOs and admin1 names with index
# there will be 1099, many more than the list of admin1s because most have an
# urban and rural location listed
site_lookup <- function(ISO){
  
  # pull site
  site <- eval(parse(text = paste0("foresite::", ISO)))
  
  # extract ISO, admin1 name, and urban / rural vars
  ISO <- site$site$iso3c
  admin <- site$site$name_1
  type <- site$site$urban_rural
  PfPR <- site$prevalence[site$prevalence$year == 2019,]$pfpr
  
  # make dataset and index by row
  site_list <- tibble(ID_0 = ISO, NAME_1 = admin, type = type, PfPR = PfPR) |>
    mutate(index = row_number())
  
  return(site_list)
  
}

all_sites <- map_dfr(unique(countries$ID_0), site_lookup)

# condensing sites - if one admin1 unit has both urban and rural locations, preferentially selecting rural
condensed_sites <- all_sites |> group_by(ID_0, NAME_1) |> slice(1)

output <- condensed_sites |> 
  left_join(admin1s |> select(ID_0, NAME_1)) |>
  st_as_sf()

p <- ggplot() + 
  geom_sf(data = countries_africa) + 
  geom_sf(data = output, aes(fill = PfPR)) + 
  geom_sf(data = borders_inner, color = "navy", size = 1) + 
  # values between 0 and 1 with the middle number as the midpoint
  scale_fill_gradientn(colors = met.brewer("VanGogh3"), 
                       values = c(0, (range - max) / range, 1)) + 
  labs(fill = expression(italic(Pf)~PR[2-10]~", 2019")) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank())

ggsave(plot = p, filename = "./03_output/PfPR_Africa.pdf", 
       width = 5, height = 5); beepr::beep(1)

