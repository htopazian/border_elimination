# Extracting population for each admin1 unit from site files.

# libraries
source("./02_code/packages_paths.R")

# pull in shapefiles from A_shapefiles_centroids.R
countries <- readRDS("./03_output/countries.rds")
admin1s <- readRDS("./03_output/admin1s.rds")

# Population -------------------------------------------------------------------
# select population for each admin1 unit from site files
pop_lookup <- function(ISO){
  
  # pull country sites
  site <- eval(parse(text = paste0("foresite::", ISO)))
  
  # extract population from year 2023
  site_pop <- site$population[site$population$year == 2023, ]
  
  # sum by admin1 across urban / rural strata
  pop <- site_pop |> group_by(iso3c, name_1) |> summarize(pop = sum(pop, na.rm = T))
  
  return(pop)
  
}

# loop through function for all malaria endemic countries
population <- map_dfr(unique(countries$ID_0), pop_lookup)

# verify that population ISO and admin1 names match with shapefiles for later use
admin1s_pop <- left_join(admin1s, population, by = c("ID_0" = "iso3c", "NAME_1" = "name_1")); summary(admin1s_pop$pop) # no NAs

# save population dataframe
saveRDS(population, "./03_output/population.rds") 
saveRDS(population, paste0(HPCpath, "./population.rds"))

