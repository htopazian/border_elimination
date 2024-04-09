# Parameterize clusters using site files data.

# libraries
source("./02_code/packages_paths.R")

# pull in shapefiles from A_shapefiles_centroids.R
countries <- readRDS('./03_output/countries.rds')
centroids_all <- readRDS("./03_output/centroids_all.rds")

# pull in net conversion values from MISC_netz.R
net_conversion <- readRDS("./03_output/net_conversion.rds")

# pull in SSA 2021 demography
africa_demog <- read.csv("./01_data/ssa_demography_2021.csv")


# Parameterize for malariasimulation using site files data ---------------------

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
  
  # make dataset and index by row
  site_list <- tibble(ID_0 = ISO, NAME_1 = admin, type = type) |>
    mutate(index = row_number())
  
  return(site_list)
  
}

all_sites <- map_dfr(unique(countries$ID_0), site_lookup)


# generate default parameters for each site
# pull in default vector species values from an example KEN site with all three
vectors <- site::single_site(foresite::"KEN", 5)$vector
vectors$prop <- rep(0, 3)
arabiensis <- vectors[1, 5:12]
funestus <- vectors[2, 5:12]
gambiae <- vectors[3, 5:12]
print(arabiensis); print(funestus); print(gambiae)



baseline_params <- function(ISO, index){
  
  year <- 365 
  
  # pull site
  site1 <- eval(parse(text = paste0("site::single_site(foresite::", ISO, ", ", index[1], ")")))
  
  if(!is.na(sum(index))){
    site2 <- eval(parse(text = paste0("site::single_site(foresite::", ISO, ", ", index[2], ")")))
  }
  
  site <- site1 # create copy for final version
  
  # combine site files if needed
  if(!is.na(sum(index))){
    
    pop1 <- site1$population$pop
    pop2 <- site2$population$pop
    popt <- site1$population$pop + site2$population$pop
    
    # prevalence
    site$prevalence$pfpr <- (site1$prevalence$pfpr * site1$population$par_pf[1:20] + site2$prevalence$pfpr * site2$population$par_pf[1:20]) / (site1$population$par_pf[1:20] + site2$population$par_pf[1:20])
    site$prevalence$pvpr <- (site1$prevalence$pvpr * site1$population$par_pv[1:20] + site2$prevalence$pvpr * site2$population$par_pv[1:20]) / (site1$population$par_pv[1:20] + site2$population$par_pv[1:20])
    
    # interventions
    site$interventions$itn_use <- (site1$interventions$itn_use * pop1[1:23] + site2$interventions$itn_use * pop2[1:23]) / popt[1:23]
    site$interventions$itn_input_dist <- (site1$interventions$itn_input_dist * pop1[1:23] + site2$interventions$itn_input_dist * pop2[1:23]) / popt[1:23]
    site$interventions$tx_cov <- (site1$interventions$tx_cov * pop1[1:23] + site2$interventions$tx_cov * pop2[1:23]) / popt[1:23]
    site$interventions$prop_act <- (site1$interventions$prop_act * pop1[1:23] + site2$interventions$prop_act * pop2[1:23]) / popt[1:23]
    site$interventions$prop_public <- (site1$interventions$prop_public * pop1[1:23] + site2$interventions$prop_public * pop2[1:23]) / popt[1:23]
    site$interventions$irs_cov <- (site1$interventions$irs_cov * pop1[1:23] + site2$interventions$irs_cov * pop2[1:23]) / popt[1:23]
    site$interventions$hh_size <- (site1$interventions$hh_size * pop1[1:23] + site2$interventions$hh_size * pop2[1:23]) / popt[1:23]
    site$interventions$smc_cov <- (site1$interventions$smc_cov * pop1[1:23] + site2$interventions$smc_cov * pop2[1:23]) / popt[1:23]
    site$interventions$pmc_cov <- (site1$interventions$pmc_cov * pop1[1:23] + site2$interventions$pmc_cov * pop2[1:23]) / popt[1:23]
    site$interventions$rtss_cov <- (site1$interventions$rtss_cov * pop1[1:23] + site2$interventions$rtss_cov * pop2[1:23]) / popt[1:23]
    site$interventions$pyrethroid_resistance <- (site1$interventions$pyrethroid_resistance * pop1[1:23] + site2$interventions$pyrethroid_resistance * pop2[1:23]) / popt[1:23]
    
    # population
    site$population$pop <- popt
    site$population$par <- site1$population$par + site2$population$par
    site$population$par_pf <- site1$population$par_pf + site2$population$par_pf
    site$population$par_pv <- site1$population$par_pv + site2$population$par_pv
    
    # pyrethroid reistance
    site$pyrethroid_resistance$pyrethroid_resistance <- (site1$pyrethroid_resistance$pyrethroid_resistance * pop1 + site2$pyrethroid_resistance$pyrethroid_resistance * pop2) / popt
    
    # EIR
    site$eir$eir[1] <- (site1$eir$eir[1] * mean(site1$population$pop[1:23], na.rm = T) + site2$eir$eir[1] * mean(site2$population$pop[1:23], na.rm = T)) / (mean(site1$population$pop[1:23], na.rm = T) + mean(site2$population$pop[1:23], na.rm = T))
    
    # other vars should be the same between urban and rural at the admin1 level
    
    # replace NAs with 0
    site$prevalence$pfpr<- tidyr::replace_na(site$prevalence$pfpr, 0)
    site$prevalence$pvpr <- tidyr::replace_na(site$prevalence$pvpr, 0)
    
    site$interventions$itn_use <- tidyr::replace_na(site$interventions$itn_use, 0)
    site$interventions$itn_input_dist <- tidyr::replace_na(site$interventions$itn_input_dist, 0)
    site$interventions$tx_cov <- tidyr::replace_na(site$interventions$tx_cov, 0)
    site$interventions$prop_act <- tidyr::replace_na(site$interventions$prop_act, 0)
    site$interventions$prop_public <- tidyr::replace_na(site$interventions$prop_public, 0)
    site$interventions$irs_cov <- tidyr::replace_na(site$interventions$irs_cov, 0)
    site$interventions$hh_size <- tidyr::replace_na(site$interventions$hh_size, 0)
    site$interventions$smc_cov <- tidyr::replace_na(site$interventions$smc_cov, 0)
    site$interventions$pmc_cov <- tidyr::replace_na(site$interventions$pmc_cov, 0)
    site$interventions$rtss_cov <- tidyr::replace_na(site$interventions$rtss_cov, 0)
    site$interventions$pyrethroid_resistance <- tidyr::replace_na(site$interventions$pyrethroid_resistance, 0)
    
    site$pyrethroid_resistance$pyrethroid_resistance <- tidyr::replace_na(site$pyrethroid_resistance$pyrethroid_resistance, 0)
    
  }
  
  # add 10 years of constant interventions
  # start by pulling the last row of the tibble, year 2022
  last_row <- slice_tail(site$interventions, n = 1)
  
  # create a new tibble with 10 copies of the last row
  new_rows <- bind_rows(replicate(10, last_row, simplify = FALSE))
  
  # replace year variable with new years
  new_rows$year <- seq(2023, 2023 + 9, 1)
  
  # replace itn_input_dist values with constant dist values from netz
  itn_use <- tibble(itn_use = new_rows$itn_use[1])
  
  # use findInterval() to match the values in itn_use with the closest value in net_conversion (the index of the closest <= itn_use + 1)
  matched_value <- findInterval(itn_use$itn_use, net_conversion$usage) + 1
  
  itn_use$dist <- net_conversion$dist[matched_value]
  
  new_rows$itn_input_dist <- rep(itn_use$dist, nrow(new_rows))
  
  # combine the original tibble and the new rows
  site$interventions <- bind_rows(site$interventions, new_rows)
  
  # check that vector table lists all three species
  if(!("arabiensis" %in% c(site$vectors$species))){
    add_row <- slice_tail(site$vectors, n = 1)
    add_row[, 5:12] <- arabiensis
    site$vectors <- bind_rows(site$vectors, add_row)
  }
  
  if(!("funestus" %in% c(site$vectors$species))){
    add_row <- slice_tail(site$vectors, n = 1)
    add_row[, 5:12] <- funestus
    site$vectors <- bind_rows(site$vectors, add_row)
  }
  
  if(!("gambiae" %in% c(site$vectors$species))){
    add_row <- slice_tail(site$vectors, n = 1)
    add_row[, 5:12] <- gambiae
    site$vectors <- bind_rows(site$vectors, add_row)
  }
  
  site$vectors <- site$vectors |> arrange(species)
  
  print(paste0(ISO, " ", index))
  
  return(site)
  
}

# generate default parameters for each site
set_params <- function(site){
  
  site <- site 
  year <- 365 
  
  # if Pf EIR == 0, then assign a small number 
  # site cannot assign parameters if EIR == 0
  if(site$eir$eir[1] == 0){
      site$eir$eir[1] = 0.00001
  }

  # create parameter inputs
  site_par <- site::site_parameters(
    interventions = site$interventions,
    demography = site$demography,
    vectors = site$vectors,
    seasonality = site$seasonality,
    eir = site$eir$eir[1],
    overrides = list(
      human_population = 100000,
      prevalence_rendering_min_ages = 2 * year,
      prevalence_rendering_max_ages = 10 * year,
      # these commands do not override at the moment, issue created on package github
      clinical_incidence_rendering_min_ages = c(0, 0) * year, 
      clinical_incidence_rendering_max_ages = c(5, 200) * year,
      severe_incidence_rendering_min_ages = c(0, 0) * year,
      severe_incidence_rendering_max_ages = c(5, 200) * year)
  )
  
  # modify demography to match SSA 2021 UN values
  ages <- round(africa_demog$age_upper * year) # top of age bracket
  deathrates <- africa_demog$mortality_rate / 365   # age-specific death rates
  
  site_par <- set_demography(
    site_par,
    agegroups = ages,
    timesteps = 0,
    deathrates = matrix(deathrates, nrow = 1))

  print(paste0("generating params for ", site$country))
  
  return(tibble(ID_0 = site$country,
                NAME_1 = site$sites$name_1,
                params = list(site_par)))
}


# create dataframe of indices to merge
selected_site <- all_sites |> dplyr::select(-type) |> group_by(ID_0, NAME_1) |>
  mutate(rownum = row_number()) |>
  pivot_wider(names_from = rownum, values_from = index, names_prefix = "R") |>
  mutate(index = map2(R1, R2, ~c(.x, .y)))


# run through function
default_params <- map_dfr(
  map2(selected_site$ID_0, 
       selected_site$index,
       baseline_params),
  set_params)


# save 
saveRDS(default_params, "./03_output/site_params.rds")
saveRDS(default_params, paste0(HPCpath, "site_params.rds"))


# should EIRs be population weighted?
# some urban / rural sites have big differences but most of these have populations mostly in one area or another
EIR_test <- function(ISO, index){
  
  # pull site
  site1 <- eval(parse(text = paste0("site::single_site(foresite::", ISO, ", ", index[1], ")")))
  
  if(!is.na(sum(index))){
    
    site2 <- eval(parse(text = paste0("site::single_site(foresite::", ISO, ", ", index[2], ")")))
    
    pop1 <- mean(site1$population$pop, na.rm = T)
    pop2 <- mean(site2$population$pop, na.rm = T)
    popt <- pop1 + pop2
    
    EIRA <- site1$eir$eir[1]
    EIRB <- site2$eir$eir[1]
    
    EIRC <- (EIRA * pop1 + EIRB * pop2) / popt
    
    popmin <- min(pop1 / popt, pop2 / popt)
    
    return(tibble(EIRA, EIRB, EIRC, popmin))
  }
  
}

EIRdiff <- map2_dfr(selected_site$ID_0, 
                 selected_site$index,
                 EIR_test)

ggplot(data = EIRdiff) + 
  geom_point(aes(x = EIRC, y = abs(EIRA - EIRB), color = popmin), alpha = 0.2) +
  labs(x = "combined EIR", y = "difference in EIRs", color = "min population") + 
  scale_color_met_c(palette_name = "Hiroshige", dir = 1) +
  theme_bw()

