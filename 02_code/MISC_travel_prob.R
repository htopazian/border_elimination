# Determine travel probabilities from DHS data

# libraries
source("./02_code/packages_paths.R")


# log-in. Prior to running script, request access to survey and GPS data from all DHS surveys in sub-Saharan Africa
rdhs::set_rdhs_config(
  email = "h.topazian@imperial.ac.uk",
  project = "The use of surveillance test-and-treat posts to reduce malaria in border regions in sub-Saharan Africa",
  config_path = "rdhs.json",
  global = FALSE) # choose option 2

# pull in shapefiles from A_shapefiles_centroids.R
admin1s <- readRDS("./03_output/admin1s.rds")
admin1s <- as(admin1s, "Spatial")

# load available surveys
countries <- dhs_countries()

# subset to SSA
countries <- subset(countries, countries$ISO3_CountryCode %in% unique(admin1s$ID_0))$DHS_CountryCode
survs <- list(
  rdhs::dhs_surveys(countryIds = countries[c(1:20)]),
  rdhs::dhs_surveys(countryIds = countries[seq(21,length(countries))])
)

# function for importing survey and GPS
get_dhs_data <- function(types, vars, year) {
  
  # select those that meet this criteria
  datasets <- rbind(
    rdhs::dhs_datasets(
      surveyIds = survs[[1]]$SurveyId,
      fileFormat = "FL", # flat file
      fileType = types # individual recode
    ),
    rdhs::dhs_datasets(
      surveyIds = survs[[2]]$SurveyId,
      fileFormat = "FL", # flat file
      fileType = types)
  )

  
  # download datasets
  broken = c("MWPR22FL.ZIP", "KMPR32FL.ZIP", "GHCR31FL.ZIP", "KEIR03FL.ZIP", "SNIR02FL.ZIP")
  datasets <- subset(datasets, !(datasets$FileName %in% broken))
  downloads <- rdhs::get_datasets(datasets$FileName)
  datasets <- subset(datasets, downloads != "Dataset is not available with your DHS login credentials")
  
  questions <- rdhs::search_variables(
    datasets$FileName,
    variables = vars
  )
  
  questions <- subset(questions, questions$dataset_filename != "ugcr41fl")
  
  # extract variables of interest and add geographic information
  extract <- rdhs::extract_dhs(questions, add_geo = TRUE)
  
  dhs <- data.frame(data.table::rbindlist(lapply(extract, haven::zap_labels), fill=TRUE, use.names=TRUE))
  
  dhs <- dhs[dhs$LATNUM!= 0,]
  dhs <- dhs[!is.na(dhs[year]),]
  dhs[year][dhs[year] < 100] <- dhs[year][dhs[year] < 100] + 1900
  coordinates(dhs) <- ~ LONGNUM + LATNUM
  proj4string(dhs) <- proj4string(admin1s)
  matches <- over(dhs, admin1s)
  dhs$iso3c <- matches$ID_0
  dhs$name_1 <- matches$NAME_1
  dhs
}

# choose variables
vars <- c(
  "caseid", # Case identification number
  "v001",   # Cluster number
  "v023",   # Strata Weight
  "v005",   # Weight
  "v007",   # Year
  "v167",   # Number of trips away from home for one or more nights in last 12 months
  "v168"    # Away for more than one month at a time in last 12 months (base = v167 <> 0)
)

mvars <- c(
  "mcaseid",
  "mv001",
  "mv023",
  "mv005",
  "mv007",  
  "mv167",  
  "mv168"   
)

# run functions for male and female datasets
fdhs <- get_dhs_data("IR", vars, "v007")
mdhs <- get_dhs_data("MR", mvars, "mv007")

# save 
saveRDS(fdhs, "03_output/fdhs.rds")
saveRDS(mdhs, "03_output/mdhs.rds")

# fdhs <- readRDS("03_output/fdhs.rds")
# mdhs <- readRDS("03_output/mdhs.rds")

# functions for calculating the probability of stay
p_travel <- function(dhs_mini, id, strata, weights, trips){
  dhs_mini <- data.frame(dhs_mini)
  dhs_mini$trips <- dhs_mini[[trips]] # assign trips variable
  # create variables for surveys with all NAs for trips
  if (all(is.na(dhs_mini$trips))) {
    return(data.frame(n_no_trips = NA, se_no_trips = NA, med_trips = NA, total = NA))
  }
  # assign 99s to NA per codebook
  if (length(which(dhs_mini$trips == 99)) != 0) {
    dhs_mini[dhs_mini$trips == 99 & !is.na(dhs_mini$trips),]$trips <- NA 
  }
  # create complex survey design
  options(survey.lonely.psu = "adjust")
  DHSdesign <- svydesign( # full data
    id = dhs_mini[id],
    strata = dhs_mini[strata],
    weights = dhs_mini[weights]/1000000,
    data = dhs_mini,
    nest = T
  )
  dhs_mini2 <- dhs_mini[dhs_mini$trips > 0 & !is.na(dhs_mini$trips),]
  DHSdesign2 <- svydesign( # data for only those who took trips
    id = dhs_mini2[id],
    strata = dhs_mini2[strata],
    weights = dhs_mini2[weights]/1000000,
    data = dhs_mini2,
    nest = T
  )
  # create variables for surveys with all 0s for trips
  if (length(unique(dhs_mini$trips)) == 1) {
    return(
      data.frame(
        n_no_trips = sum(dhs_mini$trips == 0), 
        se_no_trips = NA,
        med_trips = NA,
        total = nrow(dhs_mini)
      )
    )
  }
  # create variables for surveys with trip data
  # to get 0 trips
  m <- svytotal(~ as.factor(trips), design = DHSdesign, na.rm = TRUE)
  m <- data.frame(m)
  # to get median n trips for people who take trips
  m3 <- svyquantile(~ trips, design = DHSdesign2, na.rm = TRUE, quantiles = c(0.5))
  no_trips = as.numeric(gsub("as.factor\\(trips\\)", "", row.names(m))) == 0
  if (sum(no_trips) == 0) { # if there are all trips in the dataset
    return(
      data.frame(
        n_no_trips = 0,
        se_no_trips = NA,
        med_trips = m3$trips[1],
        total = nrow(dhs_mini)
      )
    )
  }
  data.frame( # if there are trips in the dataset
    n_no_trips = m$total[no_trips],
    se_no_trips = m$SE[no_trips],
    med_trips = m3$trips[1],
    total = sum(m$total)
  )

}

map_travels <- function(dhs, id, strata, weights, trips, year) {
  df <- data.frame(dhs)
  df <- df[!is.na(df[strata]),]
  df["year"] <- df[year]
  df |>
    dplyr::group_by(SurveyId, iso3c, name_1, year) |>
    dplyr::group_modify(~ p_travel(.x, id, strata, weights, trips))
}

# run probability of travel functions for male and female datasets
fp <- map_travels(
  fdhs,
  "v001", # cluster number
  "v023", # strata
  "v005", # sample weight
  "v167", # n trips
  "v007"  # year
)

mp <- map_travels(
  mdhs,
  "mv001",
  "mv023",
  "mv005",
  "mv167",
  "mv007"
)

fp$sex <- "female"
mp$sex <- "male"

p = rbind(fp, mp)
p$p_stay <- p$n_no_trips / p$total

# save 
saveRDS(p, "03_output/prob_travel.rds") # save


# variation
mad_t <- aggregate(p$p_stay, by = p[c("iso3c", "name_1", "sex")], FUN = mad)
mad_s <- aggregate(p$p_stay, by = p[c("iso3c", "name_1", "year")], FUN = mad)
mad_n <- aggregate(p$p_stay, by = p[c("iso3c", "year", "sex")], FUN = mad)
mad_t <- mad_t[!is.na(mad_t$x),]
mad_s <- mad_s[!is.na(mad_s$x),]
mad_n <- mad_n[!is.na(mad_n$x),]

ggplot(
  data.frame(
    value=c(mad_t$x, mad_s$x, mad_n$x),
    across=c(rep("Year", nrow(mad_t)), rep("Sex", nrow(mad_s)), rep("Sub-national", nrow(mad_n)))
  )
) +
  geom_violin(aes(x=across, y=value, fill=across)) + 
  ggtitle("MAD of p_stay") + 
  theme(legend.position="none")

# find missing values
admin1s$n_value <- apply(
  as.data.frame(admin1s[c("ID_0", "NAME_1")]),
  1,
  function(row) {sum(!is.na(p[(p$iso3c == row[["ID_0"]]) & (p$name_1 == row[["NAME_1"]]),]$p_stay))}
)
# sp::spplot(admin1s, "n_value", main="Number of relevant surveys with results")

# p_stay_1 = mean of p_stay across all surveys, sexes, and years by admin1
admin1s$p_stay_1 <- apply(
  as.data.frame(admin1s[c("ID_0", "NAME_1")]),
  1,
  function(row) {
    surveys <- p[(p$iso3c == row[["ID_0"]] & p$name_1 == row[["NAME_1"]]),]
    sum(surveys$n_no_trips, na.rm = TRUE) / sum(surveys$total, na.rm = TRUE)
    }
)

# med_trips1 = mean of med_trips across all surveys, sexes, and years by admin1
admin1s$med_trips_1 <- apply(
  as.data.frame(admin1s[c("ID_0", "NAME_1")]),
  1,
  function(row) {
    surveys <- p[(p$iso3c == row[["ID_0"]] & p$name_1 == row[["NAME_1"]]),]
    mean(surveys$med_trips, na.rm = TRUE)
  }
)

# sp::spplot(admin1s, "p_stay_1", main ="Probability that no trips were taken in a year")

# p_stay_0 = mean of p_stay across all surveys, sexes, and years by admin0
admin1s$p_stay_0 <- apply(
  as.data.frame(admin1s["ID_0"]),
  1,
  function(iso) {
    surveys <- p[(p$iso3c == iso),]
    sum(surveys$n_no_trips, na.rm = TRUE) / sum(surveys$total, na.rm = TRUE)
  }
)

# med_trips_0 = mean of p_stay across all surveys, sexes, and years by admin0
admin1s$med_trips_0 <- apply(
  as.data.frame(admin1s["ID_0"]),
  1,
  function(iso) {
    surveys <- p[(p$iso3c == iso),]
    mean(surveys$med_trips, na.rm = TRUE)
  }
)

# p_stay_SSA = mean of p_stay across all surveys, sexes, years, and admin1s (= mean of SSA)
# med_trips_SSA = mean of med_trips across all surveys, sexes, years, and admin1s (= mean of SSA)
admin1s$p_stay_SSA <- sum(p$n_no_trips, na.rm = TRUE) / sum(p$total, na.rm = TRUE)
admin1s$med_trips_SSA <- mean(p$med_trips, na.rm = TRUE)


# missing units (Botswana, Central African Republic, Djibouti, Equatorial Guinea, Eritrea, Gambia, Guinea-Bissau, Madagascar, Mauritania, Congo, Somalia, South Sudan, Sudan, Uganda)
print(data.frame(admin1s[admin1s$n_value == 0, c("ID_0", "NAME_1")]))

# countries with no results (BWA, DJI, GNB, GNQ, MDG, SSD)
n_countries <- aggregate(admin1s$n_value, list(ID_0 = admin1s$ID_0), sum)
n_countries |> filter(x == 0)

# anomalous results - low case counts
print(data.frame(p)[!is.na(p$p_stay) & (p$p_stay == 0),])
print(data.frame(p)[!is.na(p$p_stay) & (p$p_stay == 1),])


# fix missing units and missing country results
print(data.frame(admin1s[admin1s$n_value == 0 & !(admin1s$ID_0 %in% c("BWA", "DJI", "GNB", "GNQ", "MDG", "SSD")), c("ID_0", "NAME_1")]))
# CAF and GMB and MRT and COG and SOM and SDN = only one admin1 with info in p
# ERI = checked, three missing data are missing data in p
# UGA = no Lake Albert value in p
# these are all likely mis-assigned values based on DHS cluster displacement since counts are low

# check countries where no DHS survey has been conducted
countries <- dhs_countries() # load available surveys
countries <- subset(countries, countries$ISO3_CountryCode %in% unique(admin1s$ID_0))$ISO3_CountryCode
unique(admin1s[!(admin1s$ID_0 %in% countries),]$ID_0) # "DJI" "GNB" "SOM" "SSD"

# fix mismatches and assign p_stay and med_trips values to missing units
admin1s_rematch <- admin1s |> 
  st_as_sf() |>
  mutate(p_stay = case_when(ID_0 %in% c("DJI", "GNB", "SOM", "SSD") ~ p_stay_SSA, # no DHS surveys
                            ID_0 %in% c("BWA", "GNQ", "MDG") ~ p_stay_SSA, # countries with no results
                            ID_0 %in% c("CAF", "ERI", "GMB", "MRT", "COG", "SDN") ~ p_stay_SSA, # countries with mismatches (cluster displacement)
                            is.na(p_stay_1) ~ p_stay_0,
                            TRUE ~ p_stay_1),
         med_trips = case_when(ID_0 %in% c("DJI", "GNB", "SOM", "SSD") ~ med_trips_SSA, # no DHS surveys
                            ID_0 %in% c("BWA", "GNQ", "MDG") ~ med_trips_SSA, # countries with no results
                            ID_0 %in% c("CAF", "ERI", "GMB", "MRT", "COG", "SDN") ~ med_trips_SSA, # countries with mismatches (cluster displacement)
                            is.na(med_trips_1) ~ med_trips_0,
                            TRUE ~ med_trips_1)
         )


summary(admin1s_rematch$p_stay); summary(admin1s_rematch$med_trips)


# save 
saveRDS(admin1s_rematch, "./03_output/admin1s_travel.rds")
saveRDS(admin1s_rematch, paste0(HPCpath, "admin1s_travel.rds"))


ggplot() +
  geom_sf(data = admin1s_rematch, aes(fill = p_stay)) + 
  scale_fill_met_c(name = "Hiroshige") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank())

ggplot() +
  geom_sf(data = admin1s_rematch, aes(fill = med_trips)) + 
  scale_fill_met_c(name = "Hiroshige") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank())


View(admin1s_rematch |> st_drop_geometry())
