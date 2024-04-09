# This script contains the functions necessary to generate various
# iterations of mixing matrices

# libraries
source("./02_code/packages_paths.R")

# pull in distance and travel times values from MISC_refitting_gravity_model.R
trip_data <- readRDS("./01_data/Marshall et al. 2018/trip_data.rds")


# Isolated ---------------------------------------------------------------------
# diagonal matrix

IM <- diag(100)


# Marshall et al. 2018 ---------------------------------------------------------
# https://www.nature.com/articles/s41598-018-26023-1

# function to import Marshall et a. 2016 Trips data
import_xlsx <- function(sheet){
  
  readxl::read_excel("./01_data/Marshall et al. 2016/Marshall_2016_additional_file_2.xlsx", sheet = sheet) |>
    mutate(Age = is.numeric(Age)) |>
    dplyr::select(TripID, PersonID, TripType, OrigCountry, DestCountry, 
                  OrigIndex, DestIndex,
                  Distance, Duration, Purpose, Gender, Age)
  
}

# import and merge
dat <- map_dfr(.x = c("ML Trips", "BF Trips", "ZM Trips", "TZ Trips"),
               .f = import_xlsx)

# remove NAs (-1s? a little unclear)
dat <- dat |> mutate(Duration = ifelse(Duration < 0, NA, Duration))
summary(dat$Duration)

# set max trip duration to one year
dat <- dat |> mutate(Duration = ifelse(Duration > 365, 365, Duration))
summary(dat$Duration)

test <- dat |> left_join(trip_data |>
                           dplyr::select(ORIG, DEST, COUNTRY, distances, travel_times), 
                         by = c("OrigIndex" = "ORIG", 
                                           "DestIndex" = "DEST", 
                                           "OrigCountry" = "COUNTRY"), 
                         multiple = "first")

# plot trip duration by country
ggplot(data = dat,
       aes(x = Distance, y = Duration)) +
  geom_point(alpha = 0.1,  size = 0.5, color = "darkgrey") +
  geom_quantile(quantiles = 0.5, aes(color = OrigCountry, group = OrigCountry)) +
  geom_quantile(quantiles = 0.5, aes(color = "Overall fit"), color = "black") +
  # geom_smooth(method = "glm.nb", aes(color = OrigCountry, group = OrigCountry)) +
  # geom_smooth(method = "glm.nb", aes(fill = "Overall fit"), color = "black") +
  # scale_fill_manual(values = c(NA)) +
  labs(x = "Distance (km)",
       y = "Duration (days)",
       color = "",
       fill = "") +
  coord_cartesian(ylim = c(0, 365)) +
  theme_classic() +
  theme(text = element_text(size = 14))

ggsave("./03_output/Marshall_etal_trip_duration.pdf", width = 6, height = 4)


# < fit trip duration distance traveled ----
# https://malariajournal.biomedcentral.com/articles/10.1186/s12936-016-1252-3#MOESM2
# 12936_2016_1252_MOESM2_ESM.xlsx Additional file 2, Marshall et al. 2016

# fit a model to duration data vs. distance traveled
# linear
model <- glm.nb(Duration ~ Distance, data = dat)
exp(model$coefficients); model$aic

# cubic splines 2 knot
model <- glm.nb(Duration ~ splines::ns(Distance, 2), data = dat)
exp(model$coefficients); model$aic

# cubic splines 6 knot
model <- glm.nb(Duration ~ splines::ns(Distance, 6), data = dat)
exp(model$coefficients); model$aic # best AIC
saveRDS(model, file = "./03_output/dist_duration_sk6_model.rds")

# quantile regression (median) - robust to outliers
model <- rq(Duration ~ Distance, data = dat, tau = 0.5)
exp(model$coefficients); extractAIC(model) # higher AIC but better fitting model? 

model <- rq(Duration ~ poly(Distance, 2), data = dat, tau = 0.5)
exp(model$coefficients); extractAIC(model) # higher AIC but better fitting model? 
saveRDS(model, file = paste0(HPCpath, "dist_duration_quant_model.rds"))
saveRDS(model, file = paste0(HPCpath, "traveltime_duration_quant_model.rds"))


# plot by model fit
ggplot(data = dat,
       aes(x = Distance, y = Duration)) +
  geom_point(alpha = 0.1,  size = 0.5, color = "darkgrey") +
  geom_smooth(formula = y ~ x, method = "glm", aes(color = "linear")) +
  geom_smooth(formula = y ~ splines::ns(x, 2), method = "glm.nb", aes(color = "spline 2k")) +
  geom_smooth(formula = y ~ splines::ns(x, 6), method = "glm.nb", aes(color = "spline 6k")) +
  geom_smooth(method = "glm.nb", aes(color = "glm.nb")) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_color_discrete(type = c("black", "#00BFC4", "#7CAE00", "#F8766D")) +
  labs(x = "Distance (km)",
       y = "Duration (days)",
       color = "Fit") +
  theme_classic()

ggsave("./03_output/Marshall_etal_trip_duration_fit_distance.pdf", width = 6, height = 4)


# < fit trip duration travel time ----
# import travel time matrices from MISC_refitting_gravity_model.R

tt_Burkina_Faso <- read_csv(paste0(HPCpath, "./Marshall et al. 2018 code/traveltime_Burkina_Faso.csv"), col_names = F) 
tt_Mali <- read_csv(paste0(HPCpath, "./Marshall et al. 2018 code/traveltime_Mali.csv"), col_names = F) 
tt_Tanzania <- read_csv(paste0(HPCpath, "./Marshall et al. 2018 code/traveltime_Tanzania.csv"), col_names = F) 
tt_Zambia <- read_csv(paste0(HPCpath, "./Marshall et al. 2018 code/traveltime_Zambia.csv"), col_names = F) 

# link to trip data
dat <- dat |> mutate(travel_time = 0)

for (i in 1:nrow(dat)) {
  
  if (dat$OrigCountry[i] == "Burkina Faso") {tt <- tt_Burkina_Faso}
  if (dat$OrigCountry[i] == "Mali") {tt <- tt_Mali}
  if (dat$OrigCountry[i] == "Tanzania") {tt <- tt_Tanzania}
  if (dat$OrigCountry[i] == "Zambia") {tt <- tt_Zambia}
  
  orig <- dat[[i, "OrigIndex"]]
  dest <- dat[[i, "DestIndex"]]
  
  dat$travel_time[i] <- tt[[orig, dest]]
  
}

dat <- dat |> 
  mutate(travel_time = case_when(travel_time == 0 ~ Distance,
                                 TRUE ~ travel_time))


# linear
model <- glm.nb(Duration ~ travel_time, data = dat)
exp(model$coefficients); model$aic

# cubic splines 2 knot
model <- glm.nb(Duration ~ splines::ns(travel_time, 2), data = dat)
exp(model$coefficients); model$aic

# cubic splines 6 knot
model <- glm.nb(Duration ~ splines::ns(travel_time, 6), data = dat)
exp(model$coefficients); model$aic # best AIC
saveRDS(model, file = "./03_output/traveltime_duration_sk6_model.rds")

# quantile regression (median) - robust to outliers
model <- rq(Duration ~ travel_time, data = dat, tau = 0.5)
exp(model$coefficients); extractAIC(model) # higher AIC but better fitting model? 

model <- rq(Duration ~ poly(travel_time, 2), data = dat, tau = 0.5)
exp(model$coefficients); extractAIC(model) # higher AIC but better fitting model? 
saveRDS(model, file = "./03_output/traveltime_duration_quant_model.rds")
saveRDS(model, file = paste0(HPCpath, "traveltime_duration_quant_model.rds"))



# plot by model fit
ggplot(data = dat,
       aes(x = travel_time, y = Duration)) +
  geom_point(alpha = 0.1,  size = 0.5, color = "darkgrey") +
  geom_smooth(formula = y ~ x, method = "glm", aes(color = "linear")) +
  # geom_smooth(formula = y ~ splines::ns(x, 2), method = "glm.nb", aes(color = "spline 2k")) +
  # geom_smooth(formula = y ~ splines::ns(x, 6), method = "glm.nb", aes(color = "spline 6k")) +
  geom_smooth(method = "glm.nb", aes(color = "n binomial")) +
  # geom_quantile(quantiles = 0.5, aes(color = OrigCountry, group = OrigCountry)) +
  geom_quantile(quantiles = 0.5, formula = y ~ poly(x, 2), aes(color = "quantile ^ 2")) +
  geom_quantile(quantiles = 0.5, aes(color = "quantile")) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_color_discrete(type = c( "black", "#00BFC4", "#7CAE00", "#F8766D")) +
  labs(x = "Travel time (min)",
       y = "Duration (days)",
       color = "Fit") +
  theme_classic()

# plot just quantiles
ggplot(data = dat,
       aes(x = travel_time / 60, y = Duration)) +
  geom_point(alpha = 0.1,  size = 0.5, color = "darkgrey") +
  geom_quantile(quantiles = c(0.25, 0.75), formula = y ~ poly(x, 2), color = "red") +
  geom_quantile(quantiles = c(0.50), formula = y ~ poly(x, 2), color = "blue", lty = 2) +
  coord_cartesian(ylim = c(0, 100), xlim = c(0, 15)) +
  labs(x = "Travel time (hours)",
       y = "Duration (days)",
       color = "Fit") +
  theme_classic()

ggsave("./03_output/Marshall_etal_trip_duration_fit_travel_time.pdf", width = 5, height = 3)

summary(dat$Duration); summary(dat$travel_time / 60)


# plot trip duration by country
ggplot(data = dat,
       aes(x = Distance, y = Duration)) +
  geom_point(alpha = 0.1,  size = 0.5, color = "darkgrey") +
  geom_quantile(quantiles = 0.5, aes(color = OrigCountry, group = OrigCountry)) +
  geom_quantile(quantiles = 0.5, aes(fill = "Overall fit"), color = "black") +
  scale_fill_manual(values = c(NA)) +
  labs(x = "Distance (km)",
       y = "Duration (days)",
       color = "",
       fill = "") +
  coord_cartesian(ylim = c(0, 365)) +
  theme_classic() +
  theme(text = element_text(size = 14))



  

