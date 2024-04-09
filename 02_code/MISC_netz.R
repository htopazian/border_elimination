# Determine input ITN usage to reach target average population usage values

# libraries
source("./02_code/packages_paths.R")

# run through various dist values to match target ITN values

# create function
net_convert <- function(dist){
  
  # distribute nets every 1 years for 30 years
  fit <- netz::population_usage(distribution = rep(dist, 30),
                                timesteps = 30 * 365,
                                distribution_timesteps = seq(1, 30, 1) * 365,
                                half_life = 5 * 365)
  
  # print dist value and last 15-year mean
  tibble(usage = mean(tail(fit, 15 * 365)), 
         dist = dist)
  
}

# run through distribution values
dist <- seq(0, 1, 0.001)
net_conversion <- map_dfr(dist, net_convert)

# plot
ggplot() + 
  geom_line(data = net_conversion, aes(x = usage, y = dist)) + 
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  labs(x = "target usage", y = "annual distribution") +
  theme_classic() 

# save file for use in parameterization of site files
saveRDS(net_conversion, "./03_output/net_conversion.rds")
