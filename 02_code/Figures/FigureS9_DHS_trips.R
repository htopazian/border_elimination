# Figure S9. P(travel) and number of overnight trips from the DHS

# libraries
source("./02_code/packages_paths.R")

# pull in probability of staying values from MISC_travel_prob.R
admin1s_rematch <- readRDS("03_output/admin1s_travel.rds")


# plot P(travel)
pA <- ggplot() +
  geom_sf(data = admin1s_rematch, aes(fill = (1 - p_stay))) + 
  scale_fill_gradientn(colors = met.brewer("Hiroshige"),
                       values = c(0, admin1s_rematch$p_stay_SSA[1], 1),
                       breaks = seq(0, 1, .2)) +
  labs(fill = "P(travel)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank())

# plot number of overnight trips in the last 12 months
pB <- ggplot() +
  geom_sf(data = admin1s_rematch, aes(fill = med_trips)) + 
  scale_fill_gradientn(colors = met.brewer("Hokusai3"),
                       values = c(0, admin1s_rematch$med_trips_SSA[1], 1),
                       breaks = seq(1, 18, 4)) +
  labs(fill = "Number of \novernight trips") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank())

p_all <- pA + pB + plot_annotation(tag_levels = "A")

ggsave("./03_output/DHS_trips.pdf", plot = p_all, width = 8, height = 4)

