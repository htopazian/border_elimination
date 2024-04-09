# Figure S1 and S2. Population size and travel time distributions

# libraries
source("./02_code/packages_paths.R")

# pull in population estimates from B_population.R
population <- readRDS("./03_output/population.rds") 

# population -------------------------------------------------------------------
ggplot(data = population) + 
  geom_histogram(aes(x = pop), bins = 100, color = "#376795", 
                 fill = "#528fad", linewidth = 0.6) +
  labs(x = "population", y = "frequency") + 
  scale_x_continuous(trans = "log10", labels = scales::comma) + 
  theme_classic()

ggsave(filename = "./03_output/pop_dist.pdf", width = 5, height = 3)


# travel time ------------------------------------------------------------------
travel_time <- readRDS("./03_output/traveltimes.rds")

ggplot(data = travel_time) + 
  geom_histogram(aes(x = value / 60), bins = 100, color = "#e76254", 
                 fill = "#ef8a47", linewidth = 0.6) +
  labs(x = "travel time (hours)", y = "frequency") + 
  scale_x_continuous(trans = "log10", labels = scales::comma) + 
  theme_classic()

ggsave(filename = "./03_output/travel_time_dist.pdf", width = 5, height = 3)

