# script modified from Lucy Okell, September 2022
# calculating the RDT or PCR prevalence for a given microscopy prevalence. All age, all geographies averages.

library(tidyverse)


logit <- function(p) {
  log(p / (1 - p))
}

invLogit <- function(logOdds) {
  exp(logOdds) / (1 + exp(logOdds))
}

### Relationship between slide and PCR prevalence, 10.1038/NCOMMS2241
slidePrev <- seq(0.001, 0.9, 0.001) # no data above 90%
logOddsSlidePrev <- logit(slidePrev)

logOddsPcrPrev <- 0.954 + 0.868 * logOddsSlidePrev
pcrPrev <- invLogit(logOddsPcrPrev)


### Relationship between slide and RDT prevalence, 10.1038/NATURE16039
logOddsRdtPrev <- 0.108 + 0.907 * logOddsSlidePrev
rdtPrev <- invLogit(logOddsRdtPrev)


### Relationship between PCR and RDT prevalence, 10.1038/NATURE16039
logOddsRdt2Prev <- -0.968 + 1.186 * logOddsPcrPrev
rdt2Prev <- invLogit(logOddsRdt2Prev)


### Combine data
alldat <- tibble(slide = slidePrev, RDT_from_microscopy = rdtPrev, RDT_from_PCR = rdt2Prev, PCR = pcrPrev)

alldat_long <- pivot_longer(data = alldat, cols = 1:3,
                            names_to = "type", values_to = "PR")

saveRDS(alldat, "./03_output/RDT_PCR.rds")

# plot all three curves: RDT as calculated through microscopy prev, RDT:PCR and microscopy:PCR
ggplot(data = alldat_long) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  geom_line(aes(x = PCR, y = PR, color = type), size = 1) +
  labs(x = "PCR prevalence", y = "RDT or microscopy prevalence", color = "") +
  scale_color_discrete(labels = c("RDT (through microscopy)", "RDT", "microscopy")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_classic()

# plot curves from PCR
ggplot(data = alldat_long |> filter(type != "RDT_from_microscopy")) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  geom_line(aes(x = PCR, y = PR, color = type), size = 1) +
  labs(x = "PCR prevalence", y = "RDT or microscopy prevalence", color = "") +
  scale_color_discrete(labels = c("RDT", "microscopy")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_classic()
