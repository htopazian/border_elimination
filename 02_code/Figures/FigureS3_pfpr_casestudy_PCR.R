# Figure S3. PfPR case study results with PCR detection

# libraries
source("./02_code/packages_paths.R")


# PfPR case study w/ PCR -------------------------------------------------------
output <- readRDS("./03_output/PfPR_casestudy.rds") |> ungroup() |>
  filter(diagnostic == "PCR")

# tile plot
# https://r-charts.com/correlation/heat-map-ggplot2/
tilep1 <- output |>
  separate(ID, c("pfpr1", "pfpr2"), "_") |>
  filter(pfpr1 == pfpr) |>
  select(pfpr1, pfpr2, draw, case_avert_p, detect, infectious_prev)

tilep2 <- output |>
  separate(ID, c("pfpr1", "pfpr2"), "_") |>
  filter(pfpr2 == pfpr) |>
  rename(pfpr1 = pfpr2, pfpr2 = pfpr1) |>
  select(pfpr1, pfpr2, draw, case_avert_p, detect, infectious_prev)

timepall <- bind_rows(tilep1, tilep2) |> distinct() |>
  group_by(pfpr1, pfpr2) |>
  summarize(case_avert_p = median(case_avert_p),
            detect = median(detect),
            infectious_prev = median(infectious_prev))

summary(timepall$case_avert_p); summary(timepall$detect); summary(timepall$infectious_prev)

# plot cases averted
p1 <- ggplot(timepall, aes(x = factor(pfpr1), y = factor(pfpr2), fill = case_avert_p)) +
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  geom_text(aes(label = sprintf("%.0f", round(case_avert_p, 1))), color = "white", size = 2) + # or %0.1f; size = 1.5
  # values between 0 and 1 with the middle number as the midpoint
  scale_fill_gradientn(colors = met.brewer("Hiroshige"), 
                       values = c(-1, -0.08, 1),
                       breaks = seq(0, 23, 3)) + 
  labs(x = expression(italic(Pf)~PR[2-10]~" unit 1"), 
       y = expression(italic(Pf)~PR[2-10]~" unit 2"), 
       fill = "% of cases averted \nin unit 1") + 
  theme_classic()

ggsave(plot = p1, filename = "./03_output/pfpr_casestudy_per_PCR.pdf", 
       width = 7, height = 5)
