# Figure 4. PfPR and border coverage case studies results

# libraries
source("./02_code/packages_paths.R")


# PfPR case study --------------------------------------------------------------
output <- readRDS("./03_output/PfPR_casestudy.rds") |> ungroup() |>
  filter(pfpr %in% seq(0.1, 0.8, 0.1)) |>
  filter(diagnostic == "RDT") |>
  mutate(sens = detect / infectious_prev)

# tile plot
# https://r-charts.com/correlation/heat-map-ggplot2/
tilep1 <- output |>
  separate(ID, c("pfpr1", "pfpr2"), "_") |>
  filter(pfpr1 == pfpr) |>
  filter(pfpr1 %in% seq(0.1, 0.8, 0.1) & pfpr2 %in% seq(0.1, 0.8, 0.1)) |>
  dplyr::select(pfpr1, pfpr2, draw, case_avert_p, detect, infectious_prev, sens)

tilep2 <- output |>
  separate(ID, c("pfpr1", "pfpr2"), "_") |>
  filter(pfpr2 == pfpr) |>
  filter(pfpr1 %in% seq(0.1, 0.8, 0.1) & pfpr2 %in% seq(0.1, 0.8, 0.1)) |>
  rename(pfpr1 = pfpr2, pfpr2 = pfpr1) |>
  dplyr::select(pfpr1, pfpr2, draw, case_avert_p, detect, infectious_prev, sens)

timepall <- bind_rows(tilep1, tilep2) |> distinct() |>
  group_by(pfpr1, pfpr2) |>
  summarize(case_avert_p = median(case_avert_p),
            detect = median(detect),
            infectious_prev = median(infectious_prev),
            sens = median(sens))

summary(timepall$case_avert_p); summary(timepall$detect); summary(timepall$infectious_prev); summary(timepall$sens)

# plot cases averted
p1 <- ggplot(timepall, aes(x = factor(pfpr1), y = factor(pfpr2), fill = case_avert_p)) +
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  geom_text(aes(label = sprintf("%.0f", round(case_avert_p, 1))), color = "white", size = 2) + # or %0.1f; size = 1.5
  # values between 0 and 1 with the middle number as the midpoint
  scale_fill_gradientn(colors = met.brewer("Hiroshige"), 
                       values = c(-1, -0.08, 1),
                       breaks = seq(0, 14, 2)) + 
  scale_y_discrete(breaks = seq(0.1, 0.8, 0.1)) +
  scale_x_discrete(breaks = seq(0.1, 0.8, 0.1)) +
  labs(x = expression(italic(Pf)~PR[2-10]~" unit A"), 
       y = expression(italic(Pf)~PR[2-10]~" unit B"), 
       fill = "% of cases averted \nin unit A",
       title = expression("Varying "~italic(Pf)~PR[2-10])) + 
  theme_classic()

ggsave(plot = p1, filename = "./03_output/pfpr_casestudy_per.pdf", 
       width = 7, height = 5)

# plot % infections detected
p2 <- ggplot(timepall, aes(x = factor(pfpr1), y = factor(pfpr2), fill = sens * 100)) +
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  geom_text(aes(label = sprintf("%.0f", round(sens * 100, 1))), color = "white", size = 2) + # or %0.1f; size = 1.5
  # values between 0 and 1 with the middle number as the midpoint
  scale_fill_gradientn(colors = met.brewer("Hiroshige"), 
                       values = c(-1, -0.08, 1)) + 
  scale_y_discrete(breaks = seq(0.1, 0.8, 0.1)) +
  scale_x_discrete(breaks = seq(0.1, 0.8, 0.1)) +
  labs(x = expression(italic(Pf)~PR[2-10]~" unit A"), 
       y = expression(italic(Pf)~PR[2-10]~" unit B"), 
       fill = "% of infections \ndetected by \nRDT in unit A",
       title = expression("Varying "~italic(Pf)~PR[2-10])) + 
  theme_classic()

ggsave(plot = p2, filename = "./03_output/pfpr_casestudy_detect.pdf", 
       width = 7, height = 5)


# negative cases averted values are due to stochasticity; with few draws, the 
# median value will be negative, with lots of draws, it will end up slightly 
# above 0. This is due to low RDT detection (almost 0) in low PfPR settings.


# Coverage case study --------------------------------------------------------------
output <- readRDS("./03_output/coverage_casestudy.rds") |> ungroup() |>
  mutate(sens = detect / infectious_prev)

# tile plot
# https://r-charts.com/correlation/heat-map-ggplot2/
tilep1 <- output |>
  separate(ID, c("pfpr1", "pfpr2"), "_") |>
  filter(coverage %in% seq(0, 1, 0.2)) |>
  filter(pfpr1 == pfpr) |>
  dplyr::select(pfpr1, pfpr2, coverage, draw, case_avert_p, detect, infectious_prev, sens)

timepall <- tilep1 |> distinct() |> 
  group_by(pfpr1, pfpr2, coverage) |>
  summarize(case_avert_p = median(case_avert_p),
            detect = median(detect),
            infectious_prev = median(infectious_prev),
            sens = median(sens))
  
summary(timepall$case_avert_p); summary(timepall$detect); summary(timepall$infectious_prev); summary(timepall$sens)

table(timepall$pfpr1); table(timepall$pfpr2)

# look at % of infections detected
timepall |> group_by(pfpr2) |>
  summarize(median = median(sens), 
            min = min(sens),
            max = max(sens))

# plot cases averted
# round negative value to 0 for plot
timepall2 <- timepall |> mutate(case_avert_p = ifelse(case_avert_p < 0, 0, case_avert_p))

p3 <- ggplot(timepall2, aes(x = factor(pfpr2), y = coverage, fill = case_avert_p)) +
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  geom_text(aes(label = sprintf("%.0f", round(case_avert_p, 1))), color = "white", size = 2) + # or %0.1f; size = 1.5
  # values between 0 and 1 with the middle number as the midpoint
  scale_fill_gradientn(colors = met.brewer("Hiroshige"), 
                       values = c(-1, -0.08, 1)) + 
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  labs(x = expression(italic(Pf)~PR[2-10]~" unit B"), 
       y = "border post coverage", 
       fill = "% of cases averted \nin unit A",
       title = "Varying border post coverage") + 
  theme_classic()

ggsave(plot = p3, filename = "./03_output/coverage_casestudy_per.pdf", 
       width = 4, height = 5)

# plot cases detected
ggplot(timepall, aes(x = factor(pfpr2), y = coverage, fill = sens * 100)) +
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  geom_text(aes(label = sprintf("%.1f", round(sens, 1))), color = "white", size = 1.5) +
  # values between 0 and 1 with the middle number as the midpoint
  scale_fill_gradientn(colors = met.brewer("Hiroshige"), 
                       values = c(-1, -0.08, 1)) + 
  labs(x = expression(italic(Pf)~PR[2-10]~" unit 2"), 
       y = "border post coverage", 
       fill = "% of infections detected \nby RDT in unit 1") + 
  theme_classic()

ggsave(filename = "./03_output/coverage_casestudy_detect.pdf", 
       width = 4, height = 5)


# negative cases averted values are due to stochasticity; with few draws, the 
# median value will be negative, with lots of draws, it will end up slightly 
# above 0. This is due to low RDT detection (almost 0) in low PfPR settings.



# Cross-border mixing case study -----------------------------------------------
output <- readRDS("./03_output/mixing_casestudy.rds") |> ungroup() |>
  mutate(sens = detect / infectious_prev)

# tile plot
# https://r-charts.com/correlation/heat-map-ggplot2/
tilep1 <- output |>
  separate(ID, c("pfpr1", "pfpr2"), "_") |>
  filter(pfpr1 == pfpr) |>
  dplyr::select(pfpr1, pfpr2, mixing_i, draw, case_avert_p, detect, infectious_prev, sens)

tilep2 <- output |>
  separate(ID, c("pfpr1", "pfpr2"), "_") |>
  filter(pfpr2 == pfpr) |>
  rename(pfpr1 = pfpr2, pfpr2 = pfpr1) |>
  dplyr::select(pfpr1, pfpr2, mixing_i, draw, case_avert_p, detect, infectious_prev, sens)


timepall <- bind_rows(tilep1, tilep2) |> distinct() |>
  group_by(pfpr1, pfpr2, mixing_i) |>
  summarize(case_avert_p = median(case_avert_p),
            detect = median(detect),
            infectious_prev = median(infectious_prev),
            sens = median(sens))

summary(timepall$case_avert_p); summary(timepall$detect); summary(timepall$infectious_prev); summary(timepall$sens)

table(timepall$pfpr1); table(timepall$pfpr2)

# look at % of infections detected
timepall |> group_by(pfpr2) |>
  summarize(median = median(sens), 
            min = min(sens),
            max = max(sens))

summary(timepall$case_avert_p)

# plot cases averted
# round negative value to 0 for plot
timepall2 <- timepall |> mutate(case_avert_p = ifelse(case_avert_p < 0, 0, case_avert_p))

p4 <- ggplot(timepall, aes(x = factor(pfpr1), y = factor(pfpr2), fill = case_avert_p)) +
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  geom_text(aes(label = sprintf("%.0f", round(case_avert_p, 1))), color = "white", size = 2) + # or %0.1f; size = 1.5
  # values between 0 and 1 with the middle number as the midpoint
  scale_fill_gradientn(colors = met.brewer("Hiroshige"), 
                       values = c(-1, -0.08, 1),
                       breaks = seq(0, 8, 2)) + 
  scale_y_discrete(breaks = seq(0.2, 0.8, 0.2)) +
  scale_x_discrete(breaks = seq(0.2, 0.8, 0.2)) +
  facet_grid(~ mixing_i, labeller = as_labeller(c(
    `0.0066` = "25th percentile",
    `0.0127` = "median",
    `0.0223` = "75th percentile"
  ))) + 
  labs(x = expression(italic(Pf)~PR[2-10]~" unit A"), 
       y = expression(italic(Pf)~PR[2-10]~" unit B"), 
       fill = "% of cases averted \nin unit A",
       title = "Varying cross-border mixing") + 
  theme_classic()

ggsave(plot = p4, filename = "./03_output/mixing_casestudy_per.pdf", 
       width = 7, height = 5)


# negative cases averted values are due to stochasticity; with few draws, the 
# median value will be negative, with lots of draws, it will end up slightly 
# above 0. This is due to low RDT detection (almost 0) in low PfPR settings.


# combine all case studies -----------------------------------------------------
# combine plots
design <- "
1122
3444
"
plot_total <- (p1 + p2 + plot_layout(widths = c(3,3))) /
  (p3 + p4 + plot_layout(widths = c(2,4))) + 
  plot_annotation(tag_levels = "A")

plot_total

# save
ggsave(plot = plot_total, 
       filename = "./03_output/all_case_studies.pdf", 
       width = 10, height = 5)

# alternative layout
# layout <- "
# AAABBBC
# AAABBBC
# AAABBBC
# "
# 
# plot_total <- p1 + p2 + p3 + 
#   plot_layout(design = layout) + 
#   plot_annotation(tag_levels = "A")
# 
# # save
# ggsave(plot = plot_total, 
#        filename = "./03_output/pfpr_coverage_studies.pdf", 
#        width = 15, height = 4)
