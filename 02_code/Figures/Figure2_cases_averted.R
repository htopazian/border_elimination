# Figure 2. Cases averted and % change in PfPR; SSA maps

# libraries
source("./02_code/packages_paths.R")

# pull in shapefiles from A_shapefiles_centroids.R
countries_africa <- readRDS("./03_output/countries_africa.rds")
countries <- readRDS("./03_output/countries.rds")
admin1s <- readRDS("./03_output/admin1s.rds")
border_admins <- readRDS("./03_output/border_admins.rds")
centroids_all <- readRDS("./03_output/centroids_all.rds")
centroids_border <- readRDS("./03_output/centroids_border.rds")
borders <- readRDS("./03_output/borders.rds")
borders_inner <- readRDS("03_output/borders_inner.rds") 

# read in cases averted dataset from G_processing.R
avert <- readRDS("./03_output/avert_output.rds")


# map: absolute cases averted in each seed point with border intervention ------
# attach seed point indices to each admin unit
seed_point <- centroids_border |> st_drop_geometry() |>
  mutate(seed = row_number())

output <- avert |> 
  left_join(seed_point) |>
  filter(seed_point == seed) |>
  group_by(ID_0, NAME_1, pop, seed_point, draw, mix_type) |>
  summarize(pfpr_2_10 = mean(pfpr_2_10, na.rm = T),
            case_total = sum(case_total, na.rm = T),
            severe_total = sum(severe_total, na.rm = T),
            pfpr_2_10_baseline = mean(pfpr_2_10_baseline, na.rm = T),
            case_total_baseline = sum(case_total_baseline, na.rm = T),
            severe_total_baseline = sum(severe_total_baseline, na.rm = T)) |>
  rowwise() |>
  mutate(case_avert = case_total_baseline - case_total,
         case_avert_p = case_avert / case_total_baseline * 100,
         cinc_baseline = (case_total_baseline / pop),
         cinc_diff = (case_total_baseline / pop) - (case_total / pop),
         pfpr_diff = pfpr_2_10_baseline - pfpr_2_10,
         pfpr_diff_p = (pfpr_2_10_baseline - pfpr_2_10) / pfpr_2_10_baseline * 100) |>
  # find median over draws
  ungroup() |> group_by(ID_0, NAME_1, pop, seed_point, mix_type) |>
  summarize(n = n(),
            case_avert = median(case_avert),
            case_avert_p = median(case_avert_p),
            cinc_baseline = median(cinc_baseline),
            cinc_diff = median(cinc_diff),
            pfpr_diff = median(pfpr_diff),
            pfpr_diff_p = median(pfpr_diff_p))

# visualize
output_sf <- output |> 
  left_join(border_admins |> dplyr::select(ID_0, NAME_1)) |>
  st_as_sf()

summary(output_sf$case_avert); summary(output_sf$case_avert_p); summary(output_sf$cinc_diff); summary(output_sf$pfpr_diff_p)

# set those with cinc_diff or pfpr_diff < 0 to 0 due to stochasticity
output_sf <- output_sf |> 
  mutate(case_avert = ifelse(case_avert < 0, 0, case_avert),
         case_avert_p = ifelse(case_avert_p < 0, 0, case_avert_p),
         case_avert_p = ifelse(is.na(case_avert_p), 0, case_avert_p),
         pfpr_diff_p = ifelse(pfpr_diff_p < 0, 0, pfpr_diff_p),
         pfpr_diff_p = ifelse(is.na(pfpr_diff_p), 0, pfpr_diff_p))

summary(output_sf$case_avert); summary(output_sf$case_avert_p); summary(output_sf$pfpr_diff_p)
range <- diff(range(output_sf[output_sf$case_avert > -Inf,]$case_avert, na.rm = T)) 
max <- max(output_sf$case_avert, na.rm = T)

p1 <- ggplot() + 
  geom_sf(data = countries_africa) + 
  geom_sf(data = output_sf, aes(fill = case_avert)) + 
  geom_sf(data = borders_inner, color = "navy", linewidth = 0.5) + 
  # values between 0 and 1 with the middle number as the midpoint
  # scale_fill_gradientn(colors = met.brewer("Hiroshige"), 
  #                      values = c(0, (range - max) / range + 0.02, 1),
  #                      labels = scales::comma) + 
  scale_fill_stepsn(colors = met.brewer("Hiroshige", n = 8)[c(4, 5, 6, 7)],
                    breaks = c(0, 1000, 10000, 100000), 
                    values = scales::rescale(c(0, 1000, 10000, 100000)),
                    labels = scales::comma,
                    limits = c(0, max)) + 
  labs(fill = "absolute cases \naverted") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank())

ggsave(plot = p1, filename = "./03_output/case_avert_absolute.pdf", 
       width = 5, height = 3); beepr::beep(1)


# map: % cases averted in each seed point with border intervention -------------
summary(output_sf$case_avert_p)
range <- diff(range(output_sf[output_sf$case_avert_p > -Inf,]$case_avert_p, na.rm = T)) 
max <- max(output_sf$case_avert_p, na.rm = T)

p2 <- ggplot() + 
  geom_sf(data = countries_africa) + 
  geom_sf(data = output_sf, aes(fill = case_avert_p)) + 
  geom_sf(data = borders_inner, color = "navy", linewidth = 0.5) + 
  # values between 0 and 1 with the middle number as the midpoint
  # scale_fill_gradientn(colors = met.brewer("Hiroshige"), 
  #                      values = c(0, (range - max) / range + 0.02, 1),
  #                      labels = scales::comma) + 
  scale_fill_stepsn(colors = met.brewer("Hiroshige", n = 8)[c(4, 5, 6, 7)],
                    breaks = c(0, 0.5, 3, 6), 
                    values = scales::rescale(c(0.2, 0.4, 0.6, 0.8)),
                    limits = c(0, max)) + 
  labs(fill = "% cases \naverted") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank())

ggsave(plot = p2, filename = "./03_output/case_avert_percent.pdf", 
       width = 5, height = 3); beepr::beep(1)


# map: relative change in pfpr in each seed point with border intervention -----
percentiles <- output_sf |> st_drop_geometry() |> 
  mutate(q = case_when(case_avert >= 0 & case_avert < 1000 ~ 1,
                          case_avert >= 1000 & case_avert < 10000 ~ 2,
                          case_avert >= 10000 & case_avert < 100000 ~ 3,
                          case_avert >= 100000 ~ 4)) |>
  group_by(q) |>
  summarize(n = n()) |>
  mutate(t = sum(n),
         p = n / t,
         p2 = cumsum(p))

# set colors below, centering on light yellow as the midpoint
summary(output_sf$pfpr_diff_p)
range <- diff(range(output_sf[output_sf$pfpr_diff_p > -Inf,]$pfpr_diff_p, na.rm = T)) 
max <- max(output_sf$pfpr_diff_p, na.rm = T)

p3 <- ggplot() + 
  geom_sf(data = countries_africa) + 
  geom_sf(data = output_sf, aes(fill = pfpr_diff_p)) + 
  geom_sf(data = borders_inner, color = "navy", linewidth = 0.5) + 
  # values between 0 and 1 with the middle number as the midpoint
  # scale_fill_gradientn(colors = met.brewer("Hiroshige"), 
  #                      values = c(0, (range - max) / range + 0.02, 1)) + 
  scale_fill_stepsn(colors = met.brewer("Hiroshige", n = 8)[c(4, 5, 6, 7)],
                    breaks = c(0, 0.5, 1, 3), 
                    values = scales::rescale(c(0.2, 0.4, 0.6, 0.8)),
                    limits = c(0, max)) + 
  labs(fill = expression(atop("% reduction in ", italic(Pf)~PR[2-10]~            phantom(100)))) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank()) 

ggsave(plot = p3, filename = "./03_output/pfpr_diff_relative.pdf", 
       width = 5, height = 3); beepr::beep(1)


# save plots together ----------------------------------------------------------

p4 <- p1 + p2 + p3 + 
  plot_layout(guides = "keep") + plot_annotation(tag_levels = "A")

ggsave(plot = p4, filename = "./03_output/outcome_avert_ABC.pdf", 
       width = 13, height = 3); beepr::beep(1)



p5 <- p1 + p3 + 
  plot_layout(guides = "keep") + plot_annotation(tag_levels = "A")

ggsave(plot = p5, filename = "./03_output/outcome_avert_AB.pdf", 
       width = 9, height = 3); beepr::beep(1)



# zoom-ins ---------------------------------------------------------
output_zoom <- st_drop_geometry(output_sf) |> ungroup()

output_zoom |> arrange(-pfpr_diff_p) |> dplyr::select(ID_0, NAME_1, pfpr_diff_p)

# southern africa
p2 + coord_sf(xlim = c(17, 34), ylim = c(-19, -35.5), expand = FALSE) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank())

# senegal
p2 + coord_sf(xlim = c(-18, -10), ylim = c(11, 18), expand = FALSE) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank())

