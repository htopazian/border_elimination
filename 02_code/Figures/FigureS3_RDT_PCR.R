# Figure S3. RDT vs. PCR curves

# libraries
source("./02_code/packages_paths.R")

# pull in data from MISC_microscopy_RDT_PCR.R
dat <- readRDS("./03_output/RDT_PCR.rds") 

dat_long <- dat |> 
  pivot_longer(cols = 1:3,
               names_to = "type", values_to = "PR")

# plot curves from PCR
ggplot(data = dat_long |> filter(type != "RDT_from_microscopy")) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  geom_line(aes(x = PCR, y = PR, color = type), size = 1) +
  labs(x = "PCR prevalence", y = "RDT or microscopy prevalence", color = "") +
  scale_color_manual(labels = c("RDT", "microscopy"), 
                     values = met.brewer(name = "Hiroshige", 2)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_classic()

ggsave("./03_output/PCR_RDT.pdf", width = 4, height = 3)
