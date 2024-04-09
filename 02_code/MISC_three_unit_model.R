# Testing the metapopulation model on a three-unit example

# libraries
source("./02_code/packages_paths.R")

# HPC set-up -------------------------------------------------------------------
setwd(HPCpath)

# to edit HPC username and password below
# usethis::edit_r_environ()

src <- conan::conan_sources(c("github::mrc-ide/malariasimulation@feat/metapop_asym"))

ctx <- context::context_save(path = paste0(HPCpath, "contexts"),
                             sources = c(paste0(HPCpath, "function_run_metapopulation_mix_examples.R")),
                             packages = c("dplyr", "malariasimulation", "data.table"),
                             package_sources = src)

share <- didehpc::path_mapping("Malaria", "M:", "//wpia-hn/Malaria", "M:")

config <- didehpc::didehpc_config(credentials = list(
            username = Sys.getenv("DIDE_USERNAME"),
            password = Sys.getenv("DIDE_PASSWORD")),
            shares = share,
            use_rrq = FALSE,
            cores = 1,
            cluster = "wpia-hn", # "fi--dideclusthn", "fi--didemrchnb", "wpia-hn"
            # "GeneralNodes", "12Core", "16Core", 
            # "12and16Core", "20Core", "24Core", "32Core"
            # "AllNodes" for wpia-hn
            template = "AllNodes", 
            parallel = FALSE)

obj <- didehpc::queue_didehpc(ctx, config = config)

# Set up your job --------------------------------------------------------------
index <- c("isolated", "semi-mixed", "perfectly mixed")

# Run tasks --------------------------------------------------------------------
t <- obj$enqueue_bulk(index, mixing_example)
t$status()


# Process ----------------------------------------------------------------------
d1 <- readRDS(paste0(HPCpath, "3unit_mixing_example_isolated.rds"))
d2 <- readRDS(paste0(HPCpath, "3unit_mixing_example_semi-mixed.rds"))
d3 <- readRDS(paste0(HPCpath, "3unit_mixing_example_perfectly mixed.rds"))

dat <- rbind(d1, d2, d3) |>
  mutate(model = factor(model, levels = c("isolated", "semi-mixed", "perfectly mixed")))

# plot
ggplot(data = dat) +
  geom_vline(xintercept = 15 * 365, lty = 2) + 
  geom_line(aes(x = timestep, y = prev2to10, color = factor(EIR)), linewidth = 0.7) +
  facet_wrap(~ model) +
  scale_color_manual(values = met.brewer("Hiroshige", 3),
                     labels = c("low", "medium", "high")) +
  scale_y_continuous(limits = c(0, 0.7)) + 
  # scale_x_continuous(breaks = seq(14 * 365, 18 * 365, 365),
  #                    limits = c(14 * 365, 18 * 365),
  #                    labels = c(-1, 0, 1, 2, 3)) +
  scale_x_continuous(breaks = seq(14 * 365, 38 * 365, 365),
                     limits = c(14 * 365, 38 * 365),
                     labels = seq(-1, 23, 1)) +
  labs(x = "year",
       y = expression(italic(Pf) ~ PR[2-10]),
       color = "transmission") +
  theme_classic()


# ggsave(paste0(path, "03_output/figureX_mixing.pdf"), height = 3, width = 5)
ggsave(paste0(path, "03_output/figureX_mixing.pdf"), height = 3, width = 9)
