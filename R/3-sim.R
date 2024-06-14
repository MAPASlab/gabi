## library(gen3sis)
devtools::load_all("../gen3sis_rf")

sim1 <- run_simulation(
  config = "./config/gabi_config.R",
  landscape = "./landscapes/america",
  output_directory = "./out",
  verbose = 2
)
