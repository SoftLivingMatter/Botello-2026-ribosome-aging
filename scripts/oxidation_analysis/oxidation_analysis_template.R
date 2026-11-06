source("oxidation_functions.R")

cfg <- list(
  input_file    = "/path/to/summary.tsv",

  # Named vector: name = column/label, value = substring to match in sample names.
  # Order defines group_idx (1, 2, 3, ...) used as the GLM predictor.
  # Monotonicity filtering only applies when length >= 3.
  groups = c(
    T0  = "T0",
    T24 = "T24",
    T72 = "T72"
  ),

  # NULL to omit region annotation, or a data.frame with columns: name, start, end
  region_table = NULL,
  # region_table = data.frame(
  #   name  = c("18S",  "5.8S", "28S"),
  #   start = c(3657,   6623,   7935),
  #   end   = c(5500,   6778,   12968)
  # ),

  output_prefix = "my_experiment"
)

results <- run_analysis(cfg)

# Optional: human rRNA region colors — only useful when region_table is set
# region_colors <- c(
#   "18S"   = "#2a78d6",
#   "5.8S"  = "#eb6834",
#   "28S"   = "#1baf7a",
#   "other" = "#898781"
# )

make_plots(
  dat           = results$dat,
  sig_sites     = results$sig_sites,
  groups        = cfg$groups,
  output_prefix = cfg$output_prefix,
  region_colors = NULL   # swap in region_colors if using region_table
)
