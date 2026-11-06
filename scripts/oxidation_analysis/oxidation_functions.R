library(tidyverse)
library(broom)

# ── Data preparation ──────────────────────────────────────────────────────────

# groups: named character vector — names are group labels, values are substrings
# matched against the sample column. Samples matching exactly one pattern are
# assigned that group and a positional index; others are dropped.
assign_groups <- function(dat, groups) {
  patterns <- unname(groups)
  labels   <- names(groups)

  match_idx <- map_int(dat$sample, function(s) {
    m <- which(str_detect(s, fixed(patterns)))
    if (length(m) == 1L) m else NA_integer_
  })

  n_multi <- sum(map_int(dat$sample, ~ sum(str_detect(.x, fixed(patterns)))) > 1)
  if (n_multi > 0)
    warning(n_multi, " samples matched multiple group patterns and were excluded")

  dat |>
    mutate(
      group     = factor(labels[match_idx], levels = labels),
      group_idx = match_idx
    ) |>
    filter(!is.na(group))
}

filter_sites <- function(dat) {
  dat |>
    group_by(chr, pos) |>
    filter(min(depth) >= 50, sum(T_count) >= 5) |>
    ungroup()
}

# region_table: NULL to skip, or data.frame with columns: name, start, end
annotate_regions <- function(dat, region_table = NULL) {
  if (is.null(region_table)) return(dat)

  lvls <- c(region_table$name, "other")
  dat |>
    mutate(region = {
      r <- rep("other", n())
      for (i in seq_len(nrow(region_table)))
        r[pos >= region_table$start[i] & pos <= region_table$end[i]] <- region_table$name[i]
      factor(r, levels = lvls)
    })
}

# ── Analysis ──────────────────────────────────────────────────────────────────

# Monotonicity only applied when >= 3 groups; returns NA direction otherwise.
detect_monotone_sites <- function(dat) {
  if (n_distinct(dat$group_idx) < 3)
    return(distinct(dat, chr, pos) |> mutate(monotone_direction = NA_character_))

  dat |>
    group_by(chr, pos, group_idx) |>
    summarise(mean_freq = mean(G_to_T_freq), .groups = "drop") |>
    group_by(chr, pos) |>
    summarise(
      diffs         = list(diff(mean_freq[order(group_idx)])),
      monotone_up   = all(unlist(diffs) >= 0),
      monotone_down = all(unlist(diffs) <= 0),
      .groups = "drop"
    ) |>
    filter(monotone_up | monotone_down) |>
    transmute(
      chr, pos,
      monotone_direction = if_else(monotone_up, "increase", "decrease")
    )
}

.fit_trend <- function(df) {
  glm(cbind(T_count, depth - T_count) ~ group_idx, data = df, family = binomial)
}

# group_idx used as ordinal numeric predictor. p-values are identical to using
# real time values; only the slope scale changes (per group-step, not per hour).
run_glm_analysis <- function(dat) {
  nest_vars <- intersect(c("chr", "pos", "region"), names(dat))

  dat |>
    nest(data = -all_of(nest_vars)) |>
    mutate(
      model  = map(data, possibly(.fit_trend, otherwise = NULL)),
      result = map(model, ~ if (!is.null(.x)) tidy(.x) else NULL)
    ) |>
    unnest(result) |>
    filter(term == "group_idx") |>
    select(-data, -model)
}

# For >= 3 groups: inner-join monotone sites, require direction match, BH within direction.
# For 2 groups: monotone filter is trivially uninformative; BH across all sites.
get_sig_sites <- function(site_models, monotone_sites) {
  has_monotone <- any(!is.na(monotone_sites$monotone_direction))

  if (has_monotone) {
    site_models |>
      inner_join(monotone_sites, by = c("chr", "pos")) |>
      filter(
        (monotone_direction == "increase" & estimate > 0) |
        (monotone_direction == "decrease" & estimate < 0)
      ) |>
      group_by(monotone_direction) |>
      mutate(padj = p.adjust(p.value, method = "BH")) |>
      ungroup() |>
      filter(padj < 0.05)
  } else {
    site_models |>
      left_join(monotone_sites, by = c("chr", "pos")) |>
      mutate(padj = p.adjust(p.value, method = "BH")) |>
      filter(padj < 0.05)
  }
}

build_effect_summary <- function(dat) {
  dat |>
    group_by(chr, pos, group) |>
    summarise(mean_freq = mean(G_to_T_freq), .groups = "drop") |>
    pivot_wider(
      id_cols      = c(chr, pos),
      names_from   = group,
      values_from  = mean_freq,
      names_prefix = "freq_"
    )
}

build_replicate_wide <- function(dat) {
  dat |>
    select(chr, pos, sample, G_to_T_freq) |>
    pivot_wider(
      id_cols     = c(chr, pos),
      names_from  = sample,
      values_from = G_to_T_freq,
      names_prefix = "rep_"
    )
}

.fit_anova <- function(df) {
  m <- lm(G_to_T_freq ~ factor(group_idx), data = df)
  a <- anova(m)
  list(F_stat = a$`F value`[1], p_val = a$`Pr(>F)`[1])
}

.tukey_first_last <- function(df) {
  m        <- aov(G_to_T_freq ~ factor(group_idx), data = df)
  tuk      <- TukeyHSD(m)$`factor(group_idx)`
  contrast <- paste0(max(df$group_idx), "-", min(df$group_idx))
  idx      <- match(contrast, rownames(tuk))
  if (is.na(idx)) NA_real_ else tuk[idx, "p adj"]
}

run_anova_analysis <- function(dat) {
  nest_vars <- intersect(c("chr", "pos", "region"), names(dat))

  dat |>
    nest(data = -all_of(nest_vars)) |>
    mutate(
      ar               = map(data, .fit_anova),
      anova_F          = map_dbl(ar, "F_stat"),
      anova_p          = map_dbl(ar, "p_val"),
      tukey_first_last = map_dbl(data, .tukey_first_last)
    ) |>
    select(-data, -ar) |>
    mutate(anova_q = p.adjust(anova_p, method = "BH"))
}

# ── Orchestrator ──────────────────────────────────────────────────────────────

# cfg fields:
#   input_file    — path (or character vector of paths) to mpileup summary TSV(s)
#   groups        — named character vector of sample-name substrings
#   region_table  — NULL or data.frame(name, start, end)
#   output_prefix — prefix for output CSV and PDF files
run_analysis <- function(cfg) {
  col_spec <- cols(
    chr         = col_character(),
    pos         = col_integer(),
    ref         = col_character(),
    depth       = col_integer(),
    T_count     = col_integer(),
    G_to_T_freq = col_double(),
    sample      = col_character()
  )
  dat <- map(cfg$input_file, read_tsv, col_types = col_spec) |> bind_rows()

  dat <- dat |>
    mutate(G_to_T_freq = if_else(depth > 0, T_count / depth, 0)) |>
    assign_groups(cfg$groups) |>
    filter_sites() |>
    annotate_regions(cfg$region_table)

  monotone_sites <- detect_monotone_sites(dat)
  site_models    <- run_glm_analysis(dat)
  effect_summary <- build_effect_summary(dat)
  sig_sites      <- get_sig_sites(site_models, monotone_sites) |>
    left_join(effect_summary, by = c("chr", "pos"))
  anova_results  <- run_anova_analysis(dat)

  glm_all <- site_models |>
    mutate(glm_q = p.adjust(p.value, method = "BH")) |>
    left_join(monotone_sites, by = c("chr", "pos")) |>
    select(chr, pos, monotone_direction, glm_p = p.value, glm_q)

  base_cols <- intersect(c("chr", "pos", "region"), names(anova_results))
  freq_cols <- paste0("freq_", names(cfg$groups))
  rep_wide  <- build_replicate_wide(dat)

  out <- anova_results |>
    left_join(effect_summary, by = c("chr", "pos")) |>
    left_join(glm_all,        by = c("chr", "pos")) |>
    left_join(rep_wide,       by = c("chr", "pos")) |>
    select(all_of(base_cols), any_of(freq_cols),
           anova_F, anova_q, tukey_first_last,
           monotone_direction, glm_p, glm_q,
           starts_with("rep_")) |>
    arrange(chr, pos)

  fname <- paste0(cfg$output_prefix, "_sites.csv")
  write_csv(out, fname)
  message("Wrote ", fname, " (", nrow(out), " sites)")

  invisible(list(dat = dat, out = out, sig_sites = sig_sites))
}

# ── Visualization ─────────────────────────────────────────────────────────────

theme_ox <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_line(color = "#e1e0d9", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(color = "#c3c2b7", linewidth = 0.5),
      axis.ticks       = element_line(color = "#c3c2b7", linewidth = 0.4),
      text             = element_text(color = "#0b0b0b"),
      axis.text        = element_text(color = "#52514e"),
      strip.text       = element_text(color = "#0b0b0b", face = "bold"),
      legend.position  = "bottom",
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 9)
    )
}

.has_region <- function(df) "region" %in% names(df)

.color_layer <- function(has_reg, region_colors) {
  if (has_reg)
    list(aes(color = region), scale_color_manual(values = region_colors, name = "Region"))
  else
    list(aes(color = NULL))
}

plot_counts <- function(sig_sites, region_colors = NULL) {
  has_reg <- .has_region(sig_sites) && !is.null(region_colors)
  has_dir <- any(!is.na(sig_sites$monotone_direction))

  if (has_reg) {
    counts <- count(sig_sites, region, monotone_direction)
    p <- ggplot(counts, aes(x = region, y = n, fill = region)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = n), vjust = -0.5, size = 3.2, color = "#52514e") +
      scale_fill_manual(values = region_colors, guide = "none")
    if (has_dir) p <- p + facet_wrap(~monotone_direction, labeller = label_both)
  } else if (has_dir) {
    counts <- count(sig_sites, monotone_direction)
    p <- ggplot(counts, aes(x = monotone_direction, y = n)) +
      geom_col(width = 0.6, fill = "#2a78d6") +
      geom_text(aes(label = n), vjust = -0.5, size = 3.2, color = "#52514e")
  } else {
    counts <- tibble(n = nrow(sig_sites))
    p <- ggplot(counts, aes(x = "all", y = n)) +
      geom_col(width = 0.3, fill = "#2a78d6") +
      geom_text(aes(label = n), vjust = -0.5, size = 3.2, color = "#52514e")
  }

  p +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = "Significant sites", x = NULL, y = "Site count") +
    theme_ox()
}

plot_manhattan <- function(sig_sites, region_colors = NULL) {
  has_reg <- .has_region(sig_sites) && !is.null(region_colors)
  has_dir <- any(!is.na(sig_sites$monotone_direction))

  p <- sig_sites |>
    ggplot(aes(x = pos, y = -log10(padj))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "#c3c2b7", linewidth = 0.5)

  if (has_reg) {
    p <- p +
      geom_point(aes(color = region), size = 1.5, alpha = 0.75) +
      scale_color_manual(values = region_colors, name = "Region")
  } else {
    p <- p + geom_point(color = "#2a78d6", size = 1.5, alpha = 0.75)
  }

  if (has_dir)
    p <- p + facet_wrap(~monotone_direction, ncol = 2, labeller = label_both)

  p +
    labs(title = "Significant sites along the locus",
         x = "Position (nt)", y = expression(-log[10](p[adj]))) +
    theme_ox()
}

plot_volcano <- function(sig_sites, groups, region_colors = NULL) {
  has_reg <- .has_region(sig_sites) && !is.null(region_colors)
  has_dir <- any(!is.na(sig_sites$monotone_direction))

  first_col <- paste0("freq_", names(groups)[1])
  last_col  <- paste0("freq_", names(groups)[length(groups)])

  df <- sig_sites |>
    mutate(freq_change = .data[[last_col]] - .data[[first_col]])

  p <- df |> ggplot(aes(x = freq_change, y = -log10(padj)))

  if (has_reg) {
    p <- p +
      geom_point(aes(color = region), size = 1.5, alpha = 0.75) +
      scale_color_manual(values = region_colors, name = "Region")
  } else {
    p <- p + geom_point(color = "#2a78d6", size = 1.5, alpha = 0.75)
  }

  if (has_dir)
    p <- p + facet_wrap(~monotone_direction, labeller = label_both)

  p +
    labs(
      title = "Effect size vs significance",
      x     = paste0("G->T frequency change (",
                     names(groups)[length(groups)], " - ", names(groups)[1], ")"),
      y     = expression(-log[10](p[adj]))
    ) +
    theme_ox()
}

plot_timecourse <- function(dat, sig_sites, direction = NA, groups,
                            n_top = 6, region_colors = NULL) {
  has_reg <- .has_region(sig_sites) && !is.null(region_colors)

  top_sites <- (if (is.na(direction)) sig_sites
                else filter(sig_sites, monotone_direction == direction)) |>
    slice_min(order_by = padj, n = n_top, with_ties = FALSE) |>
    mutate(site_label = if (has_reg)
             paste0(pos, "\n(", region, ")")
           else
             as.character(pos))

  plot_dat <- dat |>
    inner_join(select(top_sites, chr, pos, site_label), by = c("chr", "pos")) |>
    group_by(chr, pos, site_label, group_idx) |>
    mutate(mean_freq = mean(G_to_T_freq)) |>
    ungroup()

  dir_label <- if (is.na(direction)) "significant sites"
               else paste0(direction, "ing G->T frequency")

  p <- plot_dat |>
    ggplot(aes(x = group_idx))

  if (has_reg) {
    p <- p +
      geom_point(aes(y = G_to_T_freq, color = region), size = 2, alpha = 0.6) +
      geom_line(aes(y = mean_freq, group = interaction(chr, pos), color = region),
                linewidth = 0.9) +
      scale_color_manual(values = region_colors, guide = "none")
  } else {
    p <- p +
      geom_point(aes(y = G_to_T_freq), color = "#2a78d6", size = 2, alpha = 0.6) +
      geom_line(aes(y = mean_freq, group = interaction(chr, pos)),
                color = "#2a78d6", linewidth = 0.9)
  }

  p +
    facet_wrap(~site_label, scales = "free_y") +
    scale_x_continuous(breaks = seq_along(groups), labels = names(groups)) +
    labs(
      title    = paste0("Top ", n_top, " ", dir_label),
      subtitle = "Points = individual replicates; line = per-group mean",
      x = "Group", y = "G->T frequency"
    ) +
    theme_ox() +
    theme(strip.text = element_text(size = 7.5))
}

make_plots <- function(dat, sig_sites, groups, output_prefix,
                       region_colors = NULL) {
  directions <- unique(na.omit(sig_sites$monotone_direction))

  pdf(paste0(output_prefix, "_plots.pdf"), width = 10, height = 7)
  print(plot_counts(sig_sites, region_colors))
  print(plot_manhattan(sig_sites, region_colors))
  print(plot_volcano(sig_sites, groups, region_colors))
  dev.off()

  pdf(paste0(output_prefix, "_timecourse.pdf"), width = 14, height = 5)
  if (length(directions) == 0)
    print(plot_timecourse(dat, sig_sites, NA, groups, region_colors = region_colors))
  else
    for (dir in directions)
      print(plot_timecourse(dat, sig_sites, dir, groups, region_colors = region_colors))
  dev.off()

  message("Plots saved: ", output_prefix, "_plots.pdf, ",
          output_prefix, "_timecourse.pdf")
}
