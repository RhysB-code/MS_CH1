# =============================================================
# Long Cane (Potted) Study
# =============================================================

# --- Setup ----------------------------------------------------
rm(list = ls())
setwd("C:/Users/RhysB/OneDrive/Desktop/Chill Spreadsheets")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(RColorBrewer)
  library(emmeans)
  library(multcomp)
})

# --- Constants ------------------------------------------------
p_days           <- seq(7, 70, 7)
p_percent_traits <- c("Budbreak", "Floral Laterals")

# Treatment codes (0–5) shown as chill hours on plots
p_treat_levels <- c("0","1","2","3","4","5")
p_treat_labels <- c("0","168","336","504","672","840")

# Colorblind-safe palette mapped to chill-hour labels
p_colorblind_palette <- c(
  "0"   = "#E69F00",
  "168" = "#56B4E9",
  "336" = "#009E73",
  "504" = "#F0E442",
  "672" = "#0072B2",
  "840" = "#D55E00"
)

base_text_size <- 15

# --- Load & clean --------------------------------------------
potted <- read.csv("Analysis_Spreadsheet__Potted_Final__ForR.csv")
names(potted) <- make.names(trimws(names(potted)), unique = TRUE)

# Keep reasonable rows and factors
potted <- potted %>%
  filter(!is.na(Genotype), Genotype != "", !is.na(Treatment), Treatment != "") %>%
  mutate(
    Genotype  = factor(Genotype,
                       levels = c("A-2491T","Natchez","Navaho","Ouachita","Ponca","Von")),
    Treatment = factor(Treatment)  # keep as codes for now; label in plotting
  )

# =============================================================
# Final bar charts (Day 70) — by genotype within treatment
# =============================================================
plot_final_bars <- function(df, y_col, y_label, scale_factor = 1, ylim_max = NULL) {
  avg_df <- df %>%
    group_by(Genotype, Treatment) %>%
    summarise(mean_value = mean(.data[[y_col]], na.rm = TRUE),
              se = sd(.data[[y_col]],   na.rm = TRUE) / sqrt(n()),
              .groups = "drop") %>%
    mutate(
      mean_value = mean_value * scale_factor,
      se         = se * scale_factor,
      Treatment  = factor(as.character(Treatment),
                          levels = p_treat_levels, labels = p_treat_labels)
    )

  if (is.null(ylim_max)) {
    ylim_max <- ceiling(max(avg_df$mean_value + avg_df$se, na.rm = TRUE))
  }

  plots <- lapply(levels(avg_df$Genotype), function(g) {
    ggplot(filter(avg_df, Genotype == g),
           aes(x = Treatment, y = mean_value)) +
      geom_bar(stat = "identity", fill = "black") +
      geom_errorbar(aes(ymin = mean_value - se, ymax = mean_value + se),
                    width = 0.2) +
      labs(title = g, x = NULL, y = NULL) +
      coord_cartesian(ylim = c(0, ylim_max)) +
      theme_minimal()
  })

  grid.arrange(grobs = plots, ncol = 2,
               top  = textGrob(y_label, gp = gpar(fontsize = 16)),
               left = textGrob(y_label, rot = 90, gp = gpar(fontsize = 12)))
}

plot_final_bars(potted, "X70Day_Potted_Budbreak",        "Final Budbreak (%)",       100, 50)
plot_final_bars(potted, "X70Day_Potted_FloralLaterals",  "Final Floral Laterals (%)",100, 30)
plot_final_bars(potted, "X70Day_Potted_OpenFlowers",     "Number of Open Flowers",     1,100)

# =============================================================
# Reshape weekly data, summaries, AUPC/AUHC, weekly change
# =============================================================

# Trapezoid area helper (AUPC/AUHC)
compute_auc <- function(x, y) {
  if (length(x) < 2) return(0)
  sum((head(y, -1) + tail(y, -1)) / 2 * diff(x))
}

# Robust long-format builder (handles X7Day_... and 7Day_...)
reshape_trait_potted <- function(df, trait_suffix, trait_name) {
  id_all  <- c("Genotype","Treatment","Rep","Replicate")
  id_cols <- intersect(id_all, names(df))
  if (!length(id_cols)) stop("None of the ID columns were found: ", paste(id_all, collapse = ", "))
  
  c1 <- paste0("X", p_days, "Day_Potted_", trait_suffix)
  c2 <- paste0(     p_days, "Day_Potted_", trait_suffix)
  keep <- intersect(c(c1, c2), names(df))
  if (!length(keep)) {
    warning("Potted: missing trait columns for ", trait_name)
    return(tibble::tibble(
      Genotype=character(), Treatment=character(),
      Day_Col=character(),  Value=double(),
      Day=double(),         Trait=character()
    ))
  }
  
  df %>%
    dplyr::select(dplyr::all_of(id_cols), dplyr::all_of(keep)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(keep),
                        names_to = "Day_Col", values_to = "Value") %>%
    dplyr::mutate(
      Day   = as.numeric(stringr::str_match(Day_Col, "X?(\\d+)Day_")[,2]),
      Trait = trait_name
    ) %>%
    dplyr::filter(!is.na(Value))
}


# Build long tables
p_bd <- reshape_trait_potted(potted, "Budbreak",       "Budbreak")
p_fl <- reshape_trait_potted(potted, "FloralLaterals", "Floral Laterals")
p_op <- reshape_trait_potted(potted, "OpenFlowers",    "Open Flowers")
p_all_long <- bind_rows(p_bd, p_fl, p_op)
stopifnot(nrow(p_all_long) > 0, all(p_all_long$Day %in% p_days))

# Convert proportions to % for the percent traits
p_all_long <- p_all_long %>%
  group_by(Trait) %>%
  mutate(Value_disp = if_else(Trait %in% p_percent_traits & max(Value, na.rm = TRUE) <= 1,
                              Value * 100, Value)) %>%
  ungroup()

# Mean curves (Genotype × Treatment × Trait × Day)
p_mean_by <- p_all_long %>%
  group_by(Genotype, Treatment, Trait, Day) %>%
  summarise(Mean = mean(Value_disp, na.rm = TRUE), .groups = "drop") %>%
  arrange(Genotype, Treatment, Trait, Day)

# AUPC (percent traits) / AUHC (Open Flowers)
p_aupc <- p_mean_by %>%
  arrange(Genotype, Treatment, Trait, Day) %>%
  group_by(Genotype, Treatment, Trait) %>%
  summarise(
    AUPC_or_AUHC = compute_auc(Day, Mean),
    Units = if_else(Trait %in% p_percent_traits, "%·days", "count·days"),
    .groups = "drop"
  )
# write.csv(filter(p_aupc, Trait == "Budbreak"),          "Potted_AUPC_Budbreak.csv",         row.names = FALSE)
# write.csv(filter(p_aupc, Trait == "Floral Laterals"),   "Potted_AUPC_FloralLaterals.csv",    row.names = FALSE)
# write.csv(filter(p_aupc, Trait == "Open Flowers") %>% rename(AUHC_countdays = AUPC_or_AUHC),
#           "Potted_AUHC_OpenFlowers.csv", row.names = FALSE)

# Weekly change (Δ from prior week)
p_change <- p_mean_by %>%
  arrange(Genotype, Treatment, Trait, Day) %>%
  group_by(Genotype, Treatment, Trait) %>%
  mutate(Change = Mean - dplyr::lag(Mean, 1, default = 0)) %>%
  ungroup() %>%
  mutate(Units = if_else(Trait %in% p_percent_traits, "%-points", "count"))

# =============================================================
# Plot helpers with uniform y-axis per figure
# =============================================================
compute_ymax <- function(df, value_col = "Mean", pad = 0.05) {
  maxv <- suppressWarnings(max(df[[value_col]], na.rm = TRUE))
  if (!is.finite(maxv)) maxv <- 0
  ymax <- if (maxv <= 0) 1 else maxv * (1 + pad)
  c(0, ymax)
}

plot_progress_uniform <- function(df_progress, trait,
                                  treat_levels, treat_labels, palette,
                                  base_text_size, ylim_vec) {
  d <- df_progress %>%
    filter(Trait == trait) %>%
    mutate(Treatment = factor(as.character(Treatment),
                              levels = treat_levels, labels = treat_labels))
  ylab <- if (trait %in% p_percent_traits) paste0(trait, " (%)") else "Open Flowers (count)"

  ggplot(d, aes(Day, Mean, color = Treatment, group = Treatment)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    facet_wrap(~ Genotype, scales = "free_y") +
    scale_color_manual(values = palette, name = "Chill Hours") +
    labs(title = paste0(trait, " — Progress Over Time"), x = "Day", y = ylab) +
    coord_cartesian(ylim = ylim_vec) +
    theme_minimal(base_size = base_text_size) +
    theme(
      legend.position = "bottom",
      strip.text  = element_text(face = "bold", size = base_text_size),
      axis.title  = element_text(face = "bold"),
      axis.text   = element_text(size = base_text_size - 1),
      legend.text = element_text(size = base_text_size - 1),
      legend.title= element_text(size = base_text_size)
    )
}

plot_change_uniform <- function(df_change, trait,
                                treat_levels, treat_labels, palette,
                                base_text_size, ylim_vec) {
  d <- df_change %>%
    filter(Trait == trait) %>%
    mutate(Treatment = factor(as.character(Treatment),
                              levels = treat_levels, labels = treat_labels))
  ylab <- if (trait %in% p_percent_traits) "Change (%-points)" else "Change (count)"

  ggplot(d, aes(Day, Change, color = Treatment, group = interaction(Genotype, Treatment))) +
    geom_line(alpha = 0.7, linewidth = 0.8) +
    geom_point(size = 1.6) +
    facet_wrap(~ Genotype, scales = "free_y") +
    scale_color_manual(values = palette, name = "Chill Hours") +
    labs(title = paste0(trait, " — Weekly Change"), x = "Day", y = ylab) +
    coord_cartesian(ylim = ylim_vec) +
    theme_minimal(base_size = base_text_size) +
    theme(
      legend.position = "bottom",
      strip.text  = element_text(face = "bold", size = base_text_size),
      axis.title  = element_text(face = "bold"),
      axis.text   = element_text(size = base_text_size - 1),
      legend.text = element_text(size = base_text_size - 1),
      legend.title= element_text(size = base_text_size)
    )
}

# Heatmap of final values (by genotype × chill hours)
plot_trait_heatmap <- function(df_summary, trait,
                               treat_levels, treat_labels,
                               base_text_size = 14) {
  df <- df_summary %>%
    filter(Trait == trait) %>%
    mutate(Treatment = factor(as.character(Treatment),
                              levels = treat_levels, labels = treat_labels))

  if (trait %in% p_percent_traits) {
    df$label <- sprintf("%.1f%%", df$Final)
    fill_lab <- "Final (%)"; ttl <- paste0(trait, " — Final % Heatmap")
    limits  <- c(0, 100)
  } else {
    df$label <- sprintf("%d", round(df$Final))
    fill_lab <- "Final Count"; ttl <- paste0(trait, " — Final Count Heatmap")
    maxv <- suppressWarnings(max(df$Final, na.rm = TRUE))
    if (!is.finite(maxv)) maxv <- 1
    limits <- c(0, max(1, maxv))
  }

  ggplot(df, aes(Treatment, Genotype)) +
    geom_tile(aes(fill = Final), color = "white") +
    geom_text(aes(label = label), size = base_text_size * 0.4) +
    scale_fill_gradient(limits = limits, low = "white", high = "steelblue", name = fill_lab) +
    labs(title = ttl, x = "Chill Hours", y = "Genotype") +
    theme_minimal(base_size = base_text_size) +
    theme(panel.grid = element_blank())
}

# =============================================================
# Draw progress curves, weekly change, heatmaps (uniform y)
# =============================================================
# Progress curves
p_bud_ylim  <- compute_ymax(filter(p_mean_by, Trait == "Budbreak"),         "Mean",   0.05)
p_flor_ylim <- compute_ymax(filter(p_mean_by, Trait == "Floral Laterals"),  "Mean",   0.05)
p_open_ylim <- compute_ymax(filter(p_mean_by, Trait == "Open Flowers"),     "Mean",   0.05)

plot_progress_uniform(p_mean_by, "Budbreak",        p_treat_levels, p_treat_labels, p_colorblind_palette, base_text_size, p_bud_ylim)
plot_progress_uniform(p_mean_by, "Floral Laterals", p_treat_levels, p_treat_labels, p_colorblind_palette, base_text_size, p_flor_ylim)
plot_progress_uniform(p_mean_by, "Open Flowers",    p_treat_levels, p_treat_labels, p_colorblind_palette, base_text_size, p_open_ylim)

# Weekly change
p_bud_ch_ylim  <- compute_ymax(filter(p_change, Trait == "Budbreak"),         "Change", 0.05)
p_flor_ch_ylim <- compute_ymax(filter(p_change, Trait == "Floral Laterals"),  "Change", 0.05)
p_open_ch_ylim <- compute_ymax(filter(p_change, Trait == "Open Flowers"),     "Change", 0.05)

plot_change_uniform(p_change, "Budbreak",        p_treat_levels, p_treat_labels, p_colorblind_palette, base_text_size, p_bud_ch_ylim)
plot_change_uniform(p_change, "Floral Laterals", p_treat_levels, p_treat_labels, p_colorblind_palette, base_text_size, p_flor_ch_ylim)
plot_change_uniform(p_change, "Open Flowers",    p_treat_levels, p_treat_labels, p_colorblind_palette, base_text_size, p_open_ch_ylim)

# Heatmaps (final values per Genotype × Chill hours)
p_summary_final <- p_mean_by %>%
  group_by(Genotype, Treatment, Trait) %>%
  summarise(Final = max(Mean), .groups = "drop")

plot_trait_heatmap(p_summary_final, "Budbreak",        p_treat_levels, p_treat_labels, base_text_size)
plot_trait_heatmap(p_summary_final, "Floral Laterals", p_treat_levels, p_treat_labels, base_text_size)
plot_trait_heatmap(p_summary_final, "Open Flowers",    p_treat_levels, p_treat_labels, base_text_size)

# =============================================================
# Tukey HSD mean-separation boxplots (Day 70) for 3 traits
# =============================================================

get_letters_per_treatment <- function(df_in, trait_col) {
  df_in %>%
    group_by(Treatment) %>%
    group_modify(~{
      if (nrow(.x) < 3 || length(unique(na.omit(.x$Genotype))) < 2)
        return(tibble(Treatment = unique(.x$Treatment),
                      Genotype = character(0), Letters = character(0)))
      .x$neg_trait <- - .x[[trait_col]]  # flip so highest gets "a"
      fit <- aov(neg_trait ~ Genotype, data = .x)
      emm <- emmeans(fit, ~ Genotype)
      multcomp::cld(emm, Letters = letters, adjust = "tukey", alpha = 0.05) |>
        as.data.frame() |>
        transmute(Treatment = unique(.x$Treatment),
                  Genotype, Letters = str_trim(.group))
    }) |>
    ungroup() |>
    filter(!is.na(Genotype))
}

plot_trait_box_with_letters <- function(df_in, trait_col, title_text, y_lab) {
  df_plot <- df_in %>%
    mutate(Treatment = factor(as.character(Treatment),
                              levels = p_treat_levels, labels = p_treat_labels))

  # scale to % if stored as 0–1
  if (max(df_plot[[trait_col]], na.rm = TRUE) <= 1) {
    df_plot[[trait_col]] <- df_plot[[trait_col]] * 100
  }

  letters_df <- get_letters_per_treatment(df_plot, trait_col)

  y_pos <- df_plot %>%
    group_by(Treatment, Genotype) %>%
    summarise(y = max(.data[[trait_col]], na.rm = TRUE), .groups = "drop") %>%
    group_by(Treatment) %>%
    mutate(y = y + 0.05 * max(y, na.rm = TRUE)) %>%
    ungroup()

  letters_plot_df <- left_join(letters_df, y_pos, by = c("Treatment","Genotype"))

  pal <- setNames(brewer.pal(max(3, min(8, nlevels(df_plot$Genotype))), "Set2"),
                  levels(df_plot$Genotype))

  ggplot(df_plot, aes(Genotype, .data[[trait_col]], fill = Genotype)) +
    geom_boxplot(width = 0.65, outlier.shape = 16, outlier.size = 1.8, linewidth = 0.45) +
    geom_text(data = letters_plot_df,
              aes(Genotype, y, label = Letters),
              inherit.aes = FALSE, vjust = 0, size = 4) +
    facet_wrap(~ Treatment, nrow = 2) +
    scale_fill_manual(values = pal, name = "Genotype", drop = FALSE) +
    labs(x = "Genotype", y = y_lab, title = title_text) +
    coord_cartesian(ylim = c(0, NA)) +
    theme_bw(base_size = 13) +
    theme(
      panel.border      = element_rect(color = "grey35", fill = NA, linewidth = 1.1),
      strip.background  = element_rect(fill = "grey85", color = "grey50", linewidth = 1.1),
      strip.text        = element_text(face = "bold", size = 14),
      panel.grid.minor  = element_blank(),
      panel.spacing     = unit(9, "pt"),
      axis.text.x       = element_text(angle = 35, hjust = 1),
      legend.position   = "right"
    )
}

# Day 70 columns for the three traits
plot_trait_box_with_letters(
  potted, "X70Day_Potted_Budbreak",
  "Potted: Final budbreak (%) at Day 70 by genotype and chilling duration",
  "Final (%) Budbreak"
)
plot_trait_box_with_letters(
  potted, "X70Day_Potted_FloralLaterals",
  "Potted: Final floral laterals (%) at Day 70 by genotype and chilling duration",
  "Final (%) Floral Laterals"
)
plot_trait_box_with_letters(
  potted, "X70Day_Potted_OpenFlowers",
  "Potted: Final open flowers (#) at Day 70 by genotype and chilling duration",
  "Final (#) Open Flowers"
)



# =============================================================
# Complement Study
# =============================================================
rm(list = ls())

# 1. Setup -------------------------------
setwd("C:/Users/RhysB/OneDrive/Desktop/Chill Spreadsheets")
list.files(pattern = "ForR\\.csv$")

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)

# 2. Load Data ---------------------------
data <- read.csv("Analysis_Spreadsheet__Complement_Final_ForR.csv")
data <- data %>%
  rename(Genotype = Complement.ID) %>%
  filter(Total.Buds != "no cutting") %>%
  mutate(
    Genotype = as.factor(Genotype),
    Treatment = as.factor(Treatment),
    Total.Buds = as.numeric(Total.Buds)
  )

# 3. Treatment Labels ---------------------
treatment_labels_map <- c(
  "0"   = "0",
  "226" = "226",
  "353" = "353",
  "479" = "479",
  "596" = "596",
  "712" = "712",
  "853" = "853"
)
levels(data$Treatment) <- treatment_labels_map

# 4. Reshape Trait Data -------------------

# Reshape Trait Data (Complement) — robust + fully-qualified
reshape_trait_data <- function(df, trait_suffix, trait_name,
                               day_sequence = seq(7, 42, by = 7),
                               col_prefix = "Complement_") {
  # Use whatever ID columns exist
  id_all  <- c("Genotype","Treatment","Rep","Replicate")
  id_cols <- intersect(id_all, names(df))
  if (!length(id_cols)) stop("None of the ID columns found: ", paste(id_all, collapse=", "))
  
  # Accept both X7Day_... and 7Day_...; handle odd 35-day Budbreak name
  col_names_expected <- paste0("X", day_sequence, "Day_", col_prefix, trait_suffix)
  col_names_alt      <- paste0(     day_sequence, "Day_", col_prefix, trait_suffix)
  keep <- intersect(c(col_names_expected, col_names_alt), names(df))
  
  odd <- "X35Day_Day_Complement_Budbreak"
  if (trait_suffix == "Budbreak" && odd %in% names(df) && !("X35Day_Complement_Budbreak" %in% keep)) {
    keep <- union(keep, odd)
  }
  
  if (!length(keep)) {
    warning("Complement: missing trait columns for ", trait_name)
    return(tibble::tibble(
      Genotype=character(), Treatment=character(),
      Rep=integer(), Replicate=integer(),
      Day_Col=character(), Value=double(), Day=double(), Trait=character()
    ))
  }
  
  df %>%
    dplyr::select(dplyr::all_of(id_cols), dplyr::all_of(keep)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(keep), names_to = "Day_Col", values_to = "Value") %>%
    dplyr::mutate(
      Day   = as.numeric(stringr::str_match(Day_Col, "X?(\\d+)Day_")[,2]),
      Trait = trait_name
    ) %>%
    dplyr::filter(!is.na(Value))
}


# 4. Final Bar Charts ---------------------
plot_final_bars <- function(df, y_suffix, y_label, scale_factor = 1, ylim_max = NULL, integer_y_axis = FALSE) {
  col <- paste0("X42Day_Complement_", y_suffix)
  if (!col %in% colnames(df)) return(NULL)
  
  df_summary <- df %>%
    group_by(Genotype, Treatment) %>%
    summarise(
      mean_value = mean(.data[[col]], na.rm = TRUE),
      se = sd(.data[[col]], na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  if (is.null(ylim_max)) {
    ylim_max <- ceiling(max((df_summary$mean_value + df_summary$se) * scale_factor, na.rm = TRUE))
  }
  
  y_breaks <- if (integer_y_axis) {
    seq(0, ylim_max, by = 1)
  } else {
    waiver()
  }
  
  plots <- lapply(unique(df_summary$Genotype), function(g) {
    ggplot(df_summary %>% filter(Genotype == g), aes(x = Treatment, y = mean_value * scale_factor)) +
      geom_bar(stat = "identity", fill = "black") +
      geom_errorbar(aes(ymin = (mean_value - se) * scale_factor,
                        ymax = (mean_value + se) * scale_factor),
                    width = 0.2) +
      labs(title = g, x = NULL, y = NULL) +
      coord_cartesian(ylim = c(0, ylim_max)) +
      scale_y_continuous(breaks = y_breaks) +
      theme_minimal()
  })
  
  grid.arrange(grobs = plots, ncol = 2,
               top = textGrob(y_label, gp = gpar(fontsize = 16)),
               left = textGrob(y_label, rot = 90, gp = gpar(fontsize = 12)))
}

# Generate bar charts
plot_final_bars(data, "Budbreak", "Final Budbreak (%)", scale_factor = 100, ylim_max = 100)
plot_final_bars(data, "FloralLaterals", "Final Floral Laterals (%)", scale_factor = 100, ylim_max = 100)
plot_final_bars(data, "OpenFlowers", "Final Open Flowers (#)", scale_factor = 1, ylim_max = 5, integer_y_axis = TRUE)

# 5. Combine Reshaped Trait Data ------------------------
bd <- reshape_trait_data(data, "Budbreak", "Budbreak")
fl <- reshape_trait_data(data, "FloralLaterals", "Floral Laterals")
op <- reshape_trait_data(data, "OpenFlowers", "Open Flowers")
all_data <- bind_rows(bd, fl, op)

# =============================================================
# Complement — Reshape, AUPC CSVs (no plots), Lines & Heatmaps
# =============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(ggplot2)
})

# Use the already-loaded Complement data.frame `data`
comp <- data
names(comp) <- make.names(trimws(names(comp)), unique = TRUE)

# =============================================================
# Complement — clean reshape, summaries, uniform-axis plots, heatmaps, Tukey
# =============================================================

# Use the already-loaded Complement df
comp <- data
names(comp) <- make.names(trimws(names(comp)), unique = TRUE)

# ---- constants ----
c_days           <- seq(7, 42, 7)
c_percent_traits <- c("Budbreak", "Floral Laterals")
c_treat_levels   <- c("0","226","353","479","596","712","853")
c_treat_labels   <- c_treat_levels
c_palette <- c("0"="#E69F00","226"="#56B4E9","353"="#009E73","479"="#F0E442","596"="#0072B2","712"="#D55E00","853"="#CC79A7")
base_text_size <- 14

# ---- helpers ----
compute_auc <- function(x, y) {
  if (length(x) < 2) return(0)
  sum((head(y, -1) + tail(y, -1)) / 2 * diff(x))
}

compute_ymax <- function(df, value_col = "Mean", pad = 0.05) {
  maxv <- suppressWarnings(max(df[[value_col]], na.rm = TRUE))
  if (!is.finite(maxv)) maxv <- 0
  ymax <- if (maxv <= 0) 1 else maxv * (1 + pad)
  c(0, ymax)
}

# ---- reshape (robust: X7Day_/7Day_ + odd 35-day Budbreak) ----
reshape_trait_comp <- function(df, trait_suffix, trait_name) {
  id_all  <- c("Genotype","Treatment","Rep","Replicate")
  id_cols <- intersect(id_all, names(df))
  if (!length(id_cols)) stop("No ID columns found among: ", paste(id_all, collapse=", "))
  
  c1 <- paste0("X", c_days, "Day_Complement_", trait_suffix)
  c2 <- paste0(     c_days, "Day_Complement_", trait_suffix)
  keep <- intersect(c(c1, c2), names(df))
  
  odd <- "X35Day_Day_Complement_Budbreak"
  if (trait_suffix == "Budbreak" && odd %in% names(df) && !("X35Day_Complement_Budbreak" %in% keep)) {
    keep <- union(keep, odd)
  }
  
  if (!length(keep)) {
    warning("Complement: missing trait columns for ", trait_name)
    return(tibble::tibble(
      Genotype=character(), Treatment=character(),
      Rep=integer(), Replicate=integer(),
      Day_Col=character(), Value=double(), Day=double(), Trait=character()
    ))
  }
  
  df %>%
    dplyr::select(dplyr::all_of(id_cols), dplyr::all_of(keep)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(keep), names_to = "Day_Col", values_to = "Value") %>%
    dplyr::mutate(
      Day   = as.numeric(stringr::str_match(Day_Col, "X?(\\d+)Day_")[,2]),
      Trait = trait_name
    ) %>%
    dplyr::filter(!is.na(Value))
}

# ---- long tables & summaries ----
c_bd <- reshape_trait_comp(comp, "Budbreak",        "Budbreak")
c_fl <- reshape_trait_comp(comp, "FloralLaterals",  "Floral Laterals")
c_op <- reshape_trait_comp(comp, "OpenFlowers",     "Open Flowers")
c_all_long <- dplyr::bind_rows(c_bd, c_fl, c_op)
stopifnot(nrow(c_all_long) > 0, all(c_all_long$Day %in% c_days))

# % scale if stored 0–1
c_all_long <- c_all_long %>%
  dplyr::group_by(Trait) %>%
  dplyr::mutate(Value_disp = dplyr::case_when(
    Trait %in% c_percent_traits & max(Value, na.rm = TRUE) <= 1 ~ Value * 100,
    TRUE ~ Value
  )) %>%
  dplyr::ungroup()

# mean curves
c_mean_by <- c_all_long %>%
  dplyr::group_by(Genotype, Treatment, Trait, Day) %>%
  dplyr::summarise(Mean = mean(Value_disp, na.rm = TRUE), .groups="drop") %>%
  dplyr::arrange(Genotype, Treatment, Trait, Day)

# weekly change
c_change <- c_mean_by %>%
  dplyr::group_by(Genotype, Treatment, Trait) %>%
  dplyr::mutate(Change = Mean - dplyr::lag(Mean, 1, default = 0)) %>%
  dplyr::ungroup()

# AUPC/AUHC (and CSVs)
c_aupc <- c_mean_by %>%
  dplyr::group_by(Genotype, Treatment, Trait) %>%
  dplyr::summarise(
    AUPC_or_AUHC = compute_auc(Day, Mean),
    Units = ifelse(Trait %in% c_percent_traits, "%·days", "count·days"),
    .groups = "drop"
  )
write.csv(dplyr::filter(c_aupc, Trait == "Budbreak"),         "Complement_AUPC_Budbreak.csv",         row.names = FALSE)
write.csv(dplyr::filter(c_aupc, Trait == "Floral Laterals"),  "Complement_AUPC_FloralLaterals.csv",   row.names = FALSE)
write.csv(dplyr::filter(c_aupc, Trait == "Open Flowers") %>% dplyr::rename(AUHC_countdays = AUPC_or_AUHC),
          "Complement_AUHC_OpenFlowers.csv", row.names = FALSE)

# ---- plotting helpers (uniform y per figure) ----
plot_progress_uniform <- function(df_progress, trait, ylim_vec) {
  d <- df_progress %>%
    dplyr::filter(Trait == trait) %>%
    dplyr::mutate(Treatment = factor(as.character(Treatment),
                                     levels = c_treat_levels, labels = c_treat_labels))
  ylab <- if (trait %in% c_percent_traits) paste0(trait, " (%)") else "Open Flowers (count)"
  ggplot2::ggplot(d, ggplot2::aes(Day, Mean, color = Treatment, group = Treatment)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~ Genotype, scales = "free_y") +
    ggplot2::scale_color_manual(values = c_palette, name = "Chill Hours") +
    ggplot2::labs(title = paste0(trait, " — Progress Over Time"), x = "Day", y = ylab) +
    ggplot2::coord_cartesian(ylim = ylim_vec) +
    ggplot2::theme_minimal(base_size = base_text_size) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text  = ggplot2::element_text(face = "bold", size = base_text_size),
      axis.title  = ggplot2::element_text(face = "bold"),
      axis.text   = ggplot2::element_text(size = base_text_size - 1),
      legend.text = ggplot2::element_text(size = base_text_size - 1),
      legend.title= ggplot2::element_text(size = base_text_size)
    )
}

plot_change_uniform <- function(df_change, trait, ylim_vec) {
  d <- df_change %>%
    dplyr::filter(Trait == trait) %>%
    dplyr::mutate(Treatment = factor(as.character(Treatment),
                                     levels = c_treat_levels, labels = c_treat_labels))
  ylab <- if (trait %in% c_percent_traits) "Change (%-points)" else "Change (count)"
  ggplot2::ggplot(d, ggplot2::aes(Day, Change, color = Treatment,
                                  group = interaction(Genotype, Treatment))) +
    ggplot2::geom_line(alpha = 0.7, linewidth = 0.8) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::facet_wrap(~ Genotype, scales = "free_y") +
    ggplot2::scale_color_manual(values = c_palette, name = "Chill Hours") +
    ggplot2::labs(title = paste0(trait, " — Weekly Change"), x = "Day", y = ylab) +
    ggplot2::coord_cartesian(ylim = ylim_vec) +
    ggplot2::theme_minimal(base_size = base_text_size) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text  = ggplot2::element_text(face = "bold", size = base_text_size),
      axis.title  = ggplot2::element_text(face = "bold"),
      axis.text   = ggplot2::element_text(size = base_text_size - 1),
      legend.text = ggplot2::element_text(size = base_text_size - 1),
      legend.title= ggplot2::element_text(size = base_text_size)
    )
}

# ---- generic heatmap (restores your deleted bit) ----
plot_trait_heatmap <- function(df_summary, trait) {
  df <- df_summary %>%
    dplyr::filter(Trait == trait) %>%
    dplyr::mutate(Treatment = factor(as.character(Treatment),
                                     levels = c_treat_levels, labels = c_treat_labels))
  if (trait %in% c_percent_traits) {
    df$label <- sprintf("%.1f%%", df$Final)
    fill_lab <- "Final (%)"; ttl <- paste0(trait, " — Final % Heatmap")
  } else {
    df$label <- sprintf("%d", round(df$Final))
    fill_lab <- "Final Count"; ttl <- paste0(trait, " — Final Count Heatmap")
  }
  ggplot2::ggplot(df, ggplot2::aes(Treatment, Genotype)) +
    ggplot2::geom_tile(ggplot2::aes(fill = Final), color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = label), size = base_text_size * 0.4) +
    ggplot2::scale_fill_gradient(low = "white", high = "steelblue", name = fill_lab) +
    ggplot2::labs(title = ttl, x = "Chill Hours", y = "Genotype") +
    ggplot2::theme_minimal(base_size = base_text_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

# ---- draw progress / change (uniform y per figure) ----
c_bud_ylim  <- compute_ymax(dplyr::filter(c_mean_by, Trait == "Budbreak"),         "Mean",   0.05)
c_flor_ylim <- compute_ymax(dplyr::filter(c_mean_by, Trait == "Floral Laterals"), "Mean",   0.05)
c_open_ylim <- compute_ymax(dplyr::filter(c_mean_by, Trait == "Open Flowers"),    "Mean",   0.05)

plot_progress_uniform(c_mean_by, "Budbreak",        c_bud_ylim)
plot_progress_uniform(c_mean_by, "Floral Laterals", c_flor_ylim)
plot_progress_uniform(c_mean_by, "Open Flowers",    c_open_ylim)

c_bud_ch_ylim  <- compute_ymax(dplyr::filter(c_change, Trait == "Budbreak"),         "Change", 0.05)
c_flor_ch_ylim <- compute_ymax(dplyr::filter(c_change, Trait == "Floral Laterals"), "Change", 0.05)
c_open_ch_ylim <- compute_ymax(dplyr::filter(c_change, Trait == "Open Flowers"),    "Change", 0.05)

plot_change_uniform(c_change, "Budbreak",        c_bud_ch_ylim)
plot_change_uniform(c_change, "Floral Laterals", c_flor_ch_ylim)
plot_change_uniform(c_change, "Open Flowers",    c_open_ch_ylim)

# ---- heatmaps (final values) ----
c_summary_final <- c_mean_by %>%
  dplyr::group_by(Genotype, Treatment, Trait) %>%
  dplyr::summarise(Final = max(Mean), .groups = "drop")

plot_trait_heatmap(c_summary_final, "Budbreak")
plot_trait_heatmap(c_summary_final, "Floral Laterals")
plot_trait_heatmap(c_summary_final, "Open Flowers")

# ---- Tukey HSD mean-separation boxplots (all traits, Day 42) ----
suppressPackageStartupMessages({
  library(emmeans); library(multcomp); library(RColorBrewer)
})

make_tukey_boxletters <- function(trait_suffix, trait_label, day_col = 42, percent_trait = TRUE) {
  colname <- paste0("X", day_col, "Day_Complement_", trait_suffix)
  if (!colname %in% names(comp)) { message("Missing: ", colname); return(invisible(NULL)) }
  
  dfc <- comp %>%
    dplyr::filter(!is.na(Genotype), Genotype != "", !is.na(Treatment), Treatment != "") %>%
    dplyr::mutate(Treatment = factor(as.character(Treatment), levels = c_treat_levels, labels = c_treat_labels))
  
  if (percent_trait && max(dfc[[colname]], na.rm = TRUE) <= 1) dfc[[colname]] <- dfc[[colname]] * 100
  
  get_letters_per_treatment <- function(df_in, trait_col) {
    df_in %>%
      dplyr::group_by(Treatment) %>%
      dplyr::group_modify(~{
        if (nrow(.x) < 3 || length(unique(stats::na.omit(.x$Genotype))) < 2)
          return(tibble::tibble(Treatment = unique(.x$Treatment), Genotype = character(0), Letters = character(0)))
        .x$neg_trait <- - .x[[trait_col]]
        fit <- stats::aov(neg_trait ~ Genotype, data = .x)
        emm <- emmeans::emmeans(fit, ~ Genotype)
        multcomp::cld(emm, Letters = letters, adjust = "tukey", alpha = 0.05) |>
          as.data.frame() |>
          dplyr::transmute(Treatment = unique(.x$Treatment), Genotype, Letters = stringr::str_trim(.group))
      }) %>% dplyr::ungroup() %>% dplyr::filter(!is.na(Genotype))
  }
  
  letters_df <- get_letters_per_treatment(dfc, colname)
  
  y_pos <- dfc %>%
    dplyr::group_by(Treatment, Genotype) %>%
    dplyr::summarise(y = max(.data[[colname]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(Treatment) %>% dplyr::mutate(y = y + 0.05 * max(y, na.rm = TRUE)) %>% dplyr::ungroup()
  letters_plot_df <- dplyr::left_join(letters_df, y_pos, by = c("Treatment","Genotype"))
  
  pal <- setNames(RColorBrewer::brewer.pal(max(3, min(8, length(unique(dfc$Genotype)))), "Set2"),
                  levels(factor(dfc$Genotype)))
  
  ylab <- if (percent_trait) paste0("Final (%) ", trait_label) else paste0("Final ", trait_label, " (count)")
  ttl  <- paste0("Complement: Final ", trait_label, " at Day ", day_col, " by genotype and chilling duration")
  
  p <- ggplot2::ggplot(dfc, ggplot2::aes(Genotype, .data[[colname]], fill = Genotype)) +
    ggplot2::geom_boxplot(width = 0.65, outlier.shape = 16, outlier.size = 1.8, linewidth = 0.45) +
    ggplot2::geom_text(data = letters_plot_df, ggplot2::aes(Genotype, y, label = Letters),
                       inherit.aes = FALSE, vjust = 0, size = 4) +
    ggplot2::facet_wrap(~ Treatment, nrow = 2) +
    ggplot2::scale_fill_manual(values = pal, name = "Genotype", drop = FALSE) +
    ggplot2::labs(x = "Genotype", y = ylab, title = ttl) +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.border     = ggplot2::element_rect(color = "grey35", fill = NA, linewidth = 1.1),
      strip.background = ggplot2::element_rect(fill = "grey85", color = "grey50", linewidth = 1.1),
      strip.text       = ggplot2::element_text(face = "bold", size = 14),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing    = grid::unit(9, "pt"),
      axis.text.x      = ggplot2::element_text(angle = 35, hjust = 1),
      legend.position  = "right"
    )
  
  print(p)
  
  letters_table <- letters_df %>%
    tidyr::pivot_wider(names_from = Treatment, values_from = Letters) %>%
    dplyr::arrange(Genotype)
  print(letters_table)
  write.csv(letters_table,
            paste0("Complement_D", day_col, "_", gsub(" ", "", trait_label), "_TukeyLetters.csv"),
            row.names = FALSE)
  invisible(list(plot = p, letters = letters_table))
}

make_tukey_boxletters("Budbreak",       "Budbreak",        day_col = 42, percent_trait = TRUE)
make_tukey_boxletters("FloralLaterals", "Floral Laterals", day_col = 42, percent_trait = TRUE)
make_tukey_boxletters("OpenFlowers",    "Open Flowers",    day_col = 42, percent_trait = FALSE)
