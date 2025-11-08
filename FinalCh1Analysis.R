# =============================================================
# Long Cane (Potted) Study — reshape, summaries, AUPC/AUHC
# ANOVA + Tukey letters (a = highest); single-row boxplots (0→840)
# Outputs: C:/Users/RhysB/OneDrive/Desktop/CH1_Outputs/Potted_Results
# =============================================================

rm(list = ls())
setwd("C:/Users/RhysB/OneDrive/Desktop/Chill Spreadsheets")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(ggplot2)
  library(gridExtra); library(grid); library(RColorBrewer)
  library(emmeans); library(multcomp); library(boot); library(tidyselect)
})

# ------------------------------- Output dir
p_outdir <- file.path("C:/Users/RhysB/OneDrive/Desktop", "CH1_Outputs", "Potted_Results")
if (!dir.exists(p_outdir)) dir.create(p_outdir, recursive = TRUE, showWarnings = FALSE)

save_png <- function(p, fname, width = 18, height = 12, dpi = 600) {
  ggsave(file.path(p_outdir, fname), p, width = width, height = height, dpi = dpi, bg = "white")
}

# ------------------------------- Constants
p_days           <- seq(7, 70, 7)
p_percent_traits <- c("Budbreak", "Floral Laterals")

# Potted treatments are coded 0–5; labels are chill hours
p_treat_levels <- c("0","1","2","3","4","5")
p_treat_labels <- c("0","168","336","504","672","840")

p_palette <- c("0"="#E69F00","168"="#56B4E9","336"="#009E73",
               "504"="#F0E442","672"="#0072B2","840"="#D55E00")
base_text_size <- 14

# ------------------------------- Load & clean
potted <- read.csv("Analysis_Spreadsheet__Potted_Final__ForR.csv")
names(potted) <- make.names(trimws(names(potted)), unique = TRUE)

# keep rows; Treatment stays coded then labeled in plots
potted <- potted %>%
  filter(!is.na(Genotype), Genotype != "", !is.na(Treatment), Treatment != "") %>%
  mutate(
    Genotype  = factor(Genotype,
                       levels = c("A-2491T","Natchez","Navaho","Ouachita","Ponca","Von")),
    Treatment = factor(Treatment, levels = p_treat_levels, labels = p_treat_levels)
  )

# ------------------------------- Helpers
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

# robust long-format builder (accept X7Day_/7Day_)
reshape_trait_potted <- function(df, trait_suffix, trait_name) {
  id_all  <- c("Genotype","Treatment","Rep","Replicate","plant_ID","Plant","ID")
  id_cols <- intersect(id_all, names(df))
  if (!length(id_cols)) stop("No ID columns among: ", paste(id_all, collapse=", "))
  
  c1 <- paste0("X", p_days, "Day_Potted_", trait_suffix)
  c2 <- paste0(     p_days, "Day_Potted_", trait_suffix)
  keep <- intersect(c(c1, c2), names(df))
  if (!length(keep)) {
    warning("Potted: missing trait columns for ", trait_name)
    return(tibble(
      Genotype=character(), Treatment=character(),
      Day_Col=character(),  Value=double(),
      Day=double(),         Trait=character(), .rep=character()
    ))
  }
  
  df %>%
    dplyr::select(tidyselect::any_of(id_cols), tidyselect::any_of(keep)) %>%
    tidyr::pivot_longer(tidyselect::any_of(keep), names_to = "Day_Col", values_to = "Value") %>%
    dplyr::mutate(
      Day   = as.numeric(stringr::str_match(Day_Col, "X?(\\d+)Day_")[,2]),
      Trait = trait_name
    ) %>%
    dplyr::mutate(.rep = dplyr::coalesce(
      !!!rlang::syms(intersect(c("Replicate","Rep","plant_ID","Plant","ID"), names(.)))
    )) %>%
    dplyr::filter(!is.na(Value))
}

# ------------------------------- Long tables & summaries
p_bd <- reshape_trait_potted(potted, "Budbreak",       "Budbreak")
p_fl <- reshape_trait_potted(potted, "FloralLaterals", "Floral Laterals")
p_op <- reshape_trait_potted(potted, "OpenFlowers",    "Open Flowers")
p_all_long <- bind_rows(p_bd, p_fl, p_op)
stopifnot(nrow(p_all_long) > 0, all(p_all_long$Day %in% p_days))

# convert proportions to % (for Budbreak / Floral Laterals)
p_all_long <- p_all_long %>%
  group_by(Trait) %>%
  mutate(Value_disp = if_else(Trait %in% p_percent_traits & max(Value, na.rm = TRUE) <= 1,
                              Value * 100, Value)) %>%
  ungroup()

# mean curves (for progress/change/heatmaps)
p_mean_by <- p_all_long %>%
  group_by(Genotype, Treatment, Trait, Day) %>%
  summarise(Mean = mean(Value_disp, na.rm = TRUE), .groups = "drop") %>%
  arrange(Genotype, Treatment, Trait, Day)

# replicate-level AUPC/AUHC (for ANOVA/Tukey)
p_aupc_reps <- p_all_long %>%
  arrange(Genotype, Treatment, Trait, .rep, Day) %>%
  group_by(Genotype, Treatment, Trait, .rep) %>%
  summarise(
    AUPC_or_AUHC = compute_auc(Day, Value_disp),
    Units = if_else(first(Trait) %in% p_percent_traits, "%·days", "count·days"),
    .groups = "drop"
  )

# ===========================
# Potted — ANOVA diagnostics (QQ + Residuals vs Fitted)
# ===========================
suppressPackageStartupMessages({library(car); library(ggplot2)})

# Output dir (re-use Potted_Results, create if needed)
p_diag_outdir <- file.path("C:/Users/RhysB/OneDrive/Desktop", "CH1_Outputs", "Potted_Results")
if (!dir.exists(p_diag_outdir)) dir.create(p_diag_outdir, recursive = TRUE, showWarnings = FALSE)
.p_save <- function(p, fname, width = 12, height = 8, dpi = 600) {
  ggplot2::ggsave(file.path(p_diag_outdir, fname), p, width = width, height = height, dpi = dpi, bg = "white")
}

# Ensure required objects exist
stopifnot(exists("p_aupc_reps"))
if (!all(c("AUPC_or_AUHC","Genotype","Treatment") %in% names(p_aupc_reps))) {
  stop("`p_aupc_reps` must contain columns: AUPC_or_AUHC, Genotype, Treatment.")
}

# Fit model if not already present
if (!exists("model")) {
  model <- aov(AUPC_or_AUHC ~ Genotype * Treatment, data = p_aupc_reps)
}

# Build residual dataframe (standardized residuals)
p_res_df <- data.frame(
  Residuals_std = rstandard(model),
  Residuals     = residuals(model),
  Fitted        = fitted(model),
  Group         = interaction(p_aupc_reps$Genotype, p_aupc_reps$Treatment, drop = TRUE)
)

# --- Tests
# Shapiro–Wilk on standardized residuals
p_shap <- shapiro.test(p_res_df$Residuals_std)
p_shap_lab <- sprintf("Shapiro–Wilk  p = %.4f", p_shap$p.value)

# Levene (median center) across Genotype*Treatment groups
p_lev <- car::leveneTest(AUPC_or_AUHC ~ Genotype * Treatment, data = p_aupc_reps, center = median)
p_lev_p <- p_lev[1, "Pr(>F)"]
p_lev_lab <- sprintf("Levene (median)  p = %.4f", p_lev_p)

# --- Plot 1: Q–Q plot of standardized residuals + Shapiro p (top-right)
p_potted_qq <- ggplot(p_res_df, aes(sample = Residuals_std)) +
  stat_qq(size = 1.6) +
  stat_qq_line(color = "red", linewidth = 1) +
  labs(
    title = "Potted ANOVA Diagnostics — Q–Q plot of standardized residuals",
    x = "Theoretical normal quantiles",
    y = "Standardized residuals"
  ) +
  annotate(
    "text", x = Inf, y = Inf, label = p_shap_lab,
    hjust = 1.02, vjust = 1.2, size = 4.2, fontface = "bold", color = "darkred"
  ) +
  theme_bw(base_size = 14)

.p_save(p_potted_qq, "Potted_Diagnostics_QQ_Shapiro.png")

# --- Plot 2: Standardized residuals vs fitted + Levene p (top-right)
p_potted_rvf <- ggplot(p_res_df, aes(x = Fitted, y = Residuals_std)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.9, color = "grey35") +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Potted ANOVA Diagnostics — Standardized residuals vs fitted",
    x = "Fitted values",
    y = "Standardized residuals"
  ) +
  annotate(
    "text", x = Inf, y = Inf, label = p_lev_lab,
    hjust = 1.02, vjust = 1.2, size = 4.2, fontface = "bold", color = "darkred"
  ) +
  theme_bw(base_size = 14)

.p_save(p_potted_rvf, "Potted_Diagnostics_Residuals_vs_Fitted_Levene.png")

message("Saved Potted diagnostics to: ", normalizePath(p_diag_outdir))


# ------------------------------- Progress / change / heatmaps
plot_progress_uniform <- function(df_progress, trait, ylim_vec) {
  d <- df_progress %>%
    filter(Trait == trait) %>%
    mutate(Treatment = factor(as.character(Treatment),
                              levels = p_treat_levels, labels = p_treat_labels))
  ylab <- if (trait %in% p_percent_traits) paste0(trait, " (%)") else "Open Flowers (count)"
  ggplot(d, aes(Day, Mean, color = Treatment, group = Treatment)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    facet_wrap(~ Genotype, scales = "free_y") +
    scale_color_manual(values = p_palette, name = "Chill Hours") +
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
plot_change_uniform <- function(df_change, trait, ylim_vec) {
  d <- df_change %>%
    filter(Trait == trait) %>%
    mutate(Treatment = factor(as.character(Treatment),
                              levels = p_treat_levels, labels = p_treat_labels))
  ylab <- if (trait %in% p_percent_traits) "Change (%-points)" else "Change (count)"
  ggplot(d, aes(Day, Change, color = Treatment, group = interaction(Genotype, Treatment))) +
    geom_line(alpha = 0.7, linewidth = 0.8) +
    geom_point(size = 1.6) +
    facet_wrap(~ Genotype, scales = "free_y") +
    scale_color_manual(values = p_palette, name = "Chill Hours") +
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

p_bud_ylim  <- compute_ymax(filter(p_mean_by, Trait == "Budbreak"),        "Mean", 0.05)
p_flor_ylim <- compute_ymax(filter(p_mean_by, Trait == "Floral Laterals"), "Mean", 0.05)
p_open_ylim <- compute_ymax(filter(p_mean_by, Trait == "Open Flowers"),    "Mean", 0.05)

save_png(plot_progress_uniform(p_mean_by, "Budbreak",        p_bud_ylim),  "Potted_Progress_Budbreak.png")
save_png(plot_progress_uniform(p_mean_by, "Floral Laterals", p_flor_ylim), "Potted_Progress_FloralLaterals.png")
save_png(plot_progress_uniform(p_mean_by, "Open Flowers",    p_open_ylim), "Potted_Progress_OpenFlowers.png")

p_change <- p_mean_by %>%
  group_by(Genotype, Treatment, Trait) %>%
  mutate(Change = Mean - lag(Mean, 1, default = 0)) %>%
  ungroup()

p_bud_ch_ylim  <- compute_ymax(filter(p_change, Trait == "Budbreak"),        "Change", 0.05)
p_flor_ch_ylim <- compute_ymax(filter(p_change, Trait == "Floral Laterals"), "Change", 0.05)
p_open_ch_ylim <- compute_ymax(filter(p_change, Trait == "Open Flowers"),    "Change", 0.05)

save_png(plot_change_uniform(p_change, "Budbreak",        p_bud_ch_ylim),  "Potted_WeeklyChange_Budbreak.png")
save_png(plot_change_uniform(p_change, "Floral Laterals", p_flor_ch_ylim), "Potted_WeeklyChange_FloralLaterals.png")
save_png(plot_change_uniform(p_change, "Open Flowers",    p_open_ch_ylim), "Potted_WeeklyChange_OpenFlowers.png")

# heatmaps (final)
p_summary_final <- p_mean_by %>%
  group_by(Genotype, Treatment, Trait) %>%
  summarise(Final = max(Mean), .groups = "drop")

plot_trait_heatmap <- function(df_summary, trait) {
  df <- df_summary %>%
    filter(Trait == trait) %>%
    mutate(Treatment = factor(as.character(Treatment),
                              levels = p_treat_levels, labels = p_treat_labels))
  if (trait %in% p_percent_traits) {
    df$label <- sprintf("%.1f%%", df$Final); fill_lab <- "Final (%)"; ttl <- paste0(trait, " — Final % Heatmap")
  } else {
    df$label <- sprintf("%d", round(df$Final)); fill_lab <- "Final Count"; ttl <- paste0(trait, " — Final Count Heatmap")
  }
  ggplot(df, aes(Treatment, Genotype)) +
    geom_tile(aes(fill = Final), color = "white") +
    geom_text(aes(label = label), size = base_text_size * 0.4) +
    scale_fill_gradient(low = "white", high = "steelblue", name = fill_lab) +
    labs(title = ttl, x = "Chill Hours", y = "Genotype") +
    theme_minimal(base_size = base_text_size) +
    theme(panel.grid = element_blank())
}

save_png(plot_trait_heatmap(p_summary_final, "Budbreak"),        "Potted_Heatmap_Final_Budbreak.png")
save_png(plot_trait_heatmap(p_summary_final, "Floral Laterals"), "Potted_Heatmap_Final_FloralLaterals.png")
save_png(plot_trait_heatmap(p_summary_final, "Open Flowers"),    "Potted_Heatmap_Final_OpenFlowers.png")

# ------------------------------- Day-70 final boxplots (+ letters, a=highest)
get_letters_per_treatment <- function(df_in, trait_col) {
  df_in %>%
    group_by(Treatment) %>%
    group_modify(~{
      if (nrow(.x) < 3 || length(unique(stats::na.omit(.x$Genotype))) < 2)
        return(tibble(Treatment = unique(.x$Treatment), Genotype = character(0), Letters = character(0)))
      .x$neg_trait <- - .x[[trait_col]]                     # a = highest
      fit <- stats::aov(neg_trait ~ Genotype, data = .x)
      emm <- emmeans::emmeans(fit, ~ Genotype)
      multcomp::cld(emm, Letters = letters, adjust = "tukey", alpha = 0.05) |>
        as.data.frame() |>
        transmute(Treatment = unique(.x$Treatment), Genotype, Letters = stringr::str_trim(.group))
    }) %>% ungroup() %>% filter(!is.na(Genotype))
}

plot_day70_boxletters <- function(df_in, trait_suffix, trait_label, y_is_percent = TRUE) {
  colname <- paste0("X70Day_Potted_", trait_suffix)
  if (!colname %in% names(potted)) { message("Missing: ", colname); return(invisible(NULL)) }
  
  df_plot <- potted
  # label treatments with chill hours on facets
  df_plot$TreatmentLab <- factor(plyr::mapvalues(df_plot$Treatment, from = p_treat_levels, to = p_treat_labels),
                                 levels = p_treat_labels)
  
  if (y_is_percent && max(df_plot[[colname]], na.rm = TRUE) <= 1)
    df_plot[[colname]] <- df_plot[[colname]] * 100
  
  letters_df <- get_letters_per_treatment(df_plot, colname)
  
  # 30% global headroom
  global_max <- max(df_plot[[colname]], na.rm = TRUE)
  upper_y    <- global_max * 1.30
  
  y_pos <- df_plot %>%
    group_by(TreatmentLab, Genotype) %>%
    summarise(y = max(.data[[colname]], na.rm = TRUE), .groups = "drop") %>%
    group_by(TreatmentLab) %>% mutate(y = y + 0.20 * max(y, na.rm = TRUE)) %>% ungroup()
  
  letters_plot_df <- left_join(letters_df %>%
                                 mutate(TreatmentLab = factor(plyr::mapvalues(Treatment,
                                                                              from = p_treat_levels,
                                                                              to   = p_treat_labels),
                                                              levels = p_treat_labels)),
                               y_pos,
                               by = c("TreatmentLab","Genotype"))
  
  pal <- setNames(RColorBrewer::brewer.pal(max(3, min(8, nlevels(df_plot$Genotype))), "Set2"),
                  levels(df_plot$Genotype))
  
  ylab <- if (y_is_percent) paste0("Final (%) ", trait_label) else paste0("Final ", trait_label, " (count)")
  ttl  <- paste0("Potted: Final ", trait_label, " at Day 70 by genotype and chilling duration")
  
  p <- ggplot(df_plot, aes(Genotype, .data[[colname]], fill = Genotype)) +
    geom_boxplot(width = 0.65, outlier.shape = 16, outlier.size = 1.8, linewidth = 0.45) +
    geom_text(data = letters_plot_df, aes(Genotype, y, label = Letters),
              inherit.aes = FALSE, vjust = 0, size = 4) +
    facet_wrap(~ TreatmentLab, nrow = 1, drop = FALSE) +    # single row 0→840
    scale_fill_manual(values = pal, name = "Genotype", drop = FALSE) +
    labs(x = "Genotype", y = ylab, title = ttl) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    coord_cartesian(ylim = c(0, upper_y), expand = FALSE, clip = "off") +
    theme_bw(base_size = 13) +
    theme(
      plot.margin      = margin(t = 28, r = 18, b = 14, l = 14),
      panel.border     = element_rect(color = "grey35", fill = NA, linewidth = 1.1),
      strip.background = element_rect(fill = "grey85", color = "grey50", linewidth = 1.1),
      strip.text       = element_text(face = "bold", size = 14),
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(9, "pt"),
      axis.text.x      = element_text(angle = 35, hjust = 1),
      legend.position  = "right"
    )
  
  save_png(p, paste0("Potted_D70_", gsub(" ", "", trait_label), "_Tukey_BoxOnly.png"),
           width = 28, height = 10, dpi = 600)
  invisible(p)
}

# Day 70 boxplots
plot_day70_boxletters(potted, "Budbreak",       "Budbreak",        TRUE)
plot_day70_boxletters(potted, "FloralLaterals", "Floral Laterals",  TRUE)
plot_day70_boxletters(potted, "OpenFlowers",    "Open Flowers",     FALSE)

# ------------------------------- AUPC/AUHC ANOVA + Tukey (a = highest) + single-row boxplots
plot_aupc_boxletters <- function(df_trait, trait_label, file_stub) {
  # letters with a = highest via negation
  letters_df <- df_trait %>%
    group_by(Treatment) %>%
    group_modify(~{
      dat <- .x
      if (nrow(dat) < 3 || length(unique(dat$Genotype)) < 2)
        return(tibble(Treatment = unique(dat$Treatment), Genotype = character(0), Letters = character(0)))
      dat$neg_area <- -dat$AUPC_or_AUHC
      fit <- stats::aov(neg_area ~ Genotype, data = dat)
      emm <- emmeans::emmeans(fit, ~ Genotype)
      multcomp::cld(emm, Letters = letters, adjust = "tukey", alpha = 0.05) |>
        as.data.frame() |>
        transmute(Genotype,
                  Letters = stringr::str_trim(.group),
                  Treatment = unique(dat$Treatment))
    }) %>% ungroup()
  
  # map Treatment codes to labels for facets in order 0→840
  df_trait <- df_trait %>%
    mutate(TreatmentLab = factor(plyr::mapvalues(Treatment,
                                                 from = p_treat_levels, to = p_treat_labels),
                                 levels = p_treat_labels))
  letters_df <- letters_df %>%
    mutate(TreatmentLab = factor(plyr::mapvalues(Treatment,
                                                 from = p_treat_levels, to = p_treat_labels),
                                 levels = p_treat_labels))
  
  # 30% global headroom
  global_max <- max(df_trait$AUPC_or_AUHC, na.rm = TRUE)
  upper_y    <- global_max * 1.30
  
  y_pos <- df_trait %>%
    group_by(TreatmentLab, Genotype) %>%
    summarise(y = max(AUPC_or_AUHC, na.rm = TRUE), .groups = "drop") %>%
    group_by(TreatmentLab) %>% mutate(y = y + 0.20 * max(y, na.rm = TRUE)) %>% ungroup()
  
  lab_df <- left_join(letters_df, y_pos, by = c("TreatmentLab","Genotype")) %>%
    filter(!is.na(Genotype), !is.na(y))
  
  pal <- setNames(RColorBrewer::brewer.pal(max(3, min(8, length(unique(df_trait$Genotype)))), "Set2"),
                  levels(factor(df_trait$Genotype)))
  
  ylab <- paste0(ifelse(trait_label=="Open Flowers","AUHC","AUPC"), " (", unique(df_trait$Units), ")")
  ttl  <- paste0("Potted: ", ifelse(trait_label=="Open Flowers","AUHC","AUPC"),
                 " — ", trait_label, " by genotype and chilling duration")
  
  p <- ggplot(df_trait, aes(Genotype, AUPC_or_AUHC, fill = Genotype)) +
    geom_boxplot(width = 0.65, outlier.shape = 16, outlier.size = 1.8, linewidth = 0.45) +
    geom_text(data = lab_df, aes(Genotype, y, label = Letters),
              inherit.aes = FALSE, vjust = 0, size = 4) +
    facet_wrap(~ TreatmentLab, nrow = 1, drop = FALSE) +     # single row 0→840
    scale_fill_manual(values = pal, name = "Genotype", drop = FALSE) +
    labs(x = "Genotype", y = ylab, title = ttl) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    coord_cartesian(ylim = c(0, upper_y), expand = FALSE, clip = "off") +
    theme_bw(base_size = 13) +
    theme(
      plot.margin      = margin(t = 28, r = 18, b = 14, l = 14),
      panel.border     = element_rect(color = "grey35", fill = NA, linewidth = 1.1),
      strip.background = element_rect(fill = "grey85", color = "grey50", linewidth = 1.1),
      strip.text       = element_text(face = "bold", size = 14),
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(9, "pt"),
      axis.text.x      = element_text(angle = 35, hjust = 1),
      legend.position  = "right"
    )
  
  save_png(p, paste0(file_stub, "_Tukey_BoxOnly.png"), width = 28, height = 10, dpi = 600)
  
  letters_wide <- letters_df %>%
    select(-Treatment) %>%
    tidyr::pivot_wider(names_from = TreatmentLab, values_from = Letters) %>%
    arrange(Genotype)
  
  write.csv(letters_wide, file.path(p_outdir, paste0(file_stub, "_TukeyLetters.csv")), row.names = FALSE)
  invisible(list(plot = p, letters = letters_wide))
}

# Build AUPC/AUHC plots for each trait
make_all_potted_aupc_plots <- function(trait_label, file_base) {
  df_trait <- p_aupc_reps %>% filter(Trait == trait_label)
  if (!nrow(df_trait)) { message("No replicate AUPC rows for trait: ", trait_label); return(invisible(NULL)) }
  plot_aupc_boxletters(df_trait, trait_label, paste0(file_base))
}
make_all_potted_aupc_plots("Budbreak",        "Potted_AUPC_Budbreak")
make_all_potted_aupc_plots("Floral Laterals", "Potted_AUPC_FloralLaterals")
make_all_potted_aupc_plots("Open Flowers",    "Potted_AUHC_OpenFlowers")

message("All Potted outputs saved to: ", normalizePath(p_outdir))
print(list.files(p_outdir, full.names = TRUE))

# =============================================================
# Complement Study — reshape, summaries, AUPC/AUHC + ANOVA/Tukey
# Single-row boxplots (0→853), letters "a = highest", ~30% headroom
# Diagnostics: QQ (Shapiro) + Residuals-vs-Fitted (Levene)
# Outputs: C:/Users/RhysB/OneDrive/Desktop/CH1_Outputs/Complement_Results
# =============================================================

rm(list = ls())
setwd("C:/Users/RhysB/OneDrive/Desktop/Chill Spreadsheets")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(ggplot2)
  library(gridExtra); library(grid); library(RColorBrewer)
  library(emmeans); library(multcomp); library(boot); library(tidyselect); library(plyr)
  library(car)
})

# ------------------------------- Output dir + saver
c_outdir <- file.path("C:/Users/RhysB/OneDrive/Desktop", "CH1_Outputs", "Complement_Results")
if (!dir.exists(c_outdir)) dir.create(c_outdir, recursive = TRUE, showWarnings = FALSE)
save_png <- function(p, fname, width = 18, height = 12, dpi = 600) {
  ggplot2::ggsave(file.path(c_outdir, fname), p, width = width, height = height, dpi = dpi, bg = "white")
}

# ------------------------------- Constants
c_days           <- seq(7, 42, 7)
c_percent_traits <- c("Budbreak", "Floral Laterals")
c_treat_levels   <- c("0","226","353","479","596","712","853")
c_treat_labels   <- c_treat_levels
c_palette        <- c("0"="#E69F00","226"="#56B4E9","353"="#009E73",
                      "479"="#F0E442","596"="#0072B2","712"="#D55E00","853"="#CC79A7")
base_text_size   <- 14

# ------------------------------- Load & clean
data <- read.csv("Analysis_Spreadsheet__Complement_Final_ForR.csv")
data <- data %>%
  rename(Genotype = Complement.ID) %>%
  filter(Total.Buds != "no cutting") %>%
  mutate(
    Genotype   = as.factor(Genotype),
    Treatment  = factor(Treatment, levels = c_treat_levels, labels = c_treat_labels),
    Total.Buds = suppressWarnings(as.numeric(Total.Buds))
  )
comp <- data
names(comp) <- make.names(trimws(names(comp)), unique = TRUE)

# ------------------------------- Helpers
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

# robust reshape (accepts X7Day_/7Day_ + odd 35-day Budbreak col)
reshape_trait_comp <- function(df, trait_suffix, trait_name) {
  id_all  <- c("Genotype","Treatment","Rep","Replicate","plant_ID","Plant","ID")
  id_cols <- intersect(id_all, names(df))
  if (!length(id_cols)) stop("No ID columns among: ", paste(id_all, collapse=", "))
  
  c1 <- paste0("X", c_days, "Day_Complement_", trait_suffix)
  c2 <- paste0(     c_days, "Day_Complement_", trait_suffix)
  keep <- intersect(c(c1, c2), names(df))
  
  odd <- "X35Day_Day_Complement_Budbreak"
  if (trait_suffix == "Budbreak" && odd %in% names(df) && !("X35Day_Complement_Budbreak" %in% keep)) {
    keep <- union(keep, odd)
  }
  
  if (!length(keep)) {
    warning("Complement: missing trait columns for ", trait_name)
    return(tibble(
      Genotype=character(), Treatment=character(), Rep=integer(), Replicate=integer(),
      Day_Col=character(), Value=double(), Day=double(), Trait=character(), .rep=character()
    ))
  }
  
  df %>%
    dplyr::select(tidyselect::any_of(id_cols), tidyselect::any_of(keep)) %>%
    tidyr::pivot_longer(tidyselect::any_of(keep), names_to = "Day_Col", values_to = "Value") %>%
    dplyr::mutate(
      Day   = as.numeric(stringr::str_match(Day_Col, "X?(\\d+)Day_")[,2]),
      Trait = trait_name
    ) %>%
    dplyr::mutate(.rep = dplyr::coalesce(
      !!!rlang::syms(intersect(c("Replicate","Rep","plant_ID","Plant","ID"), names(.)))
    )) %>%
    dplyr::filter(!is.na(Value))
}

# ------------------------------- Long tables & summaries
c_bd <- reshape_trait_comp(comp, "Budbreak",       "Budbreak")
c_fl <- reshape_trait_comp(comp, "FloralLaterals", "Floral Laterals")
c_op <- reshape_trait_comp(comp, "OpenFlowers",    "Open Flowers")
c_all_long <- dplyr::bind_rows(c_bd, c_fl, c_op)
stopifnot(nrow(c_all_long) > 0, all(c_all_long$Day %in% c_days))

# scale percent traits to %
c_all_long <- c_all_long %>%
  dplyr::group_by(Trait) %>%
  dplyr::mutate(Value_disp = dplyr::if_else(Trait %in% c_percent_traits & max(Value, na.rm = TRUE) <= 1,
                                            Value * 100, Value)) %>%
  dplyr::ungroup()

# mean curves (for progress/change/heatmaps)
c_mean_by <- c_all_long %>%
  dplyr::group_by(Genotype, Treatment, Trait, Day) %>%
  dplyr::summarise(Mean = mean(Value_disp, na.rm = TRUE), .groups="drop") %>%
  dplyr::arrange(Genotype, Treatment, Trait, Day)

# replicate-level AUPC/AUHC (for ANOVA/Tukey)
c_aupc_reps <- c_all_long %>%
  dplyr::arrange(Genotype, Treatment, Trait, .rep, Day) %>%
  dplyr::group_by(Genotype, Treatment, Trait, .rep) %>%
  dplyr::summarise(
    AUPC_or_AUHC = compute_auc(Day, Value_disp),
    Units = if_else(first(Trait) %in% c_percent_traits, "%·days", "count·days"),
    .groups = "drop"
  )

# ------------------------------- Progress / change / heatmaps
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
  ggplot2::ggplot(d, ggplot2::aes(Day, Change, color = Treatment, group = interaction(Genotype, Treatment))) +
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

c_bud_ylim  <- compute_ymax(dplyr::filter(c_mean_by, Trait == "Budbreak"),        "Mean", 0.05)
c_flor_ylim <- compute_ymax(dplyr::filter(c_mean_by, Trait == "Floral Laterals"), "Mean", 0.05)
c_open_ylim <- compute_ymax(dplyr::filter(c_mean_by, Trait == "Open Flowers"),    "Mean", 0.05)

save_png(plot_progress_uniform(c_mean_by, "Budbreak",        c_bud_ylim),  "Complement_Progress_Budbreak.png")
save_png(plot_progress_uniform(c_mean_by, "Floral Laterals", c_flor_ylim), "Complement_Progress_FloralLaterals.png")
save_png(plot_progress_uniform(c_mean_by, "Open Flowers",    c_open_ylim), "Complement_Progress_OpenFlowers.png")

c_change <- c_mean_by %>%
  dplyr::group_by(Genotype, Treatment, Trait) %>%
  dplyr::mutate(Change = Mean - dplyr::lag(Mean, 1, default = 0)) %>%
  dplyr::ungroup()

c_bud_ch_ylim  <- compute_ymax(dplyr::filter(c_change, Trait == "Budbreak"),        "Change", 0.05)
c_flor_ch_ylim <- compute_ymax(dplyr::filter(c_change, Trait == "Floral Laterals"), "Change", 0.05)
c_open_ch_ylim <- compute_ymax(dplyr::filter(c_change, Trait == "Open Flowers"),    "Change", 0.05)

save_png(plot_change_uniform(c_change, "Budbreak",        c_bud_ch_ylim),  "Complement_WeeklyChange_Budbreak.png")
save_png(plot_change_uniform(c_change, "Floral Laterals", c_flor_ch_ylim), "Complement_WeeklyChange_FloralLaterals.png")
save_png(plot_change_uniform(c_change, "Open Flowers",    c_open_ch_ylim), "Complement_WeeklyChange_OpenFlowers.png")

# ------------------------------- Heatmaps (final)
c_summary_final <- c_mean_by %>%
  dplyr::group_by(Genotype, Treatment, Trait) %>%
  dplyr::summarise(Final = max(Mean), .groups = "drop")

plot_trait_heatmap <- function(df_summary, trait) {
  df <- df_summary %>%
    dplyr::filter(Trait == trait) %>%
    dplyr::mutate(Treatment = factor(as.character(Treatment),
                                     levels = c_treat_levels, labels = c_treat_labels))
  if (trait %in% c_percent_traits) {
    df$label <- sprintf("%.1f%%", df$Final); fill_lab <- "Final (%)"; ttl <- paste0(trait, " — Final % Heatmap")
  } else {
    df$label <- sprintf("%d", round(df$Final)); fill_lab <- "Final Count"; ttl <- paste0(trait, " — Final Count Heatmap")
  }
  ggplot2::ggplot(df, ggplot2::aes(Treatment, Genotype)) +
    ggplot2::geom_tile(ggplot2::aes(fill = Final), color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = label), size = base_text_size * 0.4) +
    ggplot2::scale_fill_gradient(low = "white", high = "steelblue", name = fill_lab) +
    ggplot2::labs(title = ttl, x = "Chill Hours", y = "Genotype") +
    ggplot2::theme_minimal(base_size = base_text_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}
save_png(plot_trait_heatmap(c_summary_final, "Budbreak"),        "Complement_Heatmap_Final_Budbreak.png")
save_png(plot_trait_heatmap(c_summary_final, "Floral Laterals"), "Complement_Heatmap_Final_FloralLaterals.png")
save_png(plot_trait_heatmap(c_summary_final, "Open Flowers"),    "Complement_Heatmap_Final_OpenFlowers.png")

# ------------------------------- Day-42 boxplots with Tukey (a = highest), single-row 0→853
get_letters_per_treatment <- function(df_in, trait_col) {
  df_in %>%
    dplyr::group_by(Treatment) %>%
    dplyr::group_modify(~{
      if (nrow(.x) < 3 || length(unique(stats::na.omit(.x$Genotype))) < 2)
        return(tibble::tibble(Treatment = unique(.x$Treatment), Genotype = character(0), Letters = character(0)))
      .x$neg_trait <- - .x[[trait_col]]                     # a = highest
      fit <- stats::aov(neg_trait ~ Genotype, data = .x)
      emm <- emmeans::emmeans(fit, ~ Genotype)
      multcomp::cld(emm, Letters = letters, adjust = "tukey", alpha = 0.05) |>
        as.data.frame() |>
        dplyr::transmute(Treatment = unique(.x$Treatment),
                         Genotype, Letters = stringr::str_trim(.group))
    }) %>% dplyr::ungroup() %>% dplyr::filter(!is.na(Genotype))
}

plot_day42_boxletters <- function(df_in, trait_suffix, trait_label, y_is_percent = TRUE) {
  colname <- paste0("X42Day_Complement_", trait_suffix)
  if (!colname %in% names(comp)) { message("Missing: ", colname); return(invisible(NULL)) }
  
  df_plot <- comp
  # facet labels in ascending 0→853 (single row)
  df_plot$TreatmentLab <- factor(plyr::mapvalues(as.character(df_plot$Treatment),
                                                 from = c_treat_labels, to = c_treat_labels),
                                 levels = c_treat_labels)
  
  if (y_is_percent && max(df_plot[[colname]], na.rm = TRUE) <= 1)
    df_plot[[colname]] <- df_plot[[colname]] * 100
  
  letters_df <- get_letters_per_treatment(df_plot, colname)
  
  # 30% global headroom so letters won't clip
  global_max <- max(df_plot[[colname]], na.rm = TRUE)
  upper_y    <- global_max * 1.30
  
  y_pos <- df_plot %>%
    dplyr::group_by(TreatmentLab, Genotype) %>%
    dplyr::summarise(y = max(.data[[colname]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(TreatmentLab) %>% dplyr::mutate(y = y + 0.20 * max(y, na.rm = TRUE)) %>% dplyr::ungroup()
  
  letters_plot_df <- dplyr::left_join(
    letters_df %>% dplyr::mutate(TreatmentLab = factor(plyr::mapvalues(as.character(Treatment),
                                                                       from = c_treat_levels, to = c_treat_labels),
                                                       levels = c_treat_labels)),
    y_pos, by = c("TreatmentLab","Genotype")
  ) %>% dplyr::filter(!is.na(Genotype), !is.na(y))
  
  pal <- setNames(RColorBrewer::brewer.pal(max(3, min(8, nlevels(df_plot$Genotype))), "Set2"),
                  levels(df_plot$Genotype))
  
  ylab <- if (y_is_percent) paste0("Final (%) ", trait_label) else paste0("Final ", trait_label, " (count)")
  ttl  <- paste0("Complement: Final ", trait_label, " at Day 42 by genotype and chilling duration")
  
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(Genotype, .data[[colname]], fill = Genotype)) +
    ggplot2::geom_boxplot(width = 0.65, outlier.shape = 16, outlier.size = 1.8, linewidth = 0.45) +
    ggplot2::geom_text(data = letters_plot_df, ggplot2::aes(Genotype, y, label = Letters),
                       inherit.aes = FALSE, vjust = 0, size = 4) +
    ggplot2::facet_wrap(~ TreatmentLab, nrow = 1, drop = FALSE) +  # single row 0→853
    ggplot2::scale_fill_manual(values = pal, name = "Genotype", drop = FALSE) +
    ggplot2::labs(x = "Genotype", y = ylab, title = ttl) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    ggplot2::coord_cartesian(ylim = c(0, upper_y), expand = FALSE, clip = "off") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.margin      = margin(t = 28, r = 18, b = 14, l = 14),
      panel.border     = element_rect(color = "grey35", fill = NA, linewidth = 1.1),
      strip.background = element_rect(fill = "grey85", color = "grey50", linewidth = 1.1),
      strip.text       = element_text(face = "bold", size = 14),
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(9, "pt"),
      axis.text.x      = element_text(angle = 35, hjust = 1),
      legend.position  = "right"
    )
  
  save_png(p, paste0("Complement_D42_", gsub(" ", "", trait_label), "_Tukey_BoxOnly.png"),
           width = 28, height = 10, dpi = 600)
  invisible(p)
}

# Day 42 boxplots
plot_day42_boxletters(comp, "Budbreak",       "Budbreak",        TRUE)
plot_day42_boxletters(comp, "FloralLaterals", "Floral Laterals", TRUE)
plot_day42_boxletters(comp, "OpenFlowers",    "Open Flowers",    FALSE)

# ------------------------------- AUPC/AUHC ANOVA + Tukey (a = highest) + single-row boxplots
plot_aupc_boxletters <- function(df_trait, trait_label, file_stub) {
  # Tukey letters with a = highest via negation
  letters_df <- df_trait %>%
    dplyr::group_by(Treatment) %>%
    dplyr::group_modify(~{
      dat <- .x
      if (nrow(dat) < 3 || length(unique(dat$Genotype)) < 2)
        return(tibble::tibble(Treatment = unique(dat$Treatment), Genotype = character(0), Letters = character(0)))
      dat$neg_area <- -dat$AUPC_or_AUHC
      fit <- stats::aov(neg_area ~ Genotype, data = dat)
      emm <- emmeans::emmeans(fit, ~ Genotype)
      multcomp::cld(emm, Letters = letters, adjust = "tukey", alpha = 0.05) |>
        as.data.frame() |>
        dplyr::transmute(Genotype,
                         Letters = stringr::str_trim(.group),
                         Treatment = unique(dat$Treatment))
    }) %>% dplyr::ungroup()
  
  # map codes to labels for facets 0→853
  df_trait <- df_trait %>%
    dplyr::mutate(TreatmentLab = factor(plyr::mapvalues(as.character(Treatment),
                                                        from = c_treat_levels, to = c_treat_labels),
                                        levels = c_treat_labels))
  letters_df <- letters_df %>%
    dplyr::mutate(TreatmentLab = factor(plyr::mapvalues(as.character(Treatment),
                                                        from = c_treat_levels, to = c_treat_labels),
                                        levels = c_treat_labels))
  
  # 30% global headroom
  global_max <- max(df_trait$AUPC_or_AUHC, na.rm = TRUE)
  upper_y    <- global_max * 1.30
  
  y_pos <- df_trait %>%
    dplyr::group_by(TreatmentLab, Genotype) %>%
    dplyr::summarise(y = max(AUPC_or_AUHC, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(TreatmentLab) %>% dplyr::mutate(y = y + 0.20 * max(y, na.rm = TRUE)) %>% dplyr::ungroup()
  
  lab_df <- dplyr::left_join(letters_df, y_pos, by = c("TreatmentLab","Genotype")) %>%
    dplyr::filter(!is.na(Genotype), !is.na(y))
  
  pal <- setNames(RColorBrewer::brewer.pal(max(3, min(8, length(unique(df_trait$Genotype)))), "Set2"),
                  levels(factor(df_trait$Genotype)))
  
  ylab <- paste0(ifelse(trait_label=="Open Flowers","AUHC","AUPC"), " (", unique(df_trait$Units), ")")
  ttl  <- paste0("Complement: ", ifelse(trait_label=="Open Flowers","AUHC","AUPC"),
                 " — ", trait_label, " by genotype and chilling duration")
  
  p <- ggplot2::ggplot(df_trait, ggplot2::aes(Genotype, AUPC_or_AUHC, fill = Genotype)) +
    ggplot2::geom_boxplot(width = 0.65, outlier.shape = 16, outlier.size = 1.8, linewidth = 0.45) +
    ggplot2::geom_text(data = lab_df, ggplot2::aes(Genotype, y, label = Letters),
                       inherit.aes = FALSE, vjust = 0, size = 4) +
    ggplot2::facet_wrap(~ TreatmentLab, nrow = 1, drop = FALSE) +     # single row 0→853
    ggplot2::scale_fill_manual(values = pal, name = "Genotype", drop = FALSE) +
    ggplot2::labs(x = "Genotype", y = ylab, title = ttl) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    ggplot2::coord_cartesian(ylim = c(0, upper_y), expand = FALSE, clip = "off") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.margin      = margin(t = 28, r = 18, b = 14, l = 14),
      panel.border     = element_rect(color = "grey35", fill = NA, linewidth = 1.1),
      strip.background = element_rect(fill = "grey85", color = "grey50", linewidth = 1.1),
      strip.text       = element_text(face = "bold", size = 14),
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(9, "pt"),
      axis.text.x      = element_text(angle = 35, hjust = 1),
      legend.position  = "right"
    )
  
  save_png(p, paste0(file_stub, "_Tukey_BoxOnly.png"), width = 28, height = 10, dpi = 600)
  
  letters_wide <- letters_df %>%
    dplyr::select(-Treatment) %>%
    tidyr::pivot_wider(names_from = TreatmentLab, values_from = Letters) %>%
    dplyr::arrange(Genotype)
  
  write.csv(letters_wide, file.path(c_outdir, paste0(file_stub, "_TukeyLetters.csv")), row.names = FALSE)
  invisible(list(plot = p, letters = letters_wide))
}

# Build AUPC/AUHC plots for each trait
make_all_complement_aupc_plots <- function(trait_label, file_base) {
  df_trait <- c_aupc_reps %>% dplyr::filter(Trait == trait_label)
  if (!nrow(df_trait)) { message("No replicate AUPC rows for trait: ", trait_label); return(invisible(NULL)) }
  plot_aupc_boxletters(df_trait, trait_label, paste0(file_base))
}
make_all_complement_aupc_plots("Budbreak",        "Complement_AUPC_Budbreak")
make_all_complement_aupc_plots("Floral Laterals", "Complement_AUPC_FloralLaterals")
make_all_complement_aupc_plots("Open Flowers",    "Complement_AUHC_OpenFlowers")

# ------------------------------- ANOVA diagnostics (using c_aupc_reps)
# Fit ANOVA on AUPC/AUHC with Genotype*Treatment
c_model <- aov(AUPC_or_AUHC ~ Genotype * Treatment, data = c_aupc_reps)

# standardized residuals dataframe
residuals_df <- data.frame(
  Residuals_std = rstandard(c_model),
  Residuals     = residuals(c_model),
  Fitted        = fitted(c_model),
  Group         = interaction(c_aupc_reps$Genotype, c_aupc_reps$Treatment, drop = TRUE)
)

# Shapiro–Wilk (normality of standardized residuals)
shapiro_res     <- shapiro.test(residuals_df$Residuals_std)
shapiro_p_label <- sprintf("Shapiro–Wilk  p = %.4f", shapiro_res$p.value)

# Levene (median center) for homogeneity of variance
lev_res      <- car::leveneTest(AUPC_or_AUHC ~ Genotype * Treatment, data = c_aupc_reps, center = median)
lev_p_value  <- lev_res[1, "Pr(>F)"]
levene_label <- sprintf("Levene (median)  p = %.4f", lev_p_value)

# Q–Q plot (standardized residuals) + Shapiro p
p_qq_shapiro <- ggplot2::ggplot(residuals_df, ggplot2::aes(sample = Residuals_std)) +
  ggplot2::stat_qq(size = 1.6) +
  ggplot2::stat_qq_line(color = "red", linewidth = 1) +
  ggplot2::labs(
    title = "Complement ANOVA Diagnostics — Q–Q plot of standardized residuals",
    x = "Theoretical normal quantiles",
    y = "Standardized residuals"
  ) +
  ggplot2::annotate("text", x = Inf, y = Inf, label = shapiro_p_label,
                    hjust = 1.02, vjust = 1.2, size = 4.2, fontface = "bold", color = "darkred") +
  ggplot2::theme_bw(base_size = 14)
save_png(p_qq_shapiro, "Complement_Diagnostics_QQ_Shapiro.png", width = 12, height = 8)

# Residuals vs Fitted + Levene p
p_res_vs_fit <- ggplot2::ggplot(residuals_df, ggplot2::aes(x = Fitted, y = Residuals_std)) +
  ggplot2::geom_point(alpha = 0.7) +
  ggplot2::geom_smooth(method = "loess", se = FALSE, linewidth = 0.9, color = "grey35") +
  ggplot2::geom_hline(yintercept = 0, color = "blue", linetype = "dashed", linewidth = 1) +
  ggplot2::labs(
    title = "Complement ANOVA Diagnostics — Standardized residuals vs fitted",
    x = "Fitted values",
    y = "Standardized residuals"
  ) +
  ggplot2::annotate("text", x = Inf, y = Inf, label = levene_label,
                    hjust = 1.02, vjust = 1.2, size = 4.2, fontface = "bold", color = "darkred") +
  ggplot2::theme_bw(base_size = 14)
save_png(p_res_vs_fit, "Complement_Diagnostics_Residuals_vs_Fitted_Levene.png", width = 12, height = 8)

message("All Complement outputs saved to: ", normalizePath(c_outdir))
print(list.files(c_outdir, full.names = TRUE))

