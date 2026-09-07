# ============================================================
# Figure per il paper: riepilogo feature per GRUPPO (Patients vs Controls)
# - Raincloud per variabile (per task × tipo: params / behav_indices)
# - Correlazioni (heatmap) separate per gruppo
# - Effetti standardizzati (Cohen's d) con CI bootstrap
# ============================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  readr,
  dplyr,
  tidyr,
  tibble,
  stringr,
  ggplot2,
  patchwork,
  ggridges,
  purrr,
  scales,
  forcats
)

# ----------------------------
# 0) Paths input (adatta se diverso)
# ----------------------------
paths <- list(
  wcst_params = here::here(
    "src",
    "auc_out_of_fold",
    "models_params",
    "wcst_hmm_params.csv"
  ),
  wcst_indices = here::here(
    "src",
    "auc_out_of_fold",
    "models_params",
    "wcst_behav_indices.csv"
  ),
  prl_params = here::here(
    "src",
    "auc_out_of_fold",
    "models_params",
    "prl_params.csv"
  ),
  prl_indices = here::here(
    "src",
    "auc_out_of_fold",
    "models_params",
    "prl_behav_indices.csv"
  ),
  ts_params = here::here(
    "src",
    "auc_out_of_fold",
    "models_params",
    "task_switching_params.csv"
  ),
  ts_indices = here::here(
    "src",
    "auc_out_of_fold",
    "models_params",
    "task_switching_behav_indices.csv"
  )
)

# ----------------------------
# 1) Output dirs
# ----------------------------
dir_out <- here::here("src", "auc_out_of_fold", "fig_auc")
dir.create(dir_out, showWarnings = FALSE)
dir_rain <- file.path(dir_out, "raincloud")
dir_heat <- file.path(dir_out, "corr_heatmap")
dir_d <- file.path(dir_out, "effect_size")
dir.create(dir_rain, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_heat, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_d, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 2) Utility: riconoscere colonna gruppo e ricodificare
# ----------------------------
detect_group_col <- function(df) {
  cand <- c("is_patient", "group", "group_label")
  hit <- intersect(cand, names(df))
  if (length(hit) == 0)
    stop("Colonna di gruppo non trovata (is_patient / group / group_label).")
  hit[1]
}

as_group_factor <- function(v) {
  # Supporta 0/1 o stringhe
  if (is.numeric(v)) {
    stopifnot(all(v %in% c(0, 1)))
    factor(
      ifelse(v == 1, "Patients", "Controls"),
      levels = c("Controls", "Patients")
    )
  } else {
    vv <- stringr::str_to_lower(as.character(v))
    vv <- dplyr::case_when(
      vv %in%
        c(
          "1",
          "patient",
          "patients",
          "paziente",
          "pazienti",
          "an",
          "anoressia",
          "anoressiche"
        ) ~
        "Patients",
      vv %in%
        c(
          "0",
          "control",
          "controls",
          "controllo",
          "controlli",
          "hc",
          "healthy",
          "healthy controls"
        ) ~
        "Controls",
      TRUE ~ NA_character_
    )
    factor(vv, levels = c("Controls", "Patients"))
  }
}

# ----------------------------
# 3) Lettura e messa in forma long
# ----------------------------
read_tidy <- function(path, dataset_label) {
  if (!file.exists(path)) stop("File non trovato: ", path)
  df <- readr::read_csv(path, show_col_types = FALSE)

  gcol <- detect_group_col(df)
  yfac <- as_group_factor(df[[gcol]])

  # Colonne da escludere dai predittori
  drop_cols <- c(
    gcol,
    "is_patient",
    "group",
    "group_label",
    "id",
    "ID",
    "subject",
    "subject_id",
    "user_id"
  )

  # Prendi tutte le numeriche (feature candidate)
  num_cols <- names(df)[sapply(df, is.numeric)]
  pred_cols <- setdiff(num_cols, intersect(num_cols, drop_cols))
  if (length(pred_cols) == 0)
    stop("Nessun predittore numerico trovato in: ", path)

  # Tipo = params vs indices (dalla chiave)
  type <- if (grepl("params", dataset_label, ignore.case = TRUE)) "params" else
    "indices"
  task <- case_when(
    grepl("wcst", dataset_label, ignore.case = TRUE) ~ "WCST",
    grepl("prl", dataset_label, ignore.case = TRUE) ~ "PRL",
    grepl("ts", dataset_label, ignore.case = TRUE) ~ "Task Switching",
    TRUE ~ "Task"
  )

  out <- df %>%
    mutate(.group = yfac) %>%
    select(all_of(pred_cols), .group) %>%
    tidyr::drop_na() %>%
    pivot_longer(
      cols = all_of(pred_cols),
      names_to = "feature",
      values_to = "value"
    ) %>%
    group_by(feature) %>%
    mutate(
      z_value = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(task = task, type = type, dataset = paste(task, type, sep = " · "))
  out
}

all_long <- purrr::imap_dfr(paths, ~ read_tidy(.x, .y))

# Ordina feature per coerenza estetica (facoltativo)
# all_long$feature <- forcats::fct_reorder(all_long$feature, all_long$value, median, .desc = TRUE)

# ----------------------------
# 4) Theme & palette
# ----------------------------
pal <- c(Controls = "#0072B2", Patients = "#D55E00")
theme_set(
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "top"
    )
)

# ----------------------------
# 5) RAINCLOUD per gruppo (per dataset)
#    - violin half + box + jitter
#    - valori in scala originale e, opzionale, in z (scegli "value" o "z_value")
# ----------------------------
raincloud_one <- function(df_ds, use_z = FALSE) {
  val <- if (use_z) "z_value" else "value"
  ggplot(
    df_ds,
    aes(x = .group, y = .data[[val]], fill = .group, color = .group)
  ) +
    # violin "half" usando trick di coord_flip + trim
    geom_violin(width = 0.9, alpha = 0.3, trim = TRUE) +
    geom_boxplot(
      width = 0.2,
      outlier.shape = NA,
      alpha = 0.8,
      position = position_dodge(width = 0.7)
    ) +
    geom_jitter(width = 0.12, height = 0, alpha = 0.4, size = 1) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_color_manual(values = pal, guide = "none") +
    facet_wrap(~feature, scales = "free_y", ncol = 4) +
    labs(
      x = NULL,
      y = if (use_z) "Z-score (per variabile)" else "Valore",
      title = unique(df_ds$dataset)
    )
}

# Salva raincloud per ciascun dataset
for (ds in unique(all_long$dataset)) {
  dsub <- all_long %>% filter(dataset == ds)
  # versione su scala originale
  p1 <- raincloud_one(dsub, use_z = FALSE)
  ggsave(
    file.path(dir_rain, paste0(gsub(" ", "_", ds), "_raincloud_value.png")),
    p1,
    width = 12,
    height = 8,
    dpi = 300
  )
  ggsave(
    file.path(dir_rain, paste0(gsub(" ", "_", ds), "_raincloud_value.pdf")),
    p1,
    width = 12,
    height = 8
  )

  # versione z-score (utile quando le scale sono molto diverse)
  p2 <- raincloud_one(dsub, use_z = TRUE)
  ggsave(
    file.path(dir_rain, paste0(gsub(" ", "_", ds), "_raincloud_z.png")),
    p2,
    width = 12,
    height = 8,
    dpi = 300
  )
  ggsave(
    file.path(dir_rain, paste0(gsub(" ", "_", ds), "_raincloud_z.pdf")),
    p2,
    width = 12,
    height = 8
  )
}

# ----------------------------
# 6) CORRELAZIONI per gruppo (heatmap) — per dataset
# ----------------------------
corr_mat <- function(df_wide) {
  # df_wide: colonne = feature numeriche
  cm <- suppressWarnings(cor(df_wide, use = "pairwise.complete.obs"))
  cm
}

plot_corr_heatmap <- function(C, title) {
  df <- as.data.frame(as.table(C))
  names(df) <- c("Var1", "Var2", "value")
  ggplot(df, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      oob = squish
    ) +
    coord_equal() +
    labs(x = NULL, y = NULL, title = title, fill = "r") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

for (ds in unique(all_long$dataset)) {
  dsub <- all_long %>% filter(dataset == ds)

  # Per gruppo
  for (g in levels(dsub$.group)) {
    wide <- dsub %>%
      filter(.group == g) %>%
      select(feature, value) %>%
      group_by(feature) %>%
      mutate(row = row_number()) %>%
      ungroup() %>%
      pivot_wider(names_from = feature, values_from = value)

    # rimuove colonna 'row' se c'è
    wide <- wide %>% select(where(is.numeric))

    if (ncol(wide) >= 2) {
      C <- corr_mat(wide)
      pC <- plot_corr_heatmap(C, paste0(ds, " — Correlazioni (", g, ")"))
      fname <- paste0(gsub(" ", "_", ds), "_corr_", tolower(g), ".")
      ggsave(
        file.path(dir_heat, paste0(fname, "png")),
        pC,
        width = 8,
        height = 7,
        dpi = 300
      )
      ggsave(
        file.path(dir_heat, paste0(fname, "pdf")),
        pC,
        width = 8,
        height = 7
      )
    }
  }
}

# ----------------------------
# 7) Effetto standardizzato (Cohen's d) + CI bootstrap — per dataset
# ----------------------------
# --- Cohen's d e bootstrap (mantieni come già definito) ---
cohen_d <- function(x_ctrl, x_pat) {
  n0 <- length(x_ctrl)
  n1 <- length(x_pat)
  m0 <- mean(x_ctrl)
  m1 <- mean(x_pat)
  s0 <- stats::var(x_ctrl)
  s1 <- stats::var(x_pat)
  sp2 <- ((n0 - 1) * s0 + (n1 - 1) * s1) / (n0 + n1 - 2)
  (m1 - m0) / sqrt(sp2)
}

d_boot_ci <- function(x, g, B = 4000, seed = 1L) {
  set.seed(seed)
  x_ctrl <- x[g == "Controls"]
  x_pat <- x[g == "Patients"]
  if (length(x_ctrl) < 3 || length(x_pat) < 3)
    return(c(d = NA, d_lo = NA, d_hi = NA))
  d0 <- cohen_d(x_ctrl, x_pat)
  db <- replicate(B, {
    xc <- sample(x_ctrl, length(x_ctrl), replace = TRUE)
    xp <- sample(x_pat, length(x_pat), replace = TRUE)
    cohen_d(xc, xp)
  })
  c(d = d0, d_lo = quantile(db, 0.025), d_hi = quantile(db, 0.975))
}

# --- TABELLA EFFETTI: versione robusta senza summarise/unnest_wider ---
effect_table_one_dataset <- function(df_ds) {
  feats <- sort(unique(df_ds$feature))
  purrr::map_dfr(feats, function(f) {
    x <- df_ds %>% dplyr::filter(feature == f)
    x_ctrl <- x$value[x$.group == "Controls"]
    x_pat <- x$value[x$.group == "Patients"]
    dci <- d_boot_ci(x$value, x$.group, B = 4000, seed = 42)
    tibble::tibble(
      feature = f,
      mean_ctrl = mean(x_ctrl, na.rm = TRUE),
      sd_ctrl = sd(x_ctrl, na.rm = TRUE),
      n_ctrl = length(x_ctrl),
      mean_pat = mean(x_pat, na.rm = TRUE),
      sd_pat = sd(x_pat, na.rm = TRUE),
      n_pat = length(x_pat),
      d = unname(dci["d"]),
      d_lo = unname(dci["d_lo"]),
      d_hi = unname(dci["d_hi"])
    )
  }) %>%
    dplyr::arrange(dplyr::desc(abs(d)))
}

# --- PLOT EFFETTI: usa geom_errorbar con orientation="y" ---
plot_effect_size <- function(tbl, title) {
  tbl2 <- tbl %>%
    dplyr::mutate(
      d = as.numeric(d),
      d_lo = as.numeric(d_lo),
      d_hi = as.numeric(d_hi),
      feature = forcats::fct_reorder(feature, d)
    )

  ggplot2::ggplot(tbl2, ggplot2::aes(x = d, y = feature)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, alpha = 0.6) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = d_lo, xmax = d_hi, y = feature),
      width = 0.2,
      orientation = "y"
    ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(
      x = "Cohen's d (Patients - Controls)",
      y = NULL,
      title = title
    ) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))
}

for (ds in unique(all_long$dataset)) {
  dsub <- all_long %>% filter(dataset == ds)
  eff_tbl <- effect_table_one_dataset(dsub)
  # salva tabella
  readr::write_csv(
    eff_tbl,
    file.path(dir_d, paste0(gsub(" ", "_", ds), "_effect_sizes.csv"))
  )

  # plot
  p_d <- plot_effect_size(
    eff_tbl,
    paste0(ds, " — Effetti (d) con CI bootstrap")
  )
  ggsave(
    file.path(dir_d, paste0(gsub(" ", "_", ds), "_effect_sizes.png")),
    p_d,
    width = 8,
    height = max(4, 0.35 * nrow(eff_tbl)),
    dpi = 300
  )
  ggsave(
    file.path(dir_d, paste0(gsub(" ", "_", ds), "_effect_sizes.pdf")),
    p_d,
    width = 8,
    height = max(4, 0.35 * nrow(eff_tbl))
  )
}

# ----------------------------
# 8) (Facoltativo) Figura “multipannello” riassuntiva per il paper
#    - Per ciascun dataset: top-8 feature per |d| (raincloud z + lollipop d)
# ----------------------------
make_compact_panel <- function(df_ds, top_k = 8) {
  eff_tbl <- effect_table_one_dataset(df_ds)
  top_feat <- eff_tbl %>%
    slice_max(order_by = abs(as.numeric(d)), n = min(top_k, n())) %>%
    pull(feature)

  p_left <- df_ds %>%
    filter(feature %in% top_feat) %>%
    mutate(feature = fct_relevel(feature, rev(top_feat))) %>%
    ggplot(aes(x = .group, y = z_value, fill = .group, color = .group)) +
    geom_violin(width = 0.9, alpha = 0.25, trim = TRUE) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.12, height = 0, alpha = 0.35, size = 0.9) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_color_manual(values = pal, guide = "none") +
    facet_wrap(~feature, ncol = 4) +
    labs(x = NULL, y = "Z-score", title = unique(df_ds$dataset))

  p_right <- plot_effect_size(
    eff_tbl %>% filter(feature %in% top_feat),
    title = "Top feature per |d|"
  )

  p_left + p_right + plot_layout(widths = c(3, 2))
}

dir_pan <- file.path(dir_out, "panels")
dir.create(dir_pan, showWarnings = FALSE)
for (ds in unique(all_long$dataset)) {
  dsub <- all_long %>% filter(dataset == ds)
  p <- make_compact_panel(dsub, top_k = 8)
  ggsave(
    file.path(dir_pan, paste0(gsub(" ", "_", ds), "_panel.png")),
    p,
    width = 16,
    height = 9,
    dpi = 300
  )
  ggsave(
    file.path(dir_pan, paste0(gsub(" ", "_", ds), "_panel.pdf")),
    p,
    width = 16,
    height = 9
  )
}

# ============================ FINE SCRIPT ============================
