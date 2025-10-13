# ============================================
# Estrarre valori unici: subj_id, subj_idx, group
# ============================================

# Pacchetti
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

# Task switching params --------------------------------------------------------
path_csv <- here::here(
  "data",
  "processed",
  "lookup_tables",
  "look_up_tab_task_switch.csv"
)

# Lettura
dat <- read_csv(path_csv, show_col_types = FALSE)

# Controllo colonne richieste
needed <- c("subj_id", "subj_idx", "group")
stopifnot(all(needed %in% names(dat)))

# Combinazioni uniche (una riga per soggetto)
uniq_tab <- dat |>
  distinct(subj_id, subj_idx, group) |>
  arrange(group, subj_idx, subj_id)

task_switch_params_df <- rio::import(
  here::here(
    "data",
    "processed",
    "task_switching_params.csv"
  )
)

ts_tot <- full_join(uniq_tab, task_switch_params_df, by = "subj_idx")

# Salva tabella combinazioni uniche
out_csv <- here::here(
  "data",
  "processed",
  "task_switching_params.csv"
)
write_csv(ts_tot, out_csv)
message("Salvato: ", out_csv)

# PRL params -------------------------------------------------------------------
prl_path_csv <- here::here(
  "data",
  "processed",
  "lookup_tables",
  "lookup_tab_prl.csv"
)

# Lettura
dat_prl <- read_csv(prl_path_csv, show_col_types = FALSE) |>
  dplyr::select(-`...1`)

# Controllo colonne richieste
needed <- c("subj_idx", "subj_code", "group")
stopifnot(all(needed %in% names(dat_prl)))

# Combinazioni uniche (una riga per soggetto)
uniq_tab_prl <- dat_prl |>
  distinct(subj_code, subj_idx, group) |>
  arrange(group, subj_idx, subj_code)

prl_params_df <- rio::import(
  here::here(
    "data",
    "processed",
    "prl_params.csv"
  )
) |>
  dplyr::rename(subj_idx = index)

prl_tot <- full_join(uniq_tab_prl, prl_params_df, by = "subj_idx")

# Salva tabella combinazioni uniche
prl_out_csv <- here::here(
  "data",
  "processed",
  "prl_params.csv"
)
write_csv(prl_tot, prl_out_csv)
message("Salvato: ", prl_out_csv)
