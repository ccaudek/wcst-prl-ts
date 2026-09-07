library(tidyverse)

options(max.print = 3000)


# WCSR
d_params <- rio::import(
  here::here("data", "processed", "wcst_hmm_params.csv")
) |>
  dplyr::select(-c("group", "subject"))

d_params <- d_params %>%
  mutate(
    user_id = recode(
      user_id,
      "em_or_2003_01_01_101_f" = "em_or_2003_01_02_101_f",
      "giulia_toma_1997_07_30_762_f" = "gi_to_1997_07_30_762_f"
    )
  )

d_behav <- rio::import(
  here::here("data", "processed", "wcst_behav_indices.csv")
) |>
  dplyr::rename(
    user_id = subj_name
  )

# d_params$user_id |> sort()
# d_behav$user_id |> sort()
# d <- full_join(d_params, d_behav, by = "user_id")

rio::export(
  d_params,
  here::here(
    "data",
    "processed",
    "wcst_hmm_params.csv"
  )
)


# TS ----

d_params <- rio::import(
  here::here("data", "processed", "task_switching_params.csv")
) |>
  dplyr::rename(
    user_id = subj_id
  ) |>
  dplyr::select(-c("subj_idx", "is_patient"))

d_behav <- rio::import(
  here::here("data", "processed", "task_switching_behav_indices.csv")
) |>
  dplyr::rename(
    user_id = subj_id
  )

d_behav <- d_behav %>%
  mutate(
    user_id = recode(
      user_id,
      "fe_sa_2002_05_09_08_f" = "fe_sa_2002_05_09_008_f",
      "ch_ma_2001_10_27_331_f" = "ch_ma_2001_10_27_332_f"
    )
  )


d <- full_join(d_params, d_behav, by = "user_id")

print(d |> arrange(user_id))

sort(d_behav$user_id)

rio::export(
  d_behav,
  here::here(
    "data",
    "processed",
    "task_switching_behav_indices.csv"
  )
)

rio::export(
  d_params,
  here::here("data", "processed", "task_switching_params.csv")
)


# PRL -----

d_params <- rio::import(
  here::here("data", "processed", "prl_params.csv")
) |>
  dplyr::rename(
    user_id = subj_code
  )

d_behav <- rio::import(
  here::here("data", "processed", "prl_behav_indices.csv")
) |>
  dplyr::rename(
    user_id = subj_name
  )

d <- full_join(d_params, d_behav, by = "user_id")

rio::export(
  d_params,
  here::here(
    "data",
    "processed",
    "prl_params.csv"
  )
)

rio::export(
  d_behav,
  here::here(
    "data",
    "processed",
    "prl_behav_indices.csv"
  )
)

d_params <- d_params |>
  dplyr::select(
    -c("subj_idx", "is_patient")
  )
