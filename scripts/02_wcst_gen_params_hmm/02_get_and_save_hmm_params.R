d <- readRDS(
  here::here(
    "data",
    "processed",
    "params",
    "pars_hmm_x.RDS"
  )
)

rio::export(
  d,
  here::here(
    "data",
    "processed",
    "wcst_hmm_params.csv"
  )
)
