# Overview ----------------------------------------------------------------
# Associated project: WCST, taks-switching, PRL
# Script purpose: save stan_data for further processing (models' comparison
# for the WCST task)
#
# Written by: Corrado Caudek (corrado.caudek@unifi.it)
# Version: 2025-10-01
# Last update:
# Status: final
# Notes:

# Load necessary libraries ------------------------------------------------

library(pacman)

# Load necessary libraries ------------------------------------------------

if (!requireNamespace("pacman")) install.packages("pacman")

pacman::p_load(
  here,
  tidyverse,
  cmdstanr,
  posterior,
  bayesplot,
  insight,
  pROC,
  stringr,
  loo
)


# Source functions --------------------------------------------------------

source(here::here(
  "scripts",
  "01_data_processing",
  "funs_wcst.R"
))

source(here::here(
  "scripts",
  "01_data_processing",
  "funs_model_selection_wcst.R"
))

# process_and_prepare_stan_data()
source(
  here::here(
    "scripts",
    "01_data_processing",
    "funs_input_for_stan_wcst.R"
  )
)

# Data preparation ----------------------------------------------------------

get_participant_id_and_psytoolkit_code <- function(d) {
  
  library(stringi)
  
  d$mese_c <- ifelse(
    d$`mese:1` < 10, stri_join("0", as.character(d$`mese:1`), sep=''), as.character(d$`mese:1`)
  )
  
  d$giorno_c <- ifelse(
    d$`giorno:1` < 10, 
    stri_join("0", as.character(d$`giorno:1`), sep=''), 
    as.character(d$`giorno:1`)
  )
  
  d$cellulare_c <- ifelse(
    d$`cellulare:1` < 100, 
    stri_join("0", as.character(d$`cellulare:1`), sep=''), 
    as.character(d$`cellulare:1`)
  )
  
  d$sex <- ifelse(d$`sesso:1` == 1, "f",
                  ifelse(d$`sesso:1` == 2, "m", NA))
  
  d$subj_id <- tolower(
    stri_join(d$`nome:1`, d$`cognome:1`, d$`anno:1`, 
              d$mese_c, d$giorno_c, d$cellulare_c, d$sex, 
              sep='_')
  )
  
  data.frame(
    user_id = d$subj_id,
    subj_name = d$`esperimento:1`
  )
}


# Read raw data
raw_data_patients_df <- rio::import(here::here("data", "raw", "patients", "raw_data_patients.csv"))
keys_patients_df <- rio::import(here::here("data", "raw", "patients", "data.xlsx"))
  
raw_data_controls_df <- rio::import(here::here("data", "raw", "controls", "raw_data_controls.csv"))
keys_controls_df <- rio::import(here::here("data", "raw", "controls", "data.xlsx"))

length(unique(raw_data_patients_df$subj_name))
length(unique(raw_data_controls_df$subj_name))


# Get participant's identifier
  
subj_ids_patients <- get_participant_id_and_psytoolkit_code(keys_patients_df)
subj_ids_controls <- get_participant_id_and_psytoolkit_code(keys_controls_df)

d_patients <- left_join(raw_data_patients_df, subj_ids_patients, by = "subj_name")  
d_controls <- left_join(raw_data_controls_df, subj_ids_controls, by = "subj_name")  

# Combine data and rename columns
final_data <- rbind(d_patients, d_controls) |> 
  dplyr::rename(
    code_psytoolkit = subj_name
  )

# Generate correspondence tables
patients_wcst_look_up_tbl <- gen_correspondence_table_codes("patients")
controls_wcst_look_up_tbl <- gen_correspondence_table_codes("controls")
  
# Combine correspondence tables
look_up_tbl <- rbind(patients_wcst_look_up_tbl, controls_wcst_look_up_tbl)
  
# Merge correspondence tables with raw data
raw_data_df <- left_join(raw_data_df, look_up_tbl, by = "code_psytoolkit")

# Import final data
# final_data <- rio::import(here::here("data", "processed", "wcst", "raw_data_project.csv"))

# Group subjects by their group and select unique subjects for each group
code_psytoolkit_selected_an <- unique(final_data[final_data$group == "an", ]$subj_name)
code_psytoolkit_selected_hc <- unique(final_data[final_data$group == "hc", ]$subj_name)
  # code_psytoolkit_selected_ri <- unique(final_data[final_data$group == "ri", ]$subj_name)
  
  # If DEBUGGING is TRUE, take only two subjects from each group
  if (DEBUGGING) {
    code_psytoolkit_selected_an <- head(code_psytoolkit_selected_an, 2)
    code_psytoolkit_selected_hc <- head(code_psytoolkit_selected_hc, 2)
    # code_psytoolkit_selected_ri <- head(code_psytoolkit_selected_ri, 2)
  }
  
  # Vector of subjects who have completed both tasks.
  # Select only participants of the AN and HC groups.
  keep_psytoolkit_codes <- c(
    code_psytoolkit_selected_an, code_psytoolkit_selected_hc
    # code_psytoolkit_selected_ri
  )
  
  # Filter data to keep only the subjects who completed both tasks
  raw_data_two_groups_df <-
    raw_data_df[raw_data_df$code_psytoolkit %in% keep_psytoolkit_codes, ]
  
  # Create the 'group' column
  raw_data_two_groups_df <- raw_data_two_groups_df %>%
    mutate(group = case_when(
      code_psytoolkit %in% code_psytoolkit_selected_an ~ "an",
      code_psytoolkit %in% code_psytoolkit_selected_hc ~ "hc",
      # code_psytoolkit %in% code_psytoolkit_selected_ri ~ "ri",
      TRUE ~ NA_character_
    ))
  
  # Remove unwanted columns
  raw_df <- raw_data_two_groups_df %>%
    dplyr::select(!starts_with("V"))
  
  # Generate the input for the Stan models
  stan_data <- compile_data_for_stan(raw_df)
  
  return(stan_data)
}


# Generate the input for the Stan models
stan_data <- process_and_prepare_stan_data()
str(stan_data)

# Save data -----------------------------------------------------------------

saveRDS(
  stan_data,
  here::here(
    "src",
    "wcst",
    "data",
    "wcst_stan_list.RDS"
  )
)

# eof -----
