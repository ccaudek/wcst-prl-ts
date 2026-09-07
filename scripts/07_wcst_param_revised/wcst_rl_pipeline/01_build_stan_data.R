# 01_build_stan_data.R -----------------------------------------------------
# Ricostruisce i dati WCST per i modelli gerarchici, CORREGGENDO il bug della
# pipeline precedente.
#
# Bug: in funs_input_for_stan_wcst.R, riga 22, la scelta del partecipante era
#      codificata come  resp_choice = correct_card  (colonna V2 = pila CORRETTA)
#      invece di        resp_choice = chosen_card   (colonna V11 = pila SCELTA).
#      Conseguenza: in wcst_stan_data.RDS la sequenza di scelte e' identica per
#      tutti gli 88 soggetti e sempre corretta; l'unica informazione
#      individuale residua e' il feedback (rew).
#
# Input : data/processed/interim/raw_data_project.csv  (dati grezzi + gruppo)
#         data/raw/{patients,controls}/data.xlsx       (codici -> user_id)
#         data/processed/interim/wcst_stan_data.RDS    (solo per verifica)
# Output: wcst_stan_data_hier.RDS

suppressPackageStartupMessages({ library(dplyr); library(readxl); library(stringi) })

DIR       <- "wcst/data_z/data"
PATH_PROJ <- file.path(DIR, "processed/interim/raw_data_project.csv")
PATH_OLD  <- file.path(DIR, "processed/interim/wcst_stan_data.RDS")
OUT       <- "wcst_rl/wcst_stan_data_hier.RDS"

# --- tabella di corrispondenza codice psytoolkit -> user_id (come funs_wcst.R) ---
gen_lookup <- function(GROUP) {
  d <- readxl::read_excel(file.path(DIR, "raw", GROUP, "data.xlsx"))
  # replica esatta di gen_subj_name(): un solo "0" di prefisso quando il valore
  # e' sotto la soglia (per il cellulare la soglia e' 100, non 10)
  pre0 <- function(x, thr) ifelse(as.integer(x) < thr, paste0("0", as.character(x)), as.character(x))
  d$subj_name <- tolower(stringi::stri_join(
    d$`nome:1`, d$`cognome:1`, d$`anno:1`, pre0(d$`mese:1`, 10), pre0(d$`giorno:1`, 10),
    pre0(d$`cellulare:1`, 100), ifelse(d$`sesso:1` == 1, "f", ifelse(d$`sesso:1` == 2, "m", NA)),
    sep = "_"))
  d$code_psytoolkit <- d$`esperimento:1`
  d |> select(subj_name, code_psytoolkit) |> filter(!is.na(code_psytoolkit))
}
look_up <- bind_rows(gen_lookup("patients"), gen_lookup("controls"))

# --- dati grezzi trial-per-trial ---
proj <- read.csv(PATH_PROJ, stringsAsFactors = FALSE) |>
  select(code_psytoolkit = subj_name, group, block, trial_in_a_sequence,
         name_of_task, card_shape, card_number, card_color,
         correct_card, chosen_card, is_correct) |>
  left_join(look_up, by = "code_psytoolkit")

# --- mappatura carta -> pila (verificata contro correct_card) ---
MAP_COLOR <- c(red = 1L, green = 2L, blue = 3L, yellow = 4L)
MAP_SHAPE <- c(circle = 1L, triangle = 2L, cross = 3L, star = 4L)
MAP_RULE  <- c(color = 1L, shape = 2L, number = 3L)
proj <- proj |> mutate(
  k_color  = MAP_COLOR[card_color],
  k_shape  = MAP_SHAPE[card_shape],
  k_number = as.integer(card_number),
  rule     = MAP_RULE[name_of_task])
chk <- with(proj, ifelse(rule == 1L, k_color, ifelse(rule == 2L, k_shape, k_number)))
message("mappatura carta->pila coerente con correct_card: ", round(mean(chk == proj$correct_card), 5))
stopifnot(all(chk == proj$correct_card))

# --- selezione dei soggetti: gli stessi 88 modellati in precedenza, stesso ordine ---
old <- readRDS(PATH_OLD)
keep <- old$user_id
proj <- proj |> filter(subj_name %in% keep)
n_tr <- proj |> count(subj_name)
stopifnot(all(n_tr$n == 60), nrow(n_tr) == length(keep))
proj <- proj |> arrange(match(subj_name, keep), block, trial_in_a_sequence)

N <- length(keep); Tmax <- 60L
mat <- function(col) matrix(proj[[col]], N, Tmax, byrow = TRUE)

stan_data <- list(
  N = N, T = Tmax,
  choice   = mat("chosen_card"),          # <-- LA CORREZIONE
  k_color  = mat("k_color"),
  k_shape  = mat("k_shape"),
  k_number = mat("k_number"),
  rew      = matrix(as.integer(proj$is_correct == 1L), N, Tmax, byrow = TRUE),
  rule     = mat("rule"),
  block    = mat("block"),
  pos_in_block = mat("trial_in_a_sequence"),
  user_id  = keep,
  group_label = ifelse(proj$group[match(keep, proj$subj_name)] == "an", "Patients", "Controls")
)
stan_data$grp <- as.integer(stan_data$group_label == "Patients")
# trial senza risposta (chosen_card == 0): esclusi dalla verosimiglianza
stan_data$valid <- matrix(as.integer(stan_data$choice %in% 1:4), N, Tmax)
stan_data$choice[stan_data$valid == 0L] <- 1L   # segnaposto, mai usato

# --- verifiche ---
message("trial senza risposta: ", sum(stan_data$valid == 0), " su ", N * Tmax,
        " (soggetti coinvolti: ", sum(rowSums(stan_data$valid == 0) > 0), ")")
message("feedback identico a quello della pipeline precedente: ",
        identical(stan_data$rew, matrix(as.integer(old$rew[, 1:Tmax] == 1L), N, Tmax)))
old_choice <- matrix(as.integer(old$resp_choice[, 1:Tmax]), N, Tmax)
message("scelte identiche a quelle usate prima: ", round(mean(stan_data$choice == old_choice), 4),
        "  (la vecchia matrice era costante tra soggetti: ",
        all(apply(old_choice, 2, function(x) length(unique(x)) == 1)), ")")
correct_pile <- with(stan_data, ifelse(rule == 1L, k_color, ifelse(rule == 2L, k_shape, k_number)))
ok <- (stan_data$choice == correct_pile) == (stan_data$rew == 1L)
message("coerenza scelta/feedback nei trial validi: ",
        round(mean(ok[stan_data$valid == 1L]), 4))
message("sequenze di scelte distinte: ", length(unique(apply(stan_data$choice, 1, paste, collapse = ""))),
        " su ", N)
message("accuratezza: media ", round(mean(stan_data$rew), 3),
        ", range ", paste(round(range(rowMeans(stan_data$rew)), 3), collapse = "-"))

saveRDS(stan_data, OUT)
message("salvato ", OUT, " | N=", N, " (pazienti ", sum(stan_data$grp),
        ", controlli ", N - sum(stan_data$grp), ")")
