# Analysis package — three tasks, computational parameters, AN vs HC

Everything needed to reproduce the analyses reported in the paper. All scripts
are R; paths inside the scripts are relative to **this directory**, so set the
working directory here (`setwd()` or open the folder as an RStudio project)
before sourcing anything.

## Directory layout

```
dati/                 input data (parameters and behavioural indices)
wcst_hmm/             estimation of the WCST parameters (Stan)
tre_compiti/          analyses on the three tasks jointly
risultati/            all output tables and figures, as reported
documenti/            Methods, Results, Supplementary, internal memo
```

## Order of execution

**Stage A — WCST parameters** (`wcst_hmm/`). Long: each model is an HMM fit by
HMC over ~90 participants. Requires the trial-level WCST data, which is *not*
included here (it is upstream of this package); the resulting parameter
estimates are provided in `dati/`, so Stage B runs without re-running Stage A.

| Script | Purpose |
|---|---|
| `01_build_stan_data.R` | trial-level data → Stan data list; choice-coding correction |
| `02_fit_models.R` | fits `hmm_hier`, `hmm_hier4`, `rw_dim_hier` |
| `03_loo_compare.R` | PSIS-LOO comparison, stage 1 |
| `04_reliability.R` | split-half reliability |
| `05_recovery.R`, `funs_simulate.R` | parameter recovery from simulated data |
| `06_ppc.R` | posterior predictive checks |
| `07_external.R` | convergent and cross-task validity (pre-perseveration model) |
| `08_sticky.R` | fits `hmm_sticky` (selected model); PSIS-LOO stage 2 |
| `09_external_sticky.R` | validity of the selected model → `wcst_params_final.csv` |
| `10_nogroup.R` | refit without the group predictor → `wcst_params_nogroup.csv`; quantifies group leakage |

**Stage B — joint analyses** (`tre_compiti/`). Minutes, except `04_mtmm.R`
(Stan, ~10 min) and `06_auc.R` (~10 min).

| Script | Purpose | Output |
|---|---|---|
| `01_prepara_dati.R` | harmonises the three parameter files; normalises group labels; defines the complete-case sample | `dati_armonizzati.csv`, `parametri_meta.csv` |
| `02_effect_size.R` | Cohen's $d$ with bootstrap CIs; sign-concordance test | `effect_sizes.csv`, `concordanza_segni.csv` |
| `03_correlazioni.R` | cross-task correlation structure; WCST reliability | `correlazioni_spearman.csv`, `coppie_correlazioni.csv`, `struttura_cross_task.csv`, `affidabilita_wcst.csv`, `correlazioni_h_wcst.csv` |
| `04_mtmm.R` + `stan/mtmm.stan` | trait/method variance decomposition; LOO over five nested variants | `loo_mtmm.csv`, `decomposizione_varianza.csv`, `contrasto_tratto_metodo.csv`, `varianza_tratto_per_primitiva.csv` |
| `05_eterogeneita.R` | dispersion and centroid permutation tests; per-indicator variance ratios; leave-one-out robustness | `dispersione_test.csv`, `dispersione_loo.csv`, `centroidi_test.csv`, `eterogeneita.csv` |
| `06_auc.R` | nested cross-validated AUC for seven predictor sets | `auc_nested.csv`, `auc_confronti.csv`, `auc_fold.csv` |
| `07_cfa_classica.R` | ML confirmatory factor analysis in `lavaan`, as a check on `04_mtmm.R` | `cfa_classica_fit.csv`, `cfa_classica_loadings.csv` |
| `08_validita_convergente.R` | WCST $h$ against model-free indices, group-free estimates, overall and within group | `validita_convergente_wcst.csv` |

Scripts 02–07 each depend only on `dati_armonizzati.csv`, so after running
`01_prepara_dati.R` they can be run in any order.

## Two things to know before reading the scripts

**1. Two sets of WCST estimates, deliberately.** `wcst_params_final.csv` comes
from the model *with* the group predictor; `wcst_params_nogroup.csv` from the
model *without* it. The first is correct for the group contrast and wrong for
anything individual-level, because each estimate is shrunk toward its own group's
mean (group explains 71% of the between-subject variance of $\log d$ there,
against 19% in the group-free refit). Every script in `tre_compiti/` uses the
**group-free** file. See Supplementary S5.

**2. Indicator sets differ by analysis, on purpose.** Structural and
correlational analyses use the 12 indicators with usable individual-level
reliability (11 DDM parameters + WCST $\mathrm{logit}\,h$); group-effect and
classification analyses use all 14. The reason is stated in each script's
header and in Methods: estimation noise attenuates correlations, so including
unreliable indicators would bias a structural analysis toward the conclusion we
draw, whereas for group contrasts and cross-validated prediction it is merely
conservative.

## Environment

R ≥ 4.2 with `cmdstanr` (+ CmdStan 2.3x), `loo`, `glmnet`, `lavaan`. No other
packages are loaded. Seeds are set at the top of each script; the permutation
and bootstrap results reproduce exactly, and the Stan fits reproduce to Monte
Carlo error.

## Figures

`risultati/fig1`–`fig4` were produced in Python (matplotlib) from the CSV
tables in `risultati/`; the plotting code is not part of this package, since
every plotted quantity is in the tables.

## Not included, and why

- Trial-level data for the three tasks (upstream of this package).
- PRL and task-switching estimation code: those parameters were estimated in
  prior work and are used here exactly as provided, without re-estimation.
- Clinical covariates (BMI, illness duration, subtype, comorbidity, medication)
  and test–retest data: not available for this sample.
