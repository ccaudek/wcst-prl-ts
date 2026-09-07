# Supplementary Material — Computational model of the Wisconsin Card Sorting Test

This document describes the generative model used to obtain individual-level
parameters from WCST trial-by-trial data, the model comparison that selected it,
its validation, and the reason two sets of estimates were produced from the same
model.

## S1. Data and preprocessing

Participants sorted cards for $T = 60$ trials, choosing one of four target piles
per trial. The relevant sorting rule (colour, shape, or number) changed without
warning across blocks. On each trial $t$, participant $i$ observed a stimulus
card, chose a pile $c_{it} \in \{1,2,3,4\}$, and received binary feedback
$r_{it} \in \{0,1\}$.

For each trial we recorded which pile matched the stimulus on each of the three
dimensions: $k^{\text{col}}_{it}$, $k^{\text{sha}}_{it}$, $k^{\text{num}}_{it}$.
These indices define a $4 \times 3$ indicator matrix
$\mathbf{M}_{it}$ with $[\mathbf{M}_{it}]_{jk} = 1$ if pile $j$ matches the
stimulus on dimension $k$, and $0$ otherwise. Trials without a response were
excluded from the likelihood through a validity indicator $v_{it}$.

The analysable WCST sample was $N = 88$ participants.

### S1.1 Correction of a coding error in the upstream pipeline

**This section is essential for interpreting any earlier analysis of these
data.** In the pipeline used before the present work, the participant's choice
was assigned from the column holding the *correct* pile rather than the column
holding the *chosen* pile. The consequence is severe and easy to verify: the
resulting choice sequence was identical across all participants and always
correct, so the only participant-specific information surviving in the modelling
input was the feedback vector. Any individual-difference parameter estimated
from that input cannot reflect individual choice behaviour.

The data were rebuilt from the raw trial-level records with the choice taken
from the chosen-pile column (`wcst_rl/01_build_stan_data.R`). Three checks
confirm the rebuild: the feedback matrix is bit-identical to the previous one;
the choice matrix is no longer constant across participants; and, on valid
trials, the event "chosen pile equals the pile that is correct under the active
rule" agrees with the observed feedback. All results reported in the paper use
the rebuilt data.

## S2. The selected model: hierarchical HMM rule inference with dimension-level perseveration

The participant is modelled as maintaining a belief $\mathbf{b}_{it} \in
\Delta^2$ over which of the three dimensions is currently rewarded, updating it
by Bayes' rule after each outcome, and allowing for the possibility that the
rule has changed.

**Belief initialisation.** $\mathbf{b}_{i1} = (1/3, 1/3, 1/3)^\top$.

**Choice rule.** Let $\mathbf{q}_{it} \in \{0,1\}^3$ be the indicator of which
dimensions the chosen pile matches, i.e. row $c_{it}$ of $\mathbf{M}_{it}$.
Choice probabilities combine the belief with a perseveration bonus toward the
dimension followed on the previous trial:

$$
\mathbf{p}_{it} = (1-\lambda)\,\mathrm{softmax}\!\left(\mathbf{M}_{it}\left[d_i\,\mathbf{b}_{it} + \kappa_i\,\mathbf{q}_{i,t-1}\right]\right) + \frac{\lambda}{4},
\qquad \mathbf{q}_{i0} = \mathbf{0},
$$

$$
c_{it} \sim \mathrm{Categorical}(\mathbf{p}_{it}) \quad \text{for trials with } v_{it}=1 .
$$

Here $d_i > 0$ is the **determinism** (inverse decision temperature): how
sharply belief translates into choice. $\kappa_i \in \mathbb{R}$ is
**dimension-level perseveration**: the bonus given to piles matching the
dimension just followed, *independently of belief*. $\lambda$ is a
uniform lapse rate shared across participants.

**Belief update.** Feedback is informative about the active dimension with
reliability $1-\eta$: the probability of positive feedback is $1-\eta$ if the
chosen pile matched the active dimension and $\eta$ otherwise. With
$\boldsymbol{\pi}_{it} = \mathbf{q}_{it}(1-\eta) + (\mathbf{1}-\mathbf{q}_{it})\eta$
and likelihood
$\boldsymbol{\ell}_{it} = \boldsymbol{\pi}_{it}$ if $r_{it}=1$, else
$\mathbf{1}-\boldsymbol{\pi}_{it}$:

$$
\tilde{\mathbf{b}}_{it} = \frac{\mathbf{b}_{it} \odot \boldsymbol{\ell}_{it}}{\mathbf{b}_{it}^\top \boldsymbol{\ell}_{it}},
\qquad
\mathbf{b}_{i,t+1} = (1-h_i)\,\tilde{\mathbf{b}}_{it} + h_i\,\frac{\mathbf{1}-\tilde{\mathbf{b}}_{it}}{2}.
$$

The second step is the hidden-Markov transition: with **switch rate** $h_i \in
(0,1)$ the participant entertains that the rule has moved, redistributing belief
uniformly over the other two dimensions. Large $h_i$ means belief is discarded
readily; small $h_i$ means it is retained across disconfirming feedback. This is
the parameter we interpret as *rule inference*.

**Individual-level parameterisation.** The three subject-level quantities are
estimated on unconstrained scales, which are also the scales used in all
downstream analyses:

$$
\boldsymbol{\theta}_i = (\mathrm{logit}\,h_i,\ \log d_i,\ \kappa_i)^\top,
\qquad
\boldsymbol{\theta}_i = \boldsymbol{\mu} + \boldsymbol{\beta}\,g_i + \mathbf{L}\,\mathbf{z}_i,
\qquad \mathbf{z}_i \sim \mathcal{N}(\mathbf{0}, \mathbf{I}_3),
$$

where $g_i = 1$ for patients and $0$ for controls, and $\mathbf{L} =
\mathrm{diag}(\boldsymbol{\sigma})\,\mathbf{L}_\Omega$ with $\mathbf{L}_\Omega$
the Cholesky factor of the correlation matrix of individual deviations. This is
a non-centred parameterisation.

**Priors.**

$$
\begin{aligned}
\mu_1 &\sim \mathcal{N}(-1.5,\ 1.0) & \mu_2 &\sim \mathcal{N}(1.8,\ 0.7) & \mu_3 &\sim \mathcal{N}(0,\ 1.0) \\
\boldsymbol{\beta} &\sim \mathcal{N}(\mathbf{0},\ 0.5^2 \mathbf{I}) &
\boldsymbol{\sigma} &\sim \mathcal{N}^{+}(0,\ 1) &
\mathbf{L}_\Omega &\sim \mathrm{LKJ}(2) \\
\mathrm{logit}\,\lambda &\sim \mathcal{N}(-4.0,\ 1.0) &
\mathrm{logit}(2\eta) &\sim \mathcal{N}(-2.8,\ 1.0) & &
\end{aligned}
$$

with $\lambda = \mathrm{logit}^{-1}(\cdot)$ and $\eta = \tfrac{1}{2}\,\mathrm{logit}^{-1}(\cdot) \in (0, 0.5)$.
The priors on $\mu_1$ and $\mu_2$ are weakly informative and centred on values
that generate WCST-plausible behaviour under prior predictive simulation; the
prior on $\boldsymbol{\beta}$ is deliberately conservative relative to the
effects of interest.

**Estimation.** Stan (`cmdstanr`), 4 chains, No-U-Turn Sampler.
Model file: `wcst_hmm/stan/hmm_sticky.stan`; fitting script:
`wcst_hmm/08_sticky.R`.

## S3. Model comparison

Four models were compared by expected log pointwise predictive density estimated
by Pareto-smoothed importance-sampling leave-one-out cross-validation (PSIS-LOO,
`loo` package), computed at the trial level. Chance performance is
$\log(1/4) = -1.386$ per trial.

The comparison was run in two stages. First, among the three models without
perseveration:

| Model | Description | $\Delta\mathrm{elpd}$ | SE |
|---|---|---|---|
| `hmm_hier` | HMM rule inference, three-dimensional belief | 0.0 | — |
| `hmm_hier4` | as above, four-dimensional belief | −1.8 | 0.8 |
| `rw_dim_hier` | Rescorla–Wagner learning over dimension values | −679.2 | 42.9 |

Then, the selected model against the best of these:

| Model | $\Delta\mathrm{elpd}$ | SE |
|---|---|---|
| `hmm_sticky` | 0.0 | — |
| `hmm_hier` | −63.8 | 12.2 |

The selected model achieves −0.304 elpd per trial against
−1.386 for chance, with no Pareto $\hat{k} > 0.7$. Model-free reinforcement
learning over dimension values is decisively worse and additionally shows 90
observations with $\hat{k} > 0.7$, indicating misfit rather than merely lower
predictive accuracy.

Adding dimension perseveration was motivated by systematic posterior predictive
failure of `hmm_hier`, which underpredicted win–stay (0.978 observed vs 0.956
predicted) and dimension repetition (0.866 vs 0.842) while overpredicting
lose–shift (0.671 vs 0.731) and post-switch recovery speed — four discrepancies
with a single interpretation: participants repeat the just-followed dimension
more than belief alone warrants.

## S4. Validation

**Convergence.** For the final model: no divergent transitions, $\hat{R}_{\max}
= 1.006$, minimum effective sample size 681.

**Posterior predictive checks** (`wcst_hmm/06_ppc.R`, output `ppc_sticky.csv`).
Eight behavioural signatures were simulated from the posterior. Five fall inside
the 90% predictive interval (overall accuracy, perseverative-error rate,
dimension repetition, lose–shift, late-block accuracy); three remain outside
(non-perseverative errors, win–stay, early-block accuracy), all by small
margins. Correlations between observed and predicted values *across
participants* range from 0.50 to 0.93, with 0.93 for overall accuracy and 0.82
for perseverative errors: the model tracks who is more and less impaired, not
only the group average.

**Parameter recovery** (`wcst_hmm/05_recovery.R`, output `recovery_sticky.csv`).
Data were simulated from the posterior with the design and $T$ of the real
experiment and refitted. Correlations between generating and recovered values:

| Parameter | Raw $r$ | $r$ net of group |
|---|---|---|
| $h$ | 0.83 | 0.79 |
| $\log d$ | 0.72 | 0.68 |
| $\kappa$ | 0.32 | −0.22 |

$h$ and $\log d$ recover acceptably. $\kappa$ does not recover at the individual
level: its apparent raw correlation is carried by the group difference, and
vanishes once that is removed. $\kappa$ is therefore reported as a group-level
effect only.

**Split-half reliability within session** (`wcst_hmm/04_reliability.R`, output
`reliability_sticky.csv`): $\mathrm{logit}\,h$ 0.67, $\log d$ 0.51, $\kappa$
0.46.

**Convergent validity with conventional WCST indices**
(`wcst_hmm/07_external.R`, output `convergent_validity_final.csv`). Spearman
correlations of $h$ with model-free indices: accuracy $-0.68$,
non-perseverative errors $+0.58$, dimension repetition $-0.59$, perseverative
errors $+0.46$, lose–shift $-0.38$. The switch rate behaves as a measure of
failure to maintain and exploit a rule, as intended.

Because these correlations use the group-informed estimates, they could in
principle be inflated by the group difference being present in both terms. The
check was repeated on the group-free estimates for the three indices available
as a saved table (`tre_compiti/validita_convergente_wcst.csv`, $n = 86$), both
overall and after centring within group: perseverative errors $+0.63$ / $+0.63$,
non-perseverative errors $+0.75$ / $+0.70$, perseverative responses $-0.20$ /
$-0.23$. The associations are, if anything, stronger on the group-free
estimates, and are not carried by the group contrast.

## S5. Two sets of estimates from the same model, and which to use

The hierarchical model includes the group indicator $g_i$ in the individual-level
mean. This is the correct specification for estimating group effects: it gives
the group contrast $\boldsymbol{\beta}$ its own parameter with its own
uncertainty. But it makes each participant's posterior mean shrink toward *that
participant's own group mean*, so those estimates carry the group label by
construction. Using them in a correlation, a factor model, or a classifier would
be circular.

We quantified the leakage as the proportion of between-subject variance in the
posterior means explained by group, $R^2_{\text{group}}$:

| Parameter | With $g_i$ in the model | Refit without $g_i$ |
|---|---|---|
| $\mathrm{logit}\,h$ | 0.23 | 0.175 |
| $\log d$ | 0.71 | 0.185 |
| $\kappa$ | 0.63 | 0.105 |

For $\log d$ and $\kappa$ the estimates from the group-informed model are
dominated by the label. The model was therefore refit with $\boldsymbol{\beta}$
removed (`wcst_hmm/10_nogroup.R`; no divergences, $\hat{R}_{\max} = 1.006$,
minimum ESS 681), and **all individual-level analyses in the paper use the
group-free estimates** (`wcst_params_nogroup.csv`). Group effects are reported
from the group-informed model, where the contrast is a parameter rather than a
difference between shrunken point estimates:

| Parameter | $\beta$ | 90% CrI | $P(\text{sign})$ | $\beta/\sigma$ |
|---|---|---|---|---|
| $\mathrm{logit}\,h$ | $+0.40$ | $[+0.10, +0.72]$ | 0.983 | $+0.76$ |
| $\log d$ | $-0.14$ | $[-0.26, -0.03]$ | 0.981 | $-1.64$ |
| $\kappa$ | $-0.12$ | $[-0.32, +0.10]$ | 0.814 | $-0.85$ |

## S6. Individual-level reliability of the three parameters, and its consequences

For each parameter we compared the between-subject variance of the posterior
means, $\mathrm{Var}_i(\hat\theta_i)$, with the mean posterior variance within
subject, $\overline{\mathrm{SD}^2_i}$. The ratio

$$
\rho_{\text{ind}} = \frac{\mathrm{Var}_i(\hat\theta_i)}{\mathrm{Var}_i(\hat\theta_i) + \overline{\mathrm{SD}^2_i}}
$$

is the share of observed variation attributable to true between-subject
differences (computed on the $n = 80$ participants with complete data in all
three tasks):

| Parameter | $\mathrm{Var}_i(\hat\theta_i)$ | $\overline{\mathrm{SD}^2_i}$ | $\rho_{\text{ind}}$ | Used for individual differences |
|---|---|---|---|---|
| $\mathrm{logit}\,h$ | 0.208 | 0.150 | 0.58 | yes |
| $\log d$ | 0.006 | 0.012 | 0.32 | no |
| $\kappa$ | 0.003 | 0.035 | 0.09 | no |

This is consistent with the recovery result for $\kappa$ (S4) and extends it to
$\log d$. Accordingly:

- **Structural analyses** (correlations, variance-components model) use
  $\mathrm{logit}\,h$ only, as the single WCST indicator.
- **Group-effect analyses** retain all three parameters. Per-subject estimation
  noise is centred and independent of group, so it *attenuates* a group contrast
  rather than inflating it; the reported effect sizes for $\log d$ and $\kappa$
  are conservative.
- **Classification** retains all three. A noisy predictor simply predicts less
  well, and cross-validated AUC measures exactly that; no validity assumption is
  violated.

## S7. Files

| File | Content |
|---|---|
| `wcst_hmm/01_build_stan_data.R` | rebuild of trial-level data, choice-coding correction |
| `wcst_hmm/02_fit_models.R` | fitting of the three candidate HMM/RL models |
| `wcst_hmm/03_loo_compare.R` | PSIS-LOO comparison |
| `wcst_hmm/04_reliability.R` | split-half reliability |
| `wcst_hmm/05_recovery.R`, `funs_simulate.R` | simulation and parameter recovery |
| `wcst_hmm/06_ppc.R` | posterior predictive checks |
| `wcst_hmm/07_external.R` | convergent validity with conventional indices |
| `wcst_hmm/08_sticky.R` | final model with perseveration; LOO, PPC, recovery |
| `wcst_hmm/09_external_sticky.R` | convergent and cross-task validity, final model |
| `wcst_hmm/10_nogroup.R` | refit without the group predictor; leakage quantification |
| `wcst_hmm/stan/hmm_sticky.stan` | selected model |
| `wcst_hmm/stan/hmm_hier.stan`, `hmm_hier4.stan`, `rw_dim_hier.stan` | comparison models |
