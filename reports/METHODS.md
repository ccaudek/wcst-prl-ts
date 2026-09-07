# Methods

## Participants and design

Patients with anorexia nervosa and healthy controls completed three
computerised tasks: a probabilistic reversal learning task (PRL), a cued
task-switching task, and the Wisconsin Card Sorting Test (WCST). Sample sizes
differ across tasks because not every participant completed every task: PRL
$n = 94$ (44 patients, 50 controls), task switching $n = 96$ (46, 50), WCST
$n = 84$ (39, 45). Eighty participants (35 patients, 45 controls) completed all
three; this is the *complete-case* sample used for all analyses that require
data from more than one task.

## Computational models and individual-level parameters

Each task was modelled with a generative model of trial-level behaviour, and
each participant characterised by the posterior means of the model's
subject-level parameters. All models were fit in Stan by Hamiltonian Monte
Carlo with hierarchical (partial-pooling) priors.

**Probabilistic reversal learning.** A reinforcement-learning drift-diffusion
model, yielding five subject-level parameters: boundary separation $a$, drift
rate $v$, non-decision time $t_0$, and two learning rates on the logit scale —
$\alpha$ (overall) and $\alpha^{+}$ (asymmetry for positive outcomes).

**Task switching.** A drift-diffusion model with condition-specific
parameters, yielding six subject-level parameters: $a$, $v$, and $t_0$
separately for repetition and switch trials ($a_0, a_1, v_0, v_1, t_{0,0},
t_{0,1}$).

**WCST.** A hierarchical hidden-Markov model of rule inference with
dimension-level perseveration, yielding three subject-level parameters: the
belief switch rate $\mathrm{logit}\,h$, the decision determinism $\log d$, and
the perseveration bonus $\kappa$. Model specification, selection against three
alternatives by PSIS-LOO, posterior predictive checks, parameter recovery and
convergent validity are reported in Supplementary Material S2–S4.

The PRL and task-switching parameters were estimated in prior work and are used
here as provided, without re-estimation. The WCST parameters were estimated for
the present study.

### Group-free estimates for individual-level analyses

The WCST model includes the group indicator in the subject-level mean, which is
the appropriate specification for estimating the group contrast but shrinks each
participant's estimate toward that participant's own group mean. In these data
group explained 23%, 71% and 63% of the between-subject variance of the
posterior means of $\mathrm{logit}\,h$, $\log d$ and $\kappa$ respectively.
Using such estimates in a correlation, factor model or classifier would be
circular. The model was therefore refit without the group predictor, and **all
individual-level analyses use these group-free estimates** (group explains 18%,
19% and 11% respectively). Group effects on the WCST parameters are reported
from the group-informed model, where the contrast is an explicit parameter with
its own posterior (Supplementary S5).

### Selection of indicators by individual-level reliability

For each parameter we compared the between-subject variance of the posterior
means with the mean within-subject posterior variance, giving the share of
observed variation attributable to true between-subject differences (see
Supplementary S6). All eleven PRL and task-switching parameters were treated as
reliable. Among the WCST parameters only $\mathrm{logit}\,h$ was
($\rho_{\text{ind}} = 0.58$); $\log d$ (0.32) and $\kappa$ (0.09) were not, a
result consistent with their poor parameter recovery.

Indicator sets were fixed a priori by analysis type and are stated with each
analysis below:

- **12 reliable indicators** (11 DDM parameters + $\mathrm{logit}\,h$) for
  correlational and structural analyses, where estimation noise attenuates
  associations and therefore biases the result toward the conclusion we draw.
- **all 14 parameters** for group-effect analyses and for classification, where
  centred estimation noise is conservative: it attenuates group contrasts, and a
  noisy predictor simply predicts less well, which cross-validation measures
  directly.

## Statistical analyses

All analyses were performed in R. Scripts are provided in the accompanying
package; each analysis below names its script.

### Group differences in each parameter (`02_effect_size.R`)

For each of the 14 parameters we computed Cohen's $d$ (patients − controls,
pooled SD) with 95% percentile bootstrap confidence intervals (10 000
resamples, resampling within group), and a bootstrap two-sided $p$-value.
Effect sizes are reported both for each task's maximal sample and for the
complete-case sample.

**Sign-concordance test.** Under the hypothesis that a computational primitive
is altered coherently across tasks, all indicators of the same primitive should
share the sign of their group difference. Parameters were grouped into four
primitives (threshold, drift rate, non-decision time, learning/flexibility);
within each primitive we counted indicators agreeing with the primitive's modal
sign and tested the total against chance by a one-sided binomial test. The
primary test used the 11 DDM parameters, whose scales share an orientation; a
sensitivity version added WCST $\log d$ to the threshold primitive. Because the
modal sign is read from the data, this test is mildly anti-conservative and is
reported as a descriptive summary of coherence rather than as a confirmatory
test.

### Cross-task correlation structure (`03_correlazioni.R`)

Spearman correlations among the 12 reliable indicators were computed on the
complete-case sample, after removing the group difference from each indicator
(within-group centring), so that associations reflect individual differences and
not group separation. The 66 pairs were classified into four families — same
primitive within task, same primitive across tasks, different primitives within
task, different primitives across tasks — and the distributions of $|\rho|$
compared across families. The prediction under a shared-trait account is that
*same primitive across tasks* should stand out from *different primitives across
tasks*.

Because a correlation is attenuated by the unreliability of both variables, the
cross-task correlations of $\mathrm{logit}\,h$ were additionally reported after
correction for attenuation [@spearman1904], dividing by
$\sqrt{\rho_{\text{ind}}}$ with $\rho_{\text{ind}}$ the individual-level
reliability of $h$ defined above and the DDM parameters treated as perfectly
reliable, which makes the correction a lower bound on the disattenuated value.

### Variance-components model (`04_mtmm.R`, `stan/mtmm.stan`)

The 12 reliable indicators were modelled jointly as a multitrait–multimethod
structure in Stan. Each standardised indicator $y_{ij}$ loads on a
cross-task **trait** factor for its primitive, a task-specific **method**
factor, and a **primitive-within-task** factor, plus indicator-specific
residual variance:

$$
y_{ij} = \mu_j + \beta_j g_i + \lambda^{\text{T}}_j\, \tau_{i,\,\mathrm{prim}(j)} + \lambda^{\text{M}}_j\, m_{i,\,\mathrm{task}(j)} + \lambda^{\text{W}}_j\, w_{i,\,\mathrm{prim}(j),\,\mathrm{task}(j)} + \varepsilon_{ij}
$$

with all factors standard normal and independent. Known measurement error was
incorporated for the WCST indicator by adding its per-participant posterior
variance to the diagonal of the implied covariance matrix, so its loading is not
attenuated by estimation noise. The likelihood is the implied multivariate
normal, marginalising the latent factors.

Five nested variants were compared by PSIS-LOO: the full model; *no trait*
(cross-task factors removed); *no method*; *within-task only*; and
*independence*. Variance shares (trait, method, within-task, unique) were
computed per indicator from the posterior of the loadings, and the contrast
between total trait and total method variance was summarised by its posterior
distribution.

As a check on this Bayesian formulation we also fit three classical
maximum-likelihood confirmatory factor models in `lavaan`
(`07_cfa_classica.R`): one factor per primitive (cross-task traits), one factor
per task, and a single general factor. Convergence and standard fit indices
($\chi^2$, CFI, TLI, RMSEA, SRMR) are reported.

### Group heterogeneity (`05_eterogeneita.R`)

To test whether patients are more variable — not merely different on average —
we used the 12 reliable indicators standardised on the pooled sample.
Multivariate dispersion of a group was defined as the mean Euclidean distance of
its members from their own group centroid. The observed between-group difference
in dispersion was tested by a permutation test that shuffles group labels
(10 000 permutations), recomputing centroids within each permutation.
As specificity checks, the same test was run on the 11 PRL and task-switching
indicators alone and on a 13-indicator set adding WCST $\log d$, and a
leave-one-indicator-out analysis repeated the test dropping each indicator in
turn.

Centroid separation (Euclidean distance between group centroids in the same
12-dimensional standardised space) was tested by the same permutation scheme, to
confirm that a mean difference is present alongside the dispersion difference.

Per-indicator variance ratios (patients/controls) were computed with
Brown–Forsythe tests of equality of variance and Benjamini–Hochberg correction
across the 12 indicators.

For display, the 12-dimensional space was projected onto its first two principal
components (Figure 3a); the tests were performed in the full space.

### Classification (`06_auc.R`)

Group membership was predicted from the parameters by $L^2$-penalised logistic
regression (`glmnet`, $\alpha = 0$) on the complete-case sample, with all 14
parameters available. Seven predictor sets were compared: each task alone, each
pair of tasks, and all three.

Discrimination was estimated by **nested** cross-validation: an outer
stratified 5-fold loop repeated 20 times (100 folds), with the penalty
$\lambda$ and the number of retained predictors $k$ (univariate filter,
$k \in \{3, 5, \infty\}$) selected *inside* each training fold by an inner
stratified 5-fold cross-validation. Standardisation, filtering and penalty
selection were all refit within each training fold. AUC was computed on the
pooled out-of-fold predictions of each repetition and summarised as the mean
across the 100 folds, with the 5th–95th percentile range of fold-level AUCs.
Outer folds were generated once and reused across all seven predictor sets, so
that set-to-set comparisons are paired fold by fold; comparisons are reported as
the mean fold-level difference with its 5th–95th percentile range and the
proportion of folds favouring the richer set.

To quantify the optimism introduced by selecting hyperparameters outside the
cross-validation loop, the same procedure was repeated with $\lambda$ and $k$
chosen once on the full sample; the difference between the non-nested and nested
AUC is reported.

### Convergent validity of the WCST parameters (`08_validita_convergente.R`)

The belief switch rate was correlated (Spearman) with the conventional
model-free WCST indices, using the group-free estimates, both on the full WCST
sample and after centring both variables within group — the latter being the
form of the check that cannot be inflated by the group difference being present
in both terms. A fuller set of behavioural signatures, computed from the
trial-level data against the group-informed estimates, is reported in
Supplementary S4.

### Software

R 4.x with `cmdstanr`, `loo`, `glmnet` and `lavaan`; Stan 2.x. AUC was computed directly
from the rank statistic of out-of-fold predicted probabilities. Random seeds are
set in each script.
