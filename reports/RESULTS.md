# Results

## 1. Group differences are present in every computational primitive

Patients and controls differed on parameters drawn from all four primitives, in
all three tasks (Figure 1; Table 1). In the complete-case sample the largest
effects were WCST decision determinism ($\log d$: $d = -1.08$, 95% CI
$[-1.59, -0.65]$), the WCST belief switch rate ($\mathrm{logit}\,h$:
$d = +1.02$, $[0.59, 1.53]$), PRL boundary separation ($a$: $d = +0.84$,
$[0.41, 1.34]$), task-switching drift rate on repetition trials
($v_0$: $d = -0.85$, $[-1.36, -0.42]$), WCST perseveration
($\kappa$: $d = -0.77$, $[-1.27, -0.33]$) and the PRL learning rate
($\alpha$: $d = -0.71$, $[-1.18, -0.28]$).

The direction of these effects is internally consistent. Patients set higher
decision thresholds (PRL $a$, task-switching $a_1$), accumulated evidence more
slowly (task-switching $v_0$), learned more slowly from feedback (PRL $\alpha$,
$\alpha^{+}$), and discarded a currently-held rule more readily under
disconfirming feedback (WCST $h$) while being less deterministic in translating
belief into choice (WCST $\log d$).

Within each of the four primitives, every indicator carried the same sign as the
primitive's modal direction: 11 of 11 for the DDM parameters common to PRL and
task switching (one-sided binomial $p = 0.0005$), and 11 of 12 when WCST
$\log d$ is added to the threshold primitive ($p = 0.003$). Because the modal
sign is read from the data this summary is anti-conservative, but it makes the
qualitative point: the group differences do not point in inconsistent
directions across tasks.

Effect sizes computed on each task's maximal sample were close to those in the
complete-case sample (all differences within 0.15 of $d$), so the pattern is not
an artefact of restricting to complete cases.

## 2. No individual-level trait is shared across tasks

If a group difference in, say, the drift rate reflected a single altered
process, participants with a low drift rate in one task should have a low drift
rate in the other. They do not.

After removing the group means, Spearman correlations among the 12 reliable
indicators showed the expected within-task structure and no cross-task
structure (Figure 2; Table 2). Mean $|\rho|$ by pair family:

| Pair family | $n$ pairs | mean $|\rho|$ | max $|\rho|$ | mean $\rho$ |
|---|---|---|---|---|
| same primitive, within task | 4 | 0.37 | 0.54 | +0.37 |
| different primitives, within task | 21 | 0.19 | 0.48 | −0.01 |
| same primitive, across tasks | 6 | 0.17 | 0.29 | +0.11 |
| different primitives, across tasks | 35 | 0.11 | 0.36 | 0.00 |

The critical comparison is between rows 3 and 4: pairs of indicators measuring
*the same* primitive in *different* tasks are barely more strongly correlated
(mean $|\rho| = 0.17$) than pairs measuring *different* primitives in different
tasks (0.11). By contrast, indicators within the same task are clearly
associated (0.37 for the same primitive, 0.19 even across primitives). What the
data contain is task-specific covariation, not primitive-specific covariation.

The variance-components model reaches the same conclusion (Table 3). Removing
the cross-task trait factors did not worsen predictive accuracy — the *no
trait* model was in fact marginally preferred ($\Delta\mathrm{elpd} = +1.5 \pm
2.4$ in its favour, i.e. indistinguishable) — whereas removing the
task-specific method factors cost $11.0 \pm 6.0$, keeping only the
primitive-within-task factors cost $10.6 \pm 5.9$, and removing all structure
cost $48.9 \pm 11.0$. Posterior variance shares were: trait 0.180
$[0.087, 0.279]$, method 0.143 $[0.085, 0.205]$, primitive-within-task 0.200
$[0.093, 0.295]$. The trait−method contrast was $+0.037$ $[-0.082, +0.162]$
($P(\text{trait} > \text{method}) = 0.71$). Trait variance was of similar,
modest magnitude for all four primitives (threshold 0.19, drift 0.17,
non-decision time 0.15, learning/flexibility 0.21), with wide intervals in every
case and none approaching the value a shared dimension would imply.

A classical maximum-likelihood confirmatory factor analysis of the same 12
indicators pointed the same way. The model with one factor per primitive
(cross-task traits) failed to converge, whereas the model with one factor per
task converged (though with poor absolute fit: $\chi^2(43) = 140.5$, CFI 0.51,
RMSEA 0.168, SRMR 0.137) and so did a single-factor model (CFI 0.33, RMSEA
0.189). The trait-structured solution is not merely weakly supported; with these
data it is not identifiable.

A cross-task trait is therefore not excluded, but it is not needed to describe
the data and it is no larger than the task-specific component. The
group-level coherence of Section 1 does not correspond to a coherent
individual-level dimension.

## 3. Patients are more heterogeneous, not uniformly shifted

The two results above are reconciled by the dispersion of the patient group
(Figure 3). In the 12-dimensional standardised parameter space, the mean
distance of a participant from their own group centroid was 3.55 for patients
and 2.83 for controls, a difference of 0.72 (permutation $p = 0.004$). Group
centroids were also separated (distance 1.86, permutation $p < 0.001$), so this
is added heterogeneity on top of a mean shift, not instead of it.

The effect was robust to the choice of indicators. Restricting to the 11 PRL and
task-switching indicators gave a difference of 0.64 ($p = 0.010$); adding WCST
$\log d$ gave 0.76 ($p = 0.003$); and dropping each indicator in turn left the
permutation $p$ below 0.05 in 12 of 12 cases.

All 12 indicators had a variance ratio above 1 in the direction of greater
patient variability, though no single indicator survived multiple-comparison
correction on its own (smallest Brown–Forsythe $p_{\text{BH}} = 0.152$). The
heterogeneity is a distributed property of the parameter profile rather than the
signature of one parameter.

This provides a direct account of the pattern: patients depart from the control
profile in a coherent *direction* at the group level, while differing among
themselves in *which* parameters carry the departure. Averaging over patients
recovers the coherence; correlating across individuals does not.

## 4. The three tasks contribute complementary information

If the tasks measured one underlying construct, adding tasks to a classifier
would add little. They do not behave that way (Figure 4; Table 4).

Under nested cross-validation (5 folds × 20 repetitions; hyperparameters
selected inside each training fold), area under the ROC curve for
patient/control discrimination was:

| Predictor set | $n$ parameters | AUC | 5th–95th pct of folds |
|---|---|---|---|
| all three tasks | 14 | 0.842 | 0.674–1.000 |
| PRL + WCST | 8 | 0.823 | 0.650–0.984 |
| task switching + WCST | 9 | 0.800 | 0.593–0.961 |
| PRL + task switching | 11 | 0.774 | 0.531–0.968 |
| PRL alone | 5 | 0.753 | 0.523–0.953 |
| WCST alone | 3 | 0.753 | 0.556–0.905 |
| task switching alone | 6 | 0.691 | 0.467–0.898 |

No single task exceeded 0.76; the three together reached 0.84. PRL and WCST were
tied as the best single task (both 0.753), and the paired fold-level advantage
of the full set over each was $+0.089$ ($[-0.112, +0.262]$, favoured in 79% of
folds, against PRL; $[-0.104, +0.262]$, 77%, against WCST); over the weakest
single task it was $+0.152$ ($[-0.112, +0.350]$; 91% of folds). Against the best
two-task set (PRL + WCST) the advantage was small ($+0.020$, $[-0.111, +0.143]$;
59% of folds), so most of the gain is obtained once two tasks are combined. Fold-level intervals are wide, as
expected at $n = 80$, and none of the paired comparisons excludes zero; the
consistent ranking across folds rather than any individual interval is the
evidence for complementarity.

Selecting the penalty and the number of predictors on the full sample rather
than inside each training fold inflated AUC by 0.004 to 0.025 across the seven
sets (Figure 4b). The bias is modest here, but it is one-directional, and the
values reported above are the nested ones.

## 5. Convergent validity of the WCST parameters

The WCST switch rate $h$ behaved as expected against conventional model-free
scoring. Using the group-free estimates and the $n = 86$ participants with
model-free indices available, $\mathrm{logit}\,h$ correlated $\rho = +0.75$ with
the proportion of non-perseverative errors and $+0.63$ with the proportion of
perseverative errors (both $p < 0.001$); correlations computed within group
were essentially unchanged ($+0.70$ and $+0.63$), so the association is not
carried by the group difference. The proportion of perseverative *responses* was
only weakly related to $h$ ($-0.20$, $p = 0.072$; $-0.23$ within group),
as expected for a parameter that indexes rule maintenance rather than response
repetition — the latter is captured by $\kappa$. A fuller set of behavioural
signatures, computed from the trial-level data against the group-informed
estimates, is reported in Supplementary S4.

Cross-task associations of $h$ (group-free estimates, within-group centred,
$n = 80$) were weak and consistent with Section 2: the strongest was with
task-switching drift rate on repetition trials ($\rho = -0.36$), followed by
switch-trial drift ($-0.27$), PRL drift ($-0.22$), the PRL learning rate
($-0.19$) and PRL boundary separation ($+0.17$); the remaining six associations
were below 0.09 in absolute value. Correcting for the individual-level
unreliability of $h$ ($\rho_{\text{ind}} = 0.58$) raises the largest of these to
$-0.47$ and leaves the ordering unchanged, so the weakness of the cross-task
associations is not merely a measurement-noise artefact.

---

## Figures and tables

| Item | File | Content |
|---|---|---|
| Figure 1 | `fig1_effect_size.png` | Cohen's $d$ with bootstrap CIs for all 14 parameters, grouped by primitive |
| Figure 2 | `fig2_struttura.png` | correlation matrix of the 12 reliable indicators and $|\rho|$ by pair family |
| Figure 3 | `fig3_eterogeneita.png` | principal-component projection with dispersion ellipses; per-participant distances from centroid; per-indicator variance ratios |
| Figure 4 | `fig4_auc.png` | nested cross-validated AUC by predictor set; optimism from non-nested selection |
| Table 1 | `effect_sizes.csv` | group means, SDs, $d$, bootstrap CIs and $p$ per parameter, both samples |
| Table 2 | `struttura_cross_task.csv`, `coppie_correlazioni.csv` | correlations by pair family; all 66 pairs |
| Table 3 | `loo_mtmm.csv`, `decomposizione_varianza.csv`, `contrasto_tratto_metodo.csv`, `varianza_tratto_per_primitiva.csv` | model comparison and variance shares |
| Table 4 | `auc_nested.csv`, `auc_confronti.csv`, `auc_fold.csv` | AUC by predictor set, paired comparisons, fold-level values |
| — | `dispersione_test.csv`, `dispersione_loo.csv`, `centroidi_test.csv`, `eterogeneita.csv` | dispersion and centroid tests, robustness, per-indicator variance ratios |
| — | `affidabilita_wcst.csv`, `concordanza_segni.csv`, `correlazioni_h_wcst.csv` | reliability of WCST parameters, sign concordance, cross-task correlations of $h$ |
| — | `validita_convergente_wcst.csv`, `cfa_classica_fit.csv`, `cfa_classica_loadings.csv` | convergent validity of $h$; classical CFA fit and loadings |
