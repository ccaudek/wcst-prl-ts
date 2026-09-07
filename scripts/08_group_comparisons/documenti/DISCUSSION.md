# Discussion

We estimated the same computational primitives in the same participants across
three tasks that all require adapting to changing contingencies, and asked
whether the group-level and individual-level pictures agree. They do not. At the
group level the picture is coherent and consistent with two decades of
set-shifting research in anorexia nervosa: patients differed from controls on
parameters drawn from every primitive and every task, and the directions of
those differences did not contradict each other — higher decision thresholds,
slower evidence accumulation, slower learning from feedback, less deterministic
translation of belief into choice. At the individual level that coherence
disappears. Parameters measuring the same primitive in different tasks were
barely more strongly correlated (mean $|\rho| = 0.17$) than parameters measuring
different primitives in different tasks (0.11), while parameters measured within
the same task were clearly associated (0.37). Removing the cross-task trait
factors from a variance-components model cost nothing predictively, whereas
removing the task-specific factors cost $11.0 \pm 6.0$ in expected log
predictive density, and a classical factor model with cross-task trait factors
was not even identifiable in these data. What reconciles the two levels is
dispersion: patients were both displaced from the control centroid and more
scattered around their own (mean distance 3.55 vs 2.83, $p = 0.004$), so
averaging over patients recovers a coherent direction of departure while
correlating across patients does not recover a coherent dimension.

## Group differences do not license the inference usually drawn from them

The clinical reading of the set-shifting literature — that inflexibility is a
characteristic of patients, and therefore a plausible treatment target
[@tchanturia2014; @brockmeyer2018] — is an assertion about individuals. Our
group-level results reproduce the empirical basis of that reading; our
individual-level results show that the basis does not carry the assertion.
A difference between group means constrains the covariance structure within
groups only weakly, and generalisation from one level to the other has to be
demonstrated rather than assumed [@fisher2018; @molenaar2009; @kievit2013].
Concretely: knowing that a patient is impaired on one flexibility task told us
little about whether the same patient was impaired on another, even when the two
tasks were scored with the same computational primitive rather than with
task-specific summary indices. The construct behaved as a property of the
group, not of the person.

The point is not confined to anorexia nervosa. Programmes that seek dimensional,
mechanism-level descriptions cutting across diagnostic categories
[@insel2010; @gillan2016] depend on the assumption we tested and did not
confirm: that a mechanistic measure ranks individuals consistently enough for
the dimension to exist at the level of the person. In eating disorders,
compulsivity and flexibility measures have been proposed as exactly such
transdiagnostic candidates [@godier2014]; our results suggest the candidacy
should be established with multi-task measurement before it is assumed.

This is not an argument that the group differences are uninteresting or
artefactual. Their internal consistency is striking — within each primitive,
every indicator carried the same sign — and it suggests a global rather than a
selective alteration. Higher thresholds together with lower drift rates across
two independent tasks describe a decision policy that trades speed for caution
while the quality of the evidence entering the decision is lower
[@ratcliff2008]. That combination is closer to a general change in how
decisions are made under uncertainty than to a specific deficit in shifting,
and it is a more parsimonious summary of our group-level data than
"inflexibility".

## Two readings of the missing convergence, in order of plausibility

The first reading is about the construct. Cognitive flexibility, measured at the
process level, may simply not be a unitary individual-differences dimension.
This is the expected result from outside the eating-disorders field: the
executive-function literature has repeatedly found that task-specific variance
dominates the shared component [@miyake2012; @friedman2017; @karr2018;
@snyder2015], large-scale work on self-regulation measures recovered a structure
organised more by task format than by construct [@enkavi2019; @eisenberg2019],
and the interpretation of a fitted computational parameter has been shown to
depend on the task context in which it was estimated [@eckstein2022]. Our data
add the observation that this holds within a clinical sample where the
group-level construct is well supported, and that it holds for model-based
parameters, not only for model-free scores — decomposing performance into
primitives did not create the cross-task convergence that summary scores lack.

The second reading is about measurement. Cross-task correlations are attenuated
by the unreliability of both terms, and the individual-level reliability of
computational parameters is frequently lower than the analyses built on them
require [@waltmann2022; @mkrtchian2023; @karvelis2023; @sullivantoole2022] —
partly for the reason identified by the reliability paradox, that designs
optimised for within-subject effects suppress between-subject variance
[@hedge2018; @rouder2019]. This applies directly to our WCST parameters: the
ratio of between-subject variance to posterior estimation error was 0.58 for the
belief switch rate, 0.33 for decision determinism and 0.09 for perseveration,
which is why only the switch rate entered the structural and correlational
analyses. A shared trait could in principle exist and be hidden by this noise.

We order these two readings as above, rather than treating them as equally
open, for three reasons internal to our data. First, disattenuation does not
rescue the trait: correcting the cross-task correlations of the WCST switch rate
for its own unreliability raised the largest of them to $-0.47$ and left the
ordering unchanged, so the associations are weak rather than merely noisy.
Second, and more decisively, the within-task correlations of 0.37 were computed
from the *same* indicators as the cross-task correlations of 0.17. Unreliability
attenuates both equally; it cannot produce an asymmetry between them. Whatever
is limiting the cross-task associations is not shared measurement error but the
absence of shared signal. Third, the trait factors were not merely small but
unnecessary — the model without them was marginally preferred — and the
trait–method contrast was centred near zero ($+0.037$, $[-0.082, +0.162]$).
That said, the measurement reading is not eliminated. Trait variance was
estimated at 0.180 $[0.087, 0.279]$, which is not zero; the credible intervals
on all four per-primitive trait estimates were wide; we had a single
administration and therefore no test–retest information; and the
non-identifiability of the trait-structured factor model is at least partly a
consequence of fitting 12 indicators in 80 participants. A shared dimension that
is real but small, or real but poorly measured by these tasks, remains
compatible with what we observed.

## Heterogeneity as a substantive finding

The dispersion result deserves emphasis in its own right, because it is what
makes the two levels coexist rather than contradict. Patients were not a
displaced but compact group; they were a scattered one, and no single parameter
carried the scatter — all 12 indicators showed greater patient variability, none
individually surviving correction. In other words, patients departed from the
control profile in a common direction but through different combinations of
parameters. The mean patient profile is therefore a composite that may describe
few actual patients, which is precisely the situation that motivated normative
and subtyping approaches in psychiatry more broadly [@marquand2016;
@wolfers2018; @feczko2019], with the caveat that clustering such profiles into
discrete subtypes has often proved unstable [@dinga2019].

For eating-disorders research this has a practical consequence. If the
mechanistic departure differs across patients, an intervention aimed at one
mechanism will show a small average effect even if it works well for the subset
whose departure lies in that mechanism — a possible contributor to the modest
average effects reported for cognitive remediation [@tchanturia2014;
@brockmeyer2018]. We did not test this: we have no treatment data, and the
inference from parameter dispersion to differential treatment response is a
hypothesis, not a result. But it is a testable one, and it requires exactly the
kind of individual-level measurement whose limitations we document here.

## Model-based and model-free perseveration are not the same thing

One result cuts against the standard interpretation of the WCST in this
population and should be flagged for replication. Patients showed a *higher*
belief switch rate ($d = +1.02$) and *lower* dimension-level perseveration
($\kappa$: $d = -0.77$) — that is, the model located their difficulty in failing
to maintain and exploit a rule, not in failing to abandon one. Meanwhile the
switch rate correlated positively with both conventional error counts across
individuals ($\rho = +0.75$ with non-perseverative errors, $+0.63$ with
perseverative errors), while the proportion of perseverative *responses* was
essentially unrelated to it. An elevated perseverative-error count is therefore
compatible with excessive rather than insufficient belief switching, because
both error categories increase when rule representations are unstable. Since
perseverative errors are the index on which much of the AN flexibility
literature rests [@tchanturia2004; @tchanturia2012; @roberts2007], the
possibility that they do not measure perseveration in the mechanistic sense is
consequential — and it is the kind of ambiguity that model-based scoring exists
to resolve [@bishara2010; @steinke2020]. We would not build a claim on a single
sample, but the pattern was internally consistent across the group contrast and
the individual-level correlations, and it was obtained with the group-free
estimates, so it is not an artefact of hierarchical shrinkage toward group
means.

## Complementarity, and what classification does and does not show

The three tasks were complementary as predictors: no single task exceeded
AUC 0.76 while the three together reached 0.84, and most of that gain was
already obtained with two tasks. This is the same message as the correlational
analysis, seen from the prediction side — if the tasks measured one construct,
combining them would add little. It is not a diagnostic claim. Diagnosis is not
the clinical problem in anorexia nervosa, fold-level intervals at $n = 80$ were
wide and none of the paired comparisons excluded zero, and cross-validated
performance at this sample size is both unstable and optimistic relative to what
a new sample would yield [@varoquaux2018; @poldrack2020; @chekroud2024]. We
report it as evidence about the information structure of the tasks, and we note
that selecting the penalty and predictor count outside the resampling loop
inflated AUC by up to 0.025 in our own data — a reminder of how easily such
values drift upward [@yarkoni2017].

## Limitations

Three limitations bound the individual-level conclusions. The sample of 80
complete cases supports the group contrasts and the dispersion tests, but it is
small for a 12-indicator factor model, and the failure of the trait-structured
specification to converge should be read in that light. We had one
administration per participant, so the reliability figures we report are
posterior-uncertainty ratios within a single session rather than test–retest
coefficients, and the disattenuation correction inherits that limitation; a
second administration would separate "no trait" from "no reliable measurement of
a trait" far more sharply than any analysis of a single session can. And the
sequential-sampling parameters for the reversal-learning and task-switching data
were taken from prior modelling of those datasets and used as provided, so their
estimation choices are not under our control — and the precision of
sequential-sampling parameters depends strongly on the number of trials per
participant and on which across-trial variability components are estimated
[@lerche2016; @boehm2018]. Only the WCST model was fitted
here, with model comparison, parameter recovery and posterior predictive checks
reported in the supplementary material.

Two further constraints limit interpretation. Clinical covariates — body mass
index, illness duration, diagnostic subtype, comorbidity, medication — were not
available, so we could not test whether the dispersion we observe tracks
clinical variation, which is the most obvious next question and the one that
would turn heterogeneity from a nuisance into a phenotype. Patients were
assessed in the acute state, and starvation-related effects on decision
parameters cannot be separated from trait-like differences in a
cross-sectional design [@bernardoni2020]. Finally, we tested convergence in one
sample; the asymmetry between within-task and cross-task structure is the kind
of result that needs replication in an independent sample before it constrains
theory.

## Conclusion

Decomposing three flexibility tasks into computational primitives reproduced the
group-level finding that patients with anorexia nervosa differ from controls,
and localised the difference to a globally more cautious and less efficient
decision policy rather than to a selective shifting deficit. It also showed that
this group-level coherence does not correspond to an individual-level trait:
the same primitive measured in two tasks does not identify the same people, and
patients are heterogeneous in which parameters carry their departure from the
control profile. The most plausible reading is that flexibility, at the process
level, is task-specific rather than trait-like, though limited reliability of
individual parameter estimates cannot be ruled out as a contributor. For a field
that treats cognitive inflexibility as a patient characteristic and a treatment
target, the practical implication is measurement-level and immediate: claims
about individual patients require multi-task estimation and explicit reliability
reporting, and single-task evidence supports claims about groups only.
