# Introduction

Anorexia nervosa (AN) has among the highest mortality and the most persistent
course of any psychiatric disorder, and treatment response in adults remains
modest [@treasure2020; @zipfel2015]. This has motivated a long search for
cognitive mechanisms that might explain why restrictive eating, once
established, is so hard to give up. The most influential candidate is
*cognitive inflexibility*: a difficulty in updating behaviour when the
environment or the rewarded rule changes. Two meta-analyses of set-shifting in
eating disorders converge on a small-to-moderate deficit in AN
[@roberts2007; @wu2014], the Wisconsin Card Sorting Test (WCST) and the Trail
Making Test show group differences in both acute and weight-restored patients
[@tchanturia2004; @tchanturia2012; @westwood2016], and the deficit has been
reported in adolescents as well as adults [@fitzpatrick2012]. The account has
been consequential: cognitive remediation therapy was designed specifically to
target it [@tchanturia2014; @brockmeyer2018], and inflexibility features in
neurobiological models of how restriction becomes habitual and self-maintaining
[@steinglass2016; @foerde2015].

The evidence base for that account, however, rests almost entirely on
*model-free* summary scores — perseverative errors, categories completed,
switch costs in mean reaction time. Such scores are informative about
performance but ambiguous about process. A high perseverative-error count can
arise because a participant fails to infer that the rule has changed, because
they infer it but respond with insufficient determinism, because they
perseverate on a previously rewarded dimension, or because they respond more
cautiously and therefore differently under time pressure. Computational
phenotyping addresses exactly this ambiguity by fitting a generative model of
the trial-by-trial data and using its parameters as the individual-level
measures [@huys2016; @montague2012; @patzelt2018]. For the paradigms relevant
here the models are well developed: sequential-sampling models decompose choices
and response times into decision threshold, quality of evidence accumulation
and non-decision time [@ratcliff2008; @wiecki2013]; reinforcement-learning
models separate the rate at which feedback updates value from the policy that
converts value into choice [@ahn2017]; and card-sorting performance has been
modelled both as rule inference in a hidden-Markov process and as
dimension-level reinforcement learning [@bishara2010; @steinke2020]. Applied to
AN, such models have already localised group differences to specific
components — altered feedback-learning signals [@bernardoni2018; @ritschel2017],
value-based decision-making that varies with metabolic state
[@bernardoni2020], and category-learning deficits that survive weight
restoration [@filoteo2014].

Decomposing performance into parameters, however, brings into the open an
assumption that summary scores leave implicit. Calling inflexibility a *trait*
of AN — a stable characteristic of patients that a therapy could target — asserts
something about individuals, not only about group means: it asserts that the
same computational primitive, measured in different situations, identifies the
same people. Whether cognitive measures behave that way is far from settled.
The executive-function literature has documented for two decades that
task-specific variance dominates over the shared component (the "unity and
diversity" structure), and that the shared component is not recoverable from any
single task [@miyake2012; @friedman2017; @karr2018; @snyder2015]. Tasks that
produce large, replicable experimental effects can be poor instruments for
individual differences, because the experimental designs that maximise
within-subject effects minimise between-subject variance — the *reliability
paradox* [@hedge2018; @rouder2019]. Large-scale test–retest work on
self-regulation measures found modest reliabilities and a factor structure that
followed task format more than construct [@enkavi2019; @eisenberg2019].
Reliability estimates for computational parameters specifically are heterogeneous
and often lower than the analyses built on them require
[@waltmann2022; @mkrtchian2023; @karvelis2023; @sullivantoole2022], and the
meaning of a fitted parameter has been shown to depend on the task context in
which it was estimated [@eckstein2022]. There is, in addition, no logical
guarantee that a difference between group means says anything about the
covariance structure within individuals: group-to-individual generalisability
must be demonstrated rather than assumed [@fisher2018; @molenaar2009;
@kievit2013].

These two literatures rarely meet. Studies of AN typically administer one
flexibility task and interpret a group difference as evidence about a patient
characteristic; studies of measurement structure typically use healthy samples
and model-free indices. As a result, the central claim of the inflexibility
account — that the same individuals are inflexible whichever way inflexibility is
measured — has not, to our knowledge, been tested in AN with process-level
measures. The question is not whether the group difference is real, but whether
it licenses the inference about individuals that clinical reasoning draws from
it. That distinction also matters for how heterogeneity is understood. Patients
may differ from controls not by occupying a displaced but equally compact region
of parameter space, but by being *more dispersed* within it, in which case the
group mean describes no one well — the situation normative and subtyping
approaches were developed to handle [@marquand2016; @wolfers2018; @feczko2019;
@dinga2019].

We therefore estimated the same computational primitives in the same
participants across three tasks that all require adapting to changing
contingencies but differ in what has to change: a probabilistic reversal
learning task (PRL, feedback-driven value updating), a cued task-switching
paradigm (advance reconfiguration of a stimulus–response mapping
[@monsell2003; @kiesel2010; @rogers1995]), and the WCST (inference over an
unsignalled sorting rule). This yields, for 80 participants with complete data
on all three tasks (35 with AN, 45 healthy controls), 14 individual-level
parameters. Four primitives are measured in at least two of the three tasks —
decision threshold, quality of evidence accumulation, non-decision time, and
learning or flexibility — which is what makes convergence testable; two further
parameters, decision determinism and dimension-level perseveration, are specific
to the WCST. Sequential-sampling parameters were
available from prior modelling of the PRL and task-switching data; for the WCST
we fitted a hierarchical hidden-Markov model of rule inference with
dimension-level perseveration, selected by leave-one-out cross-validation
against reinforcement-learning and non-perseverative alternatives
[@vehtari2017; @carpenter2017].

The design supports four questions that are usually asked in separate papers. First,
where do the group differences fall once performance is decomposed into
primitives? Second, does the cross-task structure of the parameters support a
shared flexibility trait — that is, do parameters measuring the same primitive in
different tasks converge, as a multitrait–multimethod analysis requires
[@campbell1959]? Third, do patients differ from controls in the dispersion of
their parameters as well as in their location? Fourth, how well can group
membership be predicted out of sample from parameters alone, and does combining
tasks help — a question that must be answered with nested cross-validation, since
selection inside the resampling loop is what separates an honest estimate from
an optimistic one [@yarkoni2017; @poldrack2020; @varoquaux2018; @chekroud2024].
Because parameters estimated hierarchically are shrunk toward the mean of
whatever population is specified, we estimated the WCST model twice — once with
and once without a group predictor — and used the group-free estimates for every
individual-level analysis, reserving the group-informed model for the group
contrast itself. Anticipating the results, the four answers do not point the same
way, and the disagreement between them is the finding we develop.
