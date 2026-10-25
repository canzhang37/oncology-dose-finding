# Design: dose escalation with 3+3, BOIN and CRM, and the size of a randomized dose-optimization stage

This plan was written before any simulation was run. Changes made after results are seen go in the Amendments section at the end, with a reason, and are labeled in the readout.

## Questions

Part A. In a first-in-human oncology trial with six dose levels and a target dose-limiting toxicity (DLT) rate of 25%, how do the 3+3 design, the Bayesian optimal interval design (BOIN) and the continual reassessment method (CRM) compare? The comparison covers how often each picks the right dose, how many patients it treats at overly toxic doses, and how many patients it uses.

Part B. Dose-optimization guidance for oncology asks sponsors to compare more than one dose before the registration trial, often by randomizing patients between the maximum tolerated dose (MTD) and a lower dose. How many patients per arm does such a two-dose comparison need to choose the lower dose when response has plateaued and the higher dose when it has not?

## Part A settings

Common to all designs: K = 6 doses, target DLT rate 0.25, cohorts of 3, start at dose 1, no dose skipping when escalating, at most 36 patients. Each simulated trial draws one latent tolerance u ~ Uniform(0, 1) for each of 36 patient slots, and the i-th patient enrolled has a DLT at dose d if u < p_d. All three designs see the same patients in the same order, so their results are paired trial by trial.

3+3. Treat 3 at the current dose. With 0 of 3 DLTs, escalate; with 1 of 3, treat 3 more and escalate if 1 of 6; with 2 or more DLTs, stop escalating. After stopping at dose d, the MTD is the next lower dose if 6 patients have been treated there; if only 3 have, treat 3 more there and declare it the MTD if at most 1 of 6 had a DLT, otherwise move down again. At the top dose, 0 of 3 or 1 of 6 leads to 3 more patients if needed so that 6 have been treated, and the top dose is the MTD if at most 1 of 6 had a DLT. If dose 1 is too toxic, no dose is selected.

BOIN. Escalation and de-escalation boundaries from `BOIN::get.boundary` with the default phi1 = 0.6 x 0.25 and phi2 = 1.4 x 0.25. After each cohort, escalate if the observed DLT rate at the current dose is at or below the escalation boundary, de-escalate if it is at or above the de-escalation boundary, otherwise stay. A dose and all higher doses are eliminated when at least 3 patients have been treated there and P(p_d > 0.25 | data) > 0.95 under a Beta(1, 1) prior; the trial stops without a selection if dose 1 is eliminated. Twelve cohorts in all. The MTD is chosen with `BOIN::select.mtd` (isotonic regression).

CRM. One-parameter power model p_d = s_d ^ exp(beta), prior beta ~ Normal(0, 1.34), skeleton from `dfcrm::getprior(halfwidth = 0.05, target = 0.25, nu = 3, nlevel = 6)`. After each cohort, the posterior mean of beta is computed and the next dose is the one whose plug-in DLT rate is closest to 0.25, with two restrictions: no escalation by more than one level above the current dose, and no escalation right after a cohort with a DLT. Safety stop: the trial stops without a selection if at least 3 patients have been treated at dose 1 and P(p_1 > 0.25 | data) > 0.90. Twelve cohorts in all. The MTD is the dose, among those tried, whose plug-in DLT rate is closest to 0.25.

## Scenarios

Fixed scenarios (true DLT rates at doses 1 to 6). The correct dose is the one closest to 0.25.

| Scenario | p_1 to p_6 | Correct dose |
|---|---|---|
| S1 | 0.25 0.40 0.50 0.60 0.70 0.80 | 1 |
| S2 | 0.10 0.25 0.40 0.50 0.60 0.70 | 2 |
| S3 | 0.05 0.12 0.25 0.40 0.50 0.60 | 3 |
| S4 | 0.03 0.08 0.15 0.25 0.40 0.50 | 4 |
| S5 | 0.02 0.05 0.10 0.15 0.25 0.40 | 5 |
| S6 | 0.01 0.03 0.06 0.10 0.15 0.25 | 6 |
| S7 (steep jump) | 0.05 0.08 0.10 0.45 0.55 0.65 | 3 |
| S8 (flat) | 0.15 0.20 0.25 0.30 0.35 0.40 | 3 |
| S9 (all too toxic) | 0.45 0.55 0.65 0.75 0.85 0.90 | none: stop |

Random scenarios. 1,000 dose-toxicity curves from this generator: the position M of the correct dose is uniform on 1 to 6; p_M ~ Uniform(0.20, 0.30); the doses next to it differ from p_M by Uniform(0.12, 0.25); further steps are Uniform(0.03, 0.12) below M and Uniform(0.05, 0.20) above it; rates are kept within 0.01 and 0.95. Each random scenario gets 200 simulated trials.

Each fixed scenario gets 4,000 simulated trials per design, so the Monte Carlo standard error of a selection percentage is at most 0.8 points. Seeds are fixed.

## Metrics

1. Correct selection: percentage of trials that select the correct dose (for S9, that stop without a selection).
2. Overdose selection: percentage of trials that select a dose with p > 0.35.
3. Patients at overdoses: average percentage of each trial's patients treated at doses with p > 0.35.
4. Patients at the correct dose: average percentage.
5. Under-dosing: percentage of trials that select a dose below the correct dose.
6. Average number of patients, and percentage of trials stopped early for toxicity.

## Checks of the code

- 3+3: selection probabilities and the expected number of patients computed exactly by enumerating every path of the design; the simulation must agree with them within 3 Monte Carlo standard errors in every fixed scenario.
- BOIN: the boundaries must equal `BOIN::get.boundary`, and the percentage selecting each dose must agree with `BOIN::get.oc` (its own random numbers) within 3 combined standard errors in every fixed scenario.
- CRM: for 200 random data sets, the posterior mean of beta must agree with `dfcrm::crm` (method "bayes", model "empiric") to 1e-4, and the recommended dose must be the same. Over 500 simulated trials, a slow version that calls `dfcrm::crm` after every cohort must reach the same selection as the fast version in every trial.

## Predictions for Part A

P1. Averaged over S1 to S8, BOIN and CRM each select the correct dose at least 10 points more often than 3+3.
P2. BOIN and CRM differ by less than 5 points in average correct selection, both over S1 to S8 and over the random scenarios.
P3. Averaged over S1 to S8, 3+3 treats a smaller percentage of patients at overdoses than BOIN and CRM, but it selects a dose below the correct one in at least 30% of trials averaged over S2 to S6.
P4. Over the random scenarios, 3+3 has the lowest average correct selection.

## Part B

Two doses, a lower dose L and a higher dose H, randomized 1:1 with n patients per arm. The endpoint is objective response. The rule selects L if its observed response rate is at least that of H minus a margin delta = 0.10, and selects H otherwise. Selection probabilities are computed exactly from the two binomial distributions; no simulation is needed. Main setting: response rate 0.40 at H and differences 0, 0.05, 0.10, 0.15, 0.20 and 0.25 between H and L, with n from 10 to 100 per arm. Targets: P(select L) >= 0.80 when the difference is 0, and P(select H) >= 0.80 when it is 0.20. Sensitivity: response rate at H of 0.25 and 0.60, and margins of 0.05 and 0.15. An exact check against 100,000 simulated comparisons at three values of n is included.

P5. In the main setting, the smallest n per arm that meets both targets lies between 30 and 50.

## Limits stated in advance

Real escalation decisions also use pharmacokinetics, lower-grade toxicities and late-onset DLTs, and dose optimization weighs safety and patient-reported outcomes as well as response. The scenarios are stylized. The comparison asks how the decision rules behave, not which design a given trial should use.

## Amendments

Amendment 1 (post hoc, after the Part A results were seen). Patients treated at overdoses are also reported as a count per trial, because 3+3 enrolls about half as many patients as the other designs. The P3 verdict uses the pre-registered percentage. For Part B, the readout also reports the smallest n from which both targets hold for every larger n, since the selection probabilities move in a sawtooth with n. Both additions are labeled in the readout.
