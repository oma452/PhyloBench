# PhyloBench second experimental study: 3x3 factorial design

## Design

This study evaluated all nine combinations of duplication rates (0.01, 0.1, 0.5) and loss rates (0.01, 0.1, 0.5) on `((A:1,B:1)AB:1,(C:1,D:1)CD:1)ROOT;`. Each cell used the same 100 seed IDs (1001-1100), for 900 total trials, with unit reconciliation costs and the configured simulator limits. The original 30-trial exploratory study is preserved under `experiments/2026-10-30-first-study/`.

The pilot ran seeds 1001-1005 in each cell (45 executions); those five seed IDs per cell were reused in the full 100-seed run. Pilot executions are not additional independent observations. The final dataset has 900 trials, not 945. Reusing seed IDs across cells makes the conditions seed-matched and can induce dependence between outcomes. Cross-condition comparisons should account for this shared-seed design; the reported factorial contrasts are descriptive and do not estimate paired uncertainty.

## Statistical summaries

For each non-extinction metric, means and sample standard deviations are computed over the applicable valid observations. Monte Carlo standard error is `sample SD / sqrt(n)`. The original summary CSV reports descriptive 95% normal-approximation intervals (`mean +/- 1.96 * MCSE`); these are not model-based confidence statements. Complete-history metrics use all 100 trials per cell. Inferred and identifiable-loss metrics exclude extinct trials, so their denominators are the number of surviving trials. Extinction frequency uses all trials and is reported with a Wilson score interval in the statistical audit.

Several bounded evaluator metrics use 1.0 when their event denominator is empty. This is a convention for an undefined ratio, not evidence that any event was recovered. It applies to exact event-location agreement, duplication-branch compatibility, and loss-branch precision or recall. Interpret their means alongside the valid-observation counts and event denominators; in particular, loss precision and recall can have empty denominators when there are no inferred losses or no identifiable losses, respectively. The archived trial CSV does not record the underlying speciation or matched-duplication denominators, so this report does not give non-empty-denominator counts for location agreement or duplication compatibility. The audit reports denominator availability for loss precision and recall, which can be checked from the retained inferred-loss and identifiable-loss counts.

## Results

The complete machine-readable results are in `summary/factorial_summary.csv`; raw trial records are in `raw/all_trials.csv`. Figures use those measured cell summaries only.

The statistical audit in `reports/statistical_audit.md` and
`summary/statistical_audit.csv` reports Wilson score intervals for extinction
frequency and deterministic percentile-bootstrap intervals for bounded
per-trial metrics. These alternatives complement the original normal
approximation; the bootstrap intervals are empirical and do not guarantee
nominal coverage. The original summary table is retained for provenance.

## Comparison with exploratory study

The diagonal cells (0.01/0.01, 0.10/0.10, 0.50/0.50) can be compared with the earlier 10-seed exploratory study in `summary/comparison_original_study.csv`. Differences are expected from independent seed sets and finite Monte Carlo variation. The factorial study adds off-diagonal rate combinations and larger per-cell samples; it does not by itself establish general biological laws or causal conclusions.

## Limitations

This remains a simulation study on one species-tree topology and one branch-length scale. Rates are constant within branches and duplication and loss costs are fixed at one. The simulator does not silently truncate histories at `max_events` or `max_lineages`: when a configured limit would be exceeded, it raises `SimulationLimitExceeded`. The archived trial records do not include a limit-exception ledger, so they do not establish whether any other attempts reached a limit. The audited intervals do not account for model uncertainty, topology uncertainty, or dependence induced by shared seed IDs. Extinction reduces the valid sample size for inferred metrics in affected cells; those counts are reported explicitly rather than imputed.

## Cell-level means

The inferred-cost mean is calculated over non-extinct trials in each cell; the
valid sample sizes are therefore 100 minus the extinct count shown below. The
same non-extinct denominator applies to the other inferred metrics.

| Duplication rate | Loss rate | Extinct | Valid inferred n | Complete duplications | Complete losses | Inferred cost |
|---:|---:|---:|---:|---:|---:|---:|
| 0.01 | 0.01 | 0/100 | 100 | 0.05 | 0.04 | 0.09 |
| 0.01 | 0.1 | 1/100 | 99 | 0.06 | 0.55 | 0.5657 |
| 0.01 | 0.5 | 27/100 | 73 | 0.01 | 1.78 | 1.438 |
| 0.1 | 0.01 | 0/100 | 100 | 0.61 | 0.1 | 0.71 |
| 0.1 | 0.1 | 2/100 | 98 | 0.53 | 0.63 | 1.01 |
| 0.1 | 0.5 | 24/100 | 76 | 0.25 | 1.97 | 1.579 |
| 0.5 | 0.01 | 0/100 | 100 | 5.5 | 0.12 | 5.52 |
| 0.5 | 0.1 | 0/100 | 100 | 5.11 | 1 | 5.23 |
| 0.5 | 0.5 | 15/100 | 85 | 3.57 | 3.41 | 3.882 |
