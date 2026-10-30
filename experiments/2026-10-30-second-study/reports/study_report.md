# PhyloBench second experimental study: 3x3 factorial design

## Design

This study evaluated all nine combinations of duplication rates (0.01, 0.1, 0.5) and loss rates (0.01, 0.1, 0.5) on `((A:1,B:1)AB:1,(C:1,D:1)CD:1)ROOT;`. Each cell used 100 fixed seeds (1001-1100), for 900 total trials, with unit reconciliation costs and the configured simulator limits. The original 30-trial exploratory study is preserved under `experiments/2026-10-30-first-study/`.

The pilot used seeds 1001-1005 in every cell (45 trials total) before the full run.

## Statistical summaries

For each metric, means and sample standard deviations are computed over the applicable valid observations. Monte Carlo standard error is `sample SD / sqrt(n)`. Reported 95% intervals use the normal approximation `mean +/- 1.96 * MCSE`; they are descriptive simulation-study intervals, not model-based confidence statements. Complete-history metrics use all trials. Inferred and identifiable-loss metrics exclude extinct trials. Extinction frequency uses all trials.

## Results

The complete machine-readable results are in `summary/factorial_summary.csv`; raw trial records are in `raw/all_trials.csv`. Figures use those measured cell summaries only.

The statistical audit in `reports/statistical_audit.md` and
`summary/statistical_audit.csv` provides the preferred interval reporting:
Wilson score intervals for extinction frequency and deterministic percentile
bootstrap intervals for bounded per-trial metrics. The original summary table
is retained for provenance.

## Comparison with exploratory study

The diagonal cells (0.01/0.01, 0.10/0.10, 0.50/0.50) can be compared with the earlier 10-seed exploratory study in `summary/comparison_original_study.csv`. Differences are expected from independent seed sets and finite Monte Carlo variation. The factorial study adds off-diagonal rate combinations and larger per-cell samples; it does not by itself establish general biological laws or causal conclusions.

## Limitations

This remains a simulation study on one species-tree topology and one branch-length scale. Rates are constant within branches, duplication and loss costs are fixed at one, and the simulator limits could truncate unusually large histories. The audited intervals do not account for model uncertainty, topology uncertainty, or dependence induced by any shared design choices. Extinction reduces the valid sample size for inferred metrics in affected cells; those counts are reported explicitly rather than imputed.

## Cell-level means

| Duplication rate | Loss rate | Extinct | Complete duplications | Complete losses | Inferred cost |
|---:|---:|---:|---:|---:|---:|
| 0.01 | 0.01 | 0/100 | 0.05 | 0.04 | 0.09 |
| 0.01 | 0.1 | 1/100 | 0.06 | 0.55 | 0.5657 |
| 0.01 | 0.5 | 27/100 | 0.01 | 1.78 | 1.438 |
| 0.1 | 0.01 | 0/100 | 0.61 | 0.1 | 0.71 |
| 0.1 | 0.1 | 2/100 | 0.53 | 0.63 | 1.01 |
| 0.1 | 0.5 | 24/100 | 0.25 | 1.97 | 1.579 |
| 0.5 | 0.01 | 0/100 | 5.5 | 0.12 | 5.52 |
| 0.5 | 0.1 | 0/100 | 5.11 | 1 | 5.23 |
| 0.5 | 0.5 | 15/100 | 3.57 | 3.41 | 3.882 |
