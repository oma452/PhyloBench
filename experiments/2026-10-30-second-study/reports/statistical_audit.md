# Statistical audit: second PhyloBench factorial study

## Scope

Raw data contain 900 trials across 9 cells. Every cell contains exactly 100 trials. This audit reads existing CSV artifacts only and does not rerun simulation.

## Statistical verification

Means and sample standard deviations match the raw CSV values to floating-point tolerance. For ordinary metrics, MCSE is sample SD divided by the square root of valid observations. The original summary tables use the normal approximation `mean +/- 1.96 * MCSE`. For extinction frequency, MCSE uses the plug-in Bernoulli standard deviation `sqrt(p * (1 - p)) / sqrt(n)`; the recommended interval is Wilson score.

| Cell | Extinct | Valid inferred observations |
|---|---:|---:|
| d=0.01, l=0.01 | 0/100 | 100 |
| d=0.01, l=0.1 | 1/100 | 99 |
| d=0.01, l=0.5 | 27/100 | 73 |
| d=0.1, l=0.01 | 0/100 | 100 |
| d=0.1, l=0.1 | 2/100 | 98 |
| d=0.1, l=0.5 | 24/100 | 76 |
| d=0.5, l=0.01 | 0/100 | 100 |
| d=0.5, l=0.1 | 0/100 | 100 |
| d=0.5, l=0.5 | 15/100 | 85 |

## Confidence-interval recommendation

Normal-approximation intervals for bounded metrics can extend beyond [0,1], and clipping them can conceal skewness and boundary behavior. The audit therefore reports Wilson score intervals for extinction frequency and deterministic percentile-bootstrap intervals for bounded per-trial metrics (location agreement, duplication compatibility, and loss precision/recall). Unbounded counts and costs retain descriptive normal-approximation intervals. Percentile-bootstrap intervals reflect the empirical distribution and do not guarantee nominal coverage; intervals can also be degenerate when the observed values lie at a boundary. The original summary CSV is preserved; `summary/statistical_audit.csv` contains the recommended intervals and method labels.

## Factorial effects

Marginal means and high-minus-low contrasts are descriptive summaries of the balanced design. They show strong duplication-rate effects on duplication counts and cost, strong loss-rate effects on loss counts and extinction, and rate-dependent changes in inferred cost. Interaction contrasts are included in `summary/factorial_effects.csv`; they are not treated as formal hypothesis tests because this study is a stochastic simulation benchmark with model-specific estimands and no preregistered inferential model.

## Figures and tables

All cell means used by the heatmaps and main-effect figures match raw-data recomputation to floating-point tolerance. Raw trial counts, valid denominators, and extinction denominators are retained in the audit CSV. The original and second-study artifacts remain separate.

## Limitations

The intervals quantify Monte Carlo uncertainty conditional on this simulator, topology, branch lengths, rates, and fixed seed set; they do not quantify biological or model uncertainty. Bootstrap intervals inherit the empirical distribution and are not guarantees of exact coverage. Per-trial bounded metrics are ratios produced by the reconciliation evaluator, so pooled binomial Wilson intervals are not appropriate without numerator/denominator event counts. Extinction makes inferred metrics conditionally estimands over surviving families, not unconditional population means.
