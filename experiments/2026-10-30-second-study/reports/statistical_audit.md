# Statistical audit: second PhyloBench factorial study

## Scope

Raw data contain 900 trials across 9 cells. Every cell contains exactly 100 trials. This audit reads existing CSV artifacts only and does not rerun simulation.

The metadata records five pilot seeds (1001-1005) in each of the nine cells and 100 full-run seeds (1001-1100) in each cell. The raw CSV has exactly 900 rows: each of the 45 pilot seed/cell combinations appears in the full-run records, so pilot executions are not additional observations and the dataset is not 945 trials. The same seed IDs occur across conditions; comparisons are seed-matched and may be dependent. The factorial contrasts below are descriptive and do not include paired-contrast uncertainty.

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

## Empty-denominator metrics

The evaluator returns 1.0 when a ratio's denominator is empty. This convention is not evidence of successful recovery. It applies to exact event-location agreement when there are no matched speciation events, duplication-branch compatibility when there are no matched duplication events, loss precision when no losses are inferred, and loss recall when there are no identifiable losses. Loss precision and recall therefore need particular care: many trials have no inferred or identifiable losses. Inferred metrics are conditional on non-extinct trials.

The raw CSV retains `inferred_losses` and `identifiable_loss_count`, so the following counts show non-extinct observations with non-empty precision and recall denominators. They should be read alongside the valid inferred sample sizes above and the metric means. The raw CSV does not retain the event counts underlying exact-location agreement or duplication compatibility; their non-empty-denominator counts cannot be recovered from the archived rows without guessing.

| Cell | Non-extinct observations with inferred losses (precision denominator) | Non-extinct observations with identifiable losses (recall denominator) |
|---|---:|---:|
| d=0.01, l=0.01 | 4/100 | 4/100 |
| d=0.01, l=0.1 | 45/99 | 42/99 |
| d=0.01, l=0.5 | 66/73 | 61/73 |
| d=0.1, l=0.01 | 10/100 | 10/100 |
| d=0.1, l=0.1 | 48/98 | 50/98 |
| d=0.1, l=0.5 | 69/76 | 65/76 |
| d=0.5, l=0.01 | 6/100 | 11/100 |
| d=0.5, l=0.1 | 43/100 | 60/100 |
| d=0.5, l=0.5 | 70/85 | 79/85 |
| **Total** | **361/831** | **382/831** |

## Confidence-interval recommendation

Normal-approximation intervals for bounded metrics can extend beyond [0,1], and clipping them can conceal skewness and boundary behavior. The audit therefore reports Wilson score intervals for extinction frequency and deterministic percentile-bootstrap intervals for bounded per-trial metrics (location agreement, duplication compatibility, and loss precision/recall). Since those metric values include the empty-denominator convention described above, their means and bootstrap intervals can reflect convention-assigned 1.0 values as well as observed event recovery. Unbounded counts and costs retain descriptive normal-approximation intervals. Percentile-bootstrap intervals reflect the empirical distribution and do not guarantee nominal coverage; intervals can also be degenerate when the observed values lie at a boundary. The original summary CSV is preserved; `summary/statistical_audit.csv` contains the recommended intervals and method labels.

## Factorial effects

Marginal means and high-minus-low contrasts are descriptive summaries of the balanced design. They show strong duplication-rate effects on duplication counts and cost, strong loss-rate effects on loss counts and extinction, and rate-dependent changes in inferred cost. Interaction contrasts are included in `summary/factorial_effects.csv`; they are not treated as formal hypothesis tests because this study is a stochastic simulation benchmark with model-specific estimands and no preregistered inferential model.

## Figures and tables

All cell means used by the heatmaps and main-effect figures match raw-data recomputation to floating-point tolerance. Raw trial counts, valid denominators, and extinction denominators are retained in the audit CSV. The original and second-study artifacts remain separate.

## Limitations

The intervals quantify Monte Carlo uncertainty conditional on this simulator, topology, branch lengths, rates, and fixed seed set; they do not quantify biological or model uncertainty. Bootstrap intervals inherit the empirical distribution and are not guarantees of exact coverage. Per-trial bounded metrics are ratios produced by the reconciliation evaluator, so pooled binomial Wilson intervals are not appropriate without numerator/denominator event counts. Extinction makes inferred metrics conditional estimands over surviving families, not unconditional population means. The simulator raises `SimulationLimitExceeded` when configured event or lineage limits would be exceeded; it does not silently truncate histories. Archived records do not include a limit-exception ledger, so this audit makes no claim about whether any other simulation attempts reached a limit.
