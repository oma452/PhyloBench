# PhyloBench first reproducible experimental study

## Design

A rooted four-species tree `((A:1,B:1)AB:1,(C:1,D:1)CD:1)ROOT;` was simulated under three equal-rate duplication/loss tiers. Each tier used ten fixed seeds (101-110), unit duplication/loss reconciliation costs, and the simulator limits recorded in `metadata/study_configuration.json`.

The pilot used seeds 101-103 for each tier, completed without extinction or limit events, and preceded the full run.

## Results

| Tier | Duplication rate | Loss rate | Trials | Extinct | Complete duplications (mean ± SD) | Complete losses (mean ± SD) | Inferred cost (mean ± SD) |
|---|---:|---:|---:|---:|---:|---:|---:|
| low | 0.01 | 0.01 | 10 | 0 | 0.1 ± 0.3162 | 0 ± 0 | 0.1 ± 0.3162 |
| medium | 0.1 | 0.1 | 10 | 0 | 0.2 ± 0.4216 | 0.4 ± 0.5164 | 0.6 ± 0.6992 |
| high | 0.5 | 0.5 | 10 | 0 | 3.8 ± 2.573 | 3.1 ± 1.449 | 3.6 ± 2.459 |

## Observed trends

Within this finite study, increasing both rates increased the observed mean number of complete duplications and losses. The high-rate tier also produced higher measured inferred reconciliation cost. These are descriptive trends for this species tree, branch lengths, rate grid, seed set, and simulation limits; they do not establish general biological scaling laws or estimator accuracy.

No trials were extinct in the full matrix, so inference metrics were available for all 30 trials. Complete-history counts and inferred metrics remain distinct: complete counts include hidden simulated events, while inferred quantities come from reconciliation of extant observed trees. Loss precision/recall use only structurally identifiable losses with extant sibling witnesses.

## Reproducibility and limitations

All configurations, seeds, rates, costs, and limits are saved in `metadata/study_configuration.json`. Raw per-trial data are in `raw/`; aggregate tables are in `summary/`; figures are in `figures/`; tier-specific reports are in `reports/`.

Limitations include the small number of trials, one species-tree topology, equal duplication and loss rates within each tier, fixed unit costs, and the simulator's finite resource limits. The experiment runner reports sample standard deviations across the ten seeds and does not provide confidence intervals or parameter inference.

