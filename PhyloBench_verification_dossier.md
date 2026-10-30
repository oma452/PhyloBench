# PhyloBench manuscript: verification dossier

Companion to `PhyloBench_manuscript.md`. Repository paths are relative to the repository root. The manuscript and this dossier were checked against the current source, archived configuration and CSV artifacts, reports, tests, and available GitHub Actions evidence. The full original study-generation scripts are not present, so this review does not claim end-to-end regeneration of the studies.

## 1. Verification performed

| Check | Evidence and method | Result |
|---|---|---|
| Base revision and CI | Current `main` at `2e4845529bcf6f3ce1af139c556ea0977a108939`; GitHub Actions run 41027204289 | All four Python 3.10–3.13 jobs succeeded. Python 3.12 job log reports 64 collected and passed test items. |
| Local tests | Local runtime is Python 3.14; pytest and Biopython are not installed | The test suite was not run locally. CI evidence above is from the base revision, before these documentation edits. |
| Resource-limit behavior | `src/phylo_bench/simulation.py`, `src/phylo_bench/model/simulation.py`, and `src/phylo_bench/experiments/runner.py` | Attempting to exceed either configured limit raises `SimulationLimitExceeded`; the runner does not convert it into a truncated trial. No archived exception ledger was found, so the archive cannot establish whether any run attempt reached a limit. |
| Pilot-record reconciliation | Both study metadata files, all pilot and full-run CSVs | Study 1 has 9 pilot executions (3 seeds × 3 tiers); Study 2 has 45 (5 seeds × 9 cells); combined total is 54. Each pilot study/cell/seed key is unique, appears exactly once in its full-run dataset, and has matching shared fields. No pilot key is missing or duplicated. Full-run archives remain 30 and 900 rows. |
| Study 2 seed design | `experiments/2026-10-30-second-study/metadata/study_configuration.json` and raw CSVs | Pilot IDs 1001–1005 per cell are reused in the full 100-seed run for all nine cells. The same seed IDs are reused across factorial conditions. Full aggregate and per-cell trial CSVs agree: 900 rows, 100 unique cell/seed keys per cell, with no duplicate cell/seed key. |
| Study 2 summaries | `raw/all_trials.csv`, `summary/factorial_summary.csv` | Raw file has 900 rows; factorial summary has 99 metric/cell rows (11 metrics × 9 cells). Recomputed means and sample SDs for all 99 rows matched the stored summary. |
| Extinction intervals | Raw trial CSV and `summary/statistical_audit.csv` | Recomputed Wilson intervals for all nine extinction rates matched the archived audit values. |
| Bootstrap intervals | `analysis/statistical_audit.py` and archived audit outputs | The archived intervals and script were inspected; bootstrap intervals were not independently recomputed for this review. |
| Empty-denominator counts | Study 2 raw trial CSV, source metric implementation | Among 831 non-extinct trials, loss precision has a non-empty denominator in 361; recall has a non-empty denominator in 382. The remaining 470 precision and 449 recall entries use the empty-denominator convention. Per-cell counts are in manuscript Table 5. |
| Figure assets | Referenced Study 2 figure paths | The extinction heatmap and inferred-cost main-effects figure referenced by the manuscript exist. |
| Reference metadata | Publisher, bibliographic, or indexing records consulted for cited literature | Bibliographic details cited in the manuscript were checked where available. DOIs remain unspecified for Goodman et al. (1979) and Page and Charleston (1997); Zmasek and Eddy (2001) is identified by a verified PMID. |

## 2. Corrections and scientific interpretation

### Simulation limits

Both event-record and lineage-allocation limits raise `SimulationLimitExceeded`. The study report must not describe these limits as silently truncating histories. Successful archived trial rows alone do not prove that no attempted execution reached a limit: there is no archived limit-exception ledger, and complete duplication/loss counts are not totals of all event records or allocated lineages.

### Pilot records and seed design

The pilot counts must be separated by study. Study 1 has 9 pilot executions: seeds 101–103 in each of three tiers. Study 2 has 45: seeds 1001–1005 in each of nine cells. The combined pilot archive therefore contains 54 executions (9 + 45); 54 is not the Study 2 count. All 9 Study 1 and all 45 Study 2 pilot study/cell/seed keys occur once in their corresponding full-run datasets, with matching shared trial fields. No pilot key is missing or duplicated. These overlapping pilot executions add no trial rows: the Study 1 and Study 2 full archives remain 30 and 900 rows, rather than 39 and 945.

Within Study 2, seed IDs 1001–1100 are reused in each of the nine factorial conditions. This creates a common-seed design that can make cross-condition outcomes dependent. The archived analysis presents cell-wise summaries and descriptive contrasts; it does not estimate uncertainty for paired cross-condition contrasts. Seed reuse must not be described as establishing statistical independence.

### Empty-denominator metrics

The evaluator returns 1.0 when the relevant denominator is empty. This convention applies to location agreement, duplication compatibility, and loss precision/recall. Such a 1.0 is not evidence that an event was successfully recovered. Loss metrics require particular care: for non-extinct trials, precision is defined on 361/831 trials with inferred losses, and recall on 382/831 trials with identifiable losses. The other trials in each metric use the convention. The raw export does not include the denominators needed to separate defined from empty-denominator cases for the location and duplication compatibility metrics; their reported 1.0 values cannot be decomposed from these CSVs.

Inferred metrics are conditional on non-extinct trials. Means should therefore be read alongside valid-trial counts and, where available, event-denominator counts. The manuscript reports the loss precision/recall non-empty counts by cell and their conditional means. It retains the archived aggregate values and does not rewrite the evaluator or experimental outputs.

### Scope of reproducibility evidence

The repository contains study configuration, raw and summary outputs, and an audit script, but not the complete scripts that generated the study CSVs and figures. Source/configuration checks and summary recomputations support the stated details, but do not constitute full end-to-end regeneration. No claim is made that the original outputs were regenerated from the recorded seeds.

## 3. Numerical and artifact checks

- Study 2 extinction: 69/900 (7.7%); 831 trials are non-extinct. Of the extinctions, 66 occurred at loss rate 0.5, three at 0.1, and none at 0.01.
- Across the nine Study 2 cells, the tabled means and sample standard deviations for the 99 metric/cell combinations matched recomputation from the raw file. All nine archived Wilson intervals also matched recomputation.
- Study 2 has 900 raw rows and 99 factorial summary rows. The archived audit contains nine cell-wise extinction intervals; these were checked against the raw trials.
- The audit script's bootstrap intervals were not independently recomputed. Do not describe the whole audit output as byte-for-byte regenerated.
- Study 1 contains 30 full-run raw rows. Its three pilot seeds per tier are included again in the full run.
- Pilot files contain 9 Study 1 rows and 45 Study 2 rows (54 combined); all 54 overlap corresponding full-run study/cell/seed records. The full-run datasets have 30 and 900 rows respectively, and each aggregate agrees with its per-cell raw files.
- Raw-data maxima for complete duplication or loss counts do not establish whether any simulation exceeded configured event or lineage limits.

## 4. Repository support and limitations

The manuscript's implementation descriptions were checked against source and tests, including:

- Tree representation and validation: `src/phylo_bench/model/tree.py`, `src/phylo_bench/validation.py`.
- Newick parsing and iterative traversals: `src/phylo_bench/io/newick.py`, `src/phylo_bench/algorithms/traversal.py`.
- LCA and reconciliation: `src/phylo_bench/algorithms/lca.py`, `src/phylo_bench/algorithms/reconciliation.py`.
- Simulation, projection, and evaluation: `src/phylo_bench/simulation.py`, `src/phylo_bench/model/simulation.py`.
- Experiment execution: `src/phylo_bench/experiments/runner.py`.
- Statistical audit: `experiments/2026-10-30-second-study/analysis/statistical_audit.py`.

The following limitations remain material and should remain explicit: one rooted four-species topology with equal branch lengths; constant duplication and loss rates; a simplified duplication–loss process without HGT, ILS, gene conversion, sequence evolution, or rate heterogeneity; extinction reduces the valid sample for inferred metrics; complete history differs from events identifiable from extant genes; branch-level loss matching is not exact event recovery; branch compatibility is not event-time recovery. These studies do not establish general biological laws or algorithmic superiority.

Additional provenance limits: original study-generation and figure scripts are absent; the original Python/package versions and code revision are not recorded in study metadata; dependency declarations are lower bounds rather than a lockfile; the repository has no LICENSE or CITATION file. Author names, affiliations, and correspondence remain to be supplied by the authors.

## 5. Bibliography and author inputs

Bibliographic records were checked for the manuscript's cited works against available publisher/index records. The manuscript includes the Wilson interval and bootstrap references. DOIs remain unverified for Goodman et al. (1979) and Page and Charleston (1997); Zmasek and Eddy (2001) is identified by its verified PMID. These identifiers were not guessed. Author identity, affiliations, correspondence, funding, acknowledgements, and author contributions are not inferable from repository evidence and require author input. Generative AI assistance was used for research assessment and prose editing; authors should provide any venue-required disclosure, including tool/model details, and review all final content.
