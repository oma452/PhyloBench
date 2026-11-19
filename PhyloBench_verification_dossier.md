# PhyloBench manuscript: verification dossier

**Author:** Omar Shibl
**Affiliation:** Independent Researcher, Egypt
**Corresponding author:** Omar Shibl

Companion to `PhyloBench_manuscript.md`. Repository paths are relative to the repository root. The manuscript and this dossier were checked against the current source, archived configuration and CSV artifacts, reports, tests, and available GitHub Actions evidence. The full original study-generation and figure-generation pipeline is not present, so this dossier does not claim end-to-end reproduction of the historical studies or byte-for-byte reproduction of their artifact bundles.

## 1. Verification performed

| Check | Evidence and method | Result |
|---|---|---|
| Historical CI revision and result | GitHub Actions run [44547413391](https://github.com/oma452/PhyloBench/actions/runs/44547413391), tested commit `e4c3d6c9c6a078ebf99ffade7567c53dee2c2c7f` | The historical run completed successfully on Python 3.10, 3.11, 3.12, and 3.13. It ran 92 tests and reported 95.24% coverage under the 90% threshold. This result validates commit `e4c3d6c`, not later documentation or packaging commits, and it does not independently reproduce the archived experiments. |
| Local validation of the integrated work | Local validation was performed on the clean local `main` history before the merge commit, including packaging commit `4e170a8`, dossier/manuscript commit `f6a8fca`, and dirty-state commit `e3f3357`. | The final local run passed 100 tests. Ruff and mypy passed. A clean-wheel build from `e3f3357` embedded the exact `git rev-parse HEAD` revision and reported `working_tree_dirty: false`. These results are distinct from historical CI and do not independently validate the scientific archive. |
| Branch integration and current repository revision | Merge of remote `origin/main` (`2f8de31abf764d582af333664178ca560141d213`) with local packaging/documentation history (`4e170a8`, `f6a8fca`, `e3f3357`) | The merge preserves both histories and combines the remote refined Figure 1/dossier wording with the local provenance and dirty-state fixes. The resulting merge commit is the revision of integrated `main`; no force-push or history rewrite was used. |
| Resource-limit behavior | `src/phylo_bench/simulation.py`, `src/phylo_bench/model/simulation.py`, and `src/phylo_bench/experiments/runner.py` | Attempting to exceed either configured limit raises `SimulationLimitExceeded`; the runner does not convert it into a truncated trial. No archived exception ledger was found, so the archive cannot establish whether any run attempt reached a limit. |
| Pilot-record reconciliation | Both study metadata files, all pilot and full-run CSVs | Study 1 has 9 pilot executions (3 seeds × 3 tiers); Study 2 has 45 (5 seeds × 9 cells); combined total is 54. Each pilot study/cell/seed key is unique, appears exactly once in its full-run dataset, and has matching shared fields. No pilot key is missing or duplicated. Full-run archives remain 30 and 900 rows. |
| Study 2 seed design | `experiments/2026-10-30-second-study/metadata/study_configuration.json` and raw CSVs | Pilot IDs 1001–1005 per cell are reused in the full 100-seed run for all nine cells. The same seed IDs are reused across factorial conditions. Full aggregate and per-cell trial CSVs agree: 900 rows, 100 unique cell/seed keys per cell, with no duplicate cell/seed key. |
| Study 2 summaries | `raw/all_trials.csv`, `summary/factorial_summary.csv` | Raw file has 900 rows; factorial summary has 99 metric/cell rows (11 metrics × 9 cells). Recomputed means and sample SDs for all 99 rows matched the stored summary. |
| Extinction intervals | Raw trial CSV and `summary/statistical_audit.csv` | Recomputed Wilson intervals for all nine extinction rates matched the archived audit values. |
| Bootstrap intervals | `analysis/statistical_audit.py` and archived audit outputs | The archived intervals and script were inspected; bootstrap intervals were not independently recomputed for this review. |
| Empty-denominator counts | Study 2 raw trial CSV, source metric implementation | Among 831 non-extinct trials, loss precision has a non-empty denominator in 361; recall has a non-empty denominator in 382. The remaining 470 precision and 449 recall entries use the empty-denominator convention. Per-cell counts are in manuscript Table 5. |
| Figure assets | Referenced Study 2 figure paths | The extinction heatmap and inferred-cost main-effects figure referenced by the manuscript exist. |
| Regeneration comparison | Current `experiments/scripts/reproduce_study.py` workflow and archived study CSVs | Study 1: 30/30 trial rows matched; Study 2: 900/900 trial rows matched. The comparison covered the seed identifier plus 11 shared measurement fields: `extinct`, `complete_duplications`, `complete_losses`, `inferred_duplications`, `inferred_losses`, `inferred_cost`, `identifiable_loss_count`, `exact_event_location_agreement`, `duplication_branch_compatibility`, `loss_branch_precision`, and `loss_branch_recall`. It did not establish byte-for-byte identity for every figure, report, CSV, or artifact. |
| Paired-seed analysis | `experiments/2026-10-30-second-study/analysis/paired_seed_analysis.py` and Study 2 raw CSV | The script produces 180 descriptive contrast records: 36 cell pairs × 5 metrics. These are paired contrasts with observation counts, mean differences, paired sample SDs, and MCSEs; no hypothesis tests, significance claims, confidence intervals, or confirmatory inference are performed. |
| Reference metadata | Publisher, bibliographic, or indexing records consulted for cited literature | Bibliographic details cited in the manuscript were checked where available. DOIs remain unspecified for Goodman et al. (1979) and Page and Charleston (1997); Zmasek and Eddy (2001) is identified by a verified PMID. |

## 2. Corrections and scientific interpretation

### Simulation limits

Both event-record and lineage-allocation limits raise `SimulationLimitExceeded`. The study report must not describe these limits as silently truncating histories. Successful archived trial rows alone do not prove that no attempted execution reached a limit: there is no archived limit-exception ledger, and complete duplication/loss counts are not totals of all event records or allocated lineages.

### Pilot records and seed design

The pilot counts must be separated by study. Study 1 has 9 pilot executions: seeds 101–103 in each of three tiers. Study 2 has 45: seeds 1001–1005 in each of nine cells. The combined pilot archive therefore contains 54 executions (9 + 45); 54 is not the Study 2 count. All 9 Study 1 and all 45 Study 2 pilot study/cell/seed keys occur once in their corresponding full-run datasets, with matching shared trial fields. No pilot key is missing or duplicated. These overlapping pilot executions add no trial rows: the Study 1 and Study 2 full archives remain 30 and 900 rows, rather than 39 and 945.

Within Study 2, seed IDs 1001–1100 are reused in each of the nine factorial conditions. This creates a common-seed design that can make cross-condition outcomes dependent. The paired-seed analysis reports descriptive paired cross-condition contrasts with paired observation counts, mean differences, paired sample SDs, and MCSEs. These quantify variation in the paired differences but are not hypothesis tests, significance claims, confidence intervals, or confirmatory inference. Seed reuse must not be described as establishing statistical independence.

### Empty-denominator metrics

The evaluator returns 1.0 when the relevant denominator is empty. This convention applies to location agreement, duplication compatibility, and loss precision/recall. Such a 1.0 is not evidence that an event was successfully recovered. Loss metrics require particular care: for non-extinct trials, precision is defined on 361/831 trials with inferred losses, and recall on 382/831 trials with identifiable losses. The other trials in each metric use the convention. The raw export does not include the denominators needed to separate defined from empty-denominator cases for the location and duplication compatibility metrics; their reported 1.0 values cannot be decomposed from these CSVs.

Inferred metrics are conditional on non-extinct trials. Means should therefore be read alongside valid-trial counts and, where available, event-denominator counts. The manuscript reports the loss precision/recall non-empty counts by cell and their conditional means. It retains the archived aggregate values and does not rewrite the evaluator or experimental outputs.

### Scope of reproducibility evidence

The repository contains study configurations, archived raw and summary outputs, reports, figures, provenance manifests, artifact hashes, atomic output writing, a resolved Python 3.12 environment record, and executable regeneration scripts. The lockfile documents a resolved Python 3.12 environment; it does not reproduce the original historical environment. The regeneration scripts run the archived configurations and seeds through the current CLI, but the complete original experiment-generation and figure-generation scripts are absent. The verified comparison therefore covers only 30/30 and 900/900 matching trial rows across the seed identifier and the 11 shared measurement fields listed above; it does not claim byte-for-byte identity for every figure, report, CSV, or artifact. Newly generated outputs must be distinguished from archived historical CSVs, and no fields absent from the historical CSV schema are inferred.

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

Additional provenance limits: original study-generation and figure scripts are absent; the original Python/package versions and code revision are not recorded in study metadata; the Python 3.12 lockfile records a resolved current environment rather than the historical one. The repository now includes an MIT `LICENSE`, `CITATION.cff`, provenance manifests, artifact hashes, atomic output writing, regeneration scripts, paired-seed analysis, independent reconciliation oracle tests, analytical simulator tests, and updated CI and quality tooling.

## 5. Bibliography and author information

Bibliographic records were checked for the manuscript's cited works against available publisher/index records. The manuscript includes the Wilson interval and bootstrap references. DOIs remain unverified for Goodman et al. (1979) and Page and Charleston (1997); Zmasek and Eddy (2001) is identified by its verified PMID. These identifiers were not guessed. The author information for this submission is Omar Shibl, Independent Researcher, Egypt; Omar Shibl is the corresponding author. No correspondence email is included because none was provided. No funding source, institutional affiliation beyond the stated independent-researcher affiliation, acknowledgements, ORCID, or other personal information is asserted.
