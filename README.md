# PhyloBench

PhyloBench is a reproducible Python framework for gene-tree/species-tree
reconciliation research. It includes typed rooted trees, iterative Newick
parsing, binary-lifting LCA queries, classical duplication-loss
reconciliation, deterministic gene-family simulation, and research-oriented
visualization and experiment reporting.

## Installation

```bash
python -m venv .venv
.venv/bin/pip install -e ".[test]"
```

## Reproducible experiments

Experiments use explicit species-tree Newick input, duplication and loss
rates, repeated integer seeds, and configurable safety limits. Results contain
per-seed records, means, standard deviations, and the complete configuration.
Only measured simulator outputs are written to CSV and Markdown reports.

```python
from phylo_bench.experiments import ExperimentConfig, run_experiment
from phylo_bench.export import write_rows_csv
from phylo_bench.report import write_research_report

config = ExperimentConfig(
    species_tree_newick="(A:1,B:1)S;",
    seeds=(1, 2, 3),
    duplication_rate=0.1,
    loss_rate=0.05,
)
rows, summary = run_experiment(config)
write_rows_csv("trials.csv", rows)
write_research_report("report.md", config, rows, summary)
```

### Command-line interface

Create a JSON configuration:

```json
{
  "species_tree_newick": "(A:1,B:1)S;",
  "seeds": [1, 2, 3],
  "duplication_rate": 0.1,
  "loss_rate": 0.05,
  "duplication_cost": 1.0,
  "loss_cost": 1.0,
  "max_events": 100000,
  "max_lineages": 10000
}
```

Validate and inspect it:

```bash
phylo-bench inspect experiment.json
```

Run the experiment and write `trials.csv`, `summary.csv`, `metadata.json`,
`run_manifest.json`, `report.md`, and measured PNG charts. The manifest
records the configuration, command line, source revision, Python/platform
information, dependency versions, and SHA-256 hashes for generated artifacts:

```bash
phylo-bench run experiment.json --output-dir results/experiment-1
```

The CLI uses the listed seeds exactly and preserves the configured simulation
limits.

Archived studies can be regenerated into a new directory without modifying
the committed research artifacts:

```bash
python experiments/scripts/reproduce_study.py \
  --study first \
  --output-dir /tmp/phylo-bench-first-reproduction
python experiments/scripts/reproduce_study.py \
  --study second \
  --output-dir /tmp/phylo-bench-second-reproduction
```

These scripts use the archived configurations and seeds through the current
CLI. They provide an executable regeneration workflow, but cannot claim
byte-for-byte identity with the archived outputs because the original study
generation scripts were not preserved.

Study 2 reuses the same seed IDs across factorial cells. For descriptive
paired common-seed contrasts, run:

```bash
python experiments/2026-10-30-second-study/analysis/paired_seed_analysis.py \
  --output /tmp/phylo-bench-paired-contrasts.csv
```

Matplotlib renderers in `phylo_bench.visualization` produce species-tree,
gene-tree, reconciliation, and measured performance figures. Visualization
does not participate in reconciliation or simulation.

## Scientific conventions

Reconciliation requires rooted binary trees and explicit extant
gene-to-species mappings. LCA mapping classifies internal gene nodes as
duplications or speciations and infers losses on species-tree branches.
Simulation retains complete history separately from the observed extant gene
tree. Exact event-location agreement is reported separately from duplication
branch compatibility, because exact duplication time is generally not
identifiable from extant leaves. Extinct trials are counted separately and
are not assigned fabricated inferred metrics.

Some evaluator ratios are assigned 1.0 when their denominator is empty. This
is a convention for an undefined ratio, not evidence of successful event
recovery. Interpret bounded metric means alongside valid-observation counts
and their event denominators, especially loss precision and recall, which may
have no inferred or identifiable losses to evaluate.

## Research reports

The [second study report](experiments/2026-10-30-second-study/reports/study_report.md)
documents a 3x3 duplication/loss-rate simulation study. Its
[statistical audit](experiments/2026-10-30-second-study/reports/statistical_audit.md)
checks the reported summaries against the archived trial records and describes
the interval methods and their limitations. These reports describe results
conditional on the study's simulator settings and species-tree design.

## Validation

```bash
.venv/bin/pytest
```

The test suite covers tree algorithms, reconciliation, simulation invariants,
reproducibility, CSV/Markdown generation, and Matplotlib figure generation.

Optional quality checks use the `quality` extra:

```bash
.venv/bin/ruff check src tests experiments/scripts
.venv/bin/mypy src/phylo_bench
.venv/bin/pytest --cov=phylo_bench --cov-report=term-missing
```

The repository includes `requirements-lock-py312.txt` as a resolved
environment record for Python 3.12 validation. The package itself supports
Python 3.10 and later.

## Continuous integration

GitHub Actions runs the complete test suite on every push and pull request
using Python 3.10, 3.11, 3.12, and 3.13. Each job installs the package with
its test extra, verifies package import and `phylo-bench --help`, and then runs
`python -m pytest`. Pip's cache is enabled per Python version using
`pyproject.toml` as the dependency key.
