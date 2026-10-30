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
`report.md`, and measured PNG charts:

```bash
phylo-bench run experiment.json --output-dir results/experiment-1
```

The CLI uses the listed seeds exactly and preserves the configured simulation
limits.

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

## Validation

```bash
.venv/bin/pytest
```

The test suite covers tree algorithms, reconciliation, simulation invariants,
reproducibility, CSV/Markdown generation, and Matplotlib figure generation.

## Continuous integration

GitHub Actions runs the complete test suite on every push and pull request
using Python 3.10, 3.11, 3.12, and 3.13. Each job installs the package with
its test extra, verifies package import and `phylo-bench --help`, and then runs
`python -m pytest`. Pip's cache is enabled per Python version using
`pyproject.toml` as the dependency key.
