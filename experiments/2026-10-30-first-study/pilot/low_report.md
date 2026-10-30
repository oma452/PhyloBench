# PhyloBench experimental report

## Methodology

Experiments use the rooted binary duplication-loss simulator and classical LCA reconciliation. Duplication and loss events are simulated with configurable continuous-time rates on species-tree branches. Extinct trials are reported separately.

## Configuration

```json
{
  "species_tree_newick": "((A:1,B:1)AB:1,(C:1,D:1)CD:1)ROOT;",
  "seeds": [
    101,
    102,
    103
  ],
  "duplication_rate": 0.01,
  "loss_rate": 0.01,
  "duplication_cost": 1.0,
  "loss_cost": 1.0,
  "max_events": 100000,
  "max_lineages": 10000
}
```

## Results

- Trials: 3
- Extinct trials: 0

Inference metrics use only non-extinct measured trials. Complete history counts use every simulated trial; unavailable inference metrics are not replaced with zeros.

Complete simulated loss counts describe the full hidden history. Identifiable loss metrics use only losses with an extant sibling witness; hidden losses are excluded from precision and recall.

| Metric | Mean | Standard deviation |
|---|---:|---:|
| complete_duplications | 0 | 0 |
| complete_losses | 0 | 0 |
| inferred_duplications | 0 | 0 |
| inferred_losses | 0 | 0 |
| inferred_cost | 0 | 0 |
| exact_event_location_agreement | 1 | 0 |
| identifiable_loss_count | 0 | 0 |
| duplication_branch_compatibility | 1 | 0 |
| loss_branch_precision | 1 | 0 |
| loss_branch_recall | 1 | 0 |

## Interpretation and limitations

Exact event-location agreement applies to inferred speciation locations. Duplication branch compatibility is a weaker measure because exact within-branch duplication time is generally not identifiable from extant leaves. Complete simulated events and inferred extant-tree events are reported separately.

## Trial records

```json
[
  {
    "seed": 101,
    "extinct": false,
    "complete_duplications": 0,
    "complete_losses": 0,
    "inferred_duplications": 0,
    "inferred_losses": 0,
    "inferred_cost": 0.0,
    "identifiable_loss_count": 0,
    "exact_event_location_agreement": 1.0,
    "duplication_branch_compatibility": 1.0,
    "loss_branch_precision": 1.0,
    "loss_branch_recall": 1.0
  },
  {
    "seed": 102,
    "extinct": false,
    "complete_duplications": 0,
    "complete_losses": 0,
    "inferred_duplications": 0,
    "inferred_losses": 0,
    "inferred_cost": 0.0,
    "identifiable_loss_count": 0,
    "exact_event_location_agreement": 1.0,
    "duplication_branch_compatibility": 1.0,
    "loss_branch_precision": 1.0,
    "loss_branch_recall": 1.0
  },
  {
    "seed": 103,
    "extinct": false,
    "complete_duplications": 0,
    "complete_losses": 0,
    "inferred_duplications": 0,
    "inferred_losses": 0,
    "inferred_cost": 0.0,
    "identifiable_loss_count": 0,
    "exact_event_location_agreement": 1.0,
    "duplication_branch_compatibility": 1.0,
    "loss_branch_precision": 1.0,
    "loss_branch_recall": 1.0
  }
]
```
