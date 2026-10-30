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
    103,
    104,
    105,
    106,
    107,
    108,
    109,
    110
  ],
  "duplication_rate": 0.1,
  "loss_rate": 0.1,
  "duplication_cost": 1.0,
  "loss_cost": 1.0,
  "max_events": 100000,
  "max_lineages": 10000
}
```

## Results

- Trials: 10
- Extinct trials: 0

Inference metrics use only non-extinct measured trials. Complete history counts use every simulated trial; unavailable inference metrics are not replaced with zeros.

Complete simulated loss counts describe the full hidden history. Identifiable loss metrics use only losses with an extant sibling witness; hidden losses are excluded from precision and recall.

| Metric | Mean | Standard deviation |
|---|---:|---:|
| complete_duplications | 0.2 | 0.421637 |
| complete_losses | 0.4 | 0.516398 |
| inferred_duplications | 0.2 | 0.421637 |
| inferred_losses | 0.4 | 0.516398 |
| inferred_cost | 0.6 | 0.699206 |
| exact_event_location_agreement | 1 | 0 |
| identifiable_loss_count | 0.4 | 0.516398 |
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
    "complete_duplications": 1,
    "complete_losses": 1,
    "inferred_duplications": 1,
    "inferred_losses": 1,
    "inferred_cost": 2.0,
    "identifiable_loss_count": 1,
    "exact_event_location_agreement": 1.0,
    "duplication_branch_compatibility": 1.0,
    "loss_branch_precision": 1.0,
    "loss_branch_recall": 1.0
  },
  {
    "seed": 103,
    "extinct": false,
    "complete_duplications": 0,
    "complete_losses": 1,
    "inferred_duplications": 0,
    "inferred_losses": 1,
    "inferred_cost": 1.0,
    "identifiable_loss_count": 1,
    "exact_event_location_agreement": 1.0,
    "duplication_branch_compatibility": 1.0,
    "loss_branch_precision": 1.0,
    "loss_branch_recall": 1.0
  },
  {
    "seed": 104,
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
    "seed": 105,
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
    "seed": 106,
    "extinct": false,
    "complete_duplications": 1,
    "complete_losses": 0,
    "inferred_duplications": 1,
    "inferred_losses": 0,
    "inferred_cost": 1.0,
    "identifiable_loss_count": 0,
    "exact_event_location_agreement": 1.0,
    "duplication_branch_compatibility": 1.0,
    "loss_branch_precision": 1.0,
    "loss_branch_recall": 1.0
  },
  {
    "seed": 107,
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
    "seed": 108,
    "extinct": false,
    "complete_duplications": 0,
    "complete_losses": 1,
    "inferred_duplications": 0,
    "inferred_losses": 1,
    "inferred_cost": 1.0,
    "identifiable_loss_count": 1,
    "exact_event_location_agreement": 1.0,
    "duplication_branch_compatibility": 1.0,
    "loss_branch_precision": 1.0,
    "loss_branch_recall": 1.0
  },
  {
    "seed": 109,
    "extinct": false,
    "complete_duplications": 0,
    "complete_losses": 1,
    "inferred_duplications": 0,
    "inferred_losses": 1,
    "inferred_cost": 1.0,
    "identifiable_loss_count": 1,
    "exact_event_location_agreement": 1.0,
    "duplication_branch_compatibility": 1.0,
    "loss_branch_precision": 1.0,
    "loss_branch_recall": 1.0
  },
  {
    "seed": 110,
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
