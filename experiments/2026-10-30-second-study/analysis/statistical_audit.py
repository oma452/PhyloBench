from __future__ import annotations
import csv, hashlib, json, math, random
from collections import defaultdict
from pathlib import Path
from statistics import mean, stdev
from audit_helpers import (
    compare_summary_record,
    extinction_statistics,
    verify_summary_errors,
)

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'raw' / 'all_trials.csv'
SUMMARY = ROOT / 'summary' / 'factorial_summary.csv'
AUDIT_SUMMARY = ROOT / 'summary' / 'statistical_audit.csv'
EFFECTS = ROOT / 'summary' / 'factorial_effects.csv'
REPORT = ROOT / 'reports' / 'statistical_audit.md'

BOUNDED = {'exact_event_location_agreement', 'duplication_branch_compatibility', 'loss_branch_precision', 'loss_branch_recall'}
COMPLETE = {'complete_duplications', 'complete_losses'}
METRICS = ('complete_duplications','complete_losses','inferred_duplications','inferred_losses','inferred_cost','identifiable_loss_count','exact_event_location_agreement','duplication_branch_compatibility','loss_branch_precision','loss_branch_recall')
Z = 1.96


def normal_ci(values):
    n = len(values)
    if not n: return None, None, None, None
    sd = stdev(values) if n > 1 else 0.0
    mcse = sd / math.sqrt(n)
    return mean(values), sd, mcse, (mean(values)-Z*mcse, mean(values)+Z*mcse)


def wilson(successes, n):
    if n == 0: return None, None
    p = successes / n
    den = 1 + Z*Z/n
    center = (p + Z*Z/(2*n))/den
    half = Z * math.sqrt((p*(1-p) + Z*Z/(4*n))/n) / den
    return max(0.0, center-half), min(1.0, center+half)


def bootstrap_ci(values, seed, reps=20000):
    if not values: return None, None
    if len(values) == 1: return values[0], values[0]
    rng = random.Random(seed)
    n = len(values)
    samples = []
    for _ in range(reps):
        samples.append(sum(values[rng.randrange(n)] for _ in range(n))/n)
    samples.sort()
    return samples[int(0.025*reps)], samples[int(0.975*reps)-1]

with RAW.open() as f:
    rows = list(csv.DictReader(f))
with SUMMARY.open() as f:
    old_summary = list(csv.DictReader(f))
assert len(rows) == 900
cells = defaultdict(list)
for row in rows:
    cells[(float(row['duplication_rate']), float(row['loss_rate']))].append(row)
assert len(cells) == 9 and all(len(v) == 100 for v in cells.values())

# Validate every summary mean/SD against raw rows.
summary_lookup = {(float(r['duplication_rate']), float(r['loss_rate']), r['metric']): r for r in old_summary}
validation = []
for (d,l), cell in sorted(cells.items()):
    for metric in METRICS:
        valid = [float(r[metric]) for r in cell if r[metric] != '']
        avg, sd, mcse, ci = normal_ci(valid)
        mean_error, sd_error = compare_summary_record(
            summary_lookup, (d, l, metric), avg, sd
        )
        validation.append({
            'duplication_rate': d, 'loss_rate': l, 'metric': metric,
            'raw_trials': len(cell), 'valid_observations': len(valid),
            'extinct_trials': sum(r['extinct']=='True' for r in cell),
            'mean': avg, 'sample_sd': sd, 'mcse': mcse,
            'normal_ci_lower': None if ci is None else ci[0],
            'normal_ci_upper': None if ci is None else ci[1],
            'summary_mean_abs_error': mean_error,
            'summary_sd_abs_error': sd_error,
        })

# Add bounded alternatives: Wilson for extinction; deterministic bootstrap for bounded per-trial metrics.
for record in validation:
    d,l,metric = record['duplication_rate'], record['loss_rate'], record['metric']
    cell = cells[(d,l)]
    values = [float(r[metric]) for r in cell if r[metric] != '']
    if metric == 'extinction_frequency':
        pass
    elif metric in BOUNDED:
        seed = int(hashlib.sha256(f'{d}:{l}:{metric}'.encode()).hexdigest()[:8], 16)
        low, high = bootstrap_ci(values, seed)
        record['recommended_interval'] = 'bootstrap_percentile'
        record['recommended_ci_lower'] = low
        record['recommended_ci_upper'] = high
    else:
        record['recommended_interval'] = 'normal_descriptive'
        record['recommended_ci_lower'] = record['normal_ci_lower']
        record['recommended_ci_upper'] = record['normal_ci_upper']
# Extinction is not in the original metric list; append it with Wilson and MCSE.
for (d,l), cell in sorted(cells.items()):
    n = len(cell)
    extinct = [r['extinct'] == 'True' for r in cell]
    successes = sum(extinct)
    p, sd, mcse = extinction_statistics(extinct)
    low, high = wilson(successes,n)
    mean_error, sd_error = compare_summary_record(
        summary_lookup, (d, l, 'extinction_frequency'), p, sd
    )
    validation.append({'duplication_rate':d,'loss_rate':l,'metric':'extinction_frequency','raw_trials':n,'valid_observations':n,'extinct_trials':successes,'mean':p,'sample_sd':sd,'mcse':mcse,'normal_ci_lower':p-Z*mcse,'normal_ci_upper':p+Z*mcse,'summary_mean_abs_error':mean_error,'summary_sd_abs_error':sd_error,'recommended_interval':'wilson_score','recommended_ci_lower':low,'recommended_ci_upper':high})

verify_summary_errors(validation)

fields=['duplication_rate','loss_rate','metric','raw_trials','valid_observations','extinct_trials','mean','sample_sd','mcse','normal_ci_lower','normal_ci_upper','summary_mean_abs_error','summary_sd_abs_error','recommended_interval','recommended_ci_lower','recommended_ci_upper']
with AUDIT_SUMMARY.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(validation)

# Factorial descriptive effects using complete-history metrics over all 900 rows and inferred cost over valid rows.
effect_rows=[]
for metric in ('complete_duplications','complete_losses','inferred_cost'):
    vals_by_cell={(d,l):[float(r[metric]) for r in cell if r[metric] != ''] for (d,l),cell in cells.items()}
    grand=mean([x for v in vals_by_cell.values() for x in v])
    dmeans={d:mean([x for (dd,l),v in vals_by_cell.items() if dd==d for x in v]) for d in sorted({d for d,l in vals_by_cell})}
    lmeans={l:mean([x for (d,ll),v in vals_by_cell.items() if ll==l for x in v]) for l in sorted({l for d,l in vals_by_cell})}
    for d,dm in dmeans.items(): effect_rows.append({'metric':metric,'effect':'duplication_rate_marginal','level':d,'mean':dm,'difference_from_grand_mean':dm-grand})
    for l,lm in lmeans.items(): effect_rows.append({'metric':metric,'effect':'loss_rate_marginal','level':l,'mean':lm,'difference_from_grand_mean':lm-grand})
    # Two-by-two interaction contrasts on the low/high corners and center/high-low contrasts.
    low, high = min(dmeans), max(dmeans)
    for l in sorted(lmeans):
        contrast=(vals_by_cell[(high,l)][0] if False else mean(vals_by_cell[(high,l)])) - mean(vals_by_cell[(low,l)])
        effect_rows.append({'metric':metric,'effect':'duplication_high_minus_low_at_loss_rate','level':l,'mean':contrast,'difference_from_grand_mean':None})
    for d in sorted(dmeans):
        contrast=mean(vals_by_cell[(d,high)]) - mean(vals_by_cell[(d,low)])
        effect_rows.append({'metric':metric,'effect':'loss_high_minus_low_at_duplication_rate','level':d,'mean':contrast,'difference_from_grand_mean':None})
    interaction=(mean(vals_by_cell[(high,high)])-mean(vals_by_cell[(high,low)]))-(mean(vals_by_cell[(low,high)])-mean(vals_by_cell[(low,low)]))
    effect_rows.append({'metric':metric,'effect':'high_high_interaction_contrast','level':'high-vs-low','mean':interaction,'difference_from_grand_mean':None})
with EFFECTS.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['metric','effect','level','mean','difference_from_grand_mean']); w.writeheader(); w.writerows(effect_rows)

# Verify figure/table means against raw-derived values.
figure_check=[]
for metric in ('complete_duplications','complete_losses','inferred_cost'):
    for (d,l),cell in sorted(cells.items()):
        vals=[float(r[metric]) for r in cell if r[metric] != '']
        figure_check.append({'metric':metric,'duplication_rate':d,'loss_rate':l,'raw_mean':mean(vals) if vals else None,'summary_mean':float(summary_lookup[(d,l,metric)]['mean']) if vals else None,'matches':abs(mean(vals)-float(summary_lookup[(d,l,metric)]['mean']))<1e-12 if vals else True})
assert all(x['matches'] for x in figure_check)

lines=['# Statistical audit: second PhyloBench factorial study','','## Scope','',f'Raw data contain {len(rows)} trials across {len(cells)} cells. Every cell contains exactly 100 trials. This audit reads existing CSV artifacts only and does not rerun simulation.','','## Statistical verification','', 'Means and sample standard deviations match the raw CSV values to floating-point tolerance. MCSE is sample SD divided by the square root of valid observations. The audit CSV reports sample SD with the n-1 denominator for extinction frequency, while its MCSE retains the plug-in Bernoulli convention sqrt(p * (1 - p)) / sqrt(n). Original summary tables use the normal approximation mean +/- 1.96 * MCSE.', '', '| Cell | Extinct | Valid inferred observations |', '|---|---:|---:|']
for (d,l),cell in sorted(cells.items()):
 extinct=sum(r['extinct']=='True' for r in cell); valid=100-extinct
 lines.append(f'| d={d:g}, l={l:g} | {extinct}/100 | {valid} |')
lines += ['', '## Confidence-interval recommendation','', 'Clipping a normal interval to [0,1] is a pragmatic display constraint, not a statistically principled interval construction; it can conceal skewness and boundary behavior. The audited reporting layer therefore uses Wilson score intervals for extinction frequency and deterministic percentile-bootstrap intervals for bounded per-trial metrics (location agreement, duplication compatibility, and loss precision/recall). Unbounded counts and costs retain descriptive normal-approximation intervals. The original summary CSV is preserved; `summary/statistical_audit.csv` contains the recommended intervals and method labels.', '', '## Factorial effects','', 'Marginal means and high-minus-low contrasts are descriptive summaries of the balanced design. They show strong duplication-rate effects on duplication counts and cost, strong loss-rate effects on loss counts and extinction, and rate-dependent changes in inferred cost. Interaction contrasts are included in `summary/factorial_effects.csv`; they are not treated as formal hypothesis tests because this study is a stochastic simulation benchmark with model-specific estimands and no preregistered inferential model.', '', '## Figures and tables','', 'All cell means used by the heatmaps and main-effect figures match raw-data recomputation to floating-point tolerance. Raw trial counts, valid denominators, and extinction denominators are retained in the audit CSV. The original and second-study artifacts remain separate.', '', '## Limitations','', 'The intervals quantify Monte Carlo uncertainty conditional on this simulator, topology, branch lengths, rates, and fixed seed set; they do not quantify biological or model uncertainty. Bootstrap intervals inherit the empirical distribution and are not guarantees of exact coverage. Per-trial bounded metrics are ratios produced by the reconciliation evaluator, so pooled binomial Wilson intervals are not appropriate without numerator/denominator event counts. Extinction makes inferred metrics conditionally estimands over surviving families, not unconditional population means.']
REPORT.write_text('\n'.join(lines)+'\n')
print('audit_complete',len(rows),len(cells),'summary_rows',len(validation),'effects_rows',len(effect_rows),'figure_checks',len(figure_check))
