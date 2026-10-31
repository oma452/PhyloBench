"""Small, side-effect-free helpers used by the Study 2 audit script."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from math import isfinite, sqrt
from statistics import mean, stdev
from typing import Any


SUMMARY_TOLERANCE = 1e-12


def extinction_statistics(extinct: Sequence[bool]) -> tuple[float, float, float]:
    """Return extinction mean, sample SD, and plug-in Bernoulli MCSE.

    The sample SD follows the same n-1 convention as the other per-trial
    summary metrics. The existing extinction MCSE convention is retained:
    sqrt(p * (1 - p)) / sqrt(n).
    """

    if not extinct:
        raise ValueError("cannot summarize extinction for zero trials")
    if any(not isinstance(value, bool) for value in extinct):
        raise ValueError("extinction observations must be booleans")
    values = [float(value) for value in extinct]
    probability = mean(values)
    sample_sd = stdev(values) if len(values) > 1 else 0.0
    bernoulli_sd = sqrt(probability * (1.0 - probability))
    mcse = bernoulli_sd / sqrt(len(values))
    return probability, sample_sd, mcse


def compare_summary_record(
    summary_lookup: Mapping[tuple[float, float, str], Mapping[str, Any]],
    key: tuple[float, float, str],
    average: float,
    sample_sd: float,
) -> tuple[float, float]:
    """Compare recomputed mean and sample SD with a stored summary record."""

    try:
        stored = summary_lookup[key]
    except KeyError as exc:
        raise ValueError(f"missing stored summary record for {key!r}") from exc
    try:
        stored_mean = float(stored["mean"])
        stored_sd = float(stored["standard_deviation"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"invalid stored summary record for {key!r}") from exc
    if not all(isfinite(value) for value in (average, sample_sd, stored_mean, stored_sd)):
        raise ValueError(f"non-finite summary value for {key!r}")
    return abs(average - stored_mean), abs(sample_sd - stored_sd)


def verify_summary_errors(
    records: Sequence[Mapping[str, Any]],
    *,
    tolerance: float = SUMMARY_TOLERANCE,
) -> None:
    """Fail with cell/metric detail if stored means or SDs disagree."""

    mismatches = []
    for record in records:
        mean_error = record.get("summary_mean_abs_error")
        sd_error = record.get("summary_sd_abs_error")
        if mean_error is None or sd_error is None:
            mismatches.append(
                f"{record.get('metric')} at "
                f"({record.get('duplication_rate')}, {record.get('loss_rate')}): "
                "summary comparison was not performed"
            )
            continue
        if not isfinite(float(mean_error)) or not isfinite(float(sd_error)):
            mismatches.append(
                f"{record.get('metric')} at "
                f"({record.get('duplication_rate')}, {record.get('loss_rate')}): "
                "summary comparison produced a non-finite error"
            )
            continue
        if mean_error is not None and float(mean_error) > tolerance:
            mismatches.append(
                f"{record['metric']} at "
                f"({record['duplication_rate']}, {record['loss_rate']}): "
                f"mean absolute error {mean_error}"
            )
        if sd_error is not None and float(sd_error) > tolerance:
            mismatches.append(
                f"{record['metric']} at "
                f"({record['duplication_rate']}, {record['loss_rate']}): "
                f"sample SD absolute error {sd_error}"
            )
    if mismatches:
        details = "; ".join(mismatches)
        raise ValueError(f"stored summary verification failed: {details}")
