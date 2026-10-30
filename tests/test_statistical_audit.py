"""Tests for the side-effect-free statistical audit calculations."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest


AUDIT_HELPERS = (
    Path(__file__).resolve().parents[1]
    / "experiments"
    / "2026-10-30-second-study"
    / "analysis"
)
sys.path.insert(0, str(AUDIT_HELPERS))

from audit_helpers import (  # noqa: E402
    compare_summary_record,
    extinction_statistics,
    verify_summary_errors,
)


def test_extinction_mean_and_sample_sd_compare_to_stored_summary():
    probability, sample_sd, mcse = extinction_statistics(
        [True, False, True, False]
    )
    key = (0.5, 0.1, "extinction_frequency")
    summary = {
        key: {
            "mean": "0.5",
            "standard_deviation": str(sample_sd),
        }
    }

    errors = compare_summary_record(summary, key, probability, sample_sd)

    assert errors == (0.0, 0.0)
    assert mcse == pytest.approx(0.25)


def test_incorrect_stored_extinction_mean_is_reported_and_fails_verification():
    key = (0.5, 0.1, "extinction_frequency")
    summary = {key: {"mean": "0.4", "standard_deviation": "0.5"}}
    mean_error, sd_error = compare_summary_record(summary, key, 0.5, 0.5)
    records = [{
        "duplication_rate": 0.5,
        "loss_rate": 0.1,
        "metric": "extinction_frequency",
        "summary_mean_abs_error": mean_error,
        "summary_sd_abs_error": sd_error,
    }]

    assert mean_error == pytest.approx(0.1)
    assert sd_error == 0.0
    with pytest.raises(ValueError, match="stored summary verification failed"):
        verify_summary_errors(records)


def test_missing_stored_summary_record_fails_clearly():
    key = (0.5, 0.1, "extinction_frequency")

    with pytest.raises(ValueError, match="missing stored summary record"):
        compare_summary_record({}, key, 0.5, 0.5)


def test_extinction_sample_sd_uses_n_minus_one_and_mcse_keeps_prior_method():
    probability, sample_sd, mcse = extinction_statistics(
        [True, False, True, False]
    )

    assert probability == 0.5
    assert math.isclose(sample_sd, math.sqrt(1.0 / 3.0))
    assert not math.isclose(sample_sd, math.sqrt(probability * (1 - probability)))
    assert mcse == 0.25


def test_zero_summary_error_is_not_fabricated_for_a_mismatch():
    key = (0.1, 0.5, "extinction_frequency")
    mean_error, sd_error = compare_summary_record(
        {key: {"mean": "0.2", "standard_deviation": "0.4"}},
        key,
        0.25,
        0.4,
    )

    assert mean_error == pytest.approx(0.05)
    assert sd_error == 0.0


def test_nonfinite_stored_summary_values_are_rejected():
    key = (0.1, 0.5, "extinction_frequency")

    with pytest.raises(ValueError, match="non-finite summary value"):
        compare_summary_record(
            {key: {"mean": "nan", "standard_deviation": "0.4"}},
            key,
            0.25,
            0.4,
        )


def test_unperformed_summary_comparison_cannot_pass():
    with pytest.raises(ValueError, match="summary comparison was not performed"):
        verify_summary_errors([{
            "duplication_rate": 0.1,
            "loss_rate": 0.5,
            "metric": "extinction_frequency",
            "summary_mean_abs_error": None,
            "summary_sd_abs_error": None,
        }])
