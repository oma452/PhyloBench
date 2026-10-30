"""CLI tests using small deterministic experiments."""

import json
from pathlib import Path

import pytest

from phylo_bench.cli import main


def _config(path: Path) -> None:
    path.write_text(
        json.dumps(
            {
                "species_tree_newick": "(A:1,B:1)S;",
                "seeds": [7, 11],
                "duplication_rate": 0.0,
                "loss_rate": 0.0,
                "max_events": 100,
                "max_lineages": 100,
            }
        ),
        encoding="utf-8",
    )


def test_cli_run_exports_reproducible_artifacts(tmp_path: Path, capsys):
    config = tmp_path / "experiment.json"
    output = tmp_path / "results"
    _config(config)

    assert main(["run", str(config), "--output-dir", str(output)]) == 0
    captured = capsys.readouterr()
    assert "completed 2 trials" in captured.out
    assert (output / "trials.csv").exists()
    assert (output / "summary.csv").exists()
    assert (output / "metadata.json").exists()
    assert (output / "report.md").exists()
    assert sorted(p.name for p in (output / "figures").glob("*.png")) == [
        "complete_duplications.png",
        "complete_losses.png",
        "inferred_cost.png",
    ]
    metadata = json.loads((output / "metadata.json").read_text())
    assert metadata["configuration"]["seeds"] == [7, 11]
    assert metadata["summary"]["trials"] == 2


def test_cli_inspect_prints_normalized_configuration(tmp_path: Path, capsys):
    config = tmp_path / "experiment.json"
    _config(config)

    assert main(["inspect", str(config)]) == 0
    output = json.loads(capsys.readouterr().out)
    assert output["seeds"] == [7, 11]
    assert output["max_events"] == 100


@pytest.mark.parametrize(
    ("newick", "message"),
    [
        ("((A:1,B:1)S;", "invalid Newick input"),
        ("(A:0,B:1)S;", "positive finite branch lengths"),
        ("A;", "at least two species"),
        ("(A:1,B:1,C:1)S;", "not binary"),
    ],
)
def test_cli_inspect_rejects_invalid_simulator_trees(
    tmp_path: Path, capsys, newick: str, message: str
):
    config = tmp_path / "experiment.json"
    _config(config)
    data = json.loads(config.read_text(encoding="utf-8"))
    data["species_tree_newick"] = newick
    config.write_text(json.dumps(data), encoding="utf-8")

    with pytest.raises(SystemExit) as error:
        main(["inspect", str(config)])

    assert error.value.code == 2
    captured = capsys.readouterr()
    assert "invalid experiment configuration" in captured.err
    assert message in captured.err


def test_cli_reports_helpful_validation_errors(tmp_path: Path):
    config = tmp_path / "invalid.json"
    config.write_text(
        json.dumps(
            {
                "species_tree_newick": "(A:1,B:1)S;",
                "seeds": [1, 1],
                "duplication_rate": 0.0,
                "loss_rate": 0.0,
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(SystemExit) as error:
        main(["inspect", str(config)])
    assert error.value.code == 2


def test_cli_help_documents_commands(capsys):
    with pytest.raises(SystemExit) as error:
        main(["--help"])
    assert error.value.code == 0
    help_text = capsys.readouterr().out
    assert "run" in help_text
    assert "inspect" in help_text
