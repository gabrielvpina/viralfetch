"""End-to-end CLI smoke tests (local VMR only, no network)."""

import json

from typer.testing import CliRunner

from viralfetch.cli import app

runner = CliRunner()  # Click >= 8.2 keeps stderr separate by default


def test_tax_json_output():
    result = runner.invoke(app, ["--json", "tax", "Coronaviridae"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload["name"] == "Coronaviridae"
    assert payload["rank"] == "family"
    assert payload["lineage"]["realm"] == "Riboviria"


def test_tax_rich_output_runs():
    result = runner.invoke(app, ["tax", "Coronaviridae"])
    assert result.exit_code == 0
    assert "Coronaviridae" in result.stdout


def test_members_json_genus():
    result = runner.invoke(app, ["--json", "members", "Coronaviridae", "--rank", "genus"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    names = {m["name"] for m in payload["members"]}
    assert "Betacoronavirus" in names


def test_members_tree_json():
    result = runner.invoke(app, ["--json", "members", "Coronaviridae", "--tree"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload["tree"]["name"] == "Coronaviridae"
    assert payload["total"] > 0
    child_ranks = {c["rank"] for c in payload["tree"]["children"]}
    assert child_ranks <= {"subfamily", "genus"}  # next populated rank(s) down


def test_members_tree_rich_runs():
    result = runner.invoke(app, ["members", "Coronaviridae", "--tree"])
    assert result.exit_code == 0
    assert "Coronaviridae" in result.stdout
    assert "descendant taxa" in result.stdout


def test_tax_not_found_exit_code_and_stderr():
    result = runner.invoke(app, ["--json", "tax", "CoronaviridaX"])
    assert result.exit_code == 1
    # error payload goes to stderr, stdout stays clean for jq
    assert result.stdout.strip() == ""
    assert "taxon_not_found" in result.stderr


def test_members_invalid_rank_exit_code():
    result = runner.invoke(app, ["members", "Coronaviridae", "--rank", "realm"])
    assert result.exit_code == 2
