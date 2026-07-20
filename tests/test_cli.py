"""End-to-end CLI smoke tests (local VMR only, no network)."""

import json
from pathlib import Path

from typer.testing import CliRunner

from viralfetch import config as config_mod
from viralfetch.cache import SEQS, Cache
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


def test_seq_taxon_aggregate_is_local_json():
    # --taxon --meta is a local aggregate: no network, no email needed.
    result = runner.invoke(app, ["--json", "seq", "--taxon", "Coronaviridae", "--meta"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload["rank"] == "family"
    assert payload["species"] > 0
    assert payload["accessions"] > 0
    assert "moltype_breakdown" in payload


def test_seq_requires_species_or_taxon():
    result = runner.invoke(app, ["seq"])
    assert result.exit_code == 2


def test_seq_rejects_both_species_and_taxon():
    result = runner.invoke(app, ["seq", "Coronaviridae", "--taxon", "Coronaviridae"])
    assert result.exit_code == 2


# -- Phase 6 utilities ----------------------------------------------------

def test_diagnose_json():
    result = runner.invoke(app, ["--json", "diagnose"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload["isolates"] > 0
    assert payload["accessions"] > 0
    assert "unparsed" in payload


def test_config_show_masks_api_key(monkeypatch, tmp_path):
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    monkeypatch.setattr(config_mod, "CONFIG_FILE", tmp_path / "config.json")
    monkeypatch.setattr(config_mod, "CACHE_DIR", tmp_path / "cache")
    result = runner.invoke(app, ["--email", "a@b.co", "--api-key", "secretKEY1234", "config"])
    assert result.exit_code == 0
    assert "a@b.co" in result.stdout
    assert "secretKEY1234" not in result.stdout  # masked
    assert "1234" in result.stdout  # last 4 shown


def test_config_store_persists(monkeypatch, tmp_path):
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    monkeypatch.setattr(config_mod, "CONFIG_DIR", tmp_path)
    monkeypatch.setattr(config_mod, "CONFIG_FILE", tmp_path / "config.json")
    monkeypatch.setattr(config_mod, "CACHE_DIR", tmp_path / "cache")
    result = runner.invoke(app, ["--json", "config", "--store-ncbi-email", "stored@x.io"])
    assert result.exit_code == 0
    saved = json.loads((tmp_path / "config.json").read_text())
    assert saved["email"] == "stored@x.io"
    assert json.loads(result.stdout)["email"] == "stored@x.io"


def test_cache_info_json(monkeypatch, tmp_path):
    monkeypatch.setattr(config_mod, "CACHE_DIR", tmp_path / "cache")
    result = runner.invoke(app, ["--json", "cache", "info"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert "seqs" in payload and "texts" in payload


def test_cache_clear_removes_entries(monkeypatch, tmp_path):
    cache_dir = tmp_path / "cache"
    monkeypatch.setattr(config_mod, "CACHE_DIR", cache_dir)
    Cache(cache_dir).set(SEQS, "k", "v")  # seed one entry
    result = runner.invoke(app, ["--json", "cache", "clear", "--seqs"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload["cleared"] == 1
    assert payload["scope"] == "seqs"
