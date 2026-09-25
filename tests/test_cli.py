"""End-to-end CLI smoke tests (local VMR only, no network)."""

import json
from pathlib import Path

from typer.testing import CliRunner

from viralfetch import compare
from viralfetch import config as config_mod
from viralfetch.cache import SEQS, Cache
from viralfetch.cli import app
from viralfetch.ncbi import NcbiLineage

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


def test_tax_not_found_exit_code_and_stderr(monkeypatch):
    # NCBI fallback finds nothing either (stubbed: no network).
    monkeypatch.setenv("NCBI_EMAIL", "test@example.com")
    monkeypatch.setattr(compare, "lineage_via_ncbi", lambda client, name: None)
    result = runner.invoke(app, ["--json", "tax", "CoronaviridaX"])
    assert result.exit_code == 1
    # error payload goes to stderr, stdout stays clean for jq
    assert result.stdout.strip() == ""
    assert "taxon_not_found" in result.stderr


def test_tax_falls_back_to_ncbi_when_not_in_vmr(monkeypatch):
    monkeypatch.setenv("NCBI_EMAIL", "test@example.com")
    lineage = NcbiLineage(
        taxid="12345", name="Fooviridae", rank="family",
        lineage=[("realm", "Riboviria"), ("family", "Fooviridae")],
    )
    monkeypatch.setattr(compare, "lineage_via_ncbi", lambda client, name: lineage)
    result = runner.invoke(app, ["--json", "tax", "Fooviridae"])
    assert result.exit_code == 0
    assert "12345" in result.stdout
    assert "not in the local VMR" in result.stderr


def test_tax_fallback_without_email_warns_and_not_found(monkeypatch, tmp_path):
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    monkeypatch.setattr(config_mod, "CONFIG_FILE", tmp_path / "config.json")
    result = runner.invoke(app, ["--json", "tax", "CoronaviridaX"])
    assert result.exit_code == 1
    assert "no NCBI email is configured" in result.stderr
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


def test_missing_email_warns_on_every_command(monkeypatch, tmp_path):
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    monkeypatch.setattr(config_mod, "CONFIG_FILE", tmp_path / "config.json")
    result = runner.invoke(app, ["--json", "tax", "Coronaviridae"])
    assert result.exit_code == 0
    assert "NCBI email is not configured" in result.stderr
    assert "--store-ncbi-email" in result.stderr
    assert "export NCBI_EMAIL" in result.stderr
    json.loads(result.stdout)  # stdout stays pure JSON


def test_configured_email_does_not_warn(monkeypatch):
    monkeypatch.setenv("NCBI_EMAIL", "test@example.com")
    result = runner.invoke(app, ["tax", "Coronaviridae"])
    assert result.exit_code == 0
    assert "not configured" not in result.stderr


def test_config_command_has_single_email_warning(monkeypatch, tmp_path):
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    monkeypatch.setattr(config_mod, "CONFIG_FILE", tmp_path / "config.json")
    monkeypatch.setattr(config_mod, "CACHE_DIR", tmp_path / "cache")
    result = runner.invoke(app, ["config"])
    assert result.exit_code == 0
    assert "NCBI email is not configured" not in result.stderr
    assert "No NCBI email is stored" in result.stderr


def _stub_ncbi_lineage(monkeypatch, pairs):
    monkeypatch.setenv("NCBI_EMAIL", "test@example.com")
    lineage = NcbiLineage(taxid="1", name=pairs[-1][1], rank=pairs[-1][0], lineage=pairs)
    monkeypatch.setattr(compare, "lineage_via_ncbi", lambda client, name: lineage)


def test_members_redirects_unknown_name_via_ncbi(monkeypatch):
    _stub_ncbi_lineage(monkeypatch, [
        ("family", "Retroviridae"), ("genus", "Lentivirus"), ("no rank", "Primate lentivirus group"),
    ])
    result = runner.invoke(app, ["--json", "members", "Primate lentivirus group", "--rank", "species"])
    assert result.exit_code == 0
    assert json.loads(result.stdout)["parent"] == {"name": "Lentivirus", "rank": "genus"}
    assert "NCBI places it in genus Lentivirus" in result.stderr


def test_seq_taxon_redirects_unknown_name_via_ncbi(monkeypatch):
    _stub_ncbi_lineage(monkeypatch, [
        ("family", "Retroviridae"), ("genus", "Lentivirus"), ("no rank", "Primate lentivirus group"),
    ])
    result = runner.invoke(app, ["--json", "seq", "--taxon", "Primate lentivirus group", "--meta"])
    assert result.exit_code == 0
    assert json.loads(result.stdout)["name"] == "Lentivirus"


def test_seq_species_redirected_to_genus_asks_for_taxon(monkeypatch):
    # NCBI only places it in a genus: `seq <species>` then points to --taxon.
    _stub_ncbi_lineage(monkeypatch, [("genus", "Lentivirus"), ("no rank", "Some lentivirus")])
    result = runner.invoke(app, ["seq", "Some lentivirus"])
    assert result.exit_code == 2
    assert "Use --taxon" in result.stderr


def test_members_unknown_placed_only_above_family_is_not_found(monkeypatch):
    _stub_ncbi_lineage(monkeypatch, [("realm", "Riboviria"), ("no rank", "unclassified Riboviria")])
    result = runner.invoke(app, ["--json", "members", "unclassified Riboviria"])
    assert result.exit_code == 1
    assert "taxon_not_found" in result.stderr
