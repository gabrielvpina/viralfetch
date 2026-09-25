"""Installing a newer VMR from an ICTV workbook, and reverting to the bundled one."""

import io
import json

import openpyxl
import pytest
from typer.testing import CliRunner

from viralfetch import cli
from viralfetch import config as config_mod
from viralfetch import vmr as vmr_mod
from viralfetch import vmr_install
from viralfetch.cli import app
from viralfetch.ictv import ICTVClient, VMRUpdate

from .conftest import HEADER, _ROWS

runner = CliRunner()

NEW = "VMR_MSL99.v1.20990101.xlsx"  # newer than any bundled release
NEW_URL = f"https://ictv.global/sites/default/files/VMR/{NEW}"


def _workbook(rows=_ROWS, header=HEADER) -> bytes:
    """An ICTV-like workbook: a notes sheet first, then the VMR sheet."""
    wb = openpyxl.Workbook()
    notes = wb.active
    notes.title = "Version"
    notes.append([None, "ICTV Virus Metadata Resource"])
    ws = wb.create_sheet("VMR MSL99")
    ws.append(list(header) + [None, None])  # formatting-only trailing columns
    for i, row in enumerate(rows, start=1):
        values = [row[col] or None for col in header]
        values[header.index("Species Sort")] = float(i)  # Excel numbers are floats
        ws.append(values + [None, None])
    ws.append([None] * (len(header) + 2))  # blank trailing row
    buf = io.BytesIO()
    wb.save(buf)
    return buf.getvalue()


@pytest.fixture
def data_dir(monkeypatch, tmp_path):
    monkeypatch.setattr(config_mod, "DATA_DIR", tmp_path / "data")
    vmr_mod.load.cache_clear()
    yield tmp_path / "data"
    vmr_mod.load.cache_clear()  # never leak an installed VMR into other tests


def test_xlsx_to_tsv_finds_vmr_sheet_and_normalises(tmp_path):
    dest = tmp_path / "out.tsv"
    vmr_install.xlsx_to_tsv(_workbook(), dest)
    lines = dest.read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t") == HEADER        # trailing empty columns dropped
    assert len(lines) == 1 + len(_ROWS)          # blank row skipped
    first = dict(zip(HEADER, lines[1].split("\t")))
    assert first["Species Sort"] == "1"          # not "1.0"
    assert first["Species"] == "Alphavirus one"


def test_xlsx_without_vmr_columns_is_rejected(tmp_path):
    data = _workbook(header=[h for h in HEADER if h != "Species"])
    with pytest.raises(vmr_install.VMRInstallError, match="expected VMR columns"):
        vmr_install.xlsx_to_tsv(data, tmp_path / "out.tsv")


def test_non_xlsx_is_rejected(tmp_path):
    with pytest.raises(vmr_install.VMRInstallError, match="not a readable"):
        vmr_install.xlsx_to_tsv(b"<html>not a workbook</html>", tmp_path / "out.tsv")


def test_install_becomes_active_and_reset_reverts(data_dir):
    result = vmr_install.install(_workbook(), NEW)
    assert result.filename == "VMR_MSL99.v1.20990101.tsv"
    assert result.isolates == len(_ROWS)
    assert vmr_mod.active_filename() == result.filename
    assert vmr_mod.load().find("Alphaviridae") is not None

    assert vmr_install.reset() == [result.filename]
    assert vmr_mod.active_filename() == vmr_mod.VMR_FILENAME
    assert vmr_install.reset() == []


def test_install_replaces_previous_install(data_dir):
    vmr_install.install(_workbook(), "VMR_MSL98.v1.20980101.xlsx")
    vmr_install.install(_workbook(), NEW)
    assert [p.name for p in vmr_mod.installed_dir().iterdir()] == ["VMR_MSL99.v1.20990101.tsv"]


def test_failed_install_keeps_previous(data_dir):
    vmr_install.install(_workbook(), NEW)
    with pytest.raises(vmr_install.VMRInstallError):
        vmr_install.install(b"garbage", "VMR_MSL100.v1.21000101.xlsx")
    assert [p.name for p in vmr_mod.installed_dir().iterdir()] == ["VMR_MSL99.v1.20990101.tsv"]


def test_install_older_than_bundled_is_ignored(data_dir):
    vmr_install.install(_workbook(), "VMR_MSL30.v1.20150101.xlsx")
    assert vmr_mod.active_filename() == vmr_mod.VMR_FILENAME


# -- CLI --------------------------------------------------------------------

@pytest.fixture
def fake_ictv(monkeypatch, data_dir):
    """Stub the network: a newer VMR is published, and downloading returns it."""
    monkeypatch.setenv("NCBI_EMAIL", "test@example.com")

    def check(self, current):
        latest_newer = vmr_mod.release_key(NEW) > vmr_mod.release_key(current)
        return VMRUpdate(current=current, latest=NEW, latest_url=NEW_URL, up_to_date=not latest_newer)

    monkeypatch.setattr(ICTVClient, "check_vmr_update", check)
    monkeypatch.setattr(ICTVClient, "download_vmr", lambda self, url: _workbook())
    return data_dir


def test_update_yes_installs(fake_ictv):
    result = runner.invoke(app, ["update", "--yes"])
    assert result.exit_code == 0, result.output
    assert "Installed VMR_MSL99.v1.20990101.tsv" in result.stdout
    assert vmr_mod.active_filename() == "VMR_MSL99.v1.20990101.tsv"
    # a second check sees the installed release as current
    again = runner.invoke(app, ["update"])
    assert "up to date (VMR_MSL99.v1.20990101.tsv)" in again.stdout


def test_update_prompt_declined_installs_nothing(fake_ictv, monkeypatch):
    monkeypatch.setattr(cli, "_interactive", lambda: True)
    result = runner.invoke(app, ["update"], input="n\n")
    assert result.exit_code == 0
    assert "Download and install it?" in result.stdout
    assert vmr_mod.installed_path() is None


def test_update_prompt_accepted_installs(fake_ictv, monkeypatch):
    monkeypatch.setattr(cli, "_interactive", lambda: True)
    result = runner.invoke(app, ["update"], input="y\n")
    assert result.exit_code == 0
    assert vmr_mod.installed_path() is not None


def test_update_non_interactive_without_yes_only_reports(fake_ictv, monkeypatch):
    monkeypatch.setattr(cli, "_interactive", lambda: False)
    result = runner.invoke(app, ["update"])
    assert result.exit_code == 0
    assert "--yes" in result.stderr
    assert vmr_mod.installed_path() is None


def test_update_json_yes_emits_single_object(fake_ictv):
    result = runner.invoke(app, ["--json", "update", "--yes"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload["up_to_date"] is False
    assert payload["installed"]["filename"] == "VMR_MSL99.v1.20990101.tsv"


def test_update_download_failure_exits_4_and_keeps_vmr(fake_ictv, monkeypatch):
    monkeypatch.setattr(ICTVClient, "download_vmr", lambda self, url: b"garbage")
    result = runner.invoke(app, ["update", "--yes"])
    assert result.exit_code == 4
    assert "current VMR is unchanged" in result.stderr
    assert vmr_mod.installed_path() is None


def test_update_reset_needs_no_email_or_network(data_dir, monkeypatch, tmp_path):
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    monkeypatch.setattr(config_mod, "CONFIG_FILE", tmp_path / "config.json")
    vmr_install.install(_workbook(), NEW)
    result = runner.invoke(app, ["update", "--reset"])
    assert result.exit_code == 0
    assert "using the bundled VMR" in result.stdout
    assert vmr_mod.installed_path() is None


def test_update_reset_and_yes_conflict(data_dir):
    result = runner.invoke(app, ["update", "--reset", "--yes"])
    assert result.exit_code == 2
