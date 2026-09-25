"""Install a newer VMR downloaded from ICTV, or revert to the embedded one.

ICTV publishes the VMR as an ``.xlsx`` workbook with several sheets; the VMR
sheet is converted to the same TSV layout as the embedded file (a conversion of
the embedded release reproduces it byte for byte) and validated by loading it
before it replaces anything. The file lands in the user data dir, never in the
installed package, so it survives reinstalls and needs no write access there.

Design note (SPEC section 4): nothing in this module prints. It returns data.
"""

from __future__ import annotations

import csv
import io
import os
from dataclasses import dataclass
from pathlib import Path

from . import vmr as vmr_mod

# Columns the loader depends on; a sheet without them is not a VMR.
REQUIRED_COLUMNS = (
    "Isolate ID", "Realm", "Family", "Genus", "Species",
    "Exemplar or additional isolate", "Virus GENBANK accession",
)
_HEADER_SEARCH_ROWS = 10  # how far down a sheet to look for the header row


class VMRInstallError(Exception):
    """The downloaded workbook could not be turned into a usable VMR."""


@dataclass
class InstalledVMR:
    filename: str
    path: Path
    isolates: int
    species: int


def _cell(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and value.is_integer():
        value = int(value)  # "Species Sort" etc. come back as floats
    return str(value)


def xlsx_to_tsv(data: bytes, dest: Path) -> None:
    """Write the VMR sheet of an ICTV workbook to ``dest`` as TSV.

    The sheet is found by its header row (containing :data:`REQUIRED_COLUMNS`),
    not by name, since the sheet title changes with each release.
    """
    import openpyxl  # heavy; only needed here

    try:
        wb = openpyxl.load_workbook(io.BytesIO(data), read_only=True, data_only=True)
    except Exception as exc:  # zipfile.BadZipFile, KeyError, InvalidFileException...
        raise VMRInstallError(f"not a readable .xlsx workbook: {exc}") from exc
    try:
        for ws in wb.worksheets:
            rows = ws.iter_rows(values_only=True)
            for _ in range(_HEADER_SEARCH_ROWS):
                header = next(rows, None)
                if header is None:
                    break
                if set(REQUIRED_COLUMNS) <= {_cell(v).strip() for v in header}:
                    _write_tsv(header, rows, dest)
                    return
        raise VMRInstallError(
            "no sheet with the expected VMR columns "
            f"({', '.join(REQUIRED_COLUMNS)}) — the ICTV format may have changed"
        )
    finally:
        wb.close()


def _write_tsv(header, rows, dest: Path) -> None:
    # Trailing formatting-only columns have an empty header: drop them.
    width = next((i for i, v in enumerate(header) if v is None), len(header))
    with dest.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow([_cell(v) for v in header[:width]])
        for row in rows:
            row = row[:width]
            if all(v is None for v in row):
                continue
            writer.writerow([_cell(v) for v in row])


def install(data: bytes, xlsx_name: str) -> InstalledVMR:
    """Convert, validate and install a downloaded VMR workbook.

    The previous install (if any) is only replaced once the new file has loaded
    successfully, so a bad download never leaves the tool without a VMR.
    """
    filename = Path(xlsx_name).stem + ".tsv"
    target_dir = vmr_mod.installed_dir()
    target_dir.mkdir(parents=True, exist_ok=True)
    tmp = target_dir / f".{filename}.partial"
    try:
        xlsx_to_tsv(data, tmp)
        loaded = vmr_mod.load_from(tmp)
        if not loaded.isolates or not loaded.species:
            raise VMRInstallError("the downloaded VMR has no isolates or species")
        final = target_dir / filename
        os.replace(tmp, final)
    finally:
        tmp.unlink(missing_ok=True)

    for old in target_dir.glob("VMR_*.tsv"):
        if old != final:
            old.unlink()
    vmr_mod.load.cache_clear()
    return InstalledVMR(
        filename=filename, path=final,
        isolates=len(loaded.isolates), species=len(loaded.species),
    )


def reset() -> list[str]:
    """Remove any installed VMR, reverting to the embedded one.

    Returns the basenames removed (empty if nothing was installed).
    """
    removed = []
    for path in vmr_mod.installed_dir().glob("VMR_*.tsv"):
        path.unlink()
        removed.append(path.name)
    vmr_mod.load.cache_clear()
    return removed
