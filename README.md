# viralfetch

A command-line tool for querying and downloading viral taxonomy, metadata and
sequences. It combines three sources:

| Source | Content | Access |
|---|---|---|
| **VMR** (ICTV Virus Metadata Resource) | taxonomic hierarchy + exemplar isolates + GenBank/RefSeq accessions | **local**, embedded TSV |
| **NCBI E-utilities** | sequence metadata, sequences (nt/aa), NCBI taxonomy | **remote**, on demand |
| **ICTV Report** | descriptive chapter text per family/order | **remote**, on demand *(planned)* |

The VMR is the local index; everything else is fetched on demand and cached.

> **Status:** the taxonomy commands (`tax`, `members`) and the NCBI sequence
> command (`seq`) are implemented. The ICTV Report text command and the
> `update` / `cache` / `config` utilities are still in progress.

---

## Installation

```bash
pip install -e .
```

This installs the `viralfetch` command. Python 3.10+ is required.

## NCBI configuration

Commands that reach NCBI (`seq`, `tax --compare-ncbi`) require a real email
address, per NCBI usage policy. There is **no default** — the command fails
with an explanation if none is set.

```bash
export NCBI_EMAIL="you@example.com"
export NCBI_API_KEY="..."   # optional; raises the rate limit from 3 to 10 req/s
```

You can also pass `--email` / `--api-key` on any command.

## Global options

These go **before** the command:

| Option | Effect |
|---|---|
| `--json` | Emit pure JSON on stdout (warnings/errors go to stderr) — ready for `jq`. |
| `--no-cache` | Ignore the cache and refetch. |
| `--verbose` | Extra diagnostics on stderr. |
| `--email`, `--api-key` | Override the NCBI credentials for this run. |

Run `viralfetch COMMAND --help` to see a command's own arguments and options.

---

## `tax` — taxonomy lineage (local)

Show the full ICTV lineage of a taxon (realm → species). Case-insensitive,
with "did you mean" suggestions on a near miss. A species also gets an isolate
summary.

```bash
viralfetch tax Coronaviridae
```

```
╭─ Coronaviridae  (family) ─────────────────────────╮
│ lineage                                           │
│ └── realm: Riboviria                              │
│     └── kingdom: Orthornavirae                    │
│         └── phylum: Pisuviricota                  │
│             └── class: Pisoniviricetes            │
│                 └── order: Nidovirales            │
│                     └── suborder: Cornidovirineae │
│                         └── family: Coronaviridae │
╰───────────────────────────────────────────────────╯
```

Same query as JSON:

```bash
viralfetch --json tax Coronaviridae
```

```json
{
  "name": "Coronaviridae",
  "rank": "family",
  "lineage": {
    "realm": "Riboviria",
    "kingdom": "Orthornavirae",
    "phylum": "Pisuviricota",
    "class": "Pisoniviricetes",
    "order": "Nidovirales",
    "suborder": "Cornidovirineae",
    "family": "Coronaviridae"
  }
}
```

### `tax --compare-ncbi` (remote)

Fetch the NCBI taxonomy lineage for a representative accession and render it
beside the ICTV lineage, highlighting divergences. Divergences are expected —
NCBI commonly lags ICTV — and are the point of the command.

```bash
viralfetch tax "Betacoronavirus pandemicum" --compare-ncbi
```

```
             ICTV vs NCBI lineage — Betacoronavirus pandemicum
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ ICTV (VMR)                          ┃ NCBI (taxid 227984)                ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ realm: Riboviria                    │ acellular root: Viruses            │
│ kingdom: Orthornavirae              │ realm: Riboviria                   │
│ …                                   │ …                                  │
│ species: Betacoronavirus pandemicum │ species: Betacoronavirus pandemi…  │
│                                     │ no rank: SARS coronavirus Tor2     │
└─────────────────────────────────────┴────────────────────────────────────┘
```

---

## `members` — child taxa (local, no network)

List the taxa below a given taxon.

### At a specific rank

```bash
viralfetch members Coronaviridae --rank genus
```

```
genera in family Coronaviridae
┏━━━━━━━━━━━━━━━━━━┳━━━━━━━━━┓
┃ genus            ┃ species ┃
┡━━━━━━━━━━━━━━━━━━╇━━━━━━━━━┩
│ Alphacoronavirus │      26 │
│ Alphaletovirus   │       1 │
│ Alphapironavirus │       1 │
│ Betacoronavirus  │      14 │
│ Deltacoronavirus │       7 │
│ Gammacoronavirus │       5 │
└──────────────────┴─────────┘
           6 genera
```

### Counts only

```bash
viralfetch members Riboviria --rank family --count
```

```
158 families in realm Riboviria
```

### Per-rank breakdown (no flags)

```bash
viralfetch members Coronaviridae
```

```
  members of family
Coronaviridae by rank
┏━━━━━━━━━━━┳━━━━━━━┓
┃ rank      ┃ count ┃
┡━━━━━━━━━━━╇━━━━━━━┩
│ subfamily │     3 │
│ genus     │     6 │
│ subgenus  │    28 │
│ species   │    54 │
└───────────┴───────┘
Tip: add --tree to list every member of Coronaviridae as a hierarchy.
```

### Full descendant tree

```bash
viralfetch members Coronaviridae --tree
```

```
Coronaviridae  (family)
├── Letovirinae  (subfamily)
│   └── Alphaletovirus  (genus)
│       └── Milecovirus  (subgenus)
│           └── Alphaletovirus microhylae  (species)
├── Orthocoronavirinae  (subfamily)
│   ├── Alphacoronavirus  (genus)
│   │   ├── Amalacovirus  (subgenus)
│   │   │   └── Alphacoronavirus almalfi  (species)
│   │   └── …
│   └── …
└── …
91 descendant taxa
```

`--rank`, `--tree`, and the breakdown all work with `--json` too.

---

## `seq` — NCBI sequence data (remote)

Accessions are resolved locally from the VMR, then metadata or records are
fetched from NCBI. Output formats are mutually exclusive; `--meta` is the
default.

| Flag | Fetches |
|---|---|
| `--meta` | metadata via `esummary` (~1 KB/accession) |
| `--fasta` | FASTA sequences via `efetch` |
| `--gb` | full GenBank records via `efetch` |

### Metadata for a species

```bash
viralfetch seq "Betacoronavirus pandemicum" --meta
```

```
             nuccore metadata — Betacoronavirus pandemicum
┏━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━┳━━━━━━━━━┳━━━━━━━━━┳━━━┈
┃ accession  ┃ organism            ┃   len ┃ moltype ┃ biomol  ┃ …
┡━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━╇━━━━━━━━━╇━━━━━━━━━╇━━━┈
│ AY274119.3 │ SARS coronavirus…   │ 29751 │ rna     │ genomic │ …
│ AY613950.1 │ SARS coronavirus…   │ 29728 │ rna     │ genomic │ …
│ MN908947.3 │ SARS-CoV-2 Wuhan-…  │ 29903 │ rna     │ genomic │ …
│ KY352407.1 │ SARS-related coro…  │ 29274 │ rna     │ genomic │ …
└────────────┴─────────────────────┴───────┴─────────┴─────────┴───┈
```

(Columns `topology`, `completeness`, `sourcedb` and `updatedate` are shown too;
trimmed here for width. Add `--json` for the full, untruncated records.)

### Download FASTA to a file

```bash
viralfetch seq "Betacoronavirus pandemicum" --fasta -o out.fa
```

```
Wrote 4 fasta record(s) for Betacoronavirus pandemicum to out.fa
```

Without `-o`, records go to stdout (pipeable), and the summary goes to stderr.

### A whole taxon

`--taxon` operates on every species beneath a taxon. With `--meta` it shows a
**local aggregate** (no network) so you can decide whether a download is worth
it:

```bash
viralfetch seq --taxon Coronaviridae --meta
```

```
╭─ Coronaviridae (family) — download estimate ─╮
│ species     54                               │
│ isolates    59                               │
│ accessions  59                               │
│ RefSeq       0                               │
╰──────────────────────────────────────────────╯
   by genome composition
┏━━━━━━━━━━━━━┳━━━━━━━━━━━━┓
┃ composition ┃ accessions ┃
┡━━━━━━━━━━━━━╇━━━━━━━━━━━━┩
│ ssRNA(+)    │         59 │
└─────────────┴────────────┘
```

```json
{
  "name": "Coronaviridae",
  "rank": "family",
  "species": 54,
  "isolates": 59,
  "accessions": 59,
  "refseq": 0,
  "moltype_breakdown": { "ssRNA(+)": 59 }
}
```

Fetches of more than 500 accessions ask for confirmation unless `--yes` is
given (in `--json` / non-interactive mode they refuse without `--yes`).

### Molecule selection

`--moltype` and `--biomol` are `db=nuccore` fields, filtered locally over the
metadata (matching is lenient — `--moltype ssRNA` matches `ss-RNA`):

```bash
viralfetch seq --taxon Filoviridae --moltype ssRNA --fasta -o filo.fa
viralfetch seq "Betacoronavirus pandemicum" --biomol genomic
```

`--protein` is a **separate path**, not a nuccore filter: proteins live in
`db=protein`, reached via `elink` (nuccore → protein).

```bash
viralfetch seq "Betacoronavirus pandemicum" --protein --fasta
```

---

## Output & exit codes

- **Rich** tables/panels/trees by default; **`--json`** for clean, `jq`-ready
  output. Colour is disabled automatically when stdout is not a TTY.
- Errors go to stderr with a non-zero exit code: `1` not found, `2` bad usage,
  `3` missing NCBI email, `4` NCBI request failed.
- Partial failures are reported, never swallowed: if you request 200 accessions
  and 197 come back, the 3 missing ones are listed.

## Caching

Immutable data (sequences, accession metadata) is cached permanently; chapter
text will use a 30-day TTL. The cache lives in the platform cache directory.
Use `--no-cache` to bypass it for a single run.

## Development

```bash
pip install -e '.[dev]'
pytest            # offline tests only
pytest -m network # live NCBI integration tests (needs NCBI_EMAIL)
```

No test makes a network request by default; NCBI and ICTV responses are frozen
as fixtures under `tests/fixtures/`.
