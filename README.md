<div align="center">

<img src="viralfetch/assets/logo.png" alt="viralfetch logo" width="540">

<br>

![Version](https://img.shields.io/badge/version-0.1.0-blue)
![Python](https://img.shields.io/badge/python-3.10%2B-3776AB?logo=python&logoColor=white)
![License](https://img.shields.io/badge/license-MIT-green)
![CLI](https://img.shields.io/badge/CLI-Typer-0b7261)
![Output](https://img.shields.io/badge/output-Rich-009688)
![Tests](https://img.shields.io/badge/tests-passing-brightgreen)

</div>

<div align="center">
A command-line tool for querying and downloading viral taxonomy, metadata and
sequences. It combines three sources:
</div>

| Source | Content | Access |
|---|---|---|
| **VMR** (ICTV Virus Metadata Resource) | taxonomic hierarchy + exemplar isolates + GenBank/RefSeq accessions | **local**, embedded TSV |
| **NCBI E-utilities** | sequence metadata, sequences (nt/aa), NCBI taxonomy | **remote**, on demand |
| **ICTV Report** | descriptive chapter text per family | **remote**, on demand |

The VMR is the local index; everything else is fetched on demand and cached.

---

## Installation

```bash
pip install viralfetch
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

## `text` — ICTV Report chapter (remote)

Fetch the ICTV Report chapter for a family, convert its main content to
Markdown, and render it with headings, subsection titles, characteristic
tables, and the **italics of scientific names** preserved. The original page
URL and the chapter's references/attribution are shown at the top, and the
content is **CC BY 4.0**. Images are omitted.

```bash
viralfetch text Coronaviridae
```

```
                          Family: Coronaviridae

Source: https://ictv.global/report/chapter/coronaviridae/coronaviridae

Patrick C.Y. Woo, Raoul J. de Groot, Bart Haagmans, Susanna K.P. Lau, …

The citation for this ICTV Report chapter is the summary published as:
Woo et al., (2023), ICTV Virus Taxonomy Profile: Coronaviridae 2023,
Journal of General Virology (2023) 104, 001843

Content is licensed CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/).
────────────────────────────────────────────────────────────────────────
Summary

Members of the family Coronaviridae, a monophyletic group of viruses in
the order Nidovirales, are enveloped, positive-sense RNA viruses …
```

### A single section

`--section` restricts the output to one section, matched by heading
(case-insensitive substring). The top attribution block is always kept.

```bash
viralfetch text Coronaviridae --section summary
```

### Raw Markdown to a file

`--raw` emits the pure Markdown (no Rich decoration), ready to redirect:

```bash
viralfetch text Geminiviridae --raw > geminiviridae.md
```

```markdown
# Family: Geminiviridae

*Source: https://ictv.global/report/chapter/geminiviridae/geminiviridae*

**Elvira Fiallo-Olivé, Jean-Michel Lett, Darren P. Martin, …**

The citation for this ICTV Report chapter is the summary published as
Fiallo-Olivé et al., ICTV Virus Taxonomy Profile: *Geminiviridae* 2021, …
```

With `--json`, the command emits `{slug, title, url, doi, markdown}` instead.

### Genera and species resolve to their family chapter

The ICTV Report is organised **by family** — genera and species have no chapter
of their own. Given one, `text` looks it up in the VMR and shows its family's
chapter, with a note on stderr (so `--raw`/`--json` stdout stays clean):

```bash
viralfetch text Betacoronavirus
# stderr: 'Betacoronavirus' is a genus; the ICTV Report has no genus chapter
#         — showing its family, Coronaviridae.
# stdout: the Coronaviridae chapter
```

An unknown name gets "did you mean" suggestions and exit code `1`.

Fetching honours `ictv.global/robots.txt`, sends a descriptive `User-Agent`
carrying your contact email, waits at least 1 second between requests, and
caches chapter HTML for 30 days.

---

## Utilities

Small helper commands, all local except `update`.

```bash
viralfetch diagnose                 # VMR parser quality (zero-accession rows)
viralfetch update                   # is a newer VMR published on ictv.global?
viralfetch config                   # show email, masked API key, cache paths
viralfetch config --store-ncbi-email you@example.com   # persist credentials
viralfetch cache info               # per-namespace entry counts and size
viralfetch cache clear --texts      # drop cached ICTV chapters (or --seqs, or all)
```

`diagnose` reports the parser's quality indicator — how many VMR rows yielded
zero accessions:

```
╭─ VMR accession-parser diagnostics ─╮
│ isolates              19271        │
│ accessions            23249        │
│ empty-accession rows  141          │
│ unparsed rows         0            │
╰────────────────────────────────────╯
```

**Shell completion:** taxon names complete from the VMR
(`viralfetch tax Corona<TAB>` → `Coronaviridae`). Install it once with
`viralfetch --install-completion`.

---

## Output & exit codes

- **Rich** tables/panels/trees by default; **`--json`** for clean, `jq`-ready
  output. Colour is disabled automatically when stdout is not a TTY.
- Errors go to stderr with a non-zero exit code: `1` not found, `2` bad usage,
  `3` missing NCBI email, `4` NCBI request failed.
- Partial failures are reported, never swallowed: if you request 200 accessions
  and 197 come back, the 3 missing ones are listed.

## Caching

Immutable data (sequences, accession metadata) is cached permanently; ICTV
chapter HTML uses a 30-day TTL. The cache lives in the platform cache
directory. Use `--no-cache` to bypass it for a single run.


