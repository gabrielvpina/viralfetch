<div align="center">

<img src="viralfetch/assets/logo.png" alt="viralfetch logo" width="540">

<br>

![Version](https://img.shields.io/badge/version-0.1.3-blue)
![Python](https://img.shields.io/badge/python-3.10%2B-3776AB?logo=python&logoColor=white)
![License](https://img.shields.io/badge/license-MIT-green)
![CLI](https://img.shields.io/badge/CLI-Typer-0b7261)
![Output](https://img.shields.io/badge/output-Rich-009688)
![Tests](https://img.shields.io/badge/tests-passing-brightgreen)

</div>

<div align="center">
A command-line tool for querying and downloading viral taxonomy, metadata,
sequences, and per-family phylogenetic trees and alignments. It combines four
sources:
</div>

| Source | Content | Access |
|---|---|---|
| **VMR** (ICTV Virus Metadata Resource) | taxonomic hierarchy + exemplar isolates + GenBank/RefSeq accessions | **local**, embedded TSV |
| **NCBI E-utilities** | sequence metadata, sequences (nt/aa), NCBI taxonomy | **remote**, on demand |
| **ICTV Report** | descriptive chapter text per family (with figures) | **remote**, on demand |
| **ICTV Report trees** | per-family Newick trees + amino-acid/nucleotide alignments | **local**, bundled |

The VMR is the local index; remote sources are fetched on demand and cached.

---

## Installation

```bash
pip install viralfetch
```

This installs the `viralfetch` command and everything it needs in one step —
including Pillow (for drawing chapter figures) and `alv`/Biopython (for the
alignment viewer). Python 3.10+ is required.

## NCBI configuration

Commands that reach NCBI (`seq`, `tax --ncbi`, `tax --compare-ncbi`, `text`,
`update`) require a real email address, per NCBI usage policy. There is **no
default** — the command fails with an explanation if none is set.

```bash
export NCBI_EMAIL="you@example.com"
export NCBI_API_KEY="..."   # optional; raises the rate limit from 3 to 10 req/s
```

You can also pass `--email` / `--api-key` on any command.

**Session vs. persisted.** An `export` (or `--email`) applies only to the
current shell session. To store the email permanently, use:

```bash
viralfetch config --store-ncbi-email you@example.com
```

This writes to the config file (see `viralfetch config`), which survives across
sessions. Running `viralfetch config` warns you when no email is persisted yet.

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

### `tax --ncbi` (remote)

Look the lineage up **directly in NCBI's taxonomy** database instead of the
local VMR — useful for a name the VMR does not carry (a non-viral host, a very
recent taxon, an NCBI synonym). The name is resolved to a taxid via `esearch`,
then its lineage is fetched. Requires an NCBI email (see above).

```bash
viralfetch tax "SARS-CoV-2" --ncbi
```

```
╭─ SARS-CoV-2  (species) ────────────────────────╮
│ lineage                                        │
│ └── realm: Riboviria                           │
│     └── … → family: Coronaviridae              │
│         └── species: SARS-CoV-2                │
╰────────────────────────────────────────────────╯
NCBI taxonomy — taxid 2697049
```

With `--json` the payload is `{source: "ncbi", taxid, name, rank, lineage}`. An
unknown name exits `1`. `--ncbi` and `--compare-ncbi` are mutually exclusive.

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
content is **CC BY 4.0**. Figures are kept as image links; `--images` also
draws them in the terminal (see below).

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

With `--json`, the command emits `{slug, title, url, doi, images, markdown}`
instead — `images` is a list of `{url, alt}` for the chapter's figures, whose
absolute URLs also appear inline in the Markdown as `![alt](url)`.

### Figures in the terminal

`--images` draws the chapter's figures **in their original place** in the text,
as truecolor Unicode block-element graphics sized to your terminal (full width
by default, for the sharpest picture). Each character cell packs a 2×2 grid of
sub-pixels with a per-cell best-fit of two colours; the image is downscaled in
linear light (LANCZOS) and error-diffusion dithered, so gradients stay smooth.
It renders the same on any truecolor terminal — no special graphics protocol
required. Each figure keeps its `alt` as a caption. It applies only to the
interactive Rich view (not `--raw`, `--json`, or when piped):

```bash
viralfetch text Coronaviridae --images
viralfetch text Coronaviridae --images --fig-width 60   # cap the size
```

Figures are fetched politely (rate-limited, `robots.txt`-honouring, from the
ICTV domain only) and cached permanently. A missing or broken figure is skipped,
never fatal.

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

### Names the VMR doesn't know fall back to NCBI

If a name is absent from the local VMR, `text` asks **NCBI taxonomy** for its
family before giving up, and — if NCBI places it in one — shows that family's
chapter (with a note on stderr). This catches recent taxa and NCBI synonyms the
bundled VMR predates:

```bash
viralfetch text "some-recent-virus"
# stderr: 'some-recent-virus' is not in the local VMR; NCBI places it in
#         family Rhabdoviridae — showing that chapter.
# stdout: the Rhabdoviridae chapter
```

The fallback is best-effort: if NCBI is unreachable or knows no family, the name
is tried verbatim and, failing that, gets "did you mean" suggestions and exit
code `1`.

Fetching honours `ictv.global/robots.txt`, sends a descriptive `User-Agent`
carrying your contact email, waits at least 1 second between requests, and
caches chapter HTML for 30 days.

---

## `tree` — phylogenetic tree (local, NCBI only as a fallback)

Show the ICTV Report phylogenetic tree for a taxon's **family**, drawn as an
indented cladogram. Trees are bundled locally (Newick + metadata + a tip→virus
table per family), so this works offline.

The name is resolved through the VMR to its family, and the tip(s) it points at
are highlighted: a **species** highlights its own tip, a **genus** its whole
clade, a **family** shows the tree with nothing highlighted. A name the VMR does
not know is searched for among every tree's members (so a virus member name that
is not an ICTV taxon still finds its tree).

Failing that, the name is looked up in **NCBI taxonomy** — which knows strains,
synonyms and taxa newer than the bundled VMR. Its lineage names the family (the
same family names ICTV uses), and the deepest rank the trees do record —
species → subgenus → genus → subfamily — highlights the query's **closest
relatives** on that tree. The note on stderr always says which rank was used:

```bash
viralfetch tree "Dengue virus 2"
# 'Dengue virus 2' is not in the local VMR; NCBI places it in family
# Flaviviridae — no tip matches 'dengue virus type 2' itself, so its subgenus
# 'Euflavivirus' is highlighted instead (closest relatives on the tree).
```

This step is the only one that touches the network, and it is best-effort: with
no NCBI email configured (`viralfetch config --store-ncbi-email you@example.com`)
or no connection it is silently skipped, leaving the offline behaviour intact.
`msa` resolves names the same way.

```bash
viralfetch tree "Betacoronavirus pandemicum"
```

```
Coronaviridae · Figure 5A · RdRp (AA) · maximum likelihood · 55 tips · 2 matches
     ┌─ Alphacoronavirus BT020
   ┌─┤
   │ └─ Scotophilus bat coronavirus 512
  ─┤          …
   │     ┌ severe acute respiratory syndrome coronavirus  ← match
   └─────┤
         └ severe acute respiratory syndrome coronavirus 2  ← match

Caption: Figure 5 Coronaviridae. Phylogenetic relationships among members …
Other trees for Coronaviridae: --tree 2 (helicase)
```

When a family has several trees, pick one with `--tree N`; the query defaults to
whichever tree contains the match. Other options:

| Option | Effect |
|---|---|
| `--tree N` | Choose a tree when a family has several (1-based). |
| `--newick` | Emit the raw Newick string to stdout (for piping to other tools). |
| `--chapter` | Show the family's bundled ICTV Report chapter text instead. |

`--newick` and `--json` keep stdout clean (the redirect note goes to stderr).
With `--json`, the command emits `{family, source, note, matched_rank,
matched_value, tree: {…, matched, newick}, other_trees}` — `source` is
`vmr`, `member` or `ncbi`, and `matched_rank`/`matched_value` say what the
highlight fell back to. An unknown name gets "did you mean" suggestions and exit
code `1`; a family with no bundled tree exits `1` with a note.

---

## `msa` — multiple sequence alignment (local, no network)

Show the alignment behind a family's tree — the aligned FASTA that sits beside
each Newick — coloured by residue and wrapped into blocks with a column ruler
(via [`alv`](https://github.com/arvestad/alv)). The name is resolved to its
family exactly as `tree` does, and the query's own records are marked `▶`.

The alignments run to **thousands of columns**, so the view defaults to a
leading column window that fits your terminal; move or widen it with `--range`.

```bash
viralfetch msa "Betacoronavirus pandemicum" --consensus
```

```
Coronaviridae · tree1 · AA · cols 1–60 of 316 · 55 seqs · 2 matches
consensus              KHFFFAQDGDAAITDYDYYRYNRPTMLDICQALFVYEVVDKYFDIYEGGCITAKEVVVTN
▶ severe acute respir  KHFFFAQDGNAAISDYDYYRYNLPTMCDIRQLLFVVEVVDKYFDCYDGGCINANQVIVNN
Alphapironavirus bona  EHFFYLQPRDCAVTDFDYYRFNRPTVLDPLQFRFVYNVVKHYFKSYSAGCLKSEFVIINN
…                      0^                 20^                 40^
```

| Option | Effect |
|---|---|
| `--tree N` | Choose a tree when a family has several (1-based). |
| `--range A:B` | Column window, 1-based inclusive (`100:180`; either bound may be omitted). |
| `--consensus` | Prepend a per-column majority-residue row. |
| `--fasta` | Emit the (windowed) alignment as FASTA to stdout. |

With `--json`, the command emits `{family, tree_id, molecule, total_cols, start,
n_cols, matched, consensus, rows: [{name, seq, matched}]}`. A family whose tree
has no bundled alignment exits `1` with a note.

---

## Utilities

Small helper commands, all local except `update`.

```bash
viralfetch diagnose                 # VMR parser quality (zero-accession rows)
viralfetch update                   # is a newer VMR published on ictv.global?
viralfetch config                   # show email, masked API key, cache paths
                                    # (warns if no NCBI email is persisted yet)
viralfetch config --store-ncbi-email you@example.com    # persist email
viralfetch config --store-ncbi-apikey KEY               # persist API key
viralfetch cache info               # per-namespace entry counts and size
viralfetch cache clear --texts      # drop cached ICTV chapters (or --seqs, --images, or all)
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

Immutable data (sequences, accession metadata, and chapter figures) is cached
permanently; ICTV chapter HTML uses a 30-day TTL. The cache lives in the
platform cache directory. Use `--no-cache` to bypass it for a single run.


