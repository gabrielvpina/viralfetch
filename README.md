<div align="center">

<img src="viralfetch/assets/logo.png" alt="viralfetch logo" width="540">

<br>

![Version](https://img.shields.io/badge/version-0.2.0-blue)
![Python](https://img.shields.io/badge/python-3.10%2B-3776AB?logo=python&logoColor=white)
![License](https://img.shields.io/badge/license-MIT-green)
![CLI](https://img.shields.io/badge/CLI-Typer-0b7261)
![Output](https://img.shields.io/badge/output-Rich-009688)
![Tests](https://img.shields.io/badge/tests-passing-brightgreen)

**Viral taxonomy, sequences, trees and curated knowledge — from one command.**

</div>

viralfetch brings together the ICTV taxonomy, NCBI, the ICTV Report and
ViralZone in the terminal:

| Source | What it provides | Access |
|---|---|---|
| **ICTV VMR** | taxonomy, exemplar isolates, GenBank/RefSeq accessions | local (bundled) |
| **ICTV Report** | family chapters (text + figures), phylogenetic trees and alignments | chapters on demand · trees bundled |
| **NCBI** (E-utilities, Datasets) | taxonomy, sequences, metadata, genome counts | on demand, cached |
| **ViralZone** (SIB) | curated factsheets: host, geography, diseases, transmission, vaccines, antivirals | local (bundled) |

## Install

```bash
pip install viralfetch
viralfetch config --store-ncbi-email you@example.com   # required by NCBI for remote queries
```

Python 3.10+. An NCBI API key (`--store-ncbi-apikey KEY`) is optional and raises the
rate limit from 3 to 10 requests/s. Without an email, local commands still work.

## Quick start: profile any virus

Type a name — family, genus, species or a common virus name:

```bash
viralfetch Flaviviridae
viralfetch Orthoflavivirus denguei     # quotes are optional
viralfetch SARS-CoV-2                  # names only NCBI knows work too
```

```
╭──────────────────────────────────────────────────────────────────────────────────────────────────╮
│ Orthoflavivirus   genus · Flaviviridae                                                           │
│ ssRNA(+) · Group IV · NCBI taxid 3044782 · ✓ ICTV-accepted name            ✓ ICTV and NCBI agree │
│ Riboviria › Orthornavirae › Kitrinoviricota › Flasuviricetes › Amarillovirales › Flaviviridae ›  │
│ Orthoflavivirus                                                                                  │
╰──────────────────────────────────────────────────────────────────────────────────────────────────╯
╭───────────────────────────────────────┬──────────────────────────────────────────────────────────╮
│ ICTV (VMR)                            │ NCBI                                                     │
├───────────────────────────────────────┼──────────────────────────────────────────────────────────┤
│ 57 species                            │ 106,498 nucleotide records                               │
│ 82 isolates (57 exemplar)             │  10,792 complete genomes (by title) █░░░░░░░░░░░ 10%     │
│ ssRNA(+)                              │      85 RefSeq sequences                                 │
│ hosts: invertebrates, vertebrates     │     944 genome assemblies                                │
╰───────────────────────────────────────┴──────────────────────────────────────────────────────────╯
╭─ ViralZone · Orthoflavivirus (genus) ────────────────────────────────────────────────────────────╮
│ ● human  ● invertebrate                                                              factsheet ↗ │
│ diseases      YFV: Hemorrhagic fever, encephalitis                                               │
│               DENV: Dengue fever                                                                 │
│               ZIKV: Zika fever                                                                   │
│ transmission  Zoonosis, arthropod bite, mainly mosquitoes                                        │
│ …                                                                                                │
╰────────────────────────────────────────── ViralZone (SIB) · CC BY 4.0 · doi:10.1093/nar/gkad946 ─╯
```

The profile compares the ICTV and NCBI lineages, counts what NCBI holds and shows
the closest ViralZone factsheet. Useful variants:

```bash
viralfetch Orthoflavivirus --full                 # + lineage table, RefSeq list, virion/genome biology
viralfetch Coronaviridae --local                  # offline: VMR + ViralZone only
viralfetch --json Coronaviridae | jq .viralzone.factsheet.diseases
```

## Commands

### Taxonomy — `tax`, `members`

```bash
viralfetch tax Coronaviridae                          # ICTV lineage
viralfetch tax "Betacoronavirus pandemicum" --compare-ncbi   # ICTV vs NCBI, side by side
viralfetch tax SARS-CoV-2 --ncbi                      # NCBI taxonomy directly

viralfetch members Coronaviridae --rank genus         # genera, with species counts
viralfetch members Coronaviridae --tree               # full descendant hierarchy
viralfetch members Riboviria --rank family --count    # "158 families in realm Riboviria"
```

### Sequences — `seq`

Accessions come from the VMR; records come from NCBI.

```bash
viralfetch seq "Betacoronavirus pandemicum"                   # metadata table
viralfetch seq "Betacoronavirus pandemicum" --fasta -o sars.fa
viralfetch seq "Betacoronavirus pandemicum" --protein --fasta # proteins via elink
viralfetch seq --taxon Coronaviridae --meta                   # local download estimate
viralfetch seq --taxon Filoviridae --moltype ssRNA --gb -o filo.gb
```

Downloads above 500 accessions ask for confirmation (`--yes` to skip).

### ICTV Report — `text`

```bash
viralfetch text Coronaviridae                    # chapter rendered as Markdown
viralfetch text Coronaviridae --section summary
viralfetch text Coronaviridae --images           # figures drawn in the terminal
viralfetch text Geminiviridae --raw > gemini.md
viralfetch text Betacoronavirus                  # genera/species resolve to their family
```

### Trees and alignments — `tree`, `msa`

The ICTV Report trees and alignments are bundled, so these work offline.

```bash
viralfetch tree "Betacoronavirus pandemicum"     # family tree, the species highlighted
viralfetch tree Coronaviridae --tree 2           # another figure of the same family
viralfetch tree Geminiviridae --newick > gemini.nwk

viralfetch msa "Betacoronavirus pandemicum" --consensus
viralfetch msa Coronaviridae --range 100:180
viralfetch msa Coronaviridae --fasta > coronaviridae_msa.fasta
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
```

### Names the VMR doesn't know

Every command falls back to NCBI taxonomy for strains, synonyms and very recent taxa,
and says so on stderr:

```bash
viralfetch tree "Dengue virus 2"
# 'Dengue virus 2' is not in the local VMR; NCBI places it in family Flaviviridae —
# its subgenus 'Euflavivirus' is highlighted instead (closest relatives on the tree).
```

## Keeping data up to date

```bash
viralfetch update                  # install a newer ICTV VMR, if one was published
viralfetch update --trees          # rebuild trees/alignments from the ICTV Report (~15 min)
viralfetch update --viralzone      # rebuild the ViralZone factsheets (~20 min)
viralfetch update --reset          # back to the bundled VMR (--trees / --viralzone too)
```

Updates go to your user data directory and only replace the current data once they
complete successfully.

## Options, output and cache

| | |
|---|---|
| `--json` | clean JSON on stdout, for `jq` and scripts (goes **before** the name/command) |
| `--no-cache` | refetch instead of using the cache |
| `--email`, `--api-key` | NCBI credentials for this run only |
| exit codes | `1` not found · `2` bad usage · `3` no NCBI email · `4` remote request failed |
| cache | sequences and figures kept permanently; ICTV chapters 30 days; NCBI counts 7 days — `viralfetch cache info` / `cache clear` |

Run `viralfetch COMMAND --help` for every option.

## Credits

- **ICTV** — taxonomy (VMR), Report chapters, trees and alignments. Chapter content is
  licensed [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); `text` shows each
  chapter's citation.
- **NCBI** — E-utilities and Datasets APIs, used per the NCBI usage policies.
- **ViralZone**, SIB Swiss Institute of Bioinformatics — factsheet content (images
  CC BY 4.0; data for research and education, not for medical use). Please cite:

  > De Castro E, Hulo C, Masson P, Auchincloss A, Bridge A, Le Mercier P. ViralZone 2024
  > provides higher-resolution images and advanced virus-specific resources.
  > *Nucleic Acids Res.* 2024;52(D1):D817–D821.
  > doi:[10.1093/nar/gkad946](https://doi.org/10.1093/nar/gkad946)

## License

MIT
