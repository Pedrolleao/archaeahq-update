# archaeahq-update

Keep your own copy of **ArchaeaHQ** up to date with NCBI.

[ArchaeaHQ](https://vee-lab.eu/archaeahq) is a quality-controlled, systematically curated reference
database of archaeal genomes: 21,644 genomes from the four archaeal kingdoms
(Methanobacteriati/Euryarchaeota, Thermoproteati/TACK, Nanobdellati/DPANN,
Promethearchaeati/Asgard), all with completeness ≥ 70 % and contamination ≤ 10 %, de-replicated at
species level and annotated with environment, geography, quality and RNA-gene metadata.

**Download the database (v1.0, released 2026):**
https://figshare.com/articles/dataset/ArchaeaHQ_v_1_00_2026/32266599
(DOI 10.6084/m9.figshare.32266599, CC BY 4.0) — `fna.zip` (genomes, 11 GB), `faa.zip` (proteins),
`Archaea_HQ-16S.fasta` and the supplementary tables.

New archaeal assemblies are deposited at NCBI every week. This tool lets anyone **update the
database themselves**, with the same rules that produced v1.0: it finds the assemblies deposited
after your current version, downloads them, runs the curation pipeline (CheckM2 quality gate →
skani redundancy → environment classification → rRNA/tRNA counts), writes a recommendation table,
and — if you accept — builds the next database version in a new folder.

Created at Vee Lab - Radboud University.

---

## Contents

1. [Requirements and installation](#1-requirements-and-installation)
2. [The update cycle in four commands](#2-the-update-cycle-in-four-commands)
3. [How a genome is evaluated](#3-how-a-genome-is-evaluated)
4. [Outputs of a run](#4-outputs-of-a-run)
5. [Database versions and the `release` command](#5-database-versions-and-the-release-command)
6. [All commands and flags](#6-all-commands-and-flags)
7. [Folder layout](#7-folder-layout)
8. [Notes, limitations, citation](#8-notes-limitations-citation)

---

## 1. Requirements and installation

* Linux or macOS; Python ≥ 3.10 to launch the script (the heavy tools run inside a conda environment).
* **conda** or **mamba** (Miniforge/Miniconda). `setup` creates the `archaeahq-update` environment
  from `environment.yml` with CheckM2, NCBI Datasets CLI, skani and barrnap — only the tools that
  are missing from your `PATH` and existing environments are installed.
* Disk: ~3 GB for the CheckM2 database (downloaded once), ~35 GB for the unpacked v1.0 genomes,
  plus room for the new genomes of each run (a few GB per thousand genomes).
* Internet access to NCBI. An NCBI API key (`--api-key` or `$NCBI_API_KEY`) is optional and doubles
  the rate limit.

**Install** (pick one):

```bash
# a) pip, into any Python ≥ 3.10 (the external tools are installed by `setup` afterwards)
pip install git+https://github.com/Pedrolleao/archaeahq-update.git

# b) or just clone: the launcher archaeahq_update.py runs from the clone without installing
git clone https://github.com/Pedrolleao/archaeahq-update.git
cd archaeahq-update
```

Installed, the command is `archaeahq-update`; from a clone, `python3 archaeahq_update.py`
(same commands and flags — the examples below use the installed name).

```bash
archaeahq-update setup                # conda env (if needed) + CheckM2 database
archaeahq-update fetch-db             # ArchaeaHQ v1.0 from figshare → Archaea_HQ-v1.0/
```

`fetch-db` is the guided first-time setup: it shows what will be downloaded and the disk space
needed, asks for confirmation, downloads `fna.zip` (11 GB, resumable — just rerun it after an
interruption), verifies the MD5 published by figshare, unpacks the 21,644 genomes into
`Archaea_HQ-v1.0/fna/`, registers the version and prints the next commands. Add `--extras faa,16s,tables`
for the proteins, the 16S sequences and the supplementary tables, or `--from-zip fna.zip` if you
already downloaded the archive by hand.

The v1.0 information table (`data/ArchaeaHQ-Info.tsv`) and the list of the 35,993 assemblies that
were evaluated for v1.0 (`data/evaluated_accessions.tsv`) ship with the repository, so the tool
knows exactly what v1.0 contains and what it already rejected.

---

## 2. The update cycle in four commands

```bash
# 1. What is new at NCBI since my version?   (nothing is downloaded)
archaeahq-update check

# 2. Download and evaluate the new genomes
archaeahq-update run --threads 24

# 3. Look at the recommendations
less archaeahq_update_work/runs/run_<date>/summary.txt

# 4. Accept them: build the next version, Archaea_HQ-v1.1/, next to the script
archaeahq-update release
```

Every command finds the newest `Archaea_HQ-v<major>.<minor>/` folder automatically (the one
`fetch-db` created, then the ones `release` builds), so the cycle is simply `check` → `run` →
`release`, and each `release` produces v1.1, v1.2, … The version folders live in the current
folder when the tool is installed, or next to `archaeahq_update.py` when run from a clone; set
`$ARCHAEAHQ_RELEASES` or pass `--releases-dir` to keep them elsewhere. Run the commands from the
same folder each time so the work directory (`./archaeahq_update_work`) is found again.

If you unpacked `fna.zip` yourself instead of using `fetch-db`, point the first `run` and
`release` at that folder with `--db-fna /path/to/fna`; `release` then creates `Archaea_HQ-v1.1/`
from it.

Useful variations:

```bash
archaeahq-update run --kingdoms asgard,dpann --max-genomes 50   # small test run
archaeahq-update run --accessions-file my_accessions.txt         # evaluate your own list
archaeahq-update release --dry-run                               # preview only
archaeahq-update release --accept reviewed_accessions.txt        # apply a curated subset
```

An interrupted `run` resumes from the last finished stage when you repeat the same command.

---

## 3. How a genome is evaluated

The rules are those of the ArchaeaHQ v1.0 build; every cut-off is in `data/thresholds.json`.

| # | Stage | Tool | Rule |
|---|---|---|---|
| 1 | **List** | `datasets summary genome taxon` for the four kingdoms (taxids 3366610, 1783275, 1783276, 1935183) | *New* = accession not in the current table, not among the 35,993 assemblies evaluated for v1.0 and not in the ledger of your previous runs (`--reevaluate` ignores the ledger). GCA/GCF twins of database genomes are reported but not downloaded unless `--include-twins`. A **newer version** of a database genome (same record, higher version number) is always evaluated; if it passes the quality gate it is recommended as `Replace`. |
| 2 | **Download** | `datasets download genome accession` | batched, retried; failures are reported, never silently dropped. |
| 3 | **Metadata** | bundled `pyassembly` + `pytaxonkit` (NCBI Datasets API, NCBI Taxonomy) | BioProject, BioSample, isolation source, ENVO fields, coordinates, metagenome source, full lineage. |
| 4 | **Quality** | CheckM2 `predict` | keep if **completeness ≥ 70 % and contamination ≤ 10 %**. `Park Quality_Score = completeness − 5·contamination + 0.5·log10(contig N50)`. |
| 5 | **Redundancy** | skani `search` (new vs database) and `triangle` (new vs new) | **identical**: ANI ≥ 99.9 % and aligned fraction ≥ 99 %, or same numeric accession GCA_x/GCF_x; **same species**: ANI ≥ 95 % and aligned fraction ≥ 95 %. Genomes identical to or of the same species as a database genome are not added (a note flags when the new genome scores higher than the database representative). Within the batch one representative per species cluster is kept (highest Park score; ties → newest RefSeq). |
| 6 | **RNA** | barrnap ≥ 1.10 (`--kingdom arc`; Infernal/Rfam for rRNA, ARAGORN for tRNA) | `rRNA` / `tRNA` columns. |
| 7 | **Environment** | bundled `classify_environments.py` | the 17 ecological categories of the paper + Unknown, same rules and priority. |
| 8 | **Report** | — | `Add` only when every filter passes; `Replace` for a newer version of a database genome that passes quality; otherwise `Do_not_add` with the reason. |

Missing metadata follows the database convention: `0` for empty text fields, `N/A` for
`Metagenome_source`, `NA` for values that could not be computed.

---

## 4. Outputs of a run

Everything goes to `archaeahq_update_work/runs/run_<date>/`:

| File | Content |
|---|---|
| `ArchaeaHQ_update_<date>_full_table.tsv` | One row per new genome: `Recommendation` (`Add` / `Replace` / `Do_not_add`), `Reason`, then **the 25 columns of the ArchaeaHQ information table** (Name, Assembly ID, Archaeal_Kingdom, Taxonomy, Completeness, Contamination, Park Quality_Score, Bioproject, Biosample, Environment_Category-ArchaeaHQ, Metagenome_source, geo_loc_name, lat_lon, Isolation_Source, Env_Broad, Env_Local, Env_medium, Coding_Density, Contig_N50, Average_Gene_Length, Genome_Size, GC_Content, Total_Coding_Sequences, rRNA, tRNA). |
| `ArchaeaHQ_update_<date>_filters.tsv` | Assembly ID, kingdom, `Overall` (`PASS` / `REPLACE` / `FAIL`), one `PASS`/`FAIL` column per filter (`Downloaded`, `Completeness_ge70`, `Contamination_le10`, `Not_identical_to_DB`, `Not_same_species_as_DB`, `Representative_within_batch`) and `Notes`. |
| `summary.txt` | Per-kingdom counts along the filter chain and the reasons for rejection. |
| `sankey_input.tsv` | Input for the Sankey figure (`--sankey` renders it, like Fig. 1a of the paper). |
| `0N_<stage>/` | Raw outputs of every stage: NCBI listings, downloaded FASTA, CheckM2 report, skani tables, GFFs. |

The work directory also keeps `ledger.tsv` (every accession you ever evaluated, with the decision
and the run) and `cache/` (metadata cache, taxonomy dump, CheckM2 database, database sketch).

---

## 5. Database versions and the `release` command

A database version is a folder `Archaea_HQ-v<major>.<minor>/` containing:

| File | Content |
|---|---|
| `fna/` | one FASTA per genome |
| `ArchaeaHQ-Info.tsv` | the information table of that version (25 columns) |
| `release_changes.tsv` | every applied or skipped recommendation: `added`, `replaced`, `removed`, `skipped`, with the reason |
| `release.json` | version label, date, previous version, source run, counts |
| `comparison.txt` | the old-vs-new comparison printed at the end of `release` |

`release` takes the newest run (or `--run NAME` / `--table F`), applies its `Add` and `Replace`
rows to the current version and writes the next one (minor version + 1) into `--releases-dir`
(default: the folder of the script):

* `Add` rows are appended to the table, grouped by kingdom, and their FASTA files brought in.
* For a `Replace` row the old assembly version is removed from the table and the FASTA folder and
  the new version takes its place.
* Rows already in the database or without a FASTA file are skipped and listed. If nothing applies,
  no folder is written.
* FASTA files are hard-linked when the new folder is on the same file system as the current one
  (instant, no extra space), otherwise copied (`--copy` forces copies).
* Once the new folder is verified (every carried genome present with the same name and size),
  the **superseded version is deleted** — `--keep-previous` keeps it. A folder holding FASTA files
  that are not in the table is never deleted automatically, and the repository's `data/` folder is
  never touched. On the first release the `Archaea_HQ-v1.0/` folder made by `fetch-db` (or the
  unpacked folder given with `--db-fna`) is the one removed; the v1.0 table stays in `data/`.
* The comparison between the two versions (genomes per kingdom, quality and size statistics,
  environment categories) is printed and saved as `comparison.txt`.

```
                    Genomes per kingdom
Kingdom                            v1.0   added   replaced     v1.1   growth
Methanobacteriati-Euryarchaeota    9,894   1,730          1   11,624   +17.5%
Thermoproteati-TACK                6,379   1,567          0    7,946   +24.6%
Nanobdellati-DPANN                 4,891     466          0    5,357    +9.5%
Promethearchaeati-Asgard             480     177          0      657   +36.9%
Total                             21,644   3,940          1   25,584   +18.2%
```

---

## 6. All commands and flags

```
archaeahq_update.py setup      [common options]
archaeahq_update.py fetch-db   [--from-zip F] [--extras faa,16s,tables] [--keep-zip] [--no-verify] [--yes]
archaeahq_update.py check      [--reevaluate] [--accessions-file F]
archaeahq_update.py run        [--db-fna DIR | --db-sketch DIR] [--max-genomes N] [--accessions-file F]
                               [--reevaluate] [--include-twins] [--from-stage STAGE] [--run-name NAME]
                               [--no-rna] [--sankey] [--lowmem]
archaeahq_update.py sketch-db  --out DIR [--db-fna DIR]
archaeahq_update.py release    [--run NAME | --table F] [--new-fna DIR] [--accept F] [--version vX.Y]
                               [--keep-previous] [--copy] [--dry-run] [--db-fna DIR] [--out DIR]
```

**Common options** (every command)

| Flag | Default | Meaning |
|---|---|---|
| `--workdir DIR` | `./archaeahq_update_work` | runs, caches, ledger and logs |
| `--releases-dir DIR` | `$ARCHAEAHQ_RELEASES`, else the current folder (clone launcher: its own folder) | where the `Archaea_HQ-v*` version folders live; the newest one is the current database |
| `--db-table F` | table of the newest version, else `data/ArchaeaHQ-Info.tsv` (v1.0) | use another information table |
| `--kingdoms LIST` | `eury,tack,dpann,asgard` | subset of kingdoms to list/evaluate |
| `--threads N` | CPUs − 2 | threads for CheckM2, skani, barrnap |
| `--api-key KEY` | `$NCBI_API_KEY` | NCBI API key |
| `--conda-env NAME` | `archaeahq-update` | environment to create when tools are missing |
| `--checkm2-db PATH` | downloaded to the cache | folder or `.dmnd` file of the CheckM2 database |
| `--no-install` | | never create a conda environment; fail instead |
| `--plain` | (also `NO_COLOR=1`) | no colours, no progress bars |

**`fetch-db`** — download ArchaeaHQ v1.0 from figshare and set it up as `Archaea_HQ-v1.0/` in
`--releases-dir`.

| Flag | Meaning |
|---|---|
| `--from-zip F` | use an already downloaded `fna.zip` instead of downloading it |
| `--extras LIST` | also fetch `faa` (proteins, 6.3 GB zip → `faa/`), `16s` (`Archaea_HQ-16S.fasta`), `tables` (supplementary xlsx) |
| `--keep-zip` | keep the zip archives after unpacking (default: deleted to save space) |
| `--no-verify` | skip the MD5 check against figshare's checksums |
| `--yes`, `-y` | do not ask for confirmation |

Downloads are resumable: an interrupted `fetch-db` continues where it stopped when rerun.

**`check`** — list the new assemblies per kingdom; nothing is downloaded.

| Flag | Meaning |
|---|---|
| `--reevaluate` | ignore the ledger of previously evaluated genomes |
| `--accessions-file F` | list this file's accessions instead of the NCBI diff |

**`run`** — the full evaluation (eight stages).

| Flag | Meaning |
|---|---|
| `--db-fna DIR` | FASTA folder of the current version (default: `fna/` of the newest release). The first run sketches it with skani and caches the sketch. |
| `--db-sketch DIR` | use a sketch made with `sketch-db` instead of the FASTA files (e.g. on a machine without the 35 GB) |
| `--max-genomes N` | cap the number of new genomes per kingdom, newest first |
| `--accessions-file F` | evaluate this accession list instead of the NCBI diff (one per line, `#` comments) |
| `--reevaluate` | ignore the ledger |
| `--include-twins` | also download and process GCA/GCF twins of database genomes (default: reported, not downloaded) |
| `--from-stage STAGE` | recompute from `env, list, download, metadata, quality, redundancy, rna` or `report` on |
| `--run-name NAME` | run folder name (default `run_<date>`) |
| `--no-rna` | skip rRNA/tRNA counting (columns become `NA`) |
| `--sankey` | render the Sankey figure (needs plotly + kaleido) |
| `--lowmem` | pass `--lowmem` to CheckM2 |

**`sketch-db`** — build a reusable skani sketch of the database genomes (`--out DIR`; `--db-fna`
defaults to the newest version).

**`release`** — build the next database version.

| Flag | Meaning |
|---|---|
| `--run NAME` | apply this run (default: the newest run in the work directory) |
| `--table F` | apply any full result table instead |
| `--new-fna DIR` | FASTA of the new genomes (default: `<run>/02_download/fna`) |
| `--accept F` | apply only the accessions listed in this file |
| `--version vX.Y` | label of the new version (default: current + 0.1; must be newer) |
| `--out DIR` | destination (default: `<releases-dir>/Archaea_HQ-v<version>`) |
| `--keep-previous` | do not delete the superseded version |
| `--copy` | copy FASTA files instead of hard-linking |
| `--dry-run` | show the selection and the comparison; write and delete nothing |
| `--db-fna DIR` | FASTA folder of the current version (needed on the first release from v1.0) |

Every external command and its output is appended to `runs/<run>/commands.log`; the tool's own
log is `runs/<run>/archaeahq_update.log`.

---

## 7. Folder layout

```
archaeahq-update/
├── archaeahq_update.py          launcher for running from a clone (no install needed)
├── pyproject.toml               pip-installable package; entry point `archaeahq-update`
├── environment.yml              → src/archaeahq_update/data/environment.yml (conda env with the tools)
├── README.md  DEVELOPMENT.md  LICENSE
├── src/archaeahq_update/
│   ├── cli.py                   commands, stages, report assembly, release
│   ├── common.py  ui.py  envcheck.py  ncbi.py  quality.py  redundancy.py  rna.py  report.py  fetchdb.py
│   ├── classify_environments.py environment classifier (as used for the paper)
│   ├── compile_barrnap.py       GFF → RNA counts (as used for the paper)
│   ├── sankey_generic.py        Sankey figure
│   ├── pyassembly/ pytaxonkit/ pypipeline/   NCBI metadata + taxonomy toolkit
│   └── data/
│       ├── ArchaeaHQ-Info.tsv           v1.0 information table (21,644 genomes, 25 columns)
│       ├── evaluated_accessions.tsv     the 35,993 assemblies evaluated for v1.0 and their decision
│       ├── kingdoms.json                taxids, labels, colours of the four kingdoms
│       ├── thresholds.json              every cut-off
│       └── environment.yml              conda environment (checkm2, ncbi-datasets-cli, skani, barrnap …)
├── tests/                       pytest suite (`pip install -e .[test] && pytest`)
├── Archaea_HQ-v1.1/             (not in git) your current database version (`fetch-db` → v1.0, `release` → v1.1, …)
└── archaeahq_update_work/       (not in git) runs, caches, ledger
```

---

## 8. Notes, limitations, citation

* NCBI archaea outside the four kingdoms (e.g. unclassified Archaea) are not listed, exactly as in
  the original build. With `--accessions-file` they are processed and labelled `Other archaea`.
* CheckM2 is the slow step: roughly 1–3 minutes per 100 genomes on 24 threads once the models are
  loaded. `--lowmem` reduces memory at some speed cost.
* If only barrnap 0.9 is found, rRNA still comes from barrnap and tRNA from a direct ARAGORN call.
* The pipeline recommends; whether to add a genome, and whether to replace a database
  representative by a better new genome of the same species, remain curator decisions. `release`
  applies only what the table (or your `--accept` list) says.

**Citation.** ArchaeaHQ — Leão P. et al., *A quality-controlled, systematically curated reference
database of archaeal genomes* (in preparation); dataset DOI 10.6084/m9.figshare.32266599.
Tools: CheckM2 (Chklovski et al. 2023), skani (Shaw & Yu 2023), barrnap (Seemann) and ARAGORN
(Laslett & Canback 2004), NCBI Datasets.

Created at Vee Lab - Radboud University.
