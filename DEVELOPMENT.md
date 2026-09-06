# archaeahq-update — how the tool was created and how it works

*Written 2026-09-05, after the first complete build and test of the tool.*

This document is the developer-side record of `archaeahq-update`. The user-facing manual is
`README.md`. Here we explain **where the tool came from** (the reconstruction of the original
ArchaeaHQ build), **which design decisions were taken and why**, **what every file in the bundle
is**, **how the code is organised**, **what was tested**, and **which pitfalls were met** so that
the next person maintaining it does not rediscover them.

---

## 1. Starting point: reconstructing how ArchaeaHQ v1.0 was built

ArchaeaHQ v1.0 (21,644 genomes) was assembled in early 2026 from the February-2025 NCBI snapshot
of the four archaeal kingdoms. The build was not a single pipeline but a chain of separate
steps, partly on the local WSL machine (`~/Database/NCBI_Archaea_02-2025/`, toolkit in
`~/Scripts/`) and partly on the Radboud cluster (CheckM2, dRep). Before writing the updater the
chain was traced file by file:

| Step | Evidence found | What it did |
|---|---|---|
| NCBI listing + download | `Methanobacteriati-3366610_Euryarchaeota.tsv` etc. (dataformat-style columns), `Assembly_ID-*.txt` | NCBI Datasets CLI per kingdom taxid: 17,898 Eury + 10,250 TACK + 6,206 DPANN + 1,639 Asgard = 35,993 |
| Metadata + taxonomy | `*-Results_Tax_Env-Info/{assembly_metadata,taxonomy_lineage,combined}.tsv`; `~/Scripts/pipeline/pipeline.py` | `assembly_meta` (NCBI Datasets v2 REST API, SQLite cache) + `pytaxonkit` (NCBI taxdump) → 36-column `combined.tsv` |
| Quality | `CheckM2-*-quality_report.tsv` | CheckM2 completeness/contamination and assembly statistics |
| Info tables | `Polishing/Info-*.txt`, `Info-*.xlsx` | join of `combined.tsv` + CheckM2; `Quality (MIMAG)`; `High-Medium Quality (DB) = Selected` when completeness ≥ 70 and contamination ≤ 10 |
| Dereplication | `Polishing/{Ndb,Wdb}.csv`, `analyze_drep_clusters.py`, `dRep_analysis_README.md`; Asgard in `Asgard-Results/data_tables/` (fastANI) | dRep (ANI 0.95) → Union-Find post-processing: identity ANI ≥ 0.999 & cov ≥ 0.99 + GCA/GCF same-accession rule; species ANI ≥ 0.95 & cov ≥ 0.95 |
| Final selection | `00.Info_unify.xlsx` (`All_genomes` = `Initial_Dataset.txt`, `ArchaeHQ` = `ArchaeaHQ.txt`) | **in DB ⇔ `Selected` AND `Representative_Genome == Yes`** (verified: reproduces all 21,644 rows exactly) |
| Environment | `Env-info/classify_environments.py` | 17 rule-based categories + Unknown from five metadata fields |
| RNA counts | `Barrnap/Script/{run_barrnap_pipeline,compile_barrnap}.py`, GFF sources `infernal` + `aragorn` | tseemann barrnap dev 1.7 (`--kingdom arc`), tRNA via ARAGORN |
| Quality score | `Park Quality_Score` column | `completeness − 5·contamination + 0.5·log10(contig N50)` (checked on two rows to 8 decimals) |
| Sequence files | `2Recover*.txt`, `ArchaeaHQ-<kingdom>/{fna,faa,ORF,gbk}` | Prodigal 2.6.3 meta mode on the selected `.fna` |

Two facts from this reconstruction shaped the tool: the decision logic lived in spreadsheets, and
CheckM2/dRep were not installed locally (only CheckM v1.2.2 was). The Supplementary Table 1 layout
(`Defense_finder/ArchaeaHQ-Info.txt`, 25 columns) was taken as *the* information table of the database.

---

## 2. Requirements and decisions

The request: one script that (1) checks which genomes are in the database, (2) lists what NCBI now
has per kingdom, (3) downloads the new ones, (4) runs the complete pipeline, and (5) outputs a full
information table (exact database columns) with the recommendation in column 1, plus a simple
PASS/FAIL table per filter stage. Follow-up requirements: the folder must be usable by people who
did not build ArchaeaHQ (everything bundled), and the CLI must be informative and pleasant, in the
colour style of vee-lab.eu/archaeahq.

Decisions taken with the user:

| Question | Decision |
|---|---|
| CheckM2 not installed | check PATH and every conda env; if missing, create a conda env from a bundled `environment.yml` and download the CheckM2 database |
| Redundancy method | **skani** against the database (plus within the batch) instead of re-running dRep — same thresholds, minutes instead of hours, no MUMmer/Mash |
| Table layout | Supplementary Table 1 (25 columns, incl. `Environment_Category-ArchaeaHQ`, `Park Quality_Score`, `rRNA`, `tRNA`) |
| After the recommendation | recommendation only; the database is never modified |

Decisions taken by the developer (documented in README):

* GCA/GCF twins of database genomes are reported but **not downloaded** by default
  (`--include-twins` to process them): they are identical by definition; their quality values
  are copied from the database twin and `Downloaded` is `SKIPPED`.
* **Newer versions** of database genomes (same prefix and numeric accession, higher version —
  found when the first `check` showed 21,643 instead of 21,644 database genomes at NCBI:
  `GCA_049860595.1` had been superseded by `.2`, 1,620,528 → 1,613,967 bp) are a separate case
  (`twin_type = newer_version`, `twin_kind()` in `archaeahq_update.py`): they are downloaded and
  fully evaluated; if they pass the quality gate the recommendation is **`Replace`**
  (`Overall = REPLACE` in the filters table, reason lists the completeness/contamination/size change),
  otherwise `Do_not_add` with the note that the database keeps the previous version.
* "New" excludes the 35,993 accessions already evaluated for v1.0 (so rejected genomes are not
  re-evaluated forever) and everything in the tool's own ledger; `--reevaluate` overrides.
* Every downloaded genome is compared, so every filter column is PASS/FAIL (no blanks), and
  within-batch clusters are built over all genomes with the representative chosen among the
  eligible ones (quality pass, not redundant with the DB; highest Park score; tie → newest RefSeq).
* Missing metadata follows the database convention (`0`, `N/A` for `Metagenome_source`); values
  that could not be computed are `NA`.
* Kingdom is taken from the NCBI lineage taxids; fallback to the listing kingdom; last resort a
  lookup of the accession in any kingdom listing cached in the work dir (needed when NCBI still
  reports a taxid that has been deleted from the taxdump).

---

## 3. What is in the bundle and where it came from

```
archaeahq-update/
├── archaeahq_update.py      entry point (CLI, run context, stage orchestration, table assembly)
├── environment.yml          conda env: python≥3.12, checkm2≥1.1, ncbi-datasets-cli≥18, skani≥0.3,
│                            barrnap≥1.10, aragorn, rich, requests, tqdm, plotly, pillow, kaleido 0.2.1
├── README.md                user manual
├── DEVELOPMENT.md           this file
├── data/
│   ├── ArchaeaHQ-Info.tsv           copy of Defense_finder/ArchaeaHQ-Info.txt (CRLF → LF), 21,644 rows
│   ├── evaluated_accessions.tsv     derived from Paper/Initial_Dataset.txt + ArchaeaHQ.txt:
│   │                                 21,644 Add · 8,948 rejected for quality · 5,401 rejected as redundant
│   ├── kingdoms.json                taxid, label, colour of the four kingdoms
│   └── thresholds.json              every cut-off (quality, MIMAG, identity, species, download)
└── lib/
    ├── common.py       constants (DB_COLUMNS…), TSV/JSON helpers, accession helpers, Park score,
    │                   subprocess wrappers (run_cmd, run_pipe) with logging, stage markers
    ├── ui.py           rich theme in the site palette; banner, stage headers, tables, progress,
    │                   spinner, verdict/error panels; --plain fallback
    ├── envcheck.py     tool discovery (PATH → conda envs → create env), CheckM2 DB download,
    │                   barrnap capability detection
    ├── ncbi.py         datasets summary/dataformat listing, accession summaries, batched download
    ├── quality.py      CheckM2 wrapper, quality_report parsing, MIMAG label, gate functions
    ├── redundancy.py   skani sketch/search/triangle wrappers, Union-Find, evaluate()
    ├── rna.py          barrnap (≥1.0) or barrnap 0.9 + aragorn, parallel, per-genome temp dirs
    ├── report.py       full table, filters table, summary, sankey input, ledger
    ├── classify_environments.py   verbatim copy of Env-info/classify_environments.py
    ├── compile_barrnap.py         verbatim copy of Barrnap/Script/compile_barrnap.py
    ├── sankey_generic.py          copy of Paper/sankey_generic.py with one bug fixed (see §7)
    └── pyassembly/ pytaxonkit/ pypipeline/   copies of ~/Scripts/{Assembly_metadata,Taxonkit,pipeline}
                                              packages; the three sys.path shims in pypipeline that
                                              pointed to ../../Assembly_metadata were removed
```

No file in the bundle references a path outside it (checked with `grep` for `/home`, `/mnt`,
`~/Scripts`, `~/Database`).

---

## 4. How it works

### 4.1 Commands

* `setup` — resolve tools (create env if needed), ensure the CheckM2 database, write
  `environment.json`.
* `fetch-db` — guided first-time setup (`lib/fetchdb.py`, standard library only): figshare API
  for the file list and MD5s, resumable download (`Range` requests, `.part` file), MD5 check,
  streaming flat unzip of `fna/*.fna` (the archive also holds an empty `fna/@eaDir/`), copy of the
  bundled v1.0 table, `release.json` with version v1.0 → `Archaea_HQ-v1.0/`. Confirmation prompts
  default to yes when stdin is not a terminal.
* `check` — stages env (datasets only) + list; writes `checks/check_<date>/new_accessions.tsv`.
* `run` — all eight stages; needs `--db-fna` (FASTA folder) or `--db-sketch` (skani sketch).
* `sketch-db` — `skani sketch -l <list> -o <out>` of the database genomes (4.9 GB for 21,644
  genomes, ~2 min on 24 threads) for distribution to machines without the FASTA files.
* `release` — applies the `Add`/`Replace` rows of a run's full table to the current database
  version and writes the next one into `<releases-dir>/Archaea_HQ-v<major>.<minor>/` (`fna/` with
  hard links when on the same file system, else copies; `ArchaeaHQ-Info.tsv`, `release_changes.tsv`,
  `release.json`, `comparison.txt`). Prints the old-vs-new comparison (per kingdom, summary
  statistics, environment categories), then deletes the superseded version folder after verifying
  every carried genome (name + size) in the new folder (`--keep-previous` to skip; folders with
  orphan FASTA files, the bundle, `data/` and the work directory are never deleted). Version
  discovery is in `common.py` (`find_releases`, `newest_release`, `release_version_of_table`); the
  `Context` resolves `--db-table`/`--db-fna` defaults from the newest release. Genome names are
  matched with `strip_fna` on both sides because 480 rows of the v1.0 table carry a `.fna` suffix
  in `Name`.

### 4.2 Run layout and resumability

`<workdir>/runs/<run_name>/NN_<stage>/` for the stage outputs, `commands.log` (every external
command with its output), `archaeahq_update.log`, `environment.json`, `state.json`. Each stage
writes a `.done` marker; rerunning the same command reuses finished stages, `--from-stage X`
clears markers from X on. Caches live in `<workdir>/cache/` (assembly metadata SQLite, taxonomy
dump, CheckM2 database, DB sketch keyed by a hash of the FASTA list).

### 4.3 Stages (function names in `archaeahq_update.py`)

| # | Stage | Function | External command | Output |
|---|---|---|---|---|
| 0 | env | `stage_env` | `conda env list --json`, `checkm2 database --download` | `environment.json` |
| 1 | list | `stage_list` | `datasets summary genome taxon <taxid> --assembly-source all --as-json-lines \| dataformat tsv genome --fields …` | `01_list/ncbi_<kingdom>.tsv`, `new_accessions.tsv` |
| 2 | download | `stage_download` → `ncbi.download_genomes` | `datasets download genome accession --inputfile … --include genome` (batches of 400, unzip in Python, retry) | `02_download/fna/*.fna`, `download_status.tsv` |
| 3 | metadata | `stage_metadata` → `pypipeline.runner.run()` | NCBI Datasets REST API (in-process) | `03_metadata/combined.tsv` |
| 4 | quality | `stage_quality` → `quality.run_checkm2` | `checkm2 predict --threads N --input fna -x fna --output-directory … --database_path … --force` | `04_quality/quality_report.tsv` |
| 5 | redundancy | `stage_redundancy` → `redundancy.*` | `skani search -d <sketch> --ql … --min-af 0`, `skani triangle -l … -E --min-af 0` | `05_redundancy/{new_vs_db,new_vs_new,redundancy}.tsv` |
| 6 | rna | `stage_rna` → `rna.run_all` | `barrnap --kingdom arc --threads 1 --quiet --rrna --trna [--no-ncrna --no-mrna --no-operon]` per genome, in a thread pool | `06_rna/gff/*.gff`, `rna_counts.tsv` |
| 7 | report | `stage_report` → `build_rows` | (`sankey_generic.py` with `--sankey`) | the two result tables, `summary.txt`, `sankey_input.tsv`, ledger |

The environment category (stage 7 in the plan) is computed inside `build_rows` by calling
`classify_environments.classify()` with the classifier's expected keys (note the historical
spelling `Metagenome_soource`).

### 4.4 Decision rules (`data/thresholds.json`, `redundancy.evaluate`, `build_rows`)

```
Downloaded                  PASS if a valid FASTA was obtained (SKIPPED for non-downloaded twins)
Completeness_ge70           CheckM2 completeness ≥ 70
Contamination_le10          CheckM2 contamination ≤ 10
Not_identical_to_DB         FAIL if skani hit ANI ≥ 99.9 and min(AF_ref, AF_query) ≥ 99,
                            or same numeric accession as a DB genome (GCA_x / GCF_x, version-agnostic)
Not_same_species_as_DB      FAIL if skani hit ANI ≥ 95 and min(AF) ≥ 95   (Notes: "higher quality score
                            than DB genome X" when the new genome's Park score is higher)
Representative_within_batch Union-Find over new genomes at ANI ≥ 95 / AF ≥ 95 (+ twins); representative =
                            min over (not eligible, redundant-with-DB, −Park, not RefSeq, −version, name)
Recommendation              Add  ⇔  every filter PASS
                            Replace  ⇔  newer version of a DB genome (same GCA/GCF record, higher version)
                                        that was downloaded and passes completeness/contamination
```

`Park Quality_Score = completeness − 5·contamination + 0.5·log10(contig N50)` — the dRep-style
score the database column was computed with.

### 4.5 The two result tables (`report.py`)

* `ArchaeaHQ_update_<date>_full_table.tsv` — `Recommendation`, `Reason`, then `common.DB_COLUMNS`
  (the 25 database columns in order). Sorted Add first, then kingdom, then accession.
* `ArchaeaHQ_update_<date>_filters.tsv` — `Assembly ID`, `Archaeal_Kingdom`, `Overall`, the six
  filter columns, `Notes`.

`report.stage_counts_from_rows` derives the nested per-kingdom counts used by the console table,
`summary.txt` and `sankey_input.tsv` (New → Downloaded → Contam ≤10% → Compl ≥70% → Not in DB →
New species → Batch rep.). `append_ledger` replaces the rows of an earlier attempt of the same run.

### 4.6 User interface (`ui.py`)

rich `Theme` built from the website CSS variables: primary `#c084fc`, secondary `#54c7fc`, lime
`#b9ef00`, tertiary `#f3ffcd`, muted `#adaaaa`/`#6f6f6f`, outline `#484847`; kingdom colours
Eury `#5fa4ee`, TACK `#858ef6`, DPANN `#3427f0`, Asgard `#914cf5`. Elements: banner panel, stage
rule with elapsed time, kingdom-coloured tables with a total row, progress bars (download, barrnap),
status spinners (CheckM2, skani, metadata), a verdict panel (Add / Do not add / reasons / files) and a
red error panel with the failing command, log tail and resume hint. `--plain` or `NO_COLOR` gives a
text-only version; log files are always plain.

---

## 5. Build chronology (what was actually done, in order)

1. Traced the v1.0 build (Section 1) and confirmed the exact final-selection rule and the Park
   score formula against the data.
2. Inventoried tools on the machine: `datasets` 18.26 (conda env `ncbi_datasets`), skani 0.3.1,
   prodigal 2.6.3, barrnap 1.7.0 dev (`~/software/barrnap`), aragorn, taxonkit; **no** CheckM2 or dRep.
3. Checked external interfaces: `dataformat tsv genome` field names (`assminfo-paired-assm-accession`,
   not `-acc`); `skani sketch` writes a consolidated database (`sketches.db`, `markers.bin`,
   `index.db`) that is queried with `skani search -d`, not `skani dist`; a dehydrated/normal
   `datasets download` zip has `ncbi_dataset/data/<acc>/<acc>_<asm>_genomic.fna`; bioconda
   `barrnap 1.10.x` bundles infernal + aragorn (same lineage as the 1.7 dev fork).
4. Scaffolded the bundle: copied the toolkit packages, classifier, compile_barrnap, sankey script,
   the 25-column table; removed the sibling-path shims in `pypipeline`; generated
   `evaluated_accessions.tsv` from `Initial_Dataset.txt`.
5. Wrote `common.py`, `ui.py`, `envcheck.py`, `ncbi.py`, `quality.py`, `redundancy.py`, `rna.py`,
   `report.py`, `archaeahq_update.py`, `environment.yml`, `README.md`.
6. Created the conda env from `environment.yml` (checkm2 1.1.0, datasets 18.36, skani 0.3.2,
   barrnap 1.10.5 resolved together without conflicts) and downloaded the CheckM2 database (3.1 GB).
7. Built the database sketch with `sketch-db` (21,644 genomes, 1 min 40 s, 4.9 GB).
8. Tested, fixed, re-tested (Sections 6 and 7).

---

## 6. Tests performed

| Test | Result |
|---|---|
| `check` against live NCBI, all kingdoms | 43,055 assemblies listed; 21,643 in DB; 14,347 evaluated before; **7,065 new** (3,077 Eury, 3,039 TACK, 614 DPANN, 335 Asgard); 14 GCA/GCF twins; 1 min 19 s |
| Bundled metadata toolkit on 16 accessions | 36-column `combined.tsv`, taxdump downloaded and indexed on first use (32 s) |
| barrnap 1.7 (PATH) and 1.10.5 (conda) on two DB genomes | rRNA/tRNA = 4/29 and 10/45, identical to the database columns |
| Control `run --accessions-file` (4 DB genomes, 2 twins, 10 new; `--db-sketch`) | DB genomes: FAIL `Not_identical_to_DB` ("same accession", skani ANI 100 / AF 99.97); twins: `Downloaded=SKIPPED`, twin rule, quality copied; GCA/GCF pair of a non-DB genome: GCF kept, GCA FAIL `Representative_within_batch`; one new DPANN genome FAIL `Not_same_species_as_DB` (ANI 99.98, AF 95.5) with a "higher quality score" note; CheckM2 values of the DB genomes reproduced (contamination 0.55 vs 0.57 on one, CheckM2 version difference); header from column 3 byte-identical to the DB table; 69 s |
| Resume | `--from-stage rna` and `--from-stage report` reused stages 1–6 / 1–7 |
| `--sankey` | `sankey.html` + `sankey.png` rendered from `sankey_input.tsv` |
| `run --kingdoms asgard --max-genomes 30` (real NCBI diff path) | 30 downloaded, 19 pass quality, none redundant → 19 Add, 11 Do_not_add (completeness); 3 deleted taxids resolved to Asgard via the listing; 1 min 42 s |
| Portability | folder copied to a fresh directory, `check --kingdoms asgard` with a fresh work dir: OK, 7 s |
| Ledger | rows of a re-run replace the previous attempt of the same run (no duplicates) |

Runtime estimate for a full update of the 7,065 currently new genomes: 4–5 h, dominated by
CheckM2 (~1 min per 30 genomes on 24 threads).

---

## 7. Pitfalls met and how they are handled

* **barrnap 1.x writes fixed-name scratch files** (`barrnap.find_operon.gff`) into the current
  directory → parallel runs collided and 12/14 genomes got no counts. Fix: each barrnap process
  runs in its own temporary directory (`rna.run_one`).
* **barrnap 1.10 vs 1.7 flags**: in 1.10 tRNA is OFF by default (`--trna` must be explicit) and
  `--no-operon` does not exist — passing it silently disabled tRNA prediction. Fix: flags are
  derived from the tool's own `--help` (`rna.barrnap_flags`).
* **skani consolidated sketch** cannot be used with `skani dist -r <folder>` / `--rl`; use
  `skani search -d <folder>`. Self-comparison gives AF ≈ 99.97, hence the identity AF threshold of
  99 (as in the original dRep post-processing) rather than 100.
* **dataformat field name** is `assminfo-paired-assm-accession`; `dataformat version` prints
  `undefined` for this build, so the version column may show `?`.
* **Deleted taxids**: NCBI Datasets still reports taxids that the current taxdump marks as deleted
  (three of 30 recent Asgard MAGs). Lineage is then empty; kingdom falls back to the listing, then to
  cached listings, else `Other archaea`.
* **CRLF**: the exported information table had Windows line endings; the bundled copy was
  converted to LF so header comparisons are exact.
* **sankey_generic.py bug**: metadata keys were upper-cased but compared with lower-case dict
  keys, so `# TITLE` / `# START_LABEL` / `# END_LABEL` were ignored. Fixed in the bundled copy
  only (the original in the Paper folder still has the bug).
* **conda env discovery**: `mamba` is found first; `conda run -p <prefix>` is used for tools
  found inside environments so their own runtime (perl for barrnap, tensorflow for CheckM2) is used.
* NCBI lists one fewer database genome than expected (21,643 of 21,644): `GCA_049860595.1` has
  status `previous`, superseded by `GCA_049860595.2` (handled as `Replace`, see §2). The
  per-kingdom differences (Eury 9,895 vs 9,894, DPANN 4,889 vs 4,891) are two SeqCode genomes
  (`GCA_034190545.1` Afararchaeum irisae, `GCA_034190395.1` Chewarchaeum aethiopicum) moved by
  NCBI from DPANN to Euryarchaeota. Harmless for the diff, which is accession-based.

---

## 8. Maintenance notes

* **After adding genomes to the database**: replace `data/ArchaeaHQ-Info.tsv` (or pass
  `--db-table`) and rebuild the sketch (`sketch-db`) or point `--db-fna` at the enlarged folder
  (the sketch cache is keyed by the file list, so it rebuilds automatically).
* **Thresholds** are only in `data/thresholds.json`; **kingdoms** only in `data/kingdoms.json`
  (`ui.KINGDOM_COLOURS` duplicates the colours for the terminal theme).
* **Environment rules** live in `lib/classify_environments.py` (`RULES` list, evaluated in order).
* **Column layout** is `common.DB_COLUMNS`; `report.FULL_COLUMNS` prepends Recommendation/Reason.
* To evaluate arbitrary genomes (e.g. a collaborator's MAGs already at NCBI) use
  `run --accessions-file`; the ledger is ignored for that list.
* Possible future work: an `--apply` mode that appends `Add` rows to the table and copies FASTA
  files; running Prodigal for the accepted genomes; a `replace` recommendation when a new genome
  of an existing species has a clearly better Park score.
