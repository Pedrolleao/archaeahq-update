"""End-to-end checks of the `release` command on a tiny synthetic database (no external tools needed)."""

import json
from pathlib import Path

import pytest

from archaeahq_update import common as C
from archaeahq_update import envcheck, proteins, rna
from archaeahq_update.cli import main
from archaeahq_update.report import FULL_COLUMNS

KINGDOM = "Methanobacteriati-Euryarchaeota"


def _row(acc: str, name: str, **kw) -> dict:
    r = {c: "0" for c in C.DB_COLUMNS}
    r.update({"Name": name, "Assembly ID": acc, "Archaeal_Kingdom": KINGDOM, "Completeness": "95",
              "Contamination": "1", "Park Quality_Score": "90", "Genome_Size": "2000000", "Contig_N50": "50000",
              "GC_Content": "0.5", "rRNA": "3", "tRNA": "40", "Environment_Category-ArchaeaHQ": "Freshwater"})
    r.update(kw)
    return r


def _fna(folder: Path, name: str, seq: str = "ACGT" * 50) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    p = folder / f"{name}.fna"
    p.write_text(f">{name}\n{seq}\n")
    return p


def _gff(folder: Path, name: str, start: int = 11, stop: int = 30) -> None:
    """barrnap GFF with one 16S hit on the single contig of _fna (contig ID = genome name)."""
    folder.mkdir(parents=True, exist_ok=True)
    (folder / f"{name}.fna-barrnap.gff").write_text(
        f"##gff-version 3\n{name}\tinfernal\trRNA\t{start}\t{stop}\t0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA\n")


@pytest.fixture(autouse=True)
def fake_tools(monkeypatch):
    """Stand-ins for Prodigal and barrnap, so the tests need no external tools. Returns the calls made."""
    calls = {"prodigal": [], "barrnap": []}

    def predict_one(tools, fna, out_dir, log_path):
        name = C.strip_fna(fna.name)
        calls["prodigal"].append(name)
        out = out_dir / f"{name}.faa"
        out.write_text(f">{name}_1\nMPREDICTED\n")
        return out

    def run_one(tools, fna, out_dir, new_barrnap, log_path):
        name = C.strip_fna(fna.name)
        calls["barrnap"].append(name)
        _gff(out_dir, name, 1, 8)
        return name, 1, 0

    monkeypatch.setattr(proteins, "predict_one", predict_one)
    monkeypatch.setattr(rna, "run_one", run_one)
    monkeypatch.setattr(envcheck, "resolve_tools", lambda required, *a, **k: {t: None for t in required})
    monkeypatch.setattr(envcheck, "barrnap_has_trna", lambda tool: True)
    return calls


@pytest.fixture
def world(tmp_path: Path):
    """A v1.0 release folder with 3 genomes, a run folder with Add/Replace/Do_not_add rows and their FASTA."""
    rel = tmp_path / "releases"
    v10 = rel / C.release_folder_name("v1.0")
    db_rows = [
        _row("GCA_000000001.1", "GCA_000000001.1_A_genomic"),
        _row("GCA_000000002.1", "GCA_000000002.1_B_genomic.fna"),           # .fna suffix inside Name, as in v1.0
        _row("GCA_000000003.1", "GCA_000000003.1_C_genomic", Archaeal_Kingdom="Promethearchaeati-Asgard"),
    ]
    C.write_tsv(v10 / C.RELEASE_TABLE, db_rows, C.DB_COLUMNS)
    C.save_json(v10 / C.RELEASE_META, {"version": "v1.0"})
    for r in db_rows:
        _fna(v10 / "fna", C.strip_fna(r["Name"]))
    # proteins and 16S of v1.0 as fetch-db unpacks them (genome B has no 16S hit)
    (v10 / "faa").mkdir()
    for r in db_rows:
        (v10 / "faa" / f"{C.strip_fna(r['Name'])}.faa").write_text(">p_1\nMOLD\n")
    (v10 / C.RELEASE_16S_FASTA).write_text(
        ">GCA_000000001.1_A_genomic#GCA_000000001.1_A_genomic#1-8#+#16S_rRNA\nACGTACGT\n"
        ">GCA_000000003.1_C_genomic#GCA_000000003.1_C_genomic#1-4#+#16S_rRNA\nACGT\n")

    work = tmp_path / "work"
    run = work / "runs" / "run_test"
    new_fna = run / "02_download" / "fna"
    rec = [
        {"Recommendation": "Replace", "Reason": "newer version of database genome GCA_000000001.1",
         **_row("GCA_000000001.2", "GCA_000000001.2_A2_genomic")},
        {"Recommendation": "Add", "Reason": "passes all filters", **_row("GCA_000000010.1", "GCA_000000010.1_N_genomic")},
        {"Recommendation": "Add", "Reason": "passes all filters",
         **_row("GCA_000000011.1", "GCA_000000011.1_M_genomic", Archaeal_Kingdom="Promethearchaeati-Asgard")},
        {"Recommendation": "Add", "Reason": "passes all filters", **_row("GCA_000000012.1", "GCA_000000012.1_missing_genomic")},
        {"Recommendation": "Do_not_add", "Reason": "Completeness_ge70 (50%)", **_row("GCA_000000013.1", "GCA_000000013.1_bad_genomic")},
    ]
    C.write_tsv(run / "ArchaeaHQ_update_test_full_table.tsv", rec, FULL_COLUMNS)
    for r in rec[:3]:
        _fna(new_fna, r["Name"], seq="TTTT" * 60)
    for r in rec[:2]:                                    # genome M (rec[2]) has no GFF: barrnap runs at release
        _gff(run / "06_rna" / "gff", r["Name"])
    return {"rel": rel, "v10": v10, "work": work, "run": run}


def _release(world, *extra: str) -> int:
    return main(["release", "--plain", "--workdir", str(world["work"]), "--releases-dir", str(world["rel"]), *extra])


def test_dry_run_writes_nothing(world, fake_tools):
    assert _release(world, "--dry-run") == 0
    assert sorted(p.name for p in world["rel"].iterdir()) == ["Archaea_HQ-v1.0"]
    assert fake_tools == {"prodigal": [], "barrnap": []}


def test_release_builds_next_version_and_removes_previous(world):
    assert _release(world) == 0
    v11 = world["rel"] / "Archaea_HQ-v1.1"
    assert v11.is_dir() and not world["v10"].exists()
    files = sorted(p.name for p in (v11 / "fna").iterdir())
    assert files == ["GCA_000000001.2_A2_genomic.fna", "GCA_000000002.1_B_genomic.fna", "GCA_000000003.1_C_genomic.fna",
                     "GCA_000000010.1_N_genomic.fna", "GCA_000000011.1_M_genomic.fna"]
    rows = C.read_tsv(v11 / C.RELEASE_TABLE)
    assert [r["Assembly ID"] for r in rows] == ["GCA_000000001.2", "GCA_000000002.1", "GCA_000000003.1",
                                                "GCA_000000010.1", "GCA_000000011.1"]
    changes = {r["Assembly ID"]: r for r in C.read_tsv(v11 / "release_changes.tsv")}
    assert changes["GCA_000000001.2"]["Action"] == "replaced" and changes["GCA_000000001.2"]["Replaces"] == "GCA_000000001.1"
    assert changes["GCA_000000001.1"]["Action"] == "removed"
    assert changes["GCA_000000012.1"]["Action"] == "skipped" and "FASTA not found" in changes["GCA_000000012.1"]["Note"]
    assert "GCA_000000013.1" not in changes
    meta = json.loads((v11 / C.RELEASE_META).read_text())
    assert meta["version"] == "v1.1" and meta["previous_version"] == "v1.0" and meta["previous_removed"] is True
    assert (v11 / "comparison.txt").read_text().startswith("ArchaeaHQ v1.0 → v1.1")
    assert C.newest_release(world["rel"]) == v11


def test_release_builds_proteins_and_16s(world, fake_tools):
    assert _release(world) == 0
    v11 = world["rel"] / "Archaea_HQ-v1.1"
    # proteins: B and C carried over, the three new genomes predicted, A (replaced) dropped
    assert sorted(p.name for p in (v11 / "faa").iterdir()) == [
        "GCA_000000001.2_A2_genomic.faa", "GCA_000000002.1_B_genomic.faa", "GCA_000000003.1_C_genomic.faa",
        "GCA_000000010.1_N_genomic.faa", "GCA_000000011.1_M_genomic.faa"]
    assert sorted(fake_tools["prodigal"]) == ["GCA_000000001.2_A2_genomic", "GCA_000000010.1_N_genomic",
                                              "GCA_000000011.1_M_genomic"]
    assert (v11 / "faa" / "GCA_000000003.1_C_genomic.faa").read_text() == ">p_1\nMOLD\n"
    # 16S: C carried over, A dropped with its replacement, A2 and N cut at the run's GFF coordinates,
    # M run through barrnap because the run has no GFF for it
    assert fake_tools["barrnap"] == ["GCA_000000011.1_M_genomic"]
    assert list(rna.read_fasta(v11 / C.RELEASE_16S_FASTA)) == [
        ("GCA_000000001.2_A2_genomic#GCA_000000001.2_A2_genomic#11-30#+#16S_rRNA", "T" * 20),
        ("GCA_000000003.1_C_genomic#GCA_000000003.1_C_genomic#1-4#+#16S_rRNA", "ACGT"),
        ("GCA_000000010.1_N_genomic#GCA_000000010.1_N_genomic#11-30#+#16S_rRNA", "T" * 20),
        ("GCA_000000011.1_M_genomic#GCA_000000011.1_M_genomic#1-8#+#16S_rRNA", "T" * 8),
    ]
    table = C.read_tsv(v11 / C.RELEASE_16S_TABLE)
    assert [(r["Genome"], r["Length_bp"]) for r in table] == [
        ("GCA_000000001.2_A2_genomic", "20"), ("GCA_000000003.1_C_genomic", "4"),
        ("GCA_000000010.1_N_genomic", "20"), ("GCA_000000011.1_M_genomic", "8")]
    meta = json.loads((v11 / C.RELEASE_META).read_text())
    assert meta["proteins"] == {"carried": 2, "predicted": 3, "missing": 0, "prodigal_mode": "meta"}
    assert meta["rrna_16s"] == {"sequences": 4, "genomes": 4, "from_previous": 1}
    assert not (world["work"] / "release_staging").joinpath("Archaea_HQ-v1.1").exists()


def test_release_without_previous_16s_runs_barrnap_on_everything(world, fake_tools):
    (world["v10"] / C.RELEASE_16S_FASTA).unlink()
    assert _release(world) == 0
    assert sorted(fake_tools["barrnap"]) == ["GCA_000000002.1_B_genomic", "GCA_000000003.1_C_genomic",
                                             "GCA_000000011.1_M_genomic"]
    genomes = [r["Genome"] for r in C.read_tsv(world["rel"] / "Archaea_HQ-v1.1" / C.RELEASE_16S_TABLE)]
    assert genomes == ["GCA_000000001.2_A2_genomic", "GCA_000000002.1_B_genomic", "GCA_000000003.1_C_genomic",
                       "GCA_000000010.1_N_genomic", "GCA_000000011.1_M_genomic"]


def test_no_faa_keeps_previous_folder_with_proteins(world, fake_tools):
    assert _release(world, "--no-faa", "--no-16s") == 0
    v11 = world["rel"] / "Archaea_HQ-v1.1"
    assert not (v11 / "faa").exists() and not (v11 / C.RELEASE_16S_FASTA).exists()
    assert fake_tools == {"prodigal": [], "barrnap": []}
    assert world["v10"].exists()                         # its faa/ and 16S FASTA would otherwise be lost
    assert json.loads((v11 / C.RELEASE_META).read_text())["previous_removed"] is False


def test_unknown_file_keeps_previous_folder(world):
    (world["v10"] / "ArchaeaHQ-Supp_tables.xlsx").write_text("x")
    assert _release(world) == 0
    assert world["v10"].exists()


def test_keep_previous_and_explicit_version(world):
    assert _release(world, "--keep-previous", "--version", "v2.0") == 0
    assert world["v10"].exists()
    assert (world["rel"] / "Archaea_HQ-v2.0" / C.RELEASE_META).exists()
    assert json.loads((world["rel"] / "Archaea_HQ-v2.0" / C.RELEASE_META).read_text())["previous_removed"] is False


def test_accept_list_and_nothing_to_apply(world, tmp_path: Path):
    acc = tmp_path / "accept.txt"
    acc.write_text("GCA_000000010.1\n")
    assert _release(world, "--accept", str(acc)) == 0
    v11 = world["rel"] / "Archaea_HQ-v1.1"
    assert [r["Assembly ID"] for r in C.read_tsv(v11 / C.RELEASE_TABLE)][-1] == "GCA_000000010.1"
    assert len(C.read_tsv(v11 / C.RELEASE_TABLE)) == 4
    # applying the same table again: everything is already in → no v1.2 folder
    assert _release(world, "--accept", str(acc)) == 0
    assert not (world["rel"] / "Archaea_HQ-v1.2").exists()


def test_rejects_older_version_label(world):
    assert _release(world, "--version", "v0.9") != 0
