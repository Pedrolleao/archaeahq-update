"""Pure-function checks: accession parsing, version labels, release folder discovery."""

from pathlib import Path

import pytest

from archaeahq_update import common as C


def test_accession_helpers():
    assert C.acc_parts("GCA_000123.2") == ("GCA", "000123", 2)
    assert C.acc_parts("not-an-accession") is None
    assert C.acc_numeric("GCF_000123.2") == "000123"
    assert C.acc_version("GCA_000123.7") == 7
    assert C.name_to_acc("GCA_000007185.1_ASM718v1_genomic") == "GCA_000007185.1"
    assert C.strip_fna("x_genomic.fna") == "x_genomic"
    assert C.strip_fna("x_genomic.fna.gz") == "x_genomic"
    assert C.strip_fna("x_genomic") == "x_genomic"


def test_version_labels():
    assert C.parse_version("v1.1") == (1, 1)
    assert C.parse_version("1.10") == (1, 10)
    assert C.parse_version("Archaea_HQ-v2.0") == (2, 0)
    assert C.parse_version("fna") is None
    assert C.next_version("v1.9") == "v1.10"
    assert C.release_folder_name("v2.0") == "Archaea_HQ-v2.0"
    with pytest.raises(ValueError):
        C.next_version("db")


def test_release_discovery(tmp_path: Path):
    assert C.newest_release(tmp_path) is None
    for v in ("v1.0", "v1.2", "v1.10"):
        d = tmp_path / C.release_folder_name(v)
        d.mkdir()
        (d / C.RELEASE_TABLE).write_text("Name\tAssembly ID\n")
        C.save_json(d / C.RELEASE_META, {"version": v})
    (tmp_path / "Archaea_HQ-v9.9").mkdir()            # no release.json → ignored
    assert [d.name for d in C.find_releases(tmp_path)] == ["Archaea_HQ-v1.0", "Archaea_HQ-v1.2", "Archaea_HQ-v1.10"]
    assert C.newest_release(tmp_path).name == "Archaea_HQ-v1.10"
    assert C.release_version_of_table(tmp_path / "Archaea_HQ-v1.2" / C.RELEASE_TABLE) == "v1.2"
    assert C.release_version_of_table(C.DB_TABLE_PATH) == C.BUNDLED_VERSION


def test_bundled_data_present():
    assert C.DB_TABLE_PATH.is_file()
    assert C.EVALUATED_PATH.is_file()
    assert C.ENVIRONMENT_YML.is_file()
    assert {k["short"] for k in C.load_kingdoms()} == {"eury", "tack", "dpann", "asgard"}
    assert C.load_thresholds()["quality"]["min_completeness"] == 70
