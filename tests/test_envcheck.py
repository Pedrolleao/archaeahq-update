"""CheckM2 database lookup: reuse an existing copy instead of downloading it again."""

from pathlib import Path

import pytest

from archaeahq_update import envcheck


def _dmnd(folder: Path) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    f = folder / envcheck.CHECKM2_DB_FILE
    f.write_bytes(b"x")
    return f


@pytest.fixture
def no_download(monkeypatch):
    def boom(*a, **k):
        raise AssertionError("download must not be triggered")
    monkeypatch.setattr(envcheck, "run_cmd", boom)
    monkeypatch.delenv("CHECKM2DB", raising=False)
    monkeypatch.setattr(envcheck, "checkm2_registered_db", lambda tool: None)


TOOL = envcheck.Tool("checkm2", ["checkm2"], "path", "1.1.0")


def test_work_cache_first(tmp_path, no_download):
    f = _dmnd(tmp_path / "cache" / "CheckM2_database" / "CheckM2_database")
    assert envcheck.ensure_checkm2_db(TOOL, tmp_path / "cache" / "CheckM2_database", None, None) == f


def test_registered_with_checkm2(tmp_path, no_download, monkeypatch):
    f = _dmnd(tmp_path / "elsewhere")
    monkeypatch.setattr(envcheck, "checkm2_registered_db", lambda tool: f)
    assert envcheck.ensure_checkm2_db(TOOL, tmp_path / "empty_cache", None, None) == f


def test_other_work_folder(tmp_path, no_download):
    f = _dmnd(tmp_path / "repo" / "archaeahq_update_work" / "cache" / "CheckM2_database" / "CheckM2_database")
    got = envcheck.ensure_checkm2_db(TOOL, tmp_path / "empty_cache", None, None,
                                     extra_dirs=[tmp_path / "repo" / "archaeahq_update_work" / "cache" / "CheckM2_database"])
    assert got == f


def test_env_var_beats_registered(tmp_path, no_download, monkeypatch):
    env_f = _dmnd(tmp_path / "env")
    reg_f = _dmnd(tmp_path / "reg")
    monkeypatch.setenv("CHECKM2DB", str(tmp_path / "env"))
    monkeypatch.setattr(envcheck, "checkm2_registered_db", lambda tool: reg_f)
    assert envcheck.ensure_checkm2_db(TOOL, tmp_path / "empty_cache", None, None) == env_f


def test_explicit_missing_raises(tmp_path, no_download):
    with pytest.raises(RuntimeError):
        envcheck.ensure_checkm2_db(TOOL, tmp_path / "empty_cache", str(tmp_path / "nope"), None)


def test_registered_db_parses_checkm2_output(tmp_path, monkeypatch):
    f = _dmnd(tmp_path / "db")
    class P:
        stdout = ""
        stderr = f"[10/17/2026 09:43:37 AM] INFO: {f}\n"
    monkeypatch.setattr(envcheck.subprocess, "run", lambda *a, **k: P())
    assert envcheck.checkm2_registered_db(TOOL) == f
    f.unlink()
    assert envcheck.checkm2_registered_db(TOOL) is None
