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
    monkeypatch.setattr(envcheck, "checkm2_definition_db", lambda tool: None)
    monkeypatch.setattr(envcheck, "conda_share_dirs", lambda: [])
    monkeypatch.setattr(envcheck.Path, "home", classmethod(lambda cls: Path("/nonexistent-home")))


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


def test_registered_db_accepts_a_folder(tmp_path, monkeypatch):
    f = _dmnd(tmp_path / "db")
    class P:
        stdout = ""
        stderr = f"[10/17/2026 09:43:37 AM] INFO: {f.parent}\n"
    monkeypatch.setattr(envcheck.subprocess, "run", lambda *a, **k: P())
    assert envcheck.checkm2_registered_db(TOOL) == f
    P.stderr = "[10/17/2026 09:43:37 AM] INFO: Not Set\n"
    assert envcheck.checkm2_registered_db(TOOL) is None


def test_definition_file_of_conda_checkm2(tmp_path):
    f = _dmnd(tmp_path / "somewhere" / "CheckM2_database")
    prefix = tmp_path / "envs" / "checkm2"
    vdir = prefix / "lib" / "python3.9" / "site-packages" / "checkm2" / "version"
    vdir.mkdir(parents=True)
    (vdir / "diamond_path.json").write_text('{"Type": "DIAMONDDB", "DBPATH": "%s"}' % f)
    tool = envcheck.Tool("checkm2", ["conda", "run", "--no-capture-output", "-p", str(prefix), "checkm2"], "conda:checkm2")
    assert envcheck.checkm2_definition_db(tool) == f
    (vdir / "diamond_path.json").write_text('{"Type": "DIAMONDDB", "DBPATH": "Not Set"}')
    assert envcheck.checkm2_definition_db(tool) is None


def test_definition_file_of_checkm2_on_path(tmp_path):
    f = _dmnd(tmp_path / "db")
    prefix = tmp_path / "envs" / "ck"
    (prefix / "bin").mkdir(parents=True)
    (prefix / "bin" / "checkm2").write_text("#!/bin/sh\n")
    vdir = prefix / "lib" / "python3.10" / "site-packages" / "checkm2" / "version"
    vdir.mkdir(parents=True)
    (vdir / "diamond_path.json").write_text('{"Type": "DIAMONDDB", "DBPATH": "%s"}' % f)
    assert envcheck.checkm2_definition_db(envcheck.Tool("checkm2", [str(prefix / "bin" / "checkm2")], "path")) == f


def test_previous_setup_of_the_work_folder(tmp_path, no_download):
    f = _dmnd(tmp_path / "elsewhere")
    wd = tmp_path / "work"; wd.mkdir()
    (wd / "environment.json").write_text('{"checkm2_db": "%s"}' % f)
    assert envcheck.ensure_checkm2_db(TOOL, wd / "cache" / "CheckM2_database", None, None, workdir=wd) == f


def test_checkm2_default_download_folder(tmp_path, no_download, monkeypatch):
    monkeypatch.setattr(envcheck.Path, "home", classmethod(lambda cls: tmp_path))
    f = _dmnd(tmp_path / "databases" / "CheckM2_database")
    assert envcheck.ensure_checkm2_db(TOOL, tmp_path / "empty_cache", None, None) == f


def test_share_folder_of_a_conda_env(tmp_path, no_download, monkeypatch):
    f = _dmnd(tmp_path / "envs" / "checkm2" / "share" / "checkm2" / "db")
    monkeypatch.setattr(envcheck, "conda_share_dirs", lambda: [tmp_path / "envs" / "checkm2" / "share"])
    assert envcheck.ensure_checkm2_db(TOOL, tmp_path / "empty_cache", None, None) == f


def test_registered_env_prefixes(tmp_path, monkeypatch):
    monkeypatch.setattr(envcheck.Path, "home", classmethod(lambda cls: tmp_path))
    (tmp_path / ".conda").mkdir()
    env = tmp_path / "vol" / "myenv"; env.mkdir(parents=True)
    (tmp_path / ".conda" / "environments.txt").write_text(f"{env}\n{tmp_path}/gone\n\n")
    assert envcheck.registered_env_prefixes() == [str(env)]
    assert envcheck.all_conda_prefixes(None) == [str(env)]


class _UI:
    def __init__(self, answers):
        self.answers = list(answers); self.msgs = []
    def ask(self, q, default=""):
        return self.answers.pop(0) if self.answers else default
    def info(self, m): self.msgs.append(("info", m))
    def warn(self, m): self.msgs.append(("warn", m))
    def err(self, m): self.msgs.append(("err", m))


def test_asks_on_terminal_before_downloading(tmp_path, no_download):
    f = _dmnd(tmp_path / "mine")
    ui = _UI([str(tmp_path / "wrong"), str(tmp_path / "mine")])
    assert envcheck.ensure_checkm2_db(TOOL, tmp_path / "empty_cache", None, None, ui=ui) == f
    assert any(k == "err" for k, _ in ui.msgs)


def test_empty_answer_downloads(tmp_path, monkeypatch):
    monkeypatch.delenv("CHECKM2DB", raising=False)
    monkeypatch.setattr(envcheck, "checkm2_registered_db", lambda tool: None)
    monkeypatch.setattr(envcheck, "checkm2_definition_db", lambda tool: None)
    monkeypatch.setattr(envcheck, "conda_share_dirs", lambda: [])
    monkeypatch.setattr(envcheck.Path, "home", classmethod(lambda cls: Path("/nonexistent-home")))
    db_dir = tmp_path / "cache" / "CheckM2_database"
    calls = []
    def fake_download(cmd, log_path=None, **k):
        calls.append(cmd); _dmnd(db_dir / "CheckM2_database")
    monkeypatch.setattr(envcheck, "run_cmd", fake_download)
    ui = _UI([""])
    got = envcheck.ensure_checkm2_db(TOOL, db_dir, None, None, ui=ui)
    assert got == db_dir / "CheckM2_database" / envcheck.CHECKM2_DB_FILE
    assert calls and "--download" in calls[0]
