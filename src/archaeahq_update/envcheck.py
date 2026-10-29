"""
envcheck.py — locate the external tools, create the bundled conda environment when
something is missing, and make sure the CheckM2 database is present.

Resolution order for each tool:
  1. already on PATH
  2. inside any conda environment (prefers an env called `archaeahq-update`)
  3. create the bundled environment (environment.yml) and look again
Tools found inside a conda env are run through `conda run -p <prefix>` so that their
own runtime (perl/python/tensorflow) is used.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from .common import BUNDLE_DIR, ENVIRONMENT_YML, CommandError, log, run_cmd

DEFAULT_ENV_NAME = "archaeahq-update"
REQUIRED_TOOLS = ["datasets", "dataformat", "checkm2", "skani", "barrnap", "prodigal"]
OPTIONAL_TOOLS = ["aragorn"]
CHECKM2_DB_FILE = "uniref100.KO.1.dmnd"


@dataclass
class Tool:
    name: str
    cmd: List[str]                 # how to invoke it (may be a conda run wrapper)
    source: str                    # "path" | "conda:<prefix>"
    version: str = ""

    def __call__(self, *args) -> List[str]:
        return self.cmd + [str(a) for a in args]

    def describe(self) -> str:
        return f"{self.name} {self.version} ({self.source})"


# ── conda discovery ───────────────────────────────────────────────────────────

def find_conda() -> Optional[str]:
    """Return the conda executable (or mamba/micromamba) if any is available."""
    for name in ("mamba", "conda", "micromamba"):
        p = shutil.which(name)
        if p:
            return p
    env = os.environ.get("CONDA_EXE")
    if env and Path(env).exists():
        return env
    home = Path.home()
    for cand in ("miniforge3", "mambaforge", "miniconda3", "anaconda3", "conda"):
        for bin_name in ("conda", "mamba"):
            p = home / cand / "bin" / bin_name
            if p.exists():
                return str(p)
    return None


def _conda_base(conda: str) -> Optional[str]:
    """Return the classic `conda` executable next to a mamba binary (for `conda run`)."""
    if os.path.basename(conda) == "conda":
        return conda
    sib = Path(conda).parent / "conda"
    if sib.exists():
        return str(sib)
    p = shutil.which("conda")
    return p


def conda_env_prefixes(conda: str) -> List[str]:
    exe = _conda_base(conda) or conda
    try:
        out = subprocess.run([exe, "env", "list", "--json"], capture_output=True, text=True, timeout=120)
        if out.returncode == 0:
            return list(json.loads(out.stdout).get("envs", []))
    except (subprocess.SubprocessError, json.JSONDecodeError, OSError):
        pass
    # Fallback: scan the usual envs folder
    prefixes = []
    root = Path(exe).resolve().parent.parent if exe else None
    if root and (root / "envs").is_dir():
        prefixes = [str(p) for p in (root / "envs").iterdir() if p.is_dir()]
    return prefixes


def registered_env_prefixes() -> List[str]:
    """
    Environments listed in ~/.conda/environments.txt: every conda installation of the user writes there,
    so it also covers environments created with another conda (a cluster module, an old Miniconda) or
    with `-p` outside any envs/ folder, which `conda env list` of the conda on PATH does not show.
    """
    txt = Path.home() / ".conda" / "environments.txt"
    try:
        lines = txt.read_text().splitlines()
    except OSError:
        return []
    return [ln.strip() for ln in lines if ln.strip() and Path(ln.strip()).is_dir()]


def all_conda_prefixes(conda: Optional[str]) -> List[str]:
    """`conda env list` of the conda found, plus the registered environments, without duplicates."""
    seen: List[str] = []
    for p in (conda_env_prefixes(conda) if conda else []) + registered_env_prefixes():
        if p not in seen:
            seen.append(p)
    return seen


def env_prefix_named(conda: str, name: str) -> Optional[str]:
    for p in conda_env_prefixes(conda):
        if Path(p).name == name:
            return p
    return None


# ── tool lookup ───────────────────────────────────────────────────────────────

def _version_of(cmd: List[str], name: str) -> str:
    flags = {
        "datasets": [["--version"], ["version"]], "dataformat": [["--version"], ["version"]],
        "checkm2": [["--version"]], "skani": [["--version"]], "barrnap": [["--version"]], "aragorn": [["-h"]], "prodigal": [["-v"]],
    }
    for fl in flags.get(name, [["--version"]]):
        try:
            out = subprocess.run(cmd + fl, capture_output=True, text=True, timeout=120)
            text = (out.stdout or "") + (out.stderr or "")
            m = re.search(r"(\d+\.\d+(?:\.\d+)?)", text)
            if m:
                return m.group(1)
        except (subprocess.SubprocessError, OSError):
            continue
    return ""


def find_tool(name: str, conda: Optional[str], preferred_env: str = DEFAULT_ENV_NAME) -> Optional[Tool]:
    p = shutil.which(name)
    if p:
        t = Tool(name, [p], "path")
        t.version = _version_of(t.cmd, name)
        return t
    if not conda:
        return None
    exe = _conda_base(conda) or conda
    prefixes = all_conda_prefixes(conda)
    prefixes.sort(key=lambda x: (Path(x).name != preferred_env, x))
    for prefix in prefixes:
        if (Path(prefix) / "bin" / name).exists():
            t = Tool(name, [exe, "run", "--no-capture-output", "-p", prefix, name], f"conda:{Path(prefix).name}")
            t.version = _version_of(t.cmd, name)
            return t
    return None


def create_bundled_env(conda: str, env_name: str, log_path: Optional[Path] = None) -> str:
    """Create the bundled conda environment and return its prefix."""
    if not ENVIRONMENT_YML.exists():
        raise RuntimeError(f"environment.yml not found in bundle: {ENVIRONMENT_YML}")
    run_cmd([conda, "env", "create", "-y", "-f", str(ENVIRONMENT_YML), "-n", env_name], log_path=log_path)
    prefix = env_prefix_named(conda, env_name)
    if not prefix:
        raise RuntimeError(f"conda reported success but environment '{env_name}' was not found")
    return prefix


def resolve_tools(required: Sequence[str], optional: Sequence[str], auto_install: bool,
                  env_name: str, log_path: Optional[Path], ui=None) -> Dict[str, Tool]:
    """Find every tool; create the bundled environment once if any required tool is missing."""
    conda = find_conda()
    tools: Dict[str, Tool] = {}
    missing: List[str] = []
    for name in list(required) + list(optional):
        t = find_tool(name, conda, env_name)
        if t:
            tools[name] = t
        elif name in required:
            missing.append(name)
    if missing and auto_install:
        if not conda:
            raise RuntimeError("Missing tools (" + ", ".join(missing) +
                               ") and no conda/mamba installation was found to create the bundled environment.")
        if ui:
            ui.warn(f"missing: {', '.join(missing)} → creating conda environment '{env_name}' from environment.yml "
                    f"(this can take 10-30 min the first time)")
        create_bundled_env(conda, env_name, log_path)
        still = []
        for name in missing:
            t = find_tool(name, conda, env_name)
            if t:
                tools[name] = t
            else:
                still.append(name)
        for name in optional:
            if name not in tools:
                t = find_tool(name, conda, env_name)
                if t:
                    tools[name] = t
        missing = still
    if missing:
        raise RuntimeError("Required tools not available: " + ", ".join(missing) +
                           ". Run `archaeahq_update.py setup` or install them into a conda environment.")
    return tools


# ── CheckM2 database ──────────────────────────────────────────────────────────

def find_checkm2_db(candidates: Sequence[Optional[str]]) -> Optional[Path]:
    """Return the path of uniref100.KO.1.dmnd searching the candidate directories (recursively, 3 levels)."""
    for c in candidates:
        if not c:
            continue
        p = Path(c).expanduser()
        if p.is_file() and p.name == CHECKM2_DB_FILE:
            return p
        if p.is_dir():
            for depth in ("", "*/", "*/*/"):
                hits = list(p.glob(f"{depth}{CHECKM2_DB_FILE}"))
                if hits:
                    return hits[0]
    return None


def checkm2_install_prefix(checkm2: Tool) -> Optional[Path]:
    """The conda prefix (or install root) CheckM2 runs from."""
    if "-p" in checkm2.cmd:
        return Path(checkm2.cmd[checkm2.cmd.index("-p") + 1])
    try:
        exe = Path(checkm2.cmd[0]).resolve()
    except OSError:
        return None
    return exe.parent.parent if exe.parent.name in ("bin", "Scripts") else None


def checkm2_definition_db(checkm2: Tool) -> Optional[Path]:
    """
    The database CheckM2 itself remembers, read from its definition file
    (site-packages/checkm2/version/diamond_path.json, written by `checkm2 database --download` and
    `--setdblocation`). Reading the file avoids starting CheckM2, which loads TensorFlow and can take
    minutes on a cluster filesystem. Returns None when the file is missing or says "Not Set".
    """
    prefix = checkm2_install_prefix(checkm2)
    if not prefix:
        return None
    for f in sorted(prefix.glob("lib/python3*/site-packages/checkm2/version/diamond_path.json")):
        try:
            loc = json.loads(f.read_text()).get("DBPATH", "")
        except (OSError, json.JSONDecodeError, AttributeError):
            continue
        if loc and loc != "Not Set":
            found = find_checkm2_db([loc])
            if found:
                return found
    return None


def checkm2_registered_db(checkm2: Tool) -> Optional[Path]:
    """The database CheckM2 reports with `checkm2 database --current` (a .dmnd file or its folder)."""
    try:
        out = subprocess.run(checkm2("database", "--current"), capture_output=True, text=True, timeout=600)
    except (subprocess.SubprocessError, OSError):
        return None
    for line in reversed((out.stdout + "\n" + out.stderr).splitlines()):
        loc = line.rsplit("INFO:", 1)[-1].strip() if "INFO:" in line else line.strip()
        if loc and loc != "Not Set" and not loc.startswith("["):
            found = find_checkm2_db([loc])
            if found:
                return found
    return None


def remembered_checkm2_db(workdir: Optional[Path]) -> Optional[Path]:
    """The database a previous `setup` or run of this work folder used (environment.json)."""
    if not workdir:
        return None
    f = workdir / "environment.json"
    try:
        loc = json.loads(f.read_text()).get("checkm2_db")
    except (OSError, json.JSONDecodeError, AttributeError):
        return None
    return find_checkm2_db([loc]) if loc else None


def conda_share_dirs() -> List[Path]:
    """share/ folders of every conda environment: where a database installed next to a tool usually lives."""
    return [Path(p) / "share" for p in all_conda_prefixes(find_conda()) if (Path(p) / "share").is_dir()]


def ask_checkm2_db(ui, db_dir: Path) -> Optional[Path]:
    """
    On a terminal, ask for the path of an existing database before downloading 3 GB; an empty answer
    means download. Without a terminal `ui.ask` returns the default at once, so batch jobs never wait.
    """
    for _ in range(3):
        ans = ui.ask(f"Path of an existing CheckM2 database ({CHECKM2_DB_FILE} or its folder), "
                     f"or Enter to download it (~3 GB) into {db_dir}:", "")
        if not ans:
            return None
        found = find_checkm2_db([ans])
        if found:
            return found
        ui.err(f"{CHECKM2_DB_FILE} not found under {ans}")
    return None


def ensure_checkm2_db(checkm2: Tool, db_dir: Path, explicit: Optional[str], log_path: Optional[Path],
                      ui=None, extra_dirs: Sequence[Path] = (), workdir: Optional[Path] = None) -> Path:
    """
    Find uniref100.KO.1.dmnd without downloading it again. Search order: --checkm2-db, $CHECKM2DB, the
    database used by the previous setup/run of this work folder (environment.json), the cache of this
    work folder, the location CheckM2 itself remembers (its definition file, else
    `checkm2 database --current`), CheckM2's default download folder ~/databases, the caches of other
    work folders (`extra_dirs`), and share/ of every conda environment. When none has it, a terminal
    user is asked for a path; otherwise (or on an empty answer) it is downloaded into `db_dir`.
    """
    def say(msg):
        if ui:
            ui.info(msg)

    found = find_checkm2_db([explicit])
    if found:
        return found
    if explicit:
        raise RuntimeError(f"CheckM2 database not found under --checkm2-db {explicit}")
    found = find_checkm2_db([os.environ.get("CHECKM2DB")])
    if found:
        say(f"CheckM2 database from $CHECKM2DB: {found}")
        return found
    found = remembered_checkm2_db(workdir)
    if found:
        say(f"CheckM2 database used by the previous setup of this work folder: {found}")
        return found
    found = find_checkm2_db([db_dir])
    if found:
        return found
    found = checkm2_definition_db(checkm2) or checkm2_registered_db(checkm2)
    if found:
        say(f"CheckM2 database that checkm2 itself remembers: {found}")
        return found
    home_db = Path.home() / "databases"
    found = find_checkm2_db([home_db])
    if found:
        say(f"CheckM2 database in CheckM2's default download folder: {found}")
        return found
    for d in extra_dirs:
        found = find_checkm2_db([d])
        if found:
            say(f"CheckM2 database found in another work folder: {found}")
            return found
    for d in conda_share_dirs():
        found = find_checkm2_db([d])
        if found:
            say(f"CheckM2 database found in a conda environment: {found}")
            return found
    looked = ["--checkm2-db", "$CHECKM2DB", str(db_dir), "`checkm2 database --current`", str(home_db)] + \
             [str(d) for d in extra_dirs] + ["share/ of the conda environments"]
    if ui:
        ui.warn(f"CheckM2 database not found (looked at {', '.join(looked)})")
        found = ask_checkm2_db(ui, db_dir)
        if found:
            return found
        ui.warn(f"downloading the CheckM2 database (~3 GB) into {db_dir}; "
                f"pass --checkm2-db or set $CHECKM2DB to reuse an existing copy")
    db_dir.mkdir(parents=True, exist_ok=True)
    run_cmd(checkm2("database", "--download", "--path", str(db_dir)), log_path=log_path)
    found = find_checkm2_db([db_dir])
    if not found:
        raise RuntimeError(f"CheckM2 database download finished but {CHECKM2_DB_FILE} was not found in {db_dir}")
    return found


def barrnap_has_trna(barrnap: Tool) -> bool:
    """barrnap >= 1.0 (tseemann dev / bioconda 1.10) predicts tRNA via aragorn; 0.9 does rRNA only."""
    try:
        return int(barrnap.version.split(".")[0]) >= 1
    except (ValueError, IndexError):
        return False


def environment_summary(tools: Dict[str, Tool], checkm2_db: Optional[Path]) -> dict:
    return {
        "tools": {k: {"cmd": v.cmd, "source": v.source, "version": v.version} for k, v in tools.items()},
        "checkm2_db": str(checkm2_db) if checkm2_db else None,
        "bundle": str(BUNDLE_DIR),
    }
