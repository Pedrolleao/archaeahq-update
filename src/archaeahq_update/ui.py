"""
ui.py — terminal presentation (rich) in the ArchaeaHQ website palette.

The palette mirrors the CSS variables of vee-lab.eu/archaeahq:
  surface #0e0e0e · primary #c084fc · secondary #54c7fc · lime #b9ef00 · tertiary #f3ffcd
  muted #adaaaa / #6f6f6f · outline #484847
  kingdoms: Eury #5fa4ee · TACK #858ef6 · DPANN #3427f0 · Asgard #914cf5
"""

from __future__ import annotations

import time
from contextlib import contextmanager
from typing import Dict, Iterable, List, Optional, Sequence

from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.progress import (BarColumn, MofNCompleteColumn, Progress, SpinnerColumn,
                           TextColumn, TimeElapsedColumn, TimeRemainingColumn)
from rich.table import Table
from rich.text import Text
from rich.theme import Theme

PALETTE = {
    "primary": "#c084fc",
    "secondary": "#54c7fc",
    "lime": "#b9ef00",
    "tertiary": "#f3ffcd",
    "muted": "#adaaaa",
    "faint": "#6f6f6f",
    "outline": "#484847",
    "pink": "#f472b6",
    "danger": "#ff6b6b",
}

KINGDOM_COLOURS = {
    "Methanobacteriati-Euryarchaeota": "#5fa4ee",
    "Thermoproteati-TACK": "#858ef6",
    "Nanobdellati-DPANN": "#3427f0",
    "Promethearchaeati-Asgard": "#914cf5",
}

THEME = Theme({
    "primary": PALETTE["primary"],
    "secondary": PALETTE["secondary"],
    "lime": PALETTE["lime"],
    "tertiary": PALETTE["tertiary"],
    "muted": PALETTE["muted"],
    "faint": PALETTE["faint"],
    "ok": f"bold {PALETTE['lime']}",
    "warn": f"bold {PALETTE['tertiary']}",
    "err": f"bold {PALETTE['danger']}",
    "stage": f"bold {PALETTE['secondary']}",
    "title": f"bold {PALETTE['primary']}",
    "rule.line": PALETTE["outline"],
    "progress.percentage": PALETTE["secondary"],
    "progress.elapsed": PALETTE["faint"],
    "progress.remaining": PALETTE["faint"],
    "bar.complete": PALETTE["primary"],
    "bar.finished": PALETTE["lime"],
    "bar.pulse": PALETTE["secondary"],
})

_console: Optional[Console] = None
_plain = False


def init(plain: bool = False) -> Console:
    global _console, _plain
    _plain = plain
    _console = Console(theme=THEME, highlight=False, no_color=plain, force_terminal=None,
                       emoji=not plain)
    return _console


def console() -> Console:
    global _console
    if _console is None:
        init(False)
    return _console


def is_plain() -> bool:
    return _plain


def kingdom_style(label: str) -> str:
    return KINGDOM_COLOURS.get(label, PALETTE["muted"])


# ── Messages ──────────────────────────────────────────────────────────────────

# Shown as the last line of the banner (add the paper link here when published).
CREDIT_LINE = "Created at Vee Lab - Radboud University"


def banner(subcommand: str, version: str, subtitle: str = "") -> None:
    c = console()
    title = Text()
    title.append("Archaea", style=f"bold {PALETTE['primary']}")
    title.append("HQ", style=f"bold {PALETTE['secondary']}")
    title.append("  update", style=f"bold {PALETTE['tertiary']}")
    title.append(f"   v{version}  ·  {subcommand}", style="muted")
    body = Text(subtitle or "Quality-controlled, curated reference database of archaeal genomes", style="faint")
    credit = Text(CREDIT_LINE, style="faint")
    c.print(Panel(Text.assemble(title, "\n", body, "\n", credit), border_style=PALETTE["outline"], box=box.ROUNDED,
                  padding=(0, 2)))


def info(msg: str) -> None:
    console().print(f"[faint]·[/] {msg}")


def ok(msg: str) -> None:
    console().print(f"[ok]✓[/] {msg}" if not _plain else f"[OK] {msg}")


def warn(msg: str) -> None:
    console().print(f"[warn]![/] {msg}" if not _plain else f"[WARN] {msg}")


def err(msg: str) -> None:
    console().print(f"[err]✗[/] {msg}" if not _plain else f"[ERROR] {msg}")


def confirm(question: str, default: bool = True) -> bool:
    """Ask yes/no on the terminal; non-interactive sessions get the default."""
    import sys
    if not sys.stdin.isatty():
        return default
    hint = "[Y/n]" if default else "[y/N]"
    try:
        ans = console().input(f"[secondary]?[/] {question} {hint} ").strip().lower()
    except (EOFError, KeyboardInterrupt):
        console().print()
        return False
    if not ans:
        return default
    return ans in ("y", "yes")


def kv(key: str, value, key_width: int = 22) -> None:
    console().print(f"  [muted]{key:<{key_width}}[/] {value}")


@contextmanager
def stage(index: int, total: int, title: str):
    """Print a stage header, run the body, print the elapsed time."""
    c = console()
    c.print()
    c.rule(Text(f"[{index}/{total}]  {title}", style="stage"), align="left")
    t0 = time.time()
    try:
        yield
    finally:
        dt = time.time() - t0
        c.print(f"[faint]  finished in {fmt_duration(dt)}[/]")


def fmt_duration(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f} s"
    m, s = divmod(int(seconds), 60)
    if m < 60:
        return f"{m} min {s:02d} s"
    h, m = divmod(m, 60)
    return f"{h} h {m:02d} min"


# ── Tables ────────────────────────────────────────────────────────────────────

def table(title: str, columns: Sequence[str], rows: Iterable[Sequence], row_styles: Optional[List[str]] = None,
          justify: Optional[Dict[str, str]] = None) -> None:
    t = Table(title=title, title_style="title", box=box.SIMPLE_HEAD, header_style="secondary",
              border_style=PALETTE["outline"], pad_edge=False, show_edge=False)
    justify = justify or {}
    rows = list(rows)
    for i, col in enumerate(columns):
        if i == 0:
            width = max([len(str(col))] + [len(str(r[0])) for r in rows]) if rows else None
            t.add_column(col, justify="left", no_wrap=True, min_width=width)
        else:
            t.add_column(col, justify=justify.get(col, "right"))
    for i, r in enumerate(rows):
        style = row_styles[i] if row_styles and i < len(row_styles) else None
        t.add_row(*[str(x) for x in r], style=style)
    console().print(t)


def kingdom_table(title: str, columns: Sequence[str], per_kingdom: Dict[str, Sequence],
                  total_row: Optional[Sequence] = None) -> None:
    """Rows keyed by kingdom label, each coloured with the kingdom colour; optional total row."""
    rows, styles = [], []
    for label, vals in per_kingdom.items():
        rows.append([label] + list(vals))
        styles.append(kingdom_style(label))
    if total_row is not None:
        rows.append(["Total"] + list(total_row))
        styles.append(f"bold {PALETTE['tertiary']}")
    table(title, ["Kingdom"] + list(columns), rows, row_styles=styles)


# ── Progress ──────────────────────────────────────────────────────────────────

class _PlainProgress:
    """Minimal stand-in for rich.progress.Progress when --plain is used."""

    def __init__(self):
        self._tasks = {}
        self._n = 0

    def __enter__(self):
        return self

    def __exit__(self, *a):
        for tid, t in self._tasks.items():
            console().print(f"  {t['desc']}: {t['done']}/{t['total']}")
        return False

    def add_task(self, description, total=None, **kw):
        self._n += 1
        self._tasks[self._n] = {"desc": description, "total": total, "done": 0}
        console().print(f"  {description} (n={total})")
        return self._n

    def update(self, task_id, advance=0, **kw):
        self._tasks[task_id]["done"] += advance
        if "total" in kw:
            self._tasks[task_id]["total"] = kw["total"]

    def advance(self, task_id, n=1):
        self.update(task_id, advance=n)


def progress(transient: bool = False):
    if _plain:
        return _PlainProgress()
    return Progress(
        SpinnerColumn(style=PALETTE["secondary"]),
        TextColumn("[muted]{task.description}"),
        BarColumn(bar_width=28),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console(),
        transient=transient,
    )


@contextmanager
def spinner(message: str):
    if _plain:
        console().print(f"  {message} ...")
        yield
        return
    with console().status(f"[muted]{message}", spinner="dots", spinner_style=PALETTE["secondary"]):
        yield


# ── Panels ────────────────────────────────────────────────────────────────────

def error_panel(title: str, detail: str, tail: str = "", hint: str = "") -> None:
    body = Text()
    body.append(detail + "\n", style="bold")
    if tail:
        body.append("\n" + tail.strip()[-2500:] + "\n", style="faint")
    if hint:
        body.append("\n" + hint, style="tertiary")
    console().print(Panel(body, title=f"[err]{title}[/]", border_style=PALETTE["danger"], box=box.ROUNDED))


def verdict_panel(add: int, no_add: int, reasons: Dict[str, int], outputs: Sequence[str], replace: int = 0) -> None:
    body = Text()
    body.append(f"  Add          {add:>7}\n", style=f"bold {PALETTE['lime']}")
    if replace:
        body.append(f"  Replace      {replace:>7}   (newer version of a database genome)\n",
                    style=f"bold {PALETTE['secondary']}")
    body.append(f"  Do not add   {no_add:>7}\n", style="muted")
    if reasons:
        body.append("\n  Reasons for not adding\n", style="secondary")
        for k, v in sorted(reasons.items(), key=lambda kv_: -kv_[1]):
            body.append(f"    {v:>6}  {k}\n", style="faint")
    body.append("\n  Output files\n", style="secondary")
    for o in outputs:
        body.append(f"    {o}\n", style="tertiary")
    console().print(Panel(body, title="[title]Recommendation[/]", border_style=PALETTE["primary"], box=box.ROUNDED))
