#!/usr/bin/env python3
"""
Launcher for running archaeahq-update straight from a clone of the repository, without installing:

    python3 archaeahq_update.py <command> ...

It puts src/ on the module path and, unless $ARCHAEAHQ_RELEASES is set, keeps the database
versions (Archaea_HQ-v*/) next to this file. The installed entry point is `archaeahq-update`.
"""

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "src"))
os.environ.setdefault("ARCHAEAHQ_RELEASES", str(HERE))

from archaeahq_update.cli import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
