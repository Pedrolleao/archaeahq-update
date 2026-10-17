"""CheckM2 progress: parsing its log lines and streaming a command's output."""

import sys
from pathlib import Path

from archaeahq_update import common as C
from archaeahq_update import quality as Q


def _collect():
    events = []
    return events, (lambda label, done, total: events.append((label, done, total)))


def test_phases_and_counters(tmp_path):
    events, cb = _collect()
    t = Q.CheckM2Progress(tmp_path, cb)
    t.on_line("[10/17/2026 09:00:00 AM] INFO: Running CheckM2 version 1.1.0")
    t.on_line("[10/17/2026 09:00:01 AM] INFO: Calling genes in 600 bins with 30 threads:")
    t.on_line("    Finished processing 1 of 600 (0.17%) bins.")
    t.on_line("    Finished processing 600 of 600 (100.00%) bins.")
    t.on_line("[10/17/2026 09:05:00 AM] INFO: Calculating metadata for 600 bins with 30 threads:")
    t.on_line("    Finished processing 300 of 600 (50.00%) bins.")
    t.on_line("[10/17/2026 09:06:00 AM] INFO: Annotating input genomes with DIAMOND using 30 threads")
    t.on_line("[10/17/2026 09:40:00 AM] INFO: Processing DIAMOND output")
    t.on_line("[10/17/2026 09:41:00 AM] INFO: Predicting completeness and contamination using ML models.")
    t.on_line("[10/17/2026 09:42:00 AM] INFO: Parsing all results and constructing final output table.")
    t.on_line("[10/17/2026 09:42:01 AM] INFO: CheckM2 finished successfully.")
    labels = [e[0] for e in events]
    assert labels[0] == "starting CheckM2"
    assert ("CheckM2 1/6 calling genes (prodigal)", 0, 600) in events
    assert ("CheckM2 1/6 calling genes (prodigal)", 600, 600) in events
    assert ("CheckM2 2/6 genome metadata", 300, 600) in events
    assert ("CheckM2 3/6 DIAMOND annotation", 0, 3) in events          # ceil(600 / 250) batches
    assert ("CheckM2 4/6 processing DIAMOND output", None, None) in events
    assert labels[-1].endswith("finished")


def test_diamond_batches_from_output_folder(tmp_path):
    events, cb = _collect()
    t = Q.CheckM2Progress(tmp_path, cb)
    t.on_line("INFO: Calling genes in 1000 bins with 4 threads:")
    t.on_line("INFO: Annotating input genomes with DIAMOND using 4 threads")
    d = tmp_path / "diamond_output"
    d.mkdir()
    t.poll()
    assert events[-1] == ("CheckM2 3/6 DIAMOND annotation", 0, 4)
    (d / "DIAMOND_RESULTS_0.tsv").touch()
    t.poll()
    assert events[-1][1] == 0                      # first batch still running
    (d / "DIAMOND_RESULTS_1.tsv").touch()
    (d / "DIAMOND_RESULTS_2.tsv").touch()
    t.poll()
    assert events[-1] == ("CheckM2 3/6 DIAMOND annotation", 2, 4)


def test_phase_glued_to_counter_still_switches(tmp_path):
    events, cb = _collect()
    t = Q.CheckM2Progress(tmp_path, cb)
    t.on_line("INFO: Calling genes in 600 bins with 4 threads:")
    t.on_line("    Finished processing 600 of 600 (100.00%) bins.[10/17/2026] INFO: Calculating metadata for 600 bins with 4 threads:")
    assert events[-1] == ("CheckM2 2/6 genome metadata", 0, 600)
    t.on_line("INFO: CheckM2 finished successfully.")
    assert events[-1][0] == "CheckM2 finished"


def test_single_batch_is_labelled(tmp_path):
    events, cb = _collect()
    t = Q.CheckM2Progress(tmp_path, cb)
    t.on_line("INFO: Calling genes in 20 bins with 4 threads:")
    t.on_line("INFO: Annotating input genomes with DIAMOND using 4 threads")
    assert events[-1][2] == 1 and "one batch" in events[-1][0]


def test_run_cmd_live_splits_on_carriage_return(tmp_path):
    script = tmp_path / "noisy.py"
    script.write_text(
        "import sys, time\n"
        "sys.stdout.write('hello\\n'); sys.stdout.flush(); time.sleep(0.3)\n"
        "for i in range(1, 4): sys.stdout.write('\\r  step %d of 3' % i); sys.stdout.flush()\n"
        "sys.stderr.write('\\nINFO: done\\n')\n"
    )
    seen, polls = [], []
    log = tmp_path / "cmd.log"
    C.run_cmd_live([sys.executable, "-I", str(script)], log, seen.append, poll=lambda: polls.append(1), poll_interval=0.05)
    assert "hello" in seen and "  step 3 of 3" in seen and "INFO: done" in seen
    assert "step 2 of 3" in log.read_text()
    assert polls  # poll was called while waiting / between lines
