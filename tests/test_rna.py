"""16S rRNA extraction from barrnap GFF coordinates."""

import gzip
from pathlib import Path

from archaeahq_update import rna

GFF_16S = "Name=16S_rRNA;Alias=SSU_rRNA_archaea;Dbxref=Rfam:RF01959;product=16S ribosomal RNA"
GFF_23S = "Name=23S_rRNA;Alias=LSU_rRNA_archaea;Dbxref=Rfam:RF02540;product=23S ribosomal RNA"


def _genome(folder: Path, name: str, contigs: dict, gff_rows: list, gz: bool = False) -> Path:
    text = "".join(f">{cid} some description\n{seq[:6]}\n{seq[6:]}\n" for cid, seq in contigs.items())
    fna = folder / (f"{name}.fna.gz" if gz else f"{name}.fna")
    if gz:
        with gzip.open(fna, "wt") as fh:
            fh.write(text)
    else:
        fna.write_text(text)
    gff_dir = folder / "gff"
    gff_dir.mkdir(exist_ok=True)
    (gff_dir / f"{fna.name}-barrnap.gff").write_text("##gff-version 3\n" + "".join(
        f"{c}\tinfernal:1.1.5\trRNA\t{s}\t{e}\t1e-40\t{st}\t.\t{info}\n" for c, s, e, st, info in gff_rows))
    return fna


def test_iter_and_write_16s(tmp_path: Path):
    a = _genome(tmp_path, "GCA_000000001.1_A_genomic", {"c1": "AAAACCCCGGGGTTTT", "c2": "ACGTACGTAAGGCCTT"},
                [("c2", 3, 8, "-", GFF_16S), ("c1", 5, 12, "+", GFF_16S), ("c1", 1, 4, "+", GFF_23S)])
    b = _genome(tmp_path, "GCA_000000002.1_B_genomic", {"x": "TTTTGGGGCCCCAAAA"},
                [("x", 1, 16, "+", GFF_16S), ("missing", 1, 4, "+", GFF_16S)], gz=True)
    c = _genome(tmp_path, "GCA_000000003.1_C_genomic", {"y": "ACGTACGT"}, [("y", 1, 4, "+", GFF_23S)])
    d = tmp_path / "GCA_000000004.1_D_genomic.fna"          # no GFF at all
    d.write_text(">z\nACGT\n")

    records = list(rna.iter_rrna([d, c, b, a], tmp_path / "gff"))
    assert records == [
        ("GCA_000000001.1_A_genomic", "GCA_000000001.1_A_genomic#c2#3-8#-#16S_rRNA", "ACGTAC"),   # GTACGT reverse-complemented
        ("GCA_000000001.1_A_genomic", "GCA_000000001.1_A_genomic#c1#5-12#+#16S_rRNA", "CCCCGGGG"),
        ("GCA_000000002.1_B_genomic", "GCA_000000002.1_B_genomic#x#1-16#+#16S_rRNA", "TTTTGGGGCCCCAAAA"),
    ]

    fasta, table = tmp_path / "Archaea_HQ-16S.fasta", tmp_path / "Archaea_HQ-16S.tsv"
    counts = rna.write_rrna(records, fasta, table)
    assert counts == {"GCA_000000001.1_A_genomic": 2, "GCA_000000002.1_B_genomic": 1}
    assert list(rna.read_fasta(fasta)) == [(sid, seq) for _, sid, seq in records]
    assert table.read_text().splitlines() == [
        "Genome\tSequence_ID\tLength_bp",
        "GCA_000000001.1_A_genomic\tGCA_000000001.1_A_genomic#c2#3-8#-#16S_rRNA\t6",
        "GCA_000000001.1_A_genomic\tGCA_000000001.1_A_genomic#c1#5-12#+#16S_rRNA\t8",
        "GCA_000000002.1_B_genomic\tGCA_000000002.1_B_genomic#x#1-16#+#16S_rRNA\t16",
    ]
    assert not list(tmp_path.glob("*.part"))
    assert rna.genome_of_seq_id("GCA_000000002.1_B_genomic.fna#x#1-16#+#16S_rRNA") == "GCA_000000002.1_B_genomic"
