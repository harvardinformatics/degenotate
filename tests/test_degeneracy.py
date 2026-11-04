#############################################################################
# End-to-end tests of degeneracy annotation on the tiny fixtures in tests/data.
# Expected outputs in tests/data/expected were derived independently of degenotate;
# see tests/data/README.md.
#############################################################################

import pytest

from helpers import DATA, EXPECTED, read_bed, read_text, read_tsv, run_degenotate

#############################################################################

GENOME = DATA / "tiny.fa"
GFF = DATA / "tiny.gff3"
GTF = DATA / "tiny.gtf"
CDS = DATA / "tiny-cds.fa"

def run_genome(tmp_path, annotation=GFF, genome=GENOME, *extra):
# Runs degenotate on an annotation + genome and returns the result and output directory
    result = run_degenotate(tmp_path, "-a", annotation, "-g", genome, "-o", "out", *extra);
    return result, tmp_path / "out";

@pytest.fixture(scope="module")
def gff_run(tmp_path_factory):
    return run_genome(tmp_path_factory.mktemp("gff"));

@pytest.fixture(scope="module")
def gtf_run(tmp_path_factory):
    return run_genome(tmp_path_factory.mktemp("gtf"), GTF);

#############################################################################
# Genome + annotation input

def test_gff_run_succeeds(gff_run):
    result, outdir = gff_run;
    assert result.returncode == 0, result.stdout + result.stderr;

def test_gff_bed_matches_expected(gff_run):
    _, outdir = gff_run;
    assert read_bed(outdir / "degeneracy-all-sites.bed") == read_bed(EXPECTED / "tiny.bed");

def test_gff_transcript_counts_match_expected(gff_run):
    _, outdir = gff_run;
    observed = read_tsv(outdir / "transcript-counts.tsv");
    expected = read_tsv(EXPECTED / "tiny-counts.tsv");
    assert list(observed[0].keys()) == list(expected[0].keys());
    assert observed == expected;

def test_mrna_length_is_inclusive(gff_run):
# mrna_length is end - start + 1, since GFF coordinates are 1-based and inclusive
    _, outdir = gff_run;
    observed = [ row['mrna_length'] for row in read_tsv(outdir / "transcript-counts.tsv") ];
    expected = [ row['mrna_length'] for row in read_tsv(EXPECTED / "tiny-counts.tsv") ];
    assert observed == expected;

def test_short_and_non_coding_transcripts_are_skipped(gff_run):
# TXS is shorter than the minimum length and TXN has no CDS: both are warned about and left out
    _, outdir = gff_run;
    log = read_text(outdir / "out.log");
    assert "WARNING: transcript TXS has a length shorter than the minimum" in log;
    assert "(2 < 3)" in log;
    # TXS is 263-264, which is 2 bases
    assert "WARNING: transcript TXN has no coding exons" in log;
    transcripts = { row['transcript'] for row in read_tsv(outdir / "transcript-counts.tsv") };
    assert "TXS" not in transcripts and "TXN" not in transcripts;

def test_gtf_gives_same_output_as_gff(gff_run, gtf_run):
    (_, gff_out), (gtf_result, gtf_out) = gff_run, gtf_run;
    assert gtf_result.returncode == 0, gtf_result.stdout + gtf_result.stderr;
    assert read_bed(gtf_out / "degeneracy-all-sites.bed") == read_bed(gff_out / "degeneracy-all-sites.bed");
    assert read_tsv(gtf_out / "transcript-counts.tsv") == read_tsv(gff_out / "transcript-counts.tsv");

def test_blank_line_in_gff(tmp_path):
    lines = read_text(GFF).splitlines();
    gff_with_blank = tmp_path / "blank-line.gff3";
    gff_with_blank.write_text("\n".join(lines[:5] + [""] + lines[5:]) + "\n");
    result, outdir = run_genome(tmp_path, gff_with_blank);
    assert result.returncode == 0;
    assert read_bed(outdir / "degeneracy-all-sites.bed") == read_bed(EXPECTED / "tiny.bed");

#############################################################################
# Ambiguous bases

def genome_with_n(tmp_path):
# Puts an N in the middle codon of TX3 (GCT at chr1:61-63 -> GNT)
    seq = "".join(line.strip() for line in read_text(GENOME).splitlines()[1:]);
    seq = seq[:61] + "N" + seq[62:];
    genome = tmp_path / "tiny-with-n.fa";
    genome.write_text(">chr1\n" + seq + "\n");
    return genome;

def test_codon_with_n_has_unknown_degeneracy(tmp_path):
    result, outdir = run_genome(tmp_path, GFF, genome_with_n(tmp_path));
    assert result.returncode == 0, result.stdout + result.stderr;
    rows = { row[3] : row for row in read_bed(outdir / "degeneracy-all-sites.bed") };
    for cds_pos in [3, 4, 5]:
        assert rows["TX3:" + str(cds_pos)][4] == ".";
        assert rows["TX3:" + str(cds_pos)][7] == ".";
        assert rows["TX3:" + str(cds_pos)][8] == "";
    counts = { row['transcript'] : row for row in read_tsv(outdir / "transcript-counts.tsv") };
    assert [ counts["TX3"][f] for f in ["f0", "f2", "f3", "f4"] ] == ["4", "2", "0", "0"];
    # ATG (000) and TAA (022) only

def test_translation_with_n(tmp_path):
# TX3's GNT codon translates to X
    result, outdir = run_genome(tmp_path, GFF, genome_with_n(tmp_path), "-c", "-ca");
    assert "Traceback" not in result.stderr;
    assert result.returncode == 0;
    assert ">TX3\nMX*\n" in read_text(outdir / "cds-aa.fa");

#############################################################################
# CDS input (-s)

def test_cds_input(tmp_path):
    result = run_degenotate(tmp_path, "-s", CDS, "-o", "out");
    assert result.returncode == 0, result.stdout + result.stderr;
    outdir = tmp_path / "out";
    assert read_bed(outdir / "degeneracy-all-sites.bed") == read_bed(EXPECTED / "tiny-cds-input.bed");
    assert read_tsv(outdir / "transcript-counts.tsv") == read_tsv(EXPECTED / "tiny-cds-input-counts.tsv");
    assert "WARNING: transcript OUT1 is partial with unknown frame" in read_text(outdir / "out.log");

#############################################################################
# Sequence extraction

def test_extract_four_fold_sites(tmp_path):
    result, outdir = run_genome(tmp_path, GFF, GENOME, "-x", "4");
    assert result.returncode == 0, result.stdout + result.stderr;
    assert read_text(outdir / "cds-4-fold.fa") == read_text(EXPECTED / "tiny-4fold.fa");

def test_extract_cds_sequences(tmp_path):
    result, outdir = run_genome(tmp_path, GFF, GENOME, "-c", "-ca", "-l", "-la");
    assert result.returncode == 0, result.stdout + result.stderr;
    assert read_text(outdir / "cds-nt.fa") == read_text(EXPECTED / "tiny-cds-nt.fa");
    assert read_text(outdir / "cds-aa.fa") == read_text(EXPECTED / "tiny-cds-aa.fa");
    assert read_text(outdir / "cds-nt-longest.fa") == read_text(EXPECTED / "tiny-cds-nt-longest.fa");
    assert read_text(outdir / "cds-aa-longest.fa") == read_text(EXPECTED / "tiny-cds-aa-longest.fa");
    assert not (outdir / "degeneracy-all-sites.bed").exists();
    # Extracting sequences ends the program before degeneracy is calculated

@pytest.mark.parametrize("flag", ["-ca", "-la"])
def test_amino_acid_extraction_alone_exits_early(tmp_path, flag):
    result, outdir = run_genome(tmp_path, GFF, GENOME, flag);
    assert result.returncode == 0, result.stdout + result.stderr;
    assert not (outdir / "degeneracy-all-sites.bed").exists();

@pytest.mark.parametrize("flags, listed, not_listed", [
    ([], ["Per-site degeneracy output:", "Transcript count output:"], ["CDS dna output:"]),
    (["-c"], ["CDS dna output:"], ["Per-site degeneracy output:", "Transcript count output:"]),
    (["-ca"], ["CDS protein output:"], ["Per-site degeneracy output:", "Transcript count output:"]),
    (["-la"], ["Longest transcript protein output:"], ["Per-site degeneracy output:", "Transcript count output:"]),
])
def test_log_lists_the_outputs_written(tmp_path, flags, listed, not_listed):
# The start of the log lists only the output files the run will actually write
    result, outdir = run_genome(tmp_path, GFF, GENOME, *flags);
    assert result.returncode == 0, result.stdout + result.stderr;
    log = read_text(outdir / "out.log");
    for line in listed:
        assert line in log;
    for line in not_listed:
        assert line not in log;

#############################################################################
