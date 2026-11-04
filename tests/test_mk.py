#############################################################################
# End-to-end tests of the MK table counts and tests with a VCF file (-v).
# Every expected value here is derived by hand in tests/data/README.md.
#############################################################################

import pytest

from helpers import DATA, indexed_vcf, needs_vcf_deps, read_text, read_tsv, run_degenotate

pytestmark = needs_vcf_deps

#############################################################################

GENOME = DATA / "tiny.fa"
GFF = DATA / "tiny.gff3"
CDS = DATA / "tiny-cds.fa"

def run_mk(tmp_path, vcf_name="tiny.vcf", *extra, vcf_text=None):
# Runs degenotate with the tiny VCF (outgroup out1) and returns the result and the
# MK table as a dict keyed by transcript
    vcf = indexed_vcf(tmp_path, vcf_name, vcf_text);
    result = run_degenotate(tmp_path, "-a", GFF, "-g", GENOME, "-v", vcf, "-u", "out1", "-o", "out", *extra);
    mk_file = tmp_path / "out" / "mk.tsv";
    mk = { row['transcript'] : row for row in read_tsv(mk_file) } if mk_file.exists() else {};
    return result, mk;

def counts(row):
    return tuple(float(row[col]) for col in ["pN", "pS", "dN", "dS"]);

@pytest.fixture(scope="module")
def default_run(tmp_path_factory):
    return run_mk(tmp_path_factory.mktemp("mk"));

@pytest.fixture(scope="module")
def polarized_run(tmp_path_factory):
    return run_mk(tmp_path_factory.mktemp("mk-polarized"), "tiny-polarized.vcf");

#############################################################################
# Counts and tests

def test_mk_run_succeeds(default_run):
    result, mk = default_run;
    assert result.returncode == 0, result.stdout + result.stderr;
    assert list(next(iter(mk.values())).keys()) == ["transcript", "pN", "pS", "dN", "dS", "pval", "odds_ni", "dos"];

@pytest.mark.parametrize("transcript, expected", [
    ("TX1", (2, 3, 1, 1)),
    # pN: TGG->TCG (AC 2), GAA->GAT (AC 4); the ATT->GTT singleton is excluded by the default MAF cutoff
    # pS: GCT->GCC (AC 3 of 18), GAA->GAG (AC 2), CGA->CGG (fixed in the ingroup, counted by default)
    # dN: CGA->CAA, dS: TTA->TTG
    ("TX2", (1, 1, 1, 1)),
    # - strand: each site only gives this answer if the VCF alleles are complemented
    ("TX3", (0, 0, 0, 0)),
    ("TX4", (0, 0, 0, 0)),
    ("TX8", (0, 0, 1.5, 0.5)),
    # Two fixed differences in one codon, CAT->CGA, averaged over both shortest paths
])
def test_mk_counts(default_run, transcript, expected):
    _, mk = default_run;
    assert counts(mk[transcript]) == expected;

def test_mk_statistics(default_run):
    _, mk = default_run;
    tx1 = mk["TX1"];
    assert float(tx1['dos']) == pytest.approx(1 / 2 - 2 / 5);
    # DoS = dN / (dN + dS) - pN / (pN + pS)
    assert float(tx1['odds_ni']) == pytest.approx((2 * 1) / (3 * 1));
    # Neutrality index = (pN * dS) / (pS * dN)
    assert float(tx1['pval']) == pytest.approx(1.0);
    # Fisher's exact test on [[2, 3], [1, 1]]: every table with these margins is at most as likely
    tx2 = mk["TX2"];
    assert (float(tx2['dos']), float(tx2['odds_ni']), float(tx2['pval'])) == pytest.approx((0.0, 1.0, 1.0));

def test_mk_statistics_without_polymorphism_are_NA(default_run):
    _, mk = default_run;
    assert (mk["TX8"]['pval'], mk["TX8"]['odds_ni'], mk["TX8"]['dos']) == ("NA", "NA", "NA");

#############################################################################
# Options

def test_no_fixed_in(tmp_path):
# Excludes the CGA->CGG site where all ingroup samples are homozygous for the alternate allele
    result, mk = run_mk(tmp_path, "tiny.vcf", "--no-fixed-in");
    assert result.returncode == 0, result.stdout + result.stderr;
    assert counts(mk["TX1"]) == (2, 2, 1, 1);

def test_maf_cutoff(tmp_path):
# With -maf 0.1, only alleles with frequency above 0.1 are counted:
# pN: GAA->GAT (0.2); pS: GCT->GCC (3/18), CGA->CGG (1.0)
    result, mk = run_mk(tmp_path, "tiny.vcf", "-maf", "0.1");
    assert result.returncode == 0, result.stdout + result.stderr;
    assert counts(mk["TX1"]) == (1, 2, 1, 1);

def test_maf_cutoff_zero_includes_singletons(tmp_path):
# -maf 0 must not be replaced by the default 1/2N
    result, mk = run_mk(tmp_path, "tiny.vcf", "-maf", "0");
    assert result.returncode == 0, result.stdout + result.stderr;
    assert counts(mk["TX1"]) == (3, 3, 1, 1);

@pytest.mark.parametrize("flags, setting", [
    ([], "1 / 2N (default)"),
    (["-maf", "0.1"], "0.1"),
])
def test_maf_cutoff_in_log(tmp_path, flags, setting):
# The options summary shows the -maf setting, and the default cutoff is logged once it has
# been calculated from the number of ingroup samples in the VCF (10, so 1 / 20)
    result, _ = run_mk(tmp_path, "tiny.vcf", *flags);
    assert result.returncode == 0, result.stdout + result.stderr;
    log = read_text(tmp_path / "out" / "out.log");
    maf_lines = [ line for line in log.splitlines() if line.startswith("# -maf ") ];
    assert len(maf_lines) == 1 and maf_lines[0][len("# -maf"):].strip().startswith(setting + " ");
    assert ("the default 1 / 2N = 1 / 20 = 0.05 will be used" in log) == (not flags);

@pytest.mark.parametrize("flags, setting", [
    ([], "0.15 (default)"),
    (["-imp", "0.2"], "0.2"),
])
def test_imp_cutoff_in_log(tmp_path, flags, setting):
    result, _ = run_mk(tmp_path, "tiny-polarized.vcf", *flags);
    assert result.returncode == 0, result.stdout + result.stderr;
    imp_lines = [ line for line in read_text(tmp_path / "out" / "out.log").splitlines() if line.startswith("# -imp ") ];
    assert len(imp_lines) == 1 and imp_lines[0][len("# -imp"):].strip().startswith(setting + " ");

def test_imp_cutoff_below_maf_cutoff(tmp_path):
# -imp only splits the polymorphisms above the MAF cutoff into low and high frequency classes,
# so the counts don't change when it is below the MAF cutoff (1/20)
    result, mk = run_mk(tmp_path, "tiny.vcf", "-imp", "0.01");
    assert result.returncode == 0, result.stdout + result.stderr;
    assert counts(mk["TX1"]) == (2, 3, 1, 1);

def test_sfs_columns(tmp_path):
# -sfs lists the raw allele frequency of every polymorphism, including those below the MAF
# cutoff, in codon order, rounded to 2 decimals
    result, mk = run_mk(tmp_path, "tiny.vcf", "-sfs");
    assert result.returncode == 0, result.stdout + result.stderr;
    assert (mk["TX1"]['pn_af'], mk["TX1"]['ps_af']) == ("0.05,0.1,0.2", "0.17,1.0,0.1");
    assert (mk["TX2"]['pn_af'], mk["TX2"]['ps_af']) == ("0.1", "0.2");

#############################################################################
# Polarized VCF and the imputed MK test

def test_polarized_vcf_adds_imputed_columns(polarized_run):
    result, mk = polarized_run;
    assert result.returncode == 0, result.stdout + result.stderr;
    assert list(mk["TX1"].keys()) == ["transcript", "pN", "pS", "dN", "dS", "pval", "odds_ni", "dos", "imp.pval", "imp.odds_ni", "imp.dos"];
    assert counts(mk["TX1"]) == (2, 3, 1, 1);
    assert counts(mk["TX2"]) == (1, 1, 1, 1);

def test_imputed_dos(polarized_run):
# TX1, with low frequency defined as 1/20 < AF <= 0.15:
#   pN: low 1 (0.1), high 1 (0.2); pS: low 1 (0.1), high 2 (3/18, 1.0)
#   p_wd = p_low - p_high / p0_high * p0_low = 1 - 1/2 * 1 = 0.5
#   imp.dos = dN / (dN + dS) - (pN - p_wd) / ((pN - p_wd) + pS) = 1/2 - 1.5 / 4.5
# TX2: AA is the ALT allele for the nonsynonymous SNP, so its derived frequency flips from
# 0.1 (low) to 0.9 (high): p_wd = 0 and imp.dos = 1/2 - 1/2 = 0. Without the flip it would be 0.5
    _, mk = polarized_run;
    assert float(mk["TX1"]['imp.dos']) == pytest.approx(1 / 2 - 1.5 / 4.5);
    assert float(mk["TX2"]['imp.dos']) == pytest.approx(0.0);

def test_imputed_odds_ni_uses_exact_counts(polarized_run):
# TX1: p - p_wd = 2 - 0.5 = 1.5, so imp.odds_ni = (1.5 * dS) / (pS * dN) = 1.5 / 3. Truncating 1.5 to 1
# for Fisher's test would give 1/3
    _, mk = polarized_run;
    assert float(mk["TX1"]['imp.odds_ni']) == pytest.approx(0.5);

def test_imp_cutoff_zero(tmp_path):
# With -imp 0 there is no low frequency class, so p_wd = 0 and the imputed test is the same as
# the standard one: imp.dos = dos = 1/2 - 2/5
    result, mk = run_mk(tmp_path, "tiny-polarized.vcf", "-imp", "0");
    assert result.returncode == 0, result.stdout + result.stderr;
    assert counts(mk["TX1"]) == (2, 3, 1, 1);
    assert float(mk["TX1"]['imp.dos']) == pytest.approx(1 / 2 - 2 / 5);

def test_polarized_record_without_aa_gives_clean_error(tmp_path):
    vcf_text = read_text(DATA / "tiny-polarized.vcf").replace("\t21\t.\tT\tC\t.\tPASS\tAA=T\t", "\t21\t.\tT\tC\t.\tPASS\t.\t");
    assert "AA=T\tGT\t0/1" not in vcf_text;
    result, _ = run_mk(tmp_path, "tiny-polarized.vcf", vcf_text=vcf_text);
    assert result.returncode == 1;
    assert "Traceback" not in result.stderr;
    assert "Error VCF4:" in result.stdout;
    assert "chr1:21" in result.stdout;

def test_missing_outgroup_error(tmp_path):
    vcf = indexed_vcf(tmp_path, "tiny.vcf");
    result = run_degenotate(tmp_path, "-a", GFF, "-g", GENOME, "-v", vcf, "-u", "not-a-sample", "-o", "out");
    assert result.returncode == 1;
    assert "Error VCF3:" in result.stdout;

#############################################################################
# Variant placement edge cases (one transcript each, so they don't affect the tests above)

def test_variant_at_first_base_of_first_codon(default_run):
# TX6: ATG->GTG nonsynonymous polymorphism (AC 2) at the first CDS base
    _, mk = default_run;
    assert counts(mk["TX6"]) == (1, 0, 0, 0);

def test_minus_strand_partial_transcript(default_run):
# TX5 (- strand, 1 leading and 2 trailing bases): fixed AAA->CAA is nonsynonymous. If the CDS
# position were offset by the trailing bases instead of the leading ones, it would land on
# CAT->CAC instead, which is synonymous
    _, mk = default_run;
    assert counts(mk["TX5"]) == (0, 0, 1, 0);

def test_variant_at_first_base_of_transcript(default_run):
# TX7 (- strand) starts at its CDS; TAG->TAA synonymous polymorphism (AC 2) at that first base.
# This checks that the 1-based start is converted for pysam's 0-based fetch
    _, mk = default_run;
    assert counts(mk["TX7"]) == (0, 1, 0, 0);

def test_vcf_with_cds_input_gives_clean_error(tmp_path):
    vcf = indexed_vcf(tmp_path, "tiny.vcf");
    result = run_degenotate(tmp_path, "-s", CDS, "-v", vcf, "-u", "out1", "-o", "out");
    assert result.returncode == 1;
    assert "Traceback" not in result.stderr;
    assert "Error OP17:" in result.stdout;

#############################################################################
