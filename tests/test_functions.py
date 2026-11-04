#############################################################################
# Unit tests for standalone functions in degenotate_lib
#############################################################################

import gzip
import itertools
import math

import pytest

from helpers import REPO, needs_vcf_deps
import degenotate_lib.core as CORE
import degenotate_lib.degen as degen
import degenotate_lib.output as OUT
import degenotate_lib.params as params
import degenotate_lib.seq as SEQ

#############################################################################

DEGEN_DICT, CODON_DICT = degen.readCodonTable(REPO / "degenotate_lib" / "codon-table.csv");

def fisher_two_sided(a, b, c, d):
# An independent two-sided Fisher's exact test on [[a, b], [c, d]] from the hypergeometric
# distribution: sum the probabilities of all tables with the same margins that are no more
# likely than the observed one
    row1, col1, n = a + b, a + c, a + b + c + d;
    prob = lambda x : math.comb(col1, x) * math.comb(n - col1, row1 - x) / math.comb(n, row1);
    p_obs = prob(a);
    return sum(prob(x) for x in range(max(0, row1 + col1 - n), min(row1, col1) + 1) if prob(x) <= p_obs * (1 + 1e-7));

#############################################################################
# degen.py

def test_codonHamming():
    assert degen.codonHamming("ATG", "ATG") == 0;
    assert degen.codonHamming("ATG", "ATC") == 1;
    assert degen.codonHamming("CAT", "CGA") == 2;
    assert degen.codonHamming("AAA", "TTT") == 3;

def test_getFrame():
    assert degen.getFrame("ATGAAA") == 0;
    assert degen.getFrame("ATGAA") == 1;
    assert degen.getFrame("ATGA") == 2;

def test_frameError():
    assert degen.frameError("ATG", 0) is False;
    assert degen.frameError("ATGAA", 1) is False;
    assert degen.frameError("ATGA", 2) is False;
    assert degen.frameError("ATGA", 0) is True;
    assert degen.frameError("ATG", 1) is True;

@needs_vcf_deps
def test_codonPath_averages_over_all_shortest_paths():
# CAT (H) -> CGA (R) has two shortest paths:
#   CAT -> CGT (R) -> CGA (R): 1 nonsynonymous, 1 synonymous
#   CAT -> CAA (Q) -> CGA (R): 2 nonsynonymous
# so the averages are dS = 0.5, dN = 1.5
    globs = params.init();
    globs['codon-methods'] = ["degen", "ns"];
    _, codon_dict, codon_graph, globs = degen.readDegen(globs);
    ds, dn = degen.codonPath("CAT", "CGA", codon_graph, codon_dict, globs['shortest-paths']);
    assert (ds, dn) == (0.5, 1.5);

def test_compute_extended_MKT():
    pytest.importorskip("scipy");
    odds, alpha, dos, pval = degen.compute_extended_MKT(4, 2, 1, 3);
    # d = dN = 4, d0 = dS = 2, p_high = pN = 1, p0_high = pS = 3
    assert odds == pytest.approx((1 * 2) / (3 * 4));
    assert alpha == pytest.approx(1 - (1 / 3) * (2 / 4));
    assert dos == pytest.approx(4 / 6 - 1 / 4);
    assert pval == pytest.approx(fisher_two_sided(1, 3, 4, 2));

def test_compute_extended_MKT_no_divergence_gives_NA():
    pytest.importorskip("scipy");
    assert degen.compute_extended_MKT(0, 0, 1, 3) == ('NA', 'NA', 'NA', 'NA');

def test_compute_imputed_MKT():
    pytest.importorskip("scipy");
    odds, alpha, dos, pval = degen.compute_imputed_MKT(5, 6, 4, 2, 2, 4, 3, 2);
    # p = 5, p0 = 6, d = 4, d0 = 2, p_high = 2, p0_high = 4, p_low = 3, p0_low = 2
    # p_wd = p_low - p_high / p0_high * p0_low = 3 - 2 / 4 * 2 = 2, so p - p_wd = 3
    assert odds == pytest.approx((3 * 2) / (6 * 4));
    assert alpha == pytest.approx(1 - 3 / 6 * 2 / 4);
    assert dos == pytest.approx(4 / 6 - 3 / (3 + 6));
    assert pval == pytest.approx(fisher_two_sided(3, 6, 4, 2));

def test_compute_extended_MKT_fractional_counts():
# Two fixed differences in one codon can give fractional dN and dS. The neutrality index uses the
# exact counts, and Fisher's test, which needs integers, uses the counts rounded half-up:
# [[1, 6], [6.5, 0.5]] -> [[1, 6], [7, 1]]. Truncating them instead ([[1, 6], [6, 0]]) gives a different p
    pytest.importorskip("scipy");
    odds, alpha, dos, pval = degen.compute_extended_MKT(6.5, 0.5, 1, 6);
    assert odds == pytest.approx((1 * 0.5) / (6 * 6.5));
    assert dos == pytest.approx(6.5 / 7 - 1 / 7);
    assert pval == pytest.approx(fisher_two_sided(1, 6, 7, 1));
    assert pval != pytest.approx(fisher_two_sided(1, 6, 6, 0));

def test_compute_imputed_MKT_fractional_counts():
    pytest.importorskip("scipy");
    odds, alpha, dos, pval = degen.compute_imputed_MKT(2, 3, 1, 1, 1, 2, 1, 1);
    # p_wd = p_low - p_high / p0_high * p0_low = 1 - 1 / 2 * 1 = 0.5, so p - p_wd = 1.5
    assert odds == pytest.approx((1.5 * 1) / (3 * 1));
    assert dos == pytest.approx(1 / 2 - 1.5 / (1.5 + 3));
    assert pval == pytest.approx(fisher_two_sided(2, 3, 1, 1));

def test_neutralityIndex_matches_scipy_for_integer_counts():
# For integer counts, the neutrality index must be the same as the odds ratio from scipy's
# fisher_exact, including nan when a row or column of the table is empty and inf when pS * dN = 0,
# so MK tables with integer counts give the same output as when the odds ratio came from scipy
    scipy_stats = pytest.importorskip("scipy.stats");
    for pn, ps, dn, ds in itertools.product(range(4), repeat=4):
        expected = scipy_stats.fisher_exact([[pn, ps], [dn, ds]])[0];
        assert str(degen.neutralityIndex(pn, ps, dn, ds)) == str(expected), (pn, ps, dn, ds);

#############################################################################
# seq.py

def test_bioTranslator():
    assert SEQ.bioTranslator("ATGGCTTAA", CODON_DICT) == "MA*";
    with pytest.raises(AssertionError):
        SEQ.bioTranslator("ATGG", CODON_DICT);

def test_bioTranslator_unknown_codon_is_X():
# Codons with ambiguous bases translate to X, the standard unknown amino acid
    assert SEQ.bioTranslator("ATGGNTTAA", CODON_DICT) == "MX*";

def write_test_fasta(path, compressed):
    content = ">seq1 some description\nACGT\nAC\n>seq2\nGGGG\n";
    if compressed:
        with gzip.open(path, "wt") as fa_stream:
            fa_stream.write(content);
    else:
        with open(path, "w") as fa_stream:
            fa_stream.write(content);

@pytest.mark.parametrize("compressed", [False, True])
def test_readFasta(tmp_path, compressed):
    fa_path = tmp_path / ("test.fa.gz" if compressed else "test.fa");
    write_test_fasta(fa_path, compressed);
    compression = "gz" if compressed else "none";
    assert SEQ.readFasta(fa_path, compression, False) == { "seq1 some description" : "ACGTAC", "seq2" : "GGGG" };
    assert SEQ.readFasta(fa_path, compression, " ") == { "seq1" : "ACGTAC", "seq2" : "GGGG" };

#############################################################################
# core.py

@pytest.mark.parametrize("compressed, expected", [(False, "none"), (True, "gz")])
def test_detectCompression(tmp_path, compressed, expected):
    fa_path = tmp_path / "test.fa";
    write_test_fasta(fa_path, compressed);
    assert CORE.detectCompression(fa_path) == expected;

def test_isPosInt():
    assert CORE.isPosInt("3", minval=3) == 3;
    assert CORE.isPosInt("2", minval=3) is False;
    assert CORE.isPosInt("3.5") is False;
    assert CORE.isPosInt("abc") is False;

def test_isPosFloat():
    assert CORE.isPosFloat("0.5", maxval=1.0) == 0.5;
    assert CORE.isPosFloat("1.5", maxval=1.0) is False;
    assert CORE.isPosFloat("-0.1") is False;
    assert CORE.isPosFloat("abc") is False;

#############################################################################
# output.py

def test_compileBedLine_cds_input():
# With CDS input (-s) the coordinates are CDS coordinates and the strand is "."
    globs = params.init();
    outline = OUT.compileBedLine(globs, "tx", "tx", 2, "A", "GAA", 2, "E", "2", CODON_DICT);
    assert outline == ["tx", "2", "3", "tx:2", "2", ".", "A", "E", "T:D;C:D"];

def test_compileBedLine_genome_input():
# With genome input the coordinates are genome coordinates and the strand comes from the annotation
    globs = params.init();
    globs['gxf-file'] = "annotation.gff";
    globs['coords'] = { "tx" : { 0 : 103 } };
    globs['annotation'] = { "tx" : { 'strand' : "-" } };
    outline = OUT.compileBedLine(globs, "tx", "chr1", 0, "A", "ATG", 0, "M", "0", CODON_DICT);
    assert outline == ["chr1", "102", "103", "tx:0", "0", "-", "A", "M", "T:L;C:L;G:V"];

def test_compileBedLine_four_fold_and_out_of_frame():
    globs = params.init();
    assert OUT.compileBedLine(globs, "tx", "tx", 5, "T", "GCT", 2, "A", "4", CODON_DICT)[-1] == "";
    assert OUT.compileBedLine(globs, "tx", "tx", 0, "C", "", "", ".", ".", "") == ["tx", "0", "1", "tx:0", ".", ".", "C", ".", ""];

#############################################################################
