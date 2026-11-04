#############################################################################
# Checks the codon table that degenotate uses against an independent calculation
# of degeneracy from the standard genetic code
#############################################################################

import csv

from helpers import REPO
import degenotate_lib.degen as degen

#############################################################################

LIB = REPO / "degenotate_lib"
BASES = "ATCG"

def read_csv(filename):
    with open(LIB / filename) as csv_stream:
        return [ row for row in csv.reader(csv_stream) if row ];

GENETIC_CODE = { row[0] : row[1] for row in read_csv("genetic-code.csv") };
CODON_TABLE = { row[0] : (row[1], row[2]) for row in read_csv("codon-table.csv") };

#############################################################################

def expected_degeneracy(codon):
# For each position, count how many of the 3 possible mutations keep the same amino acid:
# 0 synonymous -> 0-fold, 1 -> 2-fold, 2 -> 3-fold, 3 -> 4-fold
    fold_by_syn_count = { 0 : "0", 1 : "2", 2 : "3", 3 : "4" };
    degen_str = "";
    for codon_pos in range(3):
        syn = 0;
        for new_base in BASES:
            if new_base == codon[codon_pos]:
                continue;
            new_codon = codon[:codon_pos] + new_base + codon[codon_pos+1:];
            if GENETIC_CODE[new_codon] == GENETIC_CODE[codon]:
                syn += 1;
        degen_str += fold_by_syn_count[syn];
    return degen_str;

#############################################################################

def test_genetic_code_has_all_64_codons():
    expected = { a + b + c for a in BASES for b in BASES for c in BASES };
    assert set(GENETIC_CODE) == expected;
    assert set(CODON_TABLE) == expected;

def test_codon_table_amino_acids_match_genetic_code():
    for codon, (aa, _) in CODON_TABLE.items():
        assert aa == GENETIC_CODE[codon], codon;

def test_codon_table_degeneracy_matches_independent_calculation():
    mismatches = { codon : (fold, expected_degeneracy(codon)) for codon, (_, fold) in CODON_TABLE.items() if fold != expected_degeneracy(codon) };
    assert mismatches == {};

def test_known_codons():
# A few hand-checked examples
    assert CODON_TABLE["ATG"] == ("M", "000");
    assert CODON_TABLE["GCT"] == ("A", "004");
    assert CODON_TABLE["TTA"] == ("L", "202");
    assert CODON_TABLE["ATT"] == ("I", "003");
    assert CODON_TABLE["TAA"] == ("*", "022");

def test_readCodonTable_returns_the_table():
    degen_dict, codon_dict = degen.readCodonTable(LIB / "codon-table.csv");
    assert degen_dict == { codon : fold for codon, (_, fold) in CODON_TABLE.items() };
    assert codon_dict == { codon : aa for codon, (aa, _) in CODON_TABLE.items() };

def test_readCodonTable_reads_the_file_it_is_given(tmp_path):
    table = tmp_path / "codon-table.csv";
    table.write_text("ATG,X,999\n");
    assert degen.readCodonTable(table) == ({ "ATG" : "999" }, { "ATG" : "X" });

#############################################################################
