#############################################################################
# Shared helpers for the degenotate tests
#############################################################################

import csv
import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

#############################################################################

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "tests" / "data"
EXPECTED = DATA / "expected"
DEGENOTATE = REPO / "degenotate.py"
# Paths to the repo, the test fixtures, and the main script

HAS_VCF_DEPS = sys.version_info >= (3, 10) and all(importlib.util.find_spec(mod) for mod in ("pysam", "networkx", "scipy"))
needs_vcf_deps = pytest.mark.skipif(not HAS_VCF_DEPS, reason="VCF/MK tests need Python 3.10+ with pysam, networkx, and scipy")
# The VCF options need extra dependencies and Python 3.10+ (for itertools.pairwise)

#############################################################################

def run_degenotate(tmp_path, *args):
# Runs degenotate.py as a subprocess from within tmp_path. Running from tmp_path matters
# because some files (e.g. the .errlog written by errors before the output directory is
# created) are written to the current working directory.

    cmd = [sys.executable, str(DEGENOTATE)] + [ str(arg) for arg in args ];
    return subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True);

#############################################################################

def read_bed(path):
# Reads the per-site degeneracy bed file into a list of lists of columns

    with open(path) as bed_stream:
        return [ line.rstrip("\n").split("\t") for line in bed_stream ];

#############################################################################

def read_tsv(path):
# Reads a tab delimited file with a header into a list of dicts

    with open(path) as tsv_stream:
        return list(csv.DictReader(tsv_stream, delimiter="\t"));

#############################################################################

def read_text(path):
# Reads a whole file into a string

    with open(path) as text_stream:
        return text_stream.read();

#############################################################################

def indexed_vcf(tmp_path, name, text=None):
# Copies a plain text VCF fixture to tmp_path, then compresses and indexes it with pysam
# so the committed fixture can stay human readable. If text is given, it is written in
# place of the fixture's contents (for tests that need a modified VCF). Returns the path
# to the .vcf.gz

    import pysam

    vcf_copy = tmp_path / name;
    if text is None:
        shutil.copy(DATA / name, vcf_copy);
    else:
        vcf_copy.write_text(text);
    return Path(pysam.tabix_index(str(vcf_copy), preset="vcf", force=True));

#############################################################################
