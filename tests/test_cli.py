#############################################################################
# Tests of the command line interface: option errors and simple flags
#############################################################################

import pytest

from helpers import DATA, REPO, run_degenotate
import degenotate_lib.params as params

#############################################################################

CDS_FILE = REPO / "test-data" / "mm10" / "ensembl" / "cds" / "test.fa"
GFF_FILE = REPO / "test-data" / "vcf" / "test.gff"
GENOME_FILE = REPO / "test-data" / "vcf" / "test.fa.gz"

#############################################################################

@pytest.mark.parametrize("args, error_code", [
    ([], "OP1"),
    # No input
    (["-a", GFF_FILE, "-s", CDS_FILE], "OP1"),
    # Both input types
    (["-a", GFF_FILE], "OP2"),
    # Annotation without a genome
    (["-a", "annotation.txt", "-g", "genome.fa"], "OP3"),
    # Annotation file extension not recognized
    (["-a", GFF_FILE, "-g", GENOME_FILE, "-v", "variants.vcf.gz"], "OP4"),
    # VCF without outgroups
    (["-s", CDS_FILE, "-v", "variants.vcf.gz", "-u", "out1"], "OP17"),
    # VCF with CDS input
    (["-s", CDS_FILE, "-m", "2"], "OP16"),
    # Minimum length below 3
    (["-s", "does-not-exist.fa"], "CORE1"),
    # Missing input file
])
def test_option_errors(tmp_path, args, error_code):
    result = run_degenotate(tmp_path, *args);
    assert result.returncode == 1;
    assert "Error " + error_code + ":" in result.stdout;

def test_genome_required_error_message(tmp_path):
    result = run_degenotate(tmp_path, "-a", GFF_FILE);
    assert "A genome fasta file must be specified with -g when an annotation file is given with -a" in result.stdout;

@pytest.mark.parametrize("flags, option", [
    (["-imp", "0.1"], "-imp"),
    (["-sfs"], "-sfs"),
])
def test_vcf_options_without_vcf_warn(tmp_path, flags, option):
    result = run_degenotate(tmp_path, "-s", DATA / "tiny-cds.fa", "-o", "out", *flags);
    assert result.returncode == 0, result.stdout + result.stderr;
    log = (tmp_path / "out" / "out.log").read_text();
    assert "(" + option + ") was specified without a VCF file (-v). This option will be ignored." in log;

def test_existing_output_directory_error(tmp_path):
    (tmp_path / "out").mkdir();
    result = run_degenotate(tmp_path, "-s", CDS_FILE, "-o", "out");
    assert result.returncode == 1;
    assert "Error OP10:" in result.stdout;

def test_cds_extraction_requires_genome_input(tmp_path):
    result = run_degenotate(tmp_path, "-s", CDS_FILE, "-c", "-o", "out");
    assert result.returncode == 1;
    assert "Error OP15:" in result.stdout;

@pytest.mark.parametrize("flag", ["--version", "-version", "--v"])
def test_version(tmp_path, flag):
    result = run_degenotate(tmp_path, flag);
    assert result.returncode == 0;
    assert "degenotate version " + params.init()['version'] in result.stdout;

#############################################################################
