#############################################################################
# This file holds some global variables for some of the input options.
# These global parameters should be read only -- they are not modified anywhere
# else in the code except when reading the input options.
#
# This dictionary will also be used to track output.
#############################################################################

import sys
import timeit
import lib.core as PC

#############################################################################

class StrictDict(dict):
# This prevents additional keys from being added to the global params dict in
# any other part of the code, just to help me limit it's scope
# https://stackoverflow.com/questions/32258706/how-to-prevent-key-creation-through-dkey-val
    def __setitem__(self, key, value):
        if key not in self:
            raise KeyError("{} is not a legal key of this StrictDict".format(repr(key)));
        dict.__setitem__(self, key, value);

#############################################################################

def init():
    globs_init = {
        'version' : 'Beta 1.0',
        'releasedate' : "October 2021",
        'authors' : "Timothy Sackton, Gregg Thomas",
        'doi' : '',
        'http' : '',
        'github' : '',
        'starttime' : timeit.default_timer(),
        'startdatetime' : PC.getOutTime(),
        # Meta info

        'pyver' :  ".".join(map(str, sys.version_info[:3])),
        # System info

        'call' : "",
        # Script call info

        'gxf-file' : False,
        'fa-file' : False,
        'gxf-type' : False,
        # Input with annotation file and genome file

        'in-seq' : False,
        'in-seq-type' : False,
        'seq-delim' : False,
        # Input by a directory with many fasta files or a single multi-fasta

        'vcf-file' : False,
        'vcf' : False,
        'vcf-ingroups' : [],
        'vcf-outgroups' : False,
        'vcf-exclude' : [],
        # Input VCF file

        'gxf-compression' : 'none',
        'seq-compression' : 'none',
        'vcf-compression' : 'none',
        # The type of compression used for input sequence files

        'outdir' : '',
        'outbed' : 'degeneracy-all-sites.bed',
        'out-transcript' : 'transcript-counts.tsv',
        'outmk'  : 'mk.tsv',
        'write-cds' : False,
        'run-name' : 'degenotate',
        'logfilename' : 'degenotate.errlog',
        'logdir' : '',
        'overwrite' : False,
        # I/O options

        'shortest-paths' : False,
        # Dependency functions

        'genome-seqs' : {},
        'cds-seqs' : {},
        'coords' : {},
        # Sequence variables

        'degeneracy' : {},
        # Degeneracy output

        'nonsyn' : {},
        #output of syn/nonsyn snp calculations

        'annotation' : {},
        'genekey' : [],
        # Annotation information

        'num-procs' : 1,
        # Number of processes to use; currently multiprocessing not implemented

        'codon-methods' : ["degen"],
        # which codon processing steps to carry out

        'bases' : ['A', 'T', 'C', 'G'],
        # List of standard nucleotides

        'leading-bases' : { 1 : 0, 2 : 1, 3 : 2 },
        # The number of bases to remove from the beginning of a seq
        # for each frame. frame:num bases

        'info' : False,
        'quiet' : False,
        # Other user options

        'warnings' : 0,
        'skip-chars' : ["-", "N"],
        'pad' : 82,
        'endprog' : False,
        'exit-code' : 0,
        'log-v' : 1,
        'progstarttime' : 0,
        'stepstarttime' : 0,
        'pids' : "",
        'psutil' : False,
        'norun' : False,
        'debug' : False,
        'nolog' : False,
        # Internal stuff
    }

    globs_init['logfilename'] = "degenotate-" + globs_init['startdatetime'] + ".errlog";
    # Add the runtime to the error log file name.

    globs = StrictDict(globs_init);
    # Restrict the dict from having keys added to it after this

    return globs;

#############################################################################