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

def init():
    globs = {
        'version' : 'Beta 1.0',
        'releasedate' : "October 2021",
        'authors' : "Timothy Sackton, Gregg Thomas, Sara Wuitchik",
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
        # Input with annotation file and genome file

        'in-seq' : False,
        'in-seq-type' : False,
        # Input by a directory with many fasta files or a single multi-fasta

        'gxf-compression' : 'none',
        'seq-compression' : 'none',
        # The type of compression used for input sequence files

        'outdir' : '',
        'run-name' : 'degenotate',
        'logfilename' : 'degenotate.errlog',
        # 'alnstatsfile' : 'phyloacc-aln-stats.csv',
        # 'scfstatsfile' : 'phyloacc-scf-stats.csv',
        # 'scftreefile' : 'phyloacc-scf.tree',
        'logdir' : '',
        'overwrite' : False,
        # I/O options

        'in-seqs' : {},
        'in-bed' : {},
        'alns' : {},
        'aln-stats' : {},
        # Sequence variables

        'num-procs' : 1,
        # Number of processes to use

        'aln-pool' : False,
        'scf-pool' : False,
        # Process pools

        'job-dir' : '',
        'job-alns' : '',
        'job-cfgs' : '',
        'job-bed' : '',
        'job-smk' : '',
        'job-out' : '',
        # Job directories

        'info' : False,
        'dryrun' : False,
        'quiet' : False,
        # Other user options

        'skip-chars' : ["-", "N"],
        'aln-stats-written' : False,
        'scf-stats-written' : False,
        'scf-tree-written' : False,
        'pad' : 82,
        'endprog' : False,
        'exit-code' : 0,
        'log-v' : 1,
        'stats' : True,
        'progstarttime' : 0,
        'stepstarttime' : 0,
        'pids' : "",
        'psutil' : False,
        'qstats' : False,
        'norun' : False,
        'debug' : False,
        'nolog' : False,
        # Internal stuff
    }

    globs['logfilename'] = "phyloacc-" + globs['startdatetime'] + ".errlog";
    # Add the runtime to the error log file name.

    return globs;