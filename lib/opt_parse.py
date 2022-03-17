#############################################################################
# Parsing and printing the options and meta-info for PhyloAcc.
# Much of the error checking is done here as well.
#############################################################################

import sys
import os
import math
import argparse
import multiprocessing as mp
import lib.core as CORE

#############################################################################

def optParse(globs):
# This function handles the command line options and prepares the output directory and files.
# Defaults are set in params.py
    try:
        import psutil
        globs['psutil'] = True;
    except:
        globs['psutil'] = False;
    # Check if psutil is installed for memory usage stats.

    parser = argparse.ArgumentParser(description="degenotate: Annotation of codon degeneracy for coding sequences");

    parser.add_argument("-a", dest="annotation_file", help="A gff or gtf file that contains the coordinates of transcripts in the provided genome file (-g). Only one of -a/-g OR -s is REQUIRED.", default=False);
    parser.add_argument("-g", dest="genome_file", help="A FASTA file containing a genome. -a must also be specified. Only one of -a/-g OR -s is REQUIRED.", default=False);

    parser.add_argument("-s", dest="in_seq", help="Either a directory containing individual, in-frame coding sequence files or a single file containing multipl in-frame coding sequences on which to calculate degeneracy. Only one of -a/-g OR -s is REQUIRED.", default=False);
    # Input

    parser.add_argument("-o", dest="out_dest", help="Desired output directory. This will be created for you if it doesn't exist. Default: degenotate-[date]-[time]", default=False);
    # Output
    
    parser.add_argument("-p", dest="num_procs", help="The total number of processes that degenotate can use. Default: 1.", type=int, default=1);
    # User params

    parser.add_argument("--overwrite", dest="ow_flag", help="Set this to overwrite existing files.", action="store_true", default=False);
    parser.add_argument("--appendlog", dest="append_log_flag", help="Set this to keep the old log file even if --overwrite is specified. New log information will instead be appended to the previous log file.", action="store_true", default=False);
    # User options

    parser.add_argument("--info", dest="info_flag", help="Print some meta information about the program and exit. No other options required.", action="store_true", default=False);
    #parser.add_argument("--dryrun", dest="dryrun", help="With all options provided, set this to run through the whole pseudo-it pipeline without executing external commands.", action="store_true", default=False);
    parser.add_argument("--version", dest="version_flag", help="Simply print the version and exit. Can also be called as '-version', '-v', or '--v'", action="store_true", default=False);
    parser.add_argument("--quiet", dest="quiet_flag", help="Set this flag to prevent degenotate from reporting detailed information about each step.", action="store_true", default=False);
    # Run options

    parser.add_argument("--norun", dest="norun", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--debug", dest="debug_opt", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--nolog", dest="nolog_opt", help=argparse.SUPPRESS, action="store_true", default=False);
    # Performance tests
    args = parser.parse_args();
    # The input options and help messages

    globs['call'] = " ".join(sys.argv);
    # Save the program call for later

    if args.info_flag:
        globs['info'] = True;
        globs['log-v'] = -1;
        startProg(globs);
        return globs;
    # Parse the --info option and call startProg early if set

    if args.norun:
        globs['norun'] = True;
        globs['log-v'] = -1;
    globs['overwrite'] = args.ow_flag;
    # Check run mode options.

    if (not args.annotation_file and not args.in_seq) or (args.annotation_file and args.in_seq):
        CORE.errorOut("OP1", "One input method must be specified: -a or -s", globs);
    # Check that only one input type is specified

    if args.annotation_file:
        if not args.genome_file:
            CORE.errorOut("OP2", "A genome fasta file be specified with -g when an annotation file is given with -g", globs);
        globs['gxf-file'] = args.annotation_file;
        globs['fa-file'] = args.genome_file;

        if any(globs['gxf-file'].endswith(gff_ext) for gff_ext in ['.gff', '.gff.gz', '.gff3', '.gff3.gz']):
            globs['gxf-type'] = 'gff';
        elif any(globs['gxf-file'].endswith(gtf_ext) for gtf_ext in ['.gtf', '.gtf.gz']):
            globs['gxf-type'] = 'gtf';
        else:
            CORE.errorOut("OP3", "Cannot guess annotation file type from extension. Make sure it ends with '.gff' or '.gtf'.", globs);
        # Guess whether the input annotation file is GFF or GTF from the file extension
        # TODO: Can probably do this better, or let the user specify an option

    elif args.in_seq:
        globs['in-seq'] = args.in_seq;
        if os.path.isfile(globs['in-seq']):
            globs['in-seq-type'] = "file";
        elif os.path.isdir(globs['in-seq']):
            globs['in-seq-type'] = "directory";
    # Save the input type as a global param

    globs = CORE.fileCheck(globs); 
    # Make sure all the input files actually exist, and get their
    # full paths

    if not args.out_dest:
        globs['outdir'] = "degenotate-out-" + globs['startdatetime'];
    else:
        globs['outdir'] = args.out_dest;

    if not globs['overwrite'] and os.path.exists(globs['outdir']):
        CORE.errorOut("OP4", "Output directory already exists: " + globs['outdir'] + ". Specify new directory name OR set --overwrite to overwrite all files in that directory.", globs);

    if not os.path.isdir(globs['outdir']) and not globs['norun'] and not globs['info']:
        os.makedirs(globs['outdir']);
    # Main output dir

    globs['run-name'] = os.path.basename(os.path.normpath(globs['outdir']));
    globs['logfilename'] = os.path.join(globs['outdir'], globs['run-name'] + ".log");
    # Log file

    if not args.append_log_flag and not globs['norun']:
        logfile = open(globs['logfilename'], "w");
        logfile.write("");
        logfile.close();
    # Prep the logfile to be overwritten if --appendlog isn't specified

    if args.quiet_flag:
        globs['quiet'] = True;
    # Check the quiet option

    #globs['aln-pool'] = mp.Pool(processes=globs['num-procs']);
    #globs['scf-pool'] = mp.Pool(processes=globs['num-procs']);
    # Create the pool of processes for sCF calculation here so we copy the memory profile of the parent process
    # before we've read any large data in

    if globs['psutil']:
        globs['pids'] = [psutil.Process(os.getpid())];
    # Get the starting process ids to calculate memory usage throughout.

    startProg(globs);
    # After all the essential options have been set, call the welcome function.

    return globs;

#############################################################################

def startProg(globs):
# A nice way to start the program.
    print("#");
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# Welcome to degenotate -- Degeneracy annotation for coding transcripts.");
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# Version " + globs['version'] + " released on " + globs['releasedate']);
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# degenotate was developed by " + globs['authors']);
    #CORE.printWrite(globs['logfilename'], globs['log-v'], "# Citation:      " + globs['doi']);
    #CORE.printWrite(globs['logfilename'], globs['log-v'], "# Website:       " + globs['http']);
    #CORE.printWrite(globs['logfilename'], globs['log-v'], "# Report issues: " + globs['github']);
    CORE.printWrite(globs['logfilename'], globs['log-v'], "#");
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# The date and time at the start is: " + CORE.getDateTime());
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# Using Python version:              " + globs['pyver'] + "\n#");
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# The program was called as:         " + globs['call'] + "\n#");

    if globs['info']:
        return;
    # If --info is set, return after printing program info

    #######################

    pad = 30;
    opt_pad = 30;
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# INPUT/OUTPUT INFO:");

    if globs['gxf-file']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Annotation file:", pad) + globs['gxf-file']);
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Annotation file type (guessed from extension):", pad) + globs['gxf-type']);
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Genome file:", pad) + globs['fa-file']);
    elif globs['in-seq']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Sequence " + globs['in-seq-type'] + ":", pad) + globs['in-seq']);

    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Output directory:", pad) + globs['outdir']);
    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Log file:", pad) + os.path.basename(globs['logfilename']));
    # Input/Output
    #######################

    CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# OPTIONS INFO:");    
    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Option", pad) + CORE.spacedOut("Current setting", opt_pad) + "Current action");      

    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Processes (-p)", pad) + 
                CORE.spacedOut(str(globs['num-procs']), opt_pad) + 
                "degenotate will use this many processes.");
    # Reporting the resource options

    if globs['overwrite']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --overwrite", pad) +
                    CORE.spacedOut("True", opt_pad) + 
                    "degenotate will OVERWRITE the existing files in the specified output directory.");
    # Reporting the overwrite option.

    if not globs['quiet']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --quiet", pad) + 
                    CORE.spacedOut("False", opt_pad) + 
                    "Time, memory, and status info will be printed to the screen while degenotate is running.");
    else:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --quiet", pad) + 
                    CORE.spacedOut("True", opt_pad) + 
                    "No further information will be printed to the screen while degenotate is running.");
        CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
        CORE.printWrite(globs['logfilename'], globs['log-v'], "# Running...");
    # Reporting the quiet option.

    # if globs['debug']:
    #     CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --debug", pad) + 
    #                 CORE.spacedOut("True", opt_pad) + 
    #                 "Printing out a bit of debug info.");
    # Reporting the debug option.

    if globs['norun']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --norun", pad) + 
                    CORE.spacedOut("True", opt_pad) + 
                    "ONLY PRINTING RUNTIME INFO.");
        CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    # Reporting the norun option.

    # Other options
    #######################

#############################################################################