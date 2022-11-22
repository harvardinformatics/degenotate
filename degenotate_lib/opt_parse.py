#############################################################################
# Parsing and printing the options and meta-info for PhyloAcc.
# Much of the error checking is done here as well.
#############################################################################

import sys
import os
import argparse
import degenotate_lib.core as CORE

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
    
    parser.add_argument("-v", dest="vcf_file", help="Optional VCF file with in and outgroups to output polymorphic and fixed differences for MK tests.", default=False);
    parser.add_argument("-u", dest="vcf_outgroups", help="A comma separated list of sample IDs in the VCF file that make up the outgroup (e.g. 'sample1,sample2') or a file with one sample per line.", default=False);
    parser.add_argument("-e", dest="vcf_exclude", help="A comma separated list of sample IDs in the VCF file to exclude (e.g. 'sample1,sample2') or a file with one sample per line.", default=False);
    # Input

    parser.add_argument("-o", dest="out_dest", help="Desired output directory. This will be created for you if it doesn't exist. Default: degenotate-[date]-[time]", default=False);
    # Output

    parser.add_argument("-d", dest="seq_delim", help="degenotate assumes the chromosome IDs in the GFF file exactly match the sequence headers in the FASTA file. If this is not the case, use this to specify a character at which the FASTA headers will be trimmed.", default=False);
    parser.add_argument("-c", dest="write_cds", help="If a file is provided, the program will extract CDS sequences from the genome and write them to the file and exit. Equivalent to '-x 0234' except this stops the program before calculating degeneracy.", default=False);
    parser.add_argument("-l", dest="write_longest", help="If a file is provided, the program will extract CDS sequences from the longest transcript for each gene and write them to the file and exit. Both -c and -l can be specified.", default=False);
    parser.add_argument("-x", dest="extract_seq", help="Extract sites of a certain degeneracy. For instance, to extract 4-fold degenerate sites enter '4'. To extract 2- and 4-fold degenerate sites enter '24' and so on.", default=False);
    #parser.add_argument("-p", dest="num_procs", help="The total number of processes that degenotate can use. Default: 1.", type=int, default=1);
    # User params

    parser.add_argument("--overwrite", dest="ow_flag", help="Set this to overwrite existing files.", action="store_true", default=False);
    parser.add_argument("--appendlog", dest="append_log_flag", help="Set this to keep the old log file even if --overwrite is specified. New log information will instead be appended to the previous log file.", action="store_true", default=False);
    # User options

    parser.add_argument("--info", dest="info_flag", help="Print some meta information about the program and exit. No other options required.", action="store_true", default=False);
    parser.add_argument("--version", dest="version_flag", help="Simply print the version and exit. Can also be called as '-version', '-v', or '--v'", action="store_true", default=False);
    parser.add_argument("--quiet", dest="quiet_flag", help="Set this flag to prevent degenotate from reporting detailed information about each step.", action="store_true", default=False);
    # Run options

    parser.add_argument("--norun", dest="norun", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--debug", dest="debug_opt", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--nolog", dest="nolog_opt", help=argparse.SUPPRESS, action="store_true", default=False);
    # Performance tests
    args = parser.parse_args();
    # The input options and help messages

    warnings = [];
    # List of warnings to print after logfile is created

    globs['call'] = " ".join(sys.argv);
    # Save the program call for later

    ####################

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

    ####################

    if (not args.annotation_file and not args.in_seq) or (args.annotation_file and args.in_seq):
        CORE.errorOut("OP1", "One input method must be specified: -a or -s", globs);
    # Check that only one input type is specified

    ####################

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

        if args.seq_delim:
            if args.seq_delim == 'space':
                globs['seq-delim'] = " ";
            else:
                globs['seq-delim'] = args.seq_delim;
    # The special case of the -delim option for a space character.

    ####################

    elif args.in_seq:
        globs['in-seq'] = args.in_seq;
        if os.path.isfile(globs['in-seq']):
            globs['in-seq-type'] = "file";
        elif os.path.isdir(globs['in-seq']):
            globs['in-seq-type'] = "directory";
    # Save the input type as a global param

    ####################

    if args.vcf_file:
        globs['vcf-file'] = args.vcf_file;
        globs['codon-methods'].append("ns");
        # If a VCF file is supplied, add the ns method to the list of methods to apply to the input

        globs['vcf-index-file'] = globs['vcf-file'] + ".tbi";
        # Default VCF index file name from tabix, will check for it below in fileCheck

        if not args.vcf_outgroups:
            CORE.errorOut("OP4", "Outgroup samples must be specified (-u) with a vcf file (-v)", globs);
        # If a VCF file is supplied, outgroup samples must be specified

        else:
            globs['vcf-outgroups'] = args.vcf_outgroups;
            if os.path.isfile(globs['vcf-outgroups']):
                outgroups = [ line.strip() for line in open(globs['vcf-outgroups']) if line.strip() and line[0] != "#" ];
                if not outgroups:
                    CORE.errorOut("OP5", "Did not read any outgroup samples from the file provided with -u.", globs);
                globs['vcf-outgroups'] = outgroups;
            # If a file is given as -u, read the outgroup samples from it

            else:
                while ", " in globs['vcf-outgroups']:
                    globs['vcf-outgroups'] = globs['vcf-outgroups'].replace(", ", ",");
                globs['vcf-outgroups'] = globs['vcf-outgroups'].split(",");
            # If outgroups are supplied directly in the command line as a comma separated list, split
            # them here
        # Parse the outgroup samples

        ##########

        if args.vcf_exclude:
            globs['vcf-exclude'] = args.vcf_exclude;
            if os.path.isfile(globs['vcf-exclude']):
                exclude = [ line.strip() for line in open(globs['vcf-exclude']) if line.strip() and line[0] != "#" ];
                if not exclude:
                    CORE.errorOut("OP6", "Did not read any samples from the file provided with -e.", globs);
                globs['vcf-exclude'] = exclude;
            # If a file is given as -e, read the samples to exclude from the VCF from it

            else:
                while ", " in globs['vcf-exclude']:
                    globs['vcf-exclude'] = globs['vcf-exclude'].replace(", ", ",");
                globs['vcf-exclude'] = globs['vcf-exclude'].split(",");
            # If samples to exclude are supplied directly in the command line as a comma separated list, split
            # them here

            globs['vcf-outgroups'] = [ outgroup for outgroup in globs['vcf-outgroups'] if outgroup not in globs['vcf-exclude'] ];
            if not globs['vcf-outgroups']:
                CORE.errorOut("OP7", "No outgroup samples left after excluding samples specified by -e.", globs);
        # Parse the samples to exclude

    elif args.vcf_outgroups or args.vcf_exclude:
        warnings.append("# WARNING: VCF outgroups (-u) or samples to exclude (-e) were provided without a VCF file (-v). They will be ignored.")
    # Check for a VCF file

    ####################

    globs = CORE.fileCheck(globs);
    # Make sure all the input files actually exist, and get their
    # full paths

    ####################

    if args.write_cds:
        globs['write-cds'] = os.path.abspath(args.write_cds);
        if os.path.isfile(globs['write-cds']):
            if not globs['gxf-file']:
                CORE.errorOut("OP8", "Extracting CDS sequences can only be done with an annotation file (-a) and a genome file (-g).", globs);    
            if not globs['overwrite']:
                CORE.errorOut("OP9", "File specified with -c to write CDS to already exists and --overwrite was not set. Please move the current file or specify to --overwrite it.", globs);    

    ####################

    if args.write_longest:
        globs['write-longest'] = os.path.abspath(args.write_longest);
        if os.path.isfile(globs['write-longest']):
            if not globs['gxf-file']:
                CORE.errorOut("OP9", "Extracting CDS sequences from longest transcripts can only be done with an annotation file (-a) and a genome file (-g).", globs);    
            if not globs['overwrite']:
                CORE.errorOut("OP10", "File specified with -l to write CDS sequences from longest transcripts to already exists and --overwrite was not set. Please move the current file or specify to --overwrite it.", globs);    

    ####################

    if not args.out_dest:
        globs['outdir'] = "degenotate-out-" + globs['startdatetime'];
    else:
        globs['outdir'] = args.out_dest;

    if not globs['overwrite'] and os.path.exists(globs['outdir']):
        CORE.errorOut("OP11", "Output directory already exists: " + globs['outdir'] + ". Specify new directory name OR set --overwrite to overwrite all files in that directory.", globs);

    if not os.path.isdir(globs['outdir']) and not globs['norun'] and not globs['info']:
        os.makedirs(globs['outdir']);
    # Main output dir

    globs['outbed'] = os.path.join(globs['outdir'], globs['outbed']);
    # Main bed file with degeneracy for all sites

    globs['outmk'] = os.path.join(globs['outdir'], globs['outmk']);
    # MK table output
     
    globs['out-transcript'] = os.path.join(globs['outdir'], globs['out-transcript']);
    # Main bed file with degeneracy for all sites

    ####################

    if args.extract_seq:
        for char in args.extract_seq:
            if char not in ["0","2","3","4"]:
                warnings.append("# WARNING: the character '" + char + "' appears in the -x input string but does not correspond to one of the accepted folds (0,2,3,4) and will be ignored.");
            else:
                globs['extract-fold'].append(char);
        globs['extract-fold'].sort();
        # Check to see that all characters input with -x correspond to a fold and if so add them to the list to extract

        globs['outseq'] = os.path.join(globs['outdir'], "cds-" + "".join(globs['extract-fold']) + "-fold.fa");
        # Sequence output file for extracted sites     

    ####################

    globs['run-name'] = os.path.basename(os.path.normpath(globs['outdir']));
    globs['logfilename'] = os.path.join(globs['outdir'], globs['run-name'] + ".log");
    # Log file

    if not args.append_log_flag and not globs['norun']:
        logfile = open(globs['logfilename'], "w");
        logfile.write("");
        logfile.close();
    # Prep the logfile to be overwritten if --appendlog isn't specified

    if warnings:
        for warning in warnings:
            CORE.printWrite(globs['logfilename'], globs['log-v'], warning);
            globs['warnings'] += 1; 
        CORE.printWrite(globs['logfilename'], globs['log-v'], "#");
    # Print any warnings here if there were any before the logfile was created   

    ####################

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
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# Website:       " + globs['http']);
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
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Annotation file type:", pad) + globs['gxf-type']);
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Genome file:", pad) + globs['fa-file']);
    elif globs['in-seq']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Sequence " + globs['in-seq-type'] + ":", pad) + globs['in-seq']);

    if globs['vcf-file']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# VCF file:", pad) + globs['vcf-file']);

    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Output directory:", pad) + globs['outdir']);
    if globs['write-cds'] or globs['write-longest']:
        if globs['write-cds']:
            CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# CDS sequence output:", pad) + globs['write-cds']);
        if globs['write-longest']:
            CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Longest transcript output:", pad) + globs['write-longest']);
    else:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Per-site degeneracy output:", pad) + globs['outbed']);
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Transcript count output:", pad) + globs['out-transcript']);

        if globs['outseq']:
            CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Extracted sequence file:", pad) + globs['outseq']);

        if "ns" in globs['codon-methods']:
            CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# MK test count output:", pad) + globs['outmk']);
        
    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Log file:", pad) + os.path.basename(globs['logfilename']));
    # Input/Output
    #######################

    CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# OPTIONS INFO:");
    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Option", pad) + CORE.spacedOut("Current setting", opt_pad) + "Current action");

    # CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Processes (-p)", pad) +
    #             CORE.spacedOut(str(globs['num-procs']), opt_pad) +
    #             "degenotate will use this many processes.");
    # Reporting the resource options

    if globs['write-cds'] or globs['write-longest']:
        if globs['write-cds']:
            CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# -c", pad) +
                        CORE.spacedOut("True", opt_pad) +
                        "CDS sequences will be extracted from the provided genome and the program will exit.");   

        if globs['write-longest']:
            CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# -l", pad) +
                        CORE.spacedOut("True", opt_pad) +
                        "CDS sequences from longest trancsripts will be extracted from the provided genome and the program will exit.");          
    # Report whether the -c or -l option is set to exit immediately after writing CDS sequences

    else:
        if globs['vcf-file']:
            CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# -u", pad) +
                        CORE.spacedOut(",".join(globs['vcf-outgroups']), opt_pad) +
                        " These samples will be used as outgroups in the VCF file and all others as ingroups.");
            # Report VCF outgroup samples

            if globs['vcf-exclude']:
                CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# -i", pad) +
                        CORE.spacedOut(",".join(globs['vcf-exclude']), opt_pad) +
                        " These samples will be excluded in the VCF file.");  
            # Report samples to exclude in the input VCF

        if globs['outseq']:
            CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# -x", pad) +
                        CORE.spacedOut(",".join(globs['extract-fold']), opt_pad) +
                        "Extracting sites of these degeneracies.");
        # Reporting the delim option

    if globs['seq-delim']:
        if globs['seq-delim'] == " ":
            delim_str = "\" \"";
        else:
            delim_str = globs['seq-delim'];
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# -d", pad) +
                    CORE.spacedOut(delim_str, opt_pad) +
                    "degenotate will split FASTA headers at this character.");
    # Reporting the delim option

    if globs['overwrite']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --overwrite", pad) +
                    CORE.spacedOut("True", opt_pad) +
                    "degenotate will OVERWRITE the existing files in the specified output directory.");
    # Reporting the overwrite option

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
    # Reporting the quiet option

    # if globs['debug']:
    #     CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --debug", pad) +
    #                 CORE.spacedOut("True", opt_pad) +
    #                 "Printing out a bit of debug info.");
    # Reporting the debug option

    if globs['norun']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --norun", pad) +
                    CORE.spacedOut("True", opt_pad) +
                    "ONLY PRINTING RUNTIME INFO.");
        CORE.printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * 125);
    # Reporting the norun option

    # Other options
    #######################

#############################################################################
