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
    
    parser.add_argument("-v", dest="vcf_file", help="Optional VCF file with in and outgroups to output polymorphic and fixed differences for MK tests. The VCF should contain SNPs only (no indels or structural variants).", default=False);
    parser.add_argument("-u", dest="vcf_outgroups", help="A comma separated list of sample IDs in the VCF file that make up the outgroup (e.g. 'sample1,sample2') or a file with one sample per line.", default=False);
    parser.add_argument("-e", dest="vcf_exclude", help="A comma separated list of sample IDs in the VCF file to exclude (e.g. 'sample1,sample2') or a file with one sample per line.", default=False);
    #parser.add_argument("-p", dest="polarized", help="Set this to specify that the provided vcf file is polarized (=has AA (ancestral allele) in the INFO field)", action='store_true', default=False)
    # Input

    parser.add_argument("-o", dest="out_dest", help="Desired output directory. This will be created for you if it doesn't exist. Default: degenotate-[date]-[time]", default=False);
    # Output

    parser.add_argument("-d", dest="seq_delim", help="degenotate assumes the chromosome IDs in the GFF file exactly match the sequence headers in the FASTA file. If this is not the case, use this to specify a character at which the FASTA headers will be trimmed.", default=False);
    parser.add_argument("-c", dest="write_cds", help="If a file is provided, the program will extract CDS sequences from the genome and write them to the file and exit. If no file is given with the option, a file with the name of 'cds-nt.fa' will be written to the output directory. Equivalent to '-x 0234' except this stops the program before calculating degeneracy.", nargs='?', const="default", default=False);
    parser.add_argument("-ca", dest="write_cds_aa", help="The same as -c, but writes translated amino acid sequences instead. Both -c and -ca can be specified. Default file name is 'cds-aa.fa'.", nargs='?', const="default", default=False);
    parser.add_argument("-l", dest="write_longest", help="If a file is provided, the program will extract CDS sequences from the longest transcript for each gene and write them to the file and exit. If no file is given with the option, a file with the name of 'cds-nt-longest.fa' will be written to the output directory. Both -c and -l can be specified.", nargs='?', const="default", default=False);
    parser.add_argument("-la", dest="write_longest_aa", help="The same as -l, but writes translated amino acid sequences instead. Both -l and -la can be specified. Default file name is 'cds-aa-longest.fa'.", nargs='?', const="default", default=False);
    parser.add_argument("-x", dest="extract_seq", help="Extract sites of a certain degeneracy. For instance, to extract 4-fold degenerate sites enter '4'. To extract 2- and 4-fold degenerate sites enter '24' and so on.", default=False);
    parser.add_argument("-m", dest="min_length", help="The minimum length of a transcript for it to be counted. Default (and global min): 3", default=False);
    parser.add_argument("-maf", dest="maf_cutoff", help="The minor allele frequency cutoff for MK tests. Sites where alternate alleles in the ingroup are below this frequency will be excluded. Default: 1 / 2N, where N is the number of ingroup samples", default=False);
    #parser.add_argument("-p", dest="num_procs", help="The total number of processes that degenotate can use. Default: 1.", type=int, default=1);
    # User params

    parser.add_argument("--no-fixed-in", dest="no_fixed_in_flag", help="Set this if you wish to exclude sites from the MK test in which all ingroup samples share the same alternate allele (only the reference differs).", action="store_true", default=False);
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
        # Check for index in fileCheck() below

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
        
        #if args.polarized:
        #    globs['vcf-polarized'] = args.polarized;
        ## Check if provided VCF is polarized (= has ancestral allele in the INFO field)

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

        ##########

        if args.no_fixed_in_flag:
            globs['count-fixed-alt-ingroups'] = False;
        # Parse the fixed ingroups option

        ##########

        if args.maf_cutoff:
            maf_cutoff = CORE.isPosFloat(args.maf_cutoff, maxval=1.0);
            if maf_cutoff or str(maf_cutoff) == "0.0":
                globs['ingroup-maf-cutoff'] = maf_cutoff;
            else:
                CORE.errorOut("OP8", "The minor allele frequency (-maf) must be a number between 0 and 1.", globs);
    #    if args.polarized:
    #        warnings.append("# WARNING: VCF was specified as polarized (-p) => will use information about ancestral allele");
    # Check for a VCF file and its associated options

    else:
        if args.vcf_outgroups:
            warnings.append("# WARNING: VCF outgroups (-u) were provided without a VCF file (-v). This option will be ignored.");
        if args.vcf_exclude:
            warnings.append("# WARNING: VCF samples to exclude (-e) were provided without a VCF file (-v). This option will be ignored.");
        if args.no_fixed_in_flag:
            warnings.append("# WARNING: --no-fixed-in was specified without a VCF file (-v). This option will be ignored.");
        if args.maf_cutoff:
            warnings.append("# WARNING: A minor allele frequency cutoff (-maf) was specified without a VCF file (-v). This option will be ignored.");
    # If some VCF options have been specified without a VCF file, throw some warnings

    ####################

    if args.min_length:
        min_len = CORE.isPosInt(args.min_length, minval=3);
        if not min_len:
            CORE.errorOut("OP9", "The minimum transcript length (-m) must be an integer 3 or larger.", globs);
        else:
            globs['min-len'] = min_len;
    # Parse the minimun transcript length option

    ####################

    globs = CORE.fileCheck(globs);
    # Make sure all the input files actually exist, and get their
    # full paths

    ####################

    if not args.out_dest:
        globs['outdir'] = "degenotate-out-" + globs['startdatetime'];
    else:
        globs['outdir'] = args.out_dest;

    if not globs['overwrite'] and os.path.exists(globs['outdir']):
        CORE.errorOut("OP10", "Output directory already exists: " + globs['outdir'] + ". Specify new directory name OR set --overwrite to overwrite all files in that directory.", globs);

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

    if args.write_cds:
        if args.write_cds == "default":
            globs['write-cds'] = os.path.join(globs['outdir'], "cds-nt.fa");
        else:
            globs['write-cds'] = os.path.abspath(args.write_cds);
        if os.path.isfile(globs['write-cds']) and not globs['overwrite']:
            CORE.errorOut("OP11", "File specified with -c to write CDS to already exists and --overwrite was not set. Please move the current file or specify to --overwrite it.", globs);    

    ####################

    if args.write_cds_aa:
        if args.write_cds_aa == "default":
            globs['write-cds-aa'] = os.path.join(globs['outdir'], "cds-aa.fa");
        else:
            globs['write-cds-aa'] = os.path.abspath(args.write_cds_aa);
        if os.path.isfile(globs['write-cds-aa']) and not globs['overwrite']:
            CORE.errorOut("OP12", "File specified with -ca to write CDS peptides to already exists and --overwrite was not set. Please move the current file or specify to --overwrite it.", globs);    

    ####################

    if args.write_longest:
        if args.write_longest == "default":
            globs['write-longest'] = os.path.join(globs['outdir'], "cds-nt-longest.fa");
        else:
            globs['write-longest'] = os.path.abspath(args.write_longest);
        if os.path.isfile(globs['write-longest']) and not globs['overwrite']:  
            CORE.errorOut("OP13", "File specified with -l to write CDS sequences from longest transcripts to already exists and --overwrite was not set. Please move the current file or specify to --overwrite it.", globs);    

    ####################

    if args.write_longest_aa:
        if args.write_longest_aa == "default":
            globs['write-longest-aa'] = os.path.join(globs['outdir'], "cds-aa-longest.fa");
        else:
            globs['write-longest-aa'] = os.path.abspath(args.write_longest_aa);
        if os.path.isfile(globs['write-longest-aa']) and not globs['overwrite']:  
            CORE.errorOut("OP14", "File specified with -la to write CDS peptides from longest transcripts to already exists and --overwrite was not set. Please move the current file or specify to --overwrite it.", globs);    

    ####################

    if any((args.write_cds, args.write_cds_aa, args.write_longest, args.write_longest_aa)) and not globs['gxf-file']:
        CORE.errorOut("OP15", "Extracting CDS sequences with -c, -ca, -l, or -la can only be done with an annotation file (-a) and a genome file (-g).", globs); 

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
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# The date and time at the start is:  " + CORE.getDateTime());
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# Using Python executable located at: " + globs['pyexe']);
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# Using Python version:               " + globs['pyver'] + "\n#");
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# The program was called as:          " + globs['call'] + "\n#");

    if globs['info']:
        return;
    # If --info is set, return after printing program info

    #######################

    pad = 40;
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
    if globs['write-cds']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# CDS dna output:", pad) + globs['write-cds']);
    if globs['write-cds-aa']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# CDS protein output:", pad) + globs['write-cds-aa']);
    if globs['write-longest']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Longest transcript dna output:", pad) + globs['write-longest']);
    if globs['write-longest-aa']:
        CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# Longest transcript protein output:", pad) + globs['write-longest-aa']);
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

    CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# -m", pad) +
                CORE.spacedOut(str(globs['min-len']), opt_pad) +
                "Transcripts shorter than this length will be ignored by degnotate.");
    # The min length (-m) options

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
                        CORE.spacedOut(",".join(globs['vcf-outgroups']) + " ", opt_pad) +
                        "These samples will be used as outgroups in the VCF file and all others as ingroups.");
            # Report VCF outgroup samples

            if globs['vcf-polarized']:
                CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# polarized", pad) +
                        CORE.spacedOut(str(globs['vcf-polarized']) + " ", opt_pad) +
            "provided VCF appears to be polarized (has ancestral allele field in the header) => will try to recalculate derived allele frequency and run imputed MKT framework.");
            # Check if VCF is polarized


            CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# -maf", pad) +
                        CORE.spacedOut(str(globs['ingroup-maf-cutoff']) + " ", opt_pad) +
                        "Alleles below this frequency in the ingroup samples will not be counted in MK tests.");
            # The minor allele frequency cutoff (-maf) option

            if globs['vcf-exclude']:
                CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# -e", pad) +
                        CORE.spacedOut(",".join(globs['vcf-exclude']), opt_pad) +
                        "These samples will be excluded in the VCF file.");  
            # Report samples to exclude in the input VCF

            fixed_str = "";
            fixed_bool_str = "False";
            if not globs['count-fixed-alt-ingroups']:
                fixed_str = " NOT ";
                fixed_bool_str = "True"
            CORE.printWrite(globs['logfilename'], globs['log-v'], CORE.spacedOut("# --no-fixed-in", pad) +
                        CORE.spacedOut(fixed_bool_str, opt_pad) +
                        "Sites that are fixed for an alternate allele in the ingroup samples will" + fixed_str  + "be counted in MK tests." );
            # Report fixed ingroup option
            

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
