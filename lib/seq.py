#############################################################################
# Functions to read sequences for degenotate
#############################################################################

import sys
import os
import gzip
import lib.core as CORE
import multiprocessing as mp
from itertools import groupby

############################################################################# 

def readFasta(filename, seq_compression, seq_delim):
# Read a FASTA formatted sequence file
# Great iterator and groupby code from: https://www.biostars.org/p/710/ 

    if seq_compression == "gz":
        file_stream = gzip.open(filename);
        fa_iter = (x[1] for x in groupby(file_stream, lambda line: line.decode()[0] == ">"));
        readstr = lambda s : s.decode().strip();
    elif seq_compression == "none":
        file_stream = open(filename); 
        fa_iter = (x[1] for x in groupby(file_stream, lambda line: line[0] == ">"));
        readstr = lambda s : s.strip();
    # Read the lines of the file depending on the compression level
    # file_stream opens the file as an iterable
    # groupby takes an iterable (file_stream) and a function that indicates the key of the group. It iterates over
    # the iterable and when it encounters a key, it groups all following items to it (perfect for FASTA files).
    # fa_iter is a generator object holding all the grouped iterators.
    # readstr is a function that changes depending on compression level -- for compressed files we also need to decode
    # each string in the iterators below.

    seqdict = {};
    # A dictionary of sequences:
    # <sequence id/header> : <sequence>

    for header_obj in fa_iter:
        header = readstr(header_obj.__next__());
        # The header object is an iterator. This gets the string.

        curkey = header[1:];
        if seq_delim:
            curkey = curkey.split(seq_delim)[0];         
        # This removes the ">" character from the header string to act as the key in seqdict
        # and splits the header based on user input from the -d option

        seq = "".join(readstr(s) for s in fa_iter.__next__());
        # The current header should correspond to the current iterator in fa_iter. This gets all those
        # lines and combines them as a string.

        #print(header, len(seq));

        seqdict[curkey] = seq;
        # Save the sequence in the dictionary

    return seqdict;

#############################################################################

def readGenome(globs):
# A function that reads an entire genome fasta file into memory

    step = "Detecting compression of genome FASTA file";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    globs['seq-compression'] = CORE.detectCompression(globs['fa-file']);
    if globs['seq-compression'] == "none":
        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: No compression detected");
    else:
        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + globs['seq-compression'] + " detected");
    # Detect the compression of the input sequence file

    step = "Reading genome FASTA file";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    globs['genome-seqs'] = readFasta(globs['fa-file'], globs['seq-compression'], globs['seq-delim']);
    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(len(globs['genome-seqs'])) + " seqs read");
    # Read the input sequence file

    #print(list(globs['genome-seqs'].keys()))

    return globs;

#############################################################################

def checkHeaders(globs):
    step = "Checking headers";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    annotation_headers = set([ globs['annotation'][t]['header'] for t in globs['annotation'] ]);
    # Extract unique headers from annotation file

    for header in annotation_headers:
        if header not in globs['genome-seqs']:
            print();
            CORE.errorOut("SEQ1", "Region in annotation file not found in genome file: " + header + ". Reminder: you can use -d to trim FASTA headers at a given character.", globs);
    # Check each header in the annotation file against those in the FASTA file and print an error if one isn't found

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
    # Status update

#############################################################################

def extractCDS(globs):
# This takes the coordiantes read from the input annotation file as well as the sequence read from the
# input genome fasta file and extracts coding sequences and coordinates for the CDS of all transcripts
# while accounting for strand

    step = "Extracting CDS";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    complement = { 'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A', 'N' : 'N',
                   'a' : 't', 'c' : 'g', 'g' : 'c', 't' : 'a', 'n' : 'n'  };
    # A dictionary with base complements for use when reverse complementing sequences

    for transcript in globs['annotation']:

        if len(globs['annotation'][transcript]['exons']) == 0:
            continue;
            #no exons means this transcript does not have a CDS, so we skip it

        cur_seq = "";
        # Initialize the sequence string for the current transcript. This will be added to the 'seqs' dict later

        globs['coords'][transcript] = {};
        globs['coords-rev'][transcript] = {};
        cds_coord = 0;
        # Initialize the coord lookup dict for this transcript and start the coord count at 0

        header = globs['annotation'][transcript]['header'];
        strand = globs['annotation'][transcript]['strand'];
        # Unpack some info about the transcript

        exons = globs['annotation'][transcript]['exons'];
        # Get the exons for the current transcript

        if not all(exons[exon]['strand'] == strand for exon in exons):
            print("\n\n");
            print(transcript, strand);
            print(exons);
            print("\n\n");
            CORE.errorOut("SEQ2", "Some exons have differing strands", globs);
        # Add check to make sure exons all have same strand as transcript?

        exon_coords = { exons[exon]['start']-1 : exons[exon]['end'] for exon in exons };
        exon_phase = { exons[exon]['start'] : exons[exon]['phase'] for exon in exons };
        # Get the coordinates of all the exons in this transcript
        # Subtract 1 from the starting coord because GXF are 1-based and python strings are 0-based

        if strand == "+":
            sorted_starts = sorted(list(exon_coords.keys()));

        elif strand == "-":
            sorted_starts = sorted(list(exon_coords.keys()), reverse=True);
         
        first_exon_genome_start = sorted_starts[0] + 1;
        first_exon_genome_end = exon_coords[sorted_starts[0]]

        if strand == "+":
            globs['annotation'][transcript]['coding-start'] = first_exon_genome_start
        elif strand == "-":
            globs['annotation'][transcript]['coding-start'] = first_exon_genome_end
        
        globs['annotation'][transcript]['start-frame'] = int(exon_phase[first_exon_genome_start])

        # Make sure the exons are sorted correctly, reversing the order if the strand is "-"

        for start in sorted_starts:
            cur_exon_seq = globs['genome-seqs'][header][start:exon_coords[start]];
            # For each exon starting coordinate, extract the sequence that corresponds to the current header and
            # start and end coordinates

            cur_exon_len = len(cur_exon_seq);
            # Get the length of the current exon to count up coordinates

            cds_coord_list = list(range(cds_coord, cds_coord+cur_exon_len));
            # The list of coordinates in the current CDS relative to the first CDS

            genome_coord = start + 1;
            # The starting genome coordinate for the current CDS

            genome_coord_list = list(range(genome_coord, genome_coord+cur_exon_len));
            # The list of genome coordinates in the current CDS

            if strand == "-": 
                cur_exon_seq = "".join(complement.get(base, base) for base in reversed(cur_exon_seq));
                # Reverse complement the sequence of the current CDS
                
                genome_coord_list.reverse(); 
                # Reverse the order of the genome coordinates

                # for i in range(cur_exon_len):
                #     globs['coords'][transcript][tcoord] = gcoord;
                #     gcoord -= 1;
                #     tcoord += 1;                    
            # If the strand is "-", get the reverse complement of the sequence and reverse the coordinates

            cur_seq += cur_exon_seq;
            # Concatenate the current exon sequence onto the overall transcript sequence

            for i in range(len(cds_coord_list)):
                globs['coords'][transcript][cds_coord_list[i]] = genome_coord_list[i];
                globs['coords-rev'][transcript][genome_coord_list[i]] = cds_coord_list[i];
            # Add the pairs of coordinates (CDS:genome) to the coords dict for this transcript

            cds_coord += cur_exon_len;
            # Increment the CDS coordinate by the length of the current CDS so the next CDS has the correct starting coord

        # End CDS loop
        ##########

        globs['cds-seqs'][transcript] = cur_seq.upper();
        # Save the current transcript sequence to the global seqs dict

        # if strand == "+":
        #     print();
        #     print(transcript);
        #     print(globs['cds-seqs'][transcript]);
        #     print(globs['coords'][transcript]);
        #     sys.exit();

        # End transcript loop
        ##########

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(len(globs['cds-seqs'])) + " CDS read");
    # Status update

    ###
    if globs['write-cds']:
        step = "Writing CDS sequences";
        step_start_time = CORE.report_step(globs, step, False, "In progress...");
        written = 0;

        #outdir = "test-data/mm10/ensembl/cds/";
        #if not os.path.isdir(outdir):
        #    os.system("mkdir " + outdir);

        #outfilename = "/n/holylfs05/LABS/informatics/Users/gthomas/spiders/genomes/tgiga/tgiga-cds.fa";

        with open(globs['write-cds'], "w") as of:
            for seq in globs['cds-seqs']:
                #outfile = os.path.join(outdir, seq + ".fa");
                #with open(outfile, "w") as of:
                of.write(">" + seq + "\n");
                of.write(globs['cds-seqs'][seq] + "\n");
                written += 1;

                # if written == 100:
                #     break;

        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(written) + " sequences written");
    # Chunk of code to write out the sequences in a concatenated file to individual files by locus -- for development
    ###

    return globs;

#############################################################################

def readCDS(globs):
    
    step = "Reading CDS FASTA file(s)";
    step_start_time = CORE.report_step(globs, step, False, "In progress...", full_update=True);

    if globs['in-seq-type'] == "directory":
        seq_files = [ f for f in os.listdir(globs['in-seq']) if any(f.endswith(fasta_ext) for fasta_ext in [".fa", ".fa.gz", ".fasta", ".fasta.gz", ".fna", ".fna.gz"]) ];
        ## TODO: Make sure these are all the plausible extensions.
        ## TODO: Add extension lists to globs so they aren't all typed out here?
        ## NOTE: Do we even want to do this check?
    else:
        seq_files = [globs['in-seq']];
    # Get a list of files from the input
    # If the input is a directory, this will be all files in that directory
    # If the input is a file, this will just be a list with only that file in it

    if seq_files == []:
        CORE.errorOut("SEQ2", "No files in the input have extensions indicating they are FASTA files.", globs);
    # Makes sure some files have been read    

    for seq_file in seq_files:
        if globs['in-seq-type'] == "directory":
            seq_file_path = os.path.join(globs['in-seq'], seq_file);
        else:
            seq_file_path = seq_file;
        # For multiple input sequence files (directory), the path will be that directory and the current file

        cur_seqs = readFasta(seq_file_path, CORE.detectCompression(seq_file_path), False);
        # Read sequences in current file
        # Not sure if it is necessary to do the compression detections for each file... seems to take a while

        if len(cur_seqs) == 0:
            CORE.printWrite(globs['logfilename'], globs['log-v'], "# WARNING: file " + seq_file + " doesn't appear to be FASTA formatted... skipping");
            globs['warnings'] += 1;
            continue;
        # Check if we have actually read any sequence

        for seq in cur_seqs:
            # if len(cur_seqs[seq]) % 3 != 0:
            #     CORE.printWrite(globs['logfilename'], globs['log-v'], "# WARNING: sequence " + seq + " in file " + seq_file + " isn't in frame 1... skipping");
            #     globs['warnings'] += 1;
            #     continue;
            # Check that the current sequence is in frame 1

            globs['cds-seqs'][seq] = cur_seqs[seq].upper();
            # Add the current sequence to the global sequence dict

    if not globs['cds-seqs']:
       CORE.errorOut("SEQ3", "No FASTA sequences were read from input. Exiting.", globs); 
    # If no sequences were read from the input, error out

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(len(globs['cds-seqs'])) + " CDS read", full_update=True);
    # Status update

    return globs;

#############################################################################