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

def readFasta(filename, globs):
# Read a FASTA formatted sequence file
# Great iterator and groupby code from: https://www.biostars.org/p/710/ 

    if globs['seq-compression'] == "gz":
        file_stream = gzip.open(filename);
        fa_iter = (x[1] for x in groupby(file_stream, lambda line: line.decode()[0] == ">"));
        readstr = lambda s : s.decode().strip();
    elif globs['seq-compression'] == "none":
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
        #curkey = header[1:header.index(" ")];
        # This removes the ">" character from the header string to act as the key in seqdict
        # TODO: Need to decide if splitting the input sequence header on " " is good for most cases, or if we should
        #       add this as a user option?

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
    globs['genome-seqs'] = readFasta(globs['fa-file'], globs);
    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(len(globs['genome-seqs'])) + " seqs read");
    # Read the input sequence file

    #print(list(globs['genome-seqs'].keys()))

    return globs;

#############################################################################

def extractCDS(globs):
# This takes the coordiantes read from the input annotation file as well as the sequence read from the
# input genome fasta file and extracts coding sequences for all transcripts while accounting for strand

    step = "Extracting CDS";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");

    complement = { 'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A', 'N' : 'N',
                   'a' : 't', 'c' : 'g', 'g' : 'c', 't' : 'a', 'n' : 'n'  };
    # A dictionary with base complements for use when reverse complementing sequences

    for transcript in globs['annotation']:

        cur_seq = "";
        # Initialize the sequence string for the current transcript. This will be added to the 'seqs' dict later

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
            CORE.errorOut("SEQ1", "Some exons have differing strands", globs);
        # Add check to make sure exons all have same strand as transcript?

        exon_coords = { exons[exon]['start']-1 : exons[exon]['end'] for exon in exons };
        # Get the coordinates of all the exons in this transcript
        # Subtract 1 from the starting coord because GXF are 1-based and python strings are 0-based

        if strand == "+":
            sorted_starts = sorted(list(exon_coords.keys()));
        elif strand == "-":
            sorted_starts = sorted(list(exon_coords.keys()), reverse=True);
        # Make sure the exons are sorted correctly, reversing the order if the strand is "-"

        for start in sorted_starts:
            cur_exon_seq = globs['genome-seqs'][header][start:exon_coords[start]];
            # For each exon starting coordinate, extract the sequence that corresponds to the current header and
            # start and end coordinates

            if strand == "-": 
                cur_exon_seq = "".join(complement.get(base, base) for base in reversed(cur_exon_seq));
            # If the strand is "-", get the reverse complement of the sequence

            cur_seq += cur_exon_seq;
            # Concatenate the current exon sequence onto the overall transcript sequence

        globs['cds-seqs'][transcript] = cur_seq;
        # Save the current transcript sequence to the global seqs dict

        # if strand == "-":
        #     print(transcript);
        #     print(globs['in-seqs'][transcript]);
        #     sys.exit();

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(len(globs['cds-seqs'])) + " CDS read");
    # Status update

    # ###
    # step = "Writing separate CDS sequences";
    # step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # written = 0;

    # outdir = "test-data/mm10/cds/";
    # if not os.path.isdir(outdir):
    #     os.system("mkdir " + outdir);

    # for seq in globs['cds-seqs']:
    #     outfile = os.path.join(outdir, seq + ".fa");
    #     with open(outfile, "w") as of:
    #         for header in globs['cds-seqs'][seq]:
    #             of.write(">" + seq + "\n");
    #             of.write(globs['cds-seqs'][seq] + "\n");
    #     written += 1;
    # step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(written) + " sequences written");
    # # Chunk of code to write out the sequences in a concatenated file to individual files by locus -- for development
    # ###

    return globs;

#############################################################################

def readCDS(globs):
    
    if globs['in-seq-type'] == "directory":
        seq_files = [ f for f in os.listdir(globs['in-seq']) if any(f.endswith(fasta_ext) for fasta_ext in [".fa", ".fa.gz", ".fasta", ".fasta.gz", ".fna", ".fna.gz"]) ];
        ## TODO: Make sure these are all the plausible extensions.
        ## TODO: Add extension lists to globs so they aren't all typed out here?
        ## TODO: Make sure seq_files isn't empty and error out if it is
    else:
        seq_files = globs['in-seq'];

    ## TODO: Loop to call readFasta() on each file and save each sequence within the file to globs['cds-seqs']
    ## TODO: Add checks to make sure each sequence is divisible by 3

#############################################################################