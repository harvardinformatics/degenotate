#############################################################################
# Functions to handle gtf and gff annotation files.
#############################################################################

import sys
import os
import gzip
import lib.core as CORE

#############################################################################

def checkIDs(l, info, id_list, step, globs):
# Each time we read feature info and look for IDs, this checks to make sure only one
# ID is found. Probably unnecessary, but easy enough to check.

    if len(id_list) != 1:
        print("\n\n");
        print(id_list);
        print(info);
        print(l);
        print("\n\n");
        CORE.errorOut("GXF1", "Invalid number of IDs found during " + step, globs);

#############################################################################

def readFeatures(globs, file_reader, line_reader, feature_list, id_format, parent_id_format, info_field_splitter):

    num_features = 0;
    for line in file_reader(globs['gxf-file']):
            line = line_reader(line);
            # Read and parse the current line of the GXF file

            if "##FASTA" in line[0]:
                break;
            # Maker GFF files sometimes include the sequence of all transcripts at the end. We need to stop reading the file
            # at that point.

            if line[0][0] == "#":
                continue;
            # Header/comment lines should be skipped. Note that this must come after the check for "##FASTA" above, or else
            # the file will keep being read into the sequences and error out.

            feature_type, seq_header, start, end, strand, feature_info = line[2], line[0], int(line[3]), int(line[4]), line[6], line[8].split(info_field_splitter);
            # Unpack the pertinent information from the current line into more readable variables.

            if feature_type in feature_list:
            # Skipping any 'unconfirmed_transcript'
                feature_id = [ info_field for info_field in feature_info if info_field.startswith(id_format) ];
                # Get the feature ID as a list of fields with the "ID=" prefix

                # print(feature_info);
                # print(feature_id);
                # sys.exit();

                checkIDs(line, feature_info, feature_id, feature_list[0] +" id parsing", globs);
                # A quick check to make sure we have read only one ID

                feature_id = feature_id[0].replace(id_format, "").replace("\"", "");
                # Unpack and parse the ID

                parent_id = [ info_field for info_field in feature_info if info_field.startswith(parent_id_format) ];
                # Get the gene ID associated with the transcript as a list of fields with the "Parent=" prefix

                checkIDs(line, feature_info, parent_id, feature_list[0] +" parent id parsing", globs);
                # A quick check to make sure we have read only one ID

                parent_id = parent_id[0].replace(parent_id_format, "").replace("\"", "");
                # Unpack and parse the gene ID

                if feature_list[0] == "transcript":
                    globs['annotation'][feature_id] = { 'header' : seq_header, 'start' : start, 'end' : end, 'strand' : strand, 'exons' : {}, "gene-id" : parent_id };
                    # Add the ID and related info to the annotation dict. This includes an empty dict for exons to be stored in a similar way

                elif feature_list[0] == "CDS":
                    num_exons = len(globs['annotation'][parent_id]['exons']);
                    exon_id = "exon-" + str(num_exons+1);
                    # Because exon IDs are not always included for CDS, or they only represent the CDS as a whole (e.g. protein ID from Ensembl), we 
                    # count the number of exons in the transcript as the ID

                    globs['annotation'][parent_id]['exons'][exon_id] = { 'header' : seq_header, 'start' : start, 'end' : end, 'strand' : strand };
                # Add the ID and related info to the annotation dict.                   

                num_features += 1;
                # Add to the number of transcripts read


    return globs['annotation'], num_features;

#############################################################################

def read(globs):

    step = "Detecting compression of annotation file";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    globs['gxf-compression'] = CORE.detectCompression(globs['gxf-file']);
    if globs['gxf-compression'] == "none":
        reader = open;
        readline = lambda l : l.strip().split("\t");
        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: No compression detected");
    else:
        reader = gzip.open;
        readline = lambda l : l.decode().strip().split("\t");
        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + globs['gxf-compression'] + " detected");
    # Detect the compression of the input annotation file

    if globs['gxf-type'] == "gff":
        field_splitter = ";";
        gene_id_format = "ID=";
        transcript_id_format = "ID=";
        exon_id_format = "ID=";
        transcript_parent_format = "Parent=";
        exon_parent_format = "Parent=";

    elif globs['gxf-type'] == "gtf":
        field_splitter = "; ";
        gene_id_format = "gene_id ";
        transcript_id_format = "transcript_id ";
        exon_id_format = "exon_id ";
        transcript_parent_format = "gene_id ";
        exon_parent_format = "transcript_id ";
    # These outline the differences between GFF and GTF

    globs['annotation'] = {};
    # The main annotation storage dict. A nested structure for genes, transcripts, and coding exons.
    # <transcript id> : { <header>, <start coord>, <end coord>, <strand>, 
    #                       { <exon id> : <exon header>, <exon start coord>, <exon end coord>, <exon strand> } } }

    ####################

    step = "Reading transcripts";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    globs['annotation'], num_transcripts = readFeatures(globs, reader, readline, ["transcript", "mRNA"], transcript_id_format, transcript_parent_format, field_splitter);
    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(num_transcripts) + " transcripts read");

    # Read transcripts
    ####################

    step = "Reading coding exons";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    globs['annotation'], num_cds_exons = readFeatures(globs, reader, readline, ["CDS"], exon_id_format, exon_parent_format, field_splitter);
    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(num_cds_exons) + " coding exons read");

    # Read coding exons
    ####################

    if num_cds_exons == 0:
        CORE.errorOut("GXF2", "No CDS exons found in input annotation file! Cannot calculate degeneracy without coding sequences.", globs);
    # Check to make sure at least one CDS sequence is found, otherwise error out

    return globs;

    #############################################################################
    





    #############################################################################
    
    # step = "Reading genes";
    # step_start_time = CORE.report_step(globs, step, False, "In progress...");

    # for line in reader(globs['gxf-file']):
    #     line = readline(line);
    #     # Read and parse the current line of the GXF file

    #     if "##FASTA" in line[0]:
    #         break;
    #     # Maker GFF files sometimes include the sequence of all transcripts at the end. We need to stop reading the file
    #     # at that point

    #     if line[0][0] == "#":
    #         continue;
    #     # Header/comment lines should be skipped. Note that this must come after the check for "##FASTA" above, or else
    #     # the file will keep being read into the sequences and error out

    #     feature_type, seq_header, start, end, strand, feature_info = line[2], line[0], int(line[3]), int(line[4]), line[6], line[8].split(field_splitter);
    #     # Unpack the pertinent information from the current line into more readable variables

    #     if feature_type == "gene":
    #         feature_id = [ info_field for info_field in feature_info if info_field.startswith(gene_id_format) ];
    #         # Get the feature ID as a list of fields with the "ID=" prefix

    #         checkIDs(line, feature_info, feature_id, "gene id parsing", globs);
    #         # A quick check to make sure we have read only one ID

    #         feature_id = feature_id[0].replace(gene_id_format, "").replace("\"", "");
    #         # Unpack and parse the ID

    #         globs['annotation'][feature_id] = { 'header' : seq_header, 'start' : start, 'end' : end, 'strand' : strand, 'transcripts' : {} };
    #         # Add the ID and related info to the annotation dict. This includes an empty dict for transcripts to be stored in a similar way

    # step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(len(globs['annotation'])) + " genes read");
    # # Status update

    # # Old code to read genes
    # ####################