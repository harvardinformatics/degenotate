#############################################################################
# Functions to compute degeneracy and syn/nonsyn subs for CDS sequences in degenotate
#############################################################################

import sys
import os
import csv
import re
import itertools
from collections import namedtuple
import lib.vcf as VCF
import lib.output as OUT
import lib.core as CORE

#############################################################################

def readDegen(globs):
# Read a codon table file into a degeneracy dict and a codon dict
# Assumes a plain text, comma-separated file with two columns
# The first column is the three lettter (DNA) codon sequence
# The second column is a three digit code for the degeneracy of the first, second, and third positions
# 0 = non-degenerate; any mutation will change the amino acid
# 2 = two nucleotides at the position code the same AA, so 1 of the three possible
#     mutations will be synonymous and 2 will be non-synonymous
# 3 = three nucleotides at the position code for the same AA, so 2 of the three possible
#     mutations will be synonymous and 1 will be non-synonymous
# 4 = four nucleotides at the position code for the same AA, so all 3 possible
#     mutations are synonymous
# The third column is the one letter AA code for that codon

    if "ns" in globs['codon-methods']:
        try:
            import networkx as nx
            # use networkx to turn codon table into a graph to allow easy computation of paths

            globs['shortest-paths'] = nx.all_shortest_paths;
        except:
            CORE.errorOut("DEGEN1", "Missing networkx dependency. Please install and try again: https://anaconda.org/conda-forge/networkx", globs);
    # Check for the networkx module and save the all_shortest_paths function to globs if present
    # Error out if not found
    # This is done here so there are no dependencies if the user doesn't want to generate the MK tables

    DEGEN_DICT = {}
    CODON_DICT = {}
    with open(os.path.join(os.path.dirname(__file__), "codon-table.csv"), "r") as dd:
        reader = csv.reader(dd)
        for row in reader:
            DEGEN_DICT[row[0]] = row[2];
            CODON_DICT[row[0]] = row[1];

    if "ns" in globs['codon-methods']:
        # compute codon graph
        CODON_GRAPH = nx.Graph()

        # add a node for every codon
        CODON_GRAPH.add_nodes_from(list(CODON_DICT.keys()))
    else:
        CODON_GRAPH = False;

    # add an edge between all codons that are 1 mutation apart
    for codon1 in CODON_DICT.keys():
        for codon2 in CODON_DICT.keys():
            ed = codonHamming(codon1,codon2)
            if ed == 1 and "ns" in globs['codon-methods']:
                CODON_GRAPH.add_edge(codon1,codon2)
            else:
                continue

    return [ DEGEN_DICT, CODON_DICT, CODON_GRAPH, globs ]

#############################################################################

def getFrame(seq):
# A function that returns the frame of a coding sequence

    seq_mod = len(seq) % 3;
    if seq_mod == 0:
        return 0;
    elif seq_mod == 2:
        return 1;
    elif seq_mod == 1:
        return 2;

#############################################################################

def frameError(seq,frame):
# Check that the sequence is the correct multiple of three for the starting frame

    seq_mod = len(seq) % 3
    if frame == 0 and seq_mod == 0:
        return False
    elif frame == 1 and seq_mod == 2:
        return False
    elif frame == 2 and seq_mod == 1:
        return False
    else:
        return True

#############################################################################

def codonPath(start_codon, end_codon, CODON_GRAPH, CODON_DICT, nx_shortest_paths):

    #function to calculate syn/nonsyn for multi-step paths
    #by default returns the average nonsyn and syn subs over all shortest paths b/w two codons
    dn=0.0
    ds=0.0

    #get all shortest paths using networkx functions
    paths = nx_shortest_paths(CODON_GRAPH, source=start_codon, target=end_codon)

    #number of possible shortest paths
    numpaths = 0

    #for each path, calculate the number of syn and nonsyn subs implied
    for path in paths:
        for pairs in itertools.pairwise(path):
            aa1 = CODON_DICT[pairs[0]]
            aa2 = CODON_DICT[pairs[1]]
            if aa1 == aa2:
                ds+=1
            if aa1 != aa2:
                dn+=1
        numpaths += 1

    #calculate average over all paths
    dn = dn/numpaths
    ds = ds/numpaths

    #return values
    return ds,dn

#############################################################################

def codonHamming(codon1,codon2):
# Calculates hamming distance between two strings (codons) 
    return sum(1 for a, b in zip(codon1, codon2) if a != b)

#############################################################################

def processCodons(globs):
# take CDS sequence and split into list of codons, computing degeneracy, ns, or both
# might need a clearer name?

    DEGEN_DICT, CODON_DICT, CODON_GRAPH, globs = readDegen(globs)
    #MKTable = namedtuple("MKTable", "pn ps dn ds")

    ####################

    num_transcripts = len(globs['cds-seqs']);

    step = "Caclulating degeneracy per transcript";
    step_start_time = CORE.report_step(globs, step, False, "Processed 0 / " + str(num_transcripts) + " transcripts...", full_update=True);
    # Status update

    ####################

    with open(globs['outbed'], "w") as bedfile, open(globs['out-transcript'], "w") as transcriptfile:
        
        OUT.initializeTranscriptSummary(transcriptfile);
        # Write the column headers to the transcript summary file

        if globs['outseq']:
            seq_stream = open(globs['outseq'], "w");
        # Open the sequence file if necessary

        if "ns" in globs['codon-methods']:
            mk_stream = OUT.initializeMKFile(globs['outmk']);
        # Open the MK file if necessary

        counter = 0;
        for transcript in globs['cds-seqs']:

            transcript_output = { 'bed' : [], 
                                  'mk' : { 'pn' : 0, 'ps' : 0, 'dn' : 0, 'ds' : 0 },
                                  'summary' : { 0 : 0, 2 : 0, 3 : 0, 4 : 0 },
                                  'seq' : "" };
            # The output lines for each transcript

            if globs['outseq']:
                transcript_output['header'] = ">" + transcript + " " + ",".join([ str(f) for f in globs['extract-fold'] ]) + "-fold degenerate sites";
            # When outputting sequences by different folds, construct the header here

            if globs['gxf-file']:
                transcript_region = globs['annotation'][transcript]['header'];
            else:
                transcript_region = transcript;
            # Get the genome region if the input was a gxf file+genome

            if globs['gxf-file']:
                frame = globs['annotation'][transcript]['start-frame']

                if frame is None:
                    CORE.printWrite(globs['logfilename'], 3, "# WARNING: transcript " + transcript + " has an unknown frame....skipping");
                    globs['warnings'] += 1;                    
                    continue;
            # Get the frame when input is a when input is a gxf+genome
            else:
                frame = getFrame(globs['cds-seqs'][transcript]);
                if frame != 0:
                    CORE.printWrite(globs['logfilename'], 3, "# WARNING: transcript " + transcript + " is partial with unknown frame....skipping");
                    globs['warnings'] += 1;                    
                    continue;
            # Get the frame when input is a dir/file of individual CDS seqs
            # In this case we just check to make sure the sequence is a multiple of 3

            extra_leading_nt = frame
            # Look up the number of leading bases given the current frame

            #if frame is not 1, need to skip the first frame-1 bases
            fasta = globs['cds-seqs'][transcript][extra_leading_nt:]

            #now check to see if there are still trailing bases
            extra_trailing_nt = len(fasta) % 3

            if extra_trailing_nt > 0:
                fasta = fasta[:-extra_trailing_nt]
 
            #make list of codons
            codons = re.findall('...', fasta)
            nt = {'A', 'T', 'G', 'C'}

            if ("degen" in globs['codon-methods']):
                degen = [ DEGEN_DICT[x] if len(set(x) - nt) == 0 else "..." for x in codons ];
                # Get the string of degeneracy integers for every codon in the current sequence (e.g. 002)

                degen = "." * extra_leading_nt + "".join(degen);
                # Convert the degeneracy string to a list, and add on dots for any leading bases that
                # were removed if the frame is not 1

                cds_coord = 0;
                # Start the CDS coord counter

                if frame != 0:
                    for out_of_frame_pos in range(extra_leading_nt):
                        outline = OUT.compileBedLine(globs, transcript, transcript_region, cds_coord, globs['cds-seqs'][transcript][cds_coord], "", "", ".", ".", "");
                        transcript_output['bed'].append(outline);
                        # Call the output function with blank values since there is no degeneracy at this position
                        # and store the line in the output dict 

                        cds_coord += 1;
                        # Increment the position in the CDS
                # If the CDS is not in frame 1, the bed output needs to be filled in for the leading bases that were removed
                ##########

                for codon in codons:
                    try:
                        aa = CODON_DICT[codon];
                    except KeyError:
                        aa = "."
                    # Look up the AA of the current codon

                    for codon_pos in [0,1,2]:
                        base = codon[codon_pos];
                        # Extract the current base from the codon string
                    
                        outline = OUT.compileBedLine(globs, transcript, transcript_region, cds_coord, base, codon, codon_pos, aa, degen[cds_coord], CODON_DICT);
                        transcript_output['bed'].append(outline);
                        # Store the output from the current position in the output dict 

                        if degen[cds_coord] != ".":
                            transcript_output['summary'][int(degen[cds_coord])] += 1;
                        # Increment the count for the current degeneracy for the transcript summary
                        # Skip positions with unknown degeneracy ('.')

                        if globs['outseq'] and degen[cds_coord] in globs['extract-fold']:
                            transcript_output['seq'] += base;
                        # Add the current base to the sequence string if extraction is specified and it is in 
                        # a fold requested for extraction

                        cds_coord += 1;
                        # Increment the position in the CDS
                    # End base loop
                    ##########
                # End codon loop
                ##########

                if extra_trailing_nt != 0:
                    for out_of_frame_pos in range(extra_trailing_nt):
                        outline = OUT.compileBedLine(globs, transcript, transcript_region, cds_coord, globs['cds-seqs'][transcript][cds_coord], "", "", ".", ".", "");
                        transcript_output['bed'].append(outline);
                        # Call the output function with blank values since there is no degeneracy at this position
                        # and store the line in the output dict 

                        cds_coord += 1;
                        # Increment the position in the CDS
                # If the CDS has extra trailing bases, the bed output needs to be filled in for the leading bases that were removed
                ##########

            ## Runtime for test chromosome without output:              6 sec
            ## Runtime for test chromosome with output without subs:    20 sec
            ## Runtime for test chromosome with output with subs:       33 sec

            ## Out of frame test seq when using -s test-data/mm10/ensembl/cds/ as input: transcript:ENSMUST00000237320

            # End degen method block
            ####################

            if "ns" in globs['codon-methods']:

                #define coordinate shift based on frame
                #transcript_position = extra_leading_nt;

                mk_codons = VCF.getVariants(globs, transcript, transcript_region, codons, extra_leading_nt, extra_trailing_nt)
                # Call get variants for this transcript: returns a dictionary with the key being the index of each codon in codons with values as follows:
                # 'poly' :       A list of codons that incorporate all SNPs in the ingroup samples, one codon
                #                per alternate allele per site. As is, this list will never have the reference codon in it,
                #                and could be an empty list if there are no SNPs in the ingroup samples.
                # 'fixed' :      A single codon string that incorporates all fixed differences in the outgroup species
                #                relative to the reference codon. A fixed difference only occurs if all the alleles in
                #                the outgroup samples 1) are not the reference allele and 2) never occur in the ingroup 
                #                samples. As currently implemented, if there are no fixed differences this will return
                #                the reference codon. If no fixed differences are present, this is just the reference
                #                codon.
                # 'fixed-flag' : A boolean that is True if fixed differences have been found and False if not.

                for codon_index in range(len(codons)):
                # Loop over each codon by index

                    codon = codons[codon_index];
                    mk_alleles = mk_codons[codon_index];
                    # Look up the codon and the results for the codon from getVariants

                    try: 
                        ref_aa = CODON_DICT[codon]
                    except KeyError:
                        continue;
                    # Look up the original amino acid to compare variant codons against

                    pn, ps, dn, ds = 0.0, 0.0, 0.0, 0.0;
                    # Initialize site counts

                    if mk_alleles['poly']:
                    # If there are polymorphisms

                        for poly_codon in mk_alleles['poly']:
                        # For in group variants, we treat each SNP as independent

                            try:
                                poly_aa = CODON_DICT[poly_codon]
                            except KeyError:
                                continue;
                            # Look up the amino acid of the polymorphic codon
                            
                            if poly_aa == ref_aa:
                                ps += 1;
                            if poly_aa != ref_aa:
                                pn += 1;
                            # If the SNP doesn't change the AA from the reference, increment ps, otherwise pn

                            print(transcript, codon_index, codon, poly_codon, ref_aa, poly_aa, sep=":")

                        # End polymorphic codon loop
                        ##########
                    # End polymorphism block
                    ##########

                    if mk_alleles['fixed-flag']:
                    # If there are fixed differences

                        diffs = codonHamming(mk_alleles['fixed'], codon);
                        # Get number of differences between the outgroup codon and reference codon

                        if diffs == 1:
                        # If there is only one difference between the outgroup codon and the reference codon, compare the AA's directly

                            try:
                                div_aa = CODON_DICT[mk_alleles['fixed']]
                            except KeyError:
                                continue;
                            # Look up the amino acid of the outgroup codon
                                
                            if div_aa == ref_aa:
                                ds += 1;
                            if div_aa != ref_aa:
                                dn += 1;
                            # If the SNP doesn't change the AA from the reference, increment ds, otherwise dn

                        if diffs >= 2:
                            ds, dn = codonPath(codon, mk_alleles['fixed'], CODON_GRAPH, CODON_DICT, globs['shortest-paths']);
                        # If there is more than one difference between the outgroup and reference codon, find the order of the SNPs
                        # to compare
                    # End fixed diff block
                    ##########

                    transcript_output['mk']['pn'] += pn;
                    transcript_output['mk']['ps'] += ps;
                    transcript_output['mk']['dn'] += dn;
                    transcript_output['mk']['ds'] += ds;
                    # Increment the counts for each site type for this transcript

                    # try:
                    #     globs['nonsyn'][transcript][transcript_position] = MKTable(pn,ps,dn,ds)
                    # except KeyError:
                    #     globs['nonsyn'].update({transcript: {transcript_position : MKTable(pn,ps,dn,ds)}})
                    # NOTE GT: do we need to add placeholders for the extra leading bases to the nonsyn dict?
                    # e.g. globs['nonsyn'][transcript] could be a list with the index being the position... not
                    # sure what is easiest here.
                    #transcript_position += 3
                # End codon loop
                ##########

            # End ns method block
            ##########

            OUT.writeBed(transcript_output['bed'], bedfile, globs['annotation'][transcript]['strand']);
            # Write the bed output for every site in this transcript

            OUT.writeTranscriptSummary(globs, transcript, transcript_output['summary'], transcriptfile);
            # Write the summary for this transcript

            if globs['outseq']:
                OUT.writeSeq(transcript_output['header'], transcript_output['seq'], seq_stream);

            if "ns" in globs['codon-methods']:
                OUT.writeMK(transcript, transcript_output['mk'], mk_stream);
            # Write the MK table for this transcript

            counter += 1;
            if counter % 100 == 0:
                cur_step_time = CORE.report_step(globs, step, step_start_time, "Processed " + str(counter) + " / " + str(num_transcripts) + " transcripts...", full_update=True);
            # A counter and a status update every 100 loci

        # End transcript loop
        ##########

    # Close bed file
    ##########

    if "ns" in globs['codon-methods']:
        mk_stream.close();
    # Close the MK file if necessary

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success", full_update=True);
    # Status update

    return globs

#############################################################################