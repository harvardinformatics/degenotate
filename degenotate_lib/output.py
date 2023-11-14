#############################################################################
# Functions to handle output of degeneracy calcs
#############################################################################

import sys
import os
import textwrap

#############################################################################

def compileBedLine(globs, transcript, transcript_region, cds_coord, base, codon, codon_pos, aa, base_degen, cdict):
# A function to compile output for the main bed file

    outline = [];
    # Initialize an empty list to add output to

    if globs['gxf-file']:
        genome_coord = globs['coords'][transcript][cds_coord];
        outline += [transcript_region, genome_coord-1, genome_coord];
    # In case the input was a gxf file and a genome, the first three columns of output
    # reference genome coordinate which are retrieved here

    else:
        outline += [transcript, cds_coord, cds_coord+1];
    # If the input was a directory of CDS sequences, the first three columns of output
    # reference the CDS coordinates

    outline += [transcript + ":" + str(cds_coord), base_degen, base, aa];
    # The rest of the output except for the substitution strings

    subs = [];
    # Initialize an empty list to store substitutions of the current base that change the AA 
    # E.g. ["T:S"] means there is one mutation here (degeneracy is 2) from the reference base
    #       to T which changes the AA from the reference AA to Serine
    # This is initialized here so it is still added as a column even when degen is 4

    if base_degen not in [4, "."]:
    # If any substitutions of the current base change the current AA, annotate them here

        for new_base in globs['bases']:
            if new_base == base:
                continue;
            # Skip if the current base is the reference base

            new_codon = list(codon);
            new_codon[codon_pos] = new_base;
            new_codon = "".join(new_codon);

            try:
                new_aa = cdict[new_codon];
                # Replace the reference base with the new base at the current codon
                # position and look up the new AA
            except KeyError:
                print("Error: ", new_codon, base_degen, codon, transcript)
                continue

            if aa != new_aa:
                subs.append(new_base + ":" + new_aa);
            # If the new AA is different, annotate the substitution as outlined above
        # End new base loop
        ##########
    # End substitutions block
    ##########

    outline += [";".join(subs)];
    outline = [ str(col) for col in outline ];

    return outline;

#############################################################################

def writeBed(line_list, bed_stream, strand):
# Writes bet output per site in a transcript

    if strand == "-":
        line_list.reverse();
    # Reverse the order of the transcript if it is on the - strand to preserve
    # ascending ordering for bed file

    for line in line_list:
        bed_stream.write("\t".join(line) + "\n");

#############################################################################

def initializeTranscriptSummary(summary_stream):
# Writes the headers for the transcript summary file

    cols = ["transcript", "gene", "cds_length", "mrna_length", "is_longest", "f0", "f2", "f3", "f4"];
    summary_stream.write("\t".join(cols) + "\n");

#############################################################################

def writeTranscriptSummary(globs, transcript, sum_dict, summary_stream):
# Writes the summary output for a transcript including counts of sites per fold

    if globs['gxf-file']:
        geneid = globs['annotation'][transcript]['gene-id']
        cdslen = str(globs['annotation'][transcript]['cdslen'])
        mrnalen = str(globs['annotation'][transcript]['len'])
        longest = str(globs['annotation'][transcript]['longest']) 
    else:
        geneid = transcript
        cdslen = str(len(globs['cds-seqs'][transcript]))
        mrnalen = "NA"
        longest = "NA"

    outline = [ transcript, geneid, cdslen,  mrnalen, longest];
    outline += [ str(sum_dict[fold]) for fold in sum_dict ];
    summary_stream.write("\t".join(outline) + "\n");
    # Compile and write the transcript summary line to the transcript outfile 

#############################################################################

def writeSeq(header, seq, seq_stream, linelen=60):
# A function to write sequences in FASTA format when -x is specified

    seq_stream.write(header + "\n");
    seq_stream.write(textwrap.fill(seq, linelen) + "\n");


#############################################################################

def initializeMKFile(mkfilename):
# Opens the MK output file and writes the headers

    mkfile = open(mkfilename, "w");
    cols = ['transcript', 'pN', 'pS', 'dN', 'dS', 'mk.raw.p.value', 'mk.odds.ni', 'dos', 'imp_pval', 'imp_dos'];
    mkfile.write("\t".join(cols) + "\n");
    return mkfile;

#############################################################################

def writeMK(transcript, outdict, mk_stream):
# A function to write out MK tables for each transcript

    outline = [ transcript,
                str(outdict['pn']),
                str(outdict['ps']),
                str(outdict['dn']),
                str(outdict['ds']),
                str(outdict['mk.pval']),
                str(outdict['mk.odds.ni']),
                str(outdict['dos']),
                str(outdict['imp_pval']),
                str(outdict['imp_dos'])
                ];
    mk_stream.write("\t".join(outline) + "\n");
        
#############################################################################
