#############################################################################
# Functions to handle output of degeneracy calcs
#############################################################################

import sys
import os
import lib.core as CORE

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
            new_aa = cdict[new_codon];
            # Replace the reference base with the new base at the current codon
            # position and look up the new AA

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