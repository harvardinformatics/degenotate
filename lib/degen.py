#############################################################################
# Functions to compute degeneracy for CDS sequences in degenotate
#############################################################################

import os
import csv
import re

#############################################################################

def readDegen():
# Read a codon table file into a degeneracy dict
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

    DEGEN_DICT = {}
    with open(os.path.join(os.path.dirname(__file__), "codon_table.csv"), "r") as dd:
        reader = csv.reader(dd)
        DEGEN_DICT = {rows[0]:rows[1] for rows in reader}

    return DEGEN_DICT


#############################################################################

def frameError(seq,frame):

#Check that the sequence is the correct multiple of three for the starting frame

    seq_mod = len(seq) % 3
    if frame == 1 and seq_mod == 0:
        return False
    elif frame == 2 and seq_mod == 2:
        return False
    elif frame == 3 and seq_mod == 1:
        return False
    else:
        return True

#############################################################################

def calcDegen(globs):
# Take CDS sequences and return a sequence of degeneracy codes

    DEGEN_DICT = readDegen()

    for transcript in globs['cds-seqs']:
        #use dict.get() to return value or a default option if key doesn't exist
        #assumes that globs['annotation'][transcript]['start-frame'] won't exist
        #unless start frame was parsed from GFF
        if globs['gxf-file']:
            frame = globs['annotation'][transcript].get('start-frame', 1)
        else:
            frame = 1
            
        #check frame
        if frameError(globs['cds-seqs'][transcript],frame):
            continue
        ## TO DO: ADD ERROR REPORTING FOR FRAME ERRORS

        #if frame is not 1, need to skip the first frame-1 bases
        fasta = globs['cds-seqs'][transcript][frame:]

        #also add . (degeneracy cannot be called) to frame-1 beginning of degeneracy string
        globs['annotation'][transcript] = "." * (frame-1) if frame > 1 else ""

        #make list of codons
        #TO DO: check to make sure we get only triples out
        codons = re.findall('...', fasta)

        #use list comprehension to make new list of degenercy codes based on codons
        #TO DO: error checking
        degen = [DEGEN_DICT[x] for x in codons]

        globs['degeneracy'][transcript] = "".join(degen)

        return globs

#############################################################################