#############################################################################
# Functions to compute degeneracy and syn/nonsyn subs for CDS sequences in degenotate
#############################################################################

import os
import csv
import re
from collections import namedtuple
from itertools import pairwise
#use networkx to turn codon table into a graph to allow easy computation of paths
import networkx as nx

#############################################################################

def readDegen():
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

    DEGEN_DICT = {}
    CODON_DICT = {}
    with open(os.path.join(os.path.dirname(__file__), "codon_table.csv"), "r") as dd:
        reader = csv.reader(dd)
        DEGEN_DICT = {rows[0]:rows[1] for rows in reader}
        CODON_DICT = {rows[0]:rows[2] for rows in reader}

    #compute codon graph
    CODON_GRAPH = nx.Graph()
    #add a node for every codon
    CODON_GRAPH.add_nodes_from([CODON_DICT.keys])
    #add an edge between all codons that are 1 mutation apart
    for codon1 in CODON_DICT.keys:
        for codon2 in CODON_DICT.keys:
            ed = codonHamming(codon1,codon2)
            if ed == 1:
                CODON_GRAPH.add_edge(codon1,codon2)
            else:
                continue

    return [ DEGEN_DICT, CODON_DICT, CODON_GRAPH ]

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
def processCodons(globs):
# take CDS sequence and split into list of codons, computing degeneracy, ns, or both
# might need a clearer name?

    DEGEN_DICT, CODON_DICT, CODON_GRAPH = readDegen()
    MKTable = namedtuple("MKTable", "pn ps dn ds")

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

        #make list of codons
        #TO DO: check to make sure we get only triples out
        codons = re.findall('...', fasta)

        if ("degen" in globs['codon-methods']):
            degen = [DEGEN_DICT[x] for x in codons]
            globs['degeneracy'][transcript] = "." * (frame-1) if frame > 1 else ""
            globs['degeneracy'][transcript] = "".join(degen)

        if ("ns" in globs['codon-methods']):

            #define coordinate shift based on frame
            coord_shift = frame-1

            #process each codon
            for i,codon in enumerate(codons):
                #DOUBLE CHECK THIS PLEASE
                transcript_position = (i*3)+1+coord_shift
                ref_aa = CODON_DICT[codon]
                ps = 0.0
                pn = 0.0
                ds = 0.0
                dn = 0.0

                #assume getVariants returns a data structure of variant codons
                #this should be a list (empty, 1, or more) for in group codons
                #but for outgroup codons, it should be a single string with fixed differences
                poly_codons,div_codon = getVariants(globs,transcript,transcript_position)

                if poly_codons:
                    #there are variants
                    for poly_codon in poly_codons:

                        #for in group variants, we treat each as independent

                        poly_aa = CODON_DICT[poly_codon]
                        ps++ if poly_aa == ref_aa
                        pn++ if poly_aa != ref_aa

                if div_codon:
                    #there are fixed differences

                    #get number of differences between the codons
                    diffs = codonHamming(div_codon,codon)

                    if diffs == 1:
                        div_aa = CODON_DICT[div_codon]
                        ds++ if div_aa == ref_aa
                        dn++ if div_aa != ref_aa
                    if diffs >= 2:
                        ds,dn = codonPath(ref_aa,div_aa)

                globs['nonsyn'][transcript][i] = MKTable(pn,ps,dn,ds)

####################################
def getVariants(globs,transcript,transcript_position):

    #TO DO - function to return a list of variant codons based on vcf
    #should return two lists: one with all ingroup codons, the other with all outgroup codons
    #because outgroup is assumed to be only fixed differences, should only ever return a single
    #outgroup codon

def codonPath(start_codon,end_codon,CODON_GRAPH):

    #function to calculate syn/nonsyn for multi-step paths
    dn=0.0
    ds=0.0
    paths = nx.all_shortest_paths(CODON_GRAPH, source=start_codon, target=end_codon)
    numpaths = len(paths)
    for path in paths:
        pathlen=len(path)
        for pairs in pairwise(path):
            aa1 = CODON_DICT[pairs[0]]
            aa2 = CODON_DICT[pairs[1]]
            if aa1 == aa2:
                ds++
            else:
                dn++
    dn = dn/numpaths
    ds = ds/numpaths
    return ds,dn

def codonHamming(codon1,codon2):
    sum(1 for a, b in zip(codon1, codon2) if a != b)
