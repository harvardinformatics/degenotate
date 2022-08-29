#############################################################################
# Functions to handle vcf files
#############################################################################

import sys
import os
import gzip
from pysam import VariantFile
import lib.core as CORE

#############################################################################

def read(globs):
# Reads a VCF file into a VariantFile object with pysam and stores object in globs
    print();
    globs['vcf'] = VariantFile(globs['vcf-file']);
    # Read the VCF

    # x = globs['vcf'].fetch("NC_046340.1", 71, 73);

    # print(x);
    # print(list(x));

    # for i in [72,73,74]:
    #     print(i);
    #     f = globs['vcf'].fetch("NC_046340.1", i, i+1);

    #     #print(len(f));

    #     if not list(f):
    #         print("HAS")

    #     for rec in f:

    #         if rec.ref:
    #             print(i);
    #             print(rec);
    #             print(rec.ref);
    #         else:
    #             print("AAHJ")
    #         #print(rec.info);

    # sys.exit();
    # for record in globs['vcf'].fetch():
    #     print(record.ref);
    #     print(record.alts);
    #     print(record.alts[0]);

    #     for sample in record.samples:
    #         print(sample);
    #         print(record.samples[sample]['GT']);
    #     sys.exit();
    ## TEST BLOCK FOR READING VCF FILES

    outgroups_missing = [ sample for sample in globs['vcf-outgroups'] if sample not in globs['vcf'].header.samples ];
    if outgroups_missing:
        CORE.errorOut("VCF1", "Not all outgroups samples specified exist in provided VCF file. Missing outgroups: " + ",".join(outgroups_missing), globs);
    # A quick check to make sure all provided outgroup species are in the VCF file

    globs['vcf-ingroups'] = [ sample for sample in globs['vcf'].header.samples if sample not in globs['vcf-outgroups'] ];

    return globs;

#############################################################################

def getVariants(globs, transcript, transcript_position):
# Gets variant codons in the in- and outgroups

    genome_region = globs['gxf'][transcript]['header'];
    genome_start_pos = globs['coords'][transcript][transcript_position]-1;
    #genome_positions = list(range(genome_start_pos, genome_start_pos+3));
    # Convert the transcript coordinates to genomic ones for the VCF
    # Pysam coordinates are 0 based, so subtract 1
    ## TODO: Make sure this is consistent with our other data structures!!

    poly_codons = { sample : "" for sample in globs['vcf-ingroups'] };
    fixed_diffs = "";
    # A list of polymorphic codons from the ingroups and a single codon string
    # for fixed differences in the outgroups

    #for rec in bcf_in.fetch('chr1', 100000, 200000):
    for rec in globs['vcf'].fetch(genome_region, genome_start_pos, genome_start_pos+3):
    # Look up any variants for the current codon
    ## +2 or +3??

        if not list(rec):
            ref_nt = globs['cds-seqs'][transcript][transcript_position];
            poly_codons = { sample : poly_codons[sample] + ref_nt for sample in poly_codons };
            fixed_diffs += ref_nt;
        # If there is no record for this position, look-up the reference allele in the sequence
        # and add it to all codon strings (in and outgroup)

        else:
            ref_nt = rec.ref;
            alt_nts = rec.alts;
            # Look up the alleles at the current position
            
            for sample in globs['vcf-ingroups']:
                if rec.samples[sample]['GT'] == (0,0):
                    poly_codons[sample] += ref_nt;
                elif rec.samples[sample]['GT'] == (0,1):
                    poly_codons[sample] += alt_nts[0];
                elif rec.samples[sample]['GT'] == (1,1):
                    poly_codons[sample] += alt_nts[0];
            # Loop over all the ingroup samples and lookup their genotypes

            if all(1 in rec.samples[sample]['GT'] for sample in globs['vcf-outgroups']):
                fixed_diffs += alt_nts[0];
            else:
                fixed_diffs += ref_nt;
            # Check if the outgroups contain a fixed difference
        # If there is a variant at this position, look up the genotypes for each sample and
        # add the appropriate nucleotide to each codon.
        ## NOTE: Currently assuming only 1 alt allele.
    ## End codon record loop


    # if len(set(poly_codons.values())) == 1:
    # else:
    ## TODO: Need to do some checks at the end here to return exactly what we want. Currently
    ## returns all codons in all samples -- ref codon if there are no variants.    


    return poly_codons, fixed_diffs;

#############################################################################