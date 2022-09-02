#############################################################################
# Functions to handle vcf files
#############################################################################

import sys
import os
import random
from collections import defaultdict
from pysam import VariantFile
import lib.core as CORE

#############################################################################

def read(globs):
# Reads a VCF file into a VariantFile object with pysam and stores object in globs

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
    ## SANDBOX FOR FIGURING OUT PYSAM

    exclude_missing = [ sample for sample in globs['vcf-exclude'] if sample not in globs['vcf'].header.samples ];
    if exclude_missing:
        CORE.printWrite(globs['logfilename'], globs['log-v'], "# WARNING: some samples specified to be exclude (-e) were not found in the VCF file: " + ",".join(exclude_missing));
        globs['warnings'] += 1;
    # Throw a warning if some of the excluded samples aren't in the VCF just in case it was a typo or something
    ## Excludes

    outgroups_missing = [ sample for sample in globs['vcf-outgroups'] if sample not in globs['vcf'].header.samples ];
    if outgroups_missing:
        CORE.errorOut("VCF1", "Not all outgroups samples specified exist in provided VCF file. Missing outgroups: " + ",".join(outgroups_missing), globs);
    # A quick check to make sure all provided outgroup species are in the VCF file
    ## Outgroups

    globs['vcf-ingroups'] = [ sample for sample in globs['vcf'].header.samples if sample not in globs['vcf-outgroups'] + globs['vcf-exclude'] ];
    if not globs['vcf-ingroups']:
        CORE.errorOut("VCF2", "No ingroup samples in VCF file -- all samples were specified as either outgroups (-u) or to be excluded (-e).", globs);
    # Check to make sure there are ingroup samples left after removing the excluded samples.
    ## Ingroups

    return globs;

#############################################################################

def getVariants(globs, transcript, transcript_position, ref_codon):
# Gets variant codons in the in- and outgroups

    genome_region = globs['gxf'][transcript]['header'];
    genome_start_pos = globs['coords'][transcript][transcript_position]-1;
    #genome_positions = list(range(genome_start_pos, genome_start_pos+3));
    # Convert the transcript coordinates to genomic ones for the VCF
    # Pysam coordinates are 0 based, so subtract 1
    ## TODO: Make sure this is consistent with our other data structures!!

    #poly_codons = { sample : "" for sample in globs['vcf-ingroups'] };
    polymorphic_codons = [];
    fixed_diff_codon = "";
    # A list of polymorphic codons from the ingroups and a single codon string
    # for fixed differences in the outgroups

    #for rec in bcf_in.fetch('chr1', 100000, 200000):
    codon_pos = 0;
    for rec in globs['vcf'].fetch(genome_region, genome_start_pos, genome_start_pos+3):
    # Look up any variants for the current codon
    ## +2 or +3??

        if not list(rec):
            fixed_diff_codon += ref_codon[codon_pos];
        # If there is no record for this position, look-up the reference allele in the sequence
        # and add it to the fixed_diff_codon

        else:
            ref_nt = rec.ref;
            alt_nts = rec.alts;
            # Look up the alleles at the current position

            in_allele_counts = defaultdict(int);
            out_allele_counts = defaultdict(int);
            # Dicts for allele counts in each group

            ####################

            for sample in globs['vcf-ingroups']:
                for allele in rec.samples[sample]['GT']:
                    in_allele_counts[allele] += 1;
            # Count alleles in the ingroups

            for allele in in_allele_counts:
                if allele == 0:
                    continue;
                # Construct a codon for each allele in the ingroups except the
                # reference allele (since that would match the reference codon)

                polymorphic_codons = ref_codon;
                # TODO: Make sure this doesn't copy just a pointer or something
                polymorphic_codons[codon_pos] = alt_nts[allele-1];
                polymorphic_codons.append(polymorphic_codons);
                # Add the polymorphic codon to the list of polymorphic codons
            ## End ingroup allele loop

            ####################

            for sample in globs['vcf-outgroups']:
                for allele in rec.samples[sample]['GT']:
                    out_allele_counts[allele] += 1;
            # Count alleles in the outgroup

            if 0 in out_allele_counts:
                fixed_diff_codon += ref_codon[codon_pos];
            # If all outgroup alleles are reference, add the reference allele to
            # the fixed_diff_codon

            if all(alleles not in in_allele_counts for alleles in out_allele_counts):
            # Sites contain fixed differences only if all the alleles in the outgroup do not
            # exist in the ingroup 
                if len(out_allele_counts) == 1:
                    fixed_diffs += alt_nts[out_allele_counts.keys()[0]-1];
                # If there is only one allele in the outrgoup (the site is not polymorphic in the outgroup)
                # add that allele to the fixed_diff_codon by looking it up in the alt_nts list
                
                else:
                    max_allele = [ allele for allele, count in out_allele_counts.items() if count == max(out_allele_counts.values()) ];
                    max_allele = random.choice(max_allele);
                    fixed_diff_codon += alt_nts[max_allele-1];
                # If there is more than one alternate allele in the outgroup (the site is polymorphic in the outgroup)
                # use the allele at the highest frequency
                # If more than one allele exists at the highest frequency, pick one randomly
            ## End outgroup block

            ####################
        ## End codon record block

        codon_pos += 1;
    ## End codon record loop

    polymorphic_codons = [ "".join(poly_codon) for poly_codon in polymorphic_codons ];
    # Convert the polymorphic_codons from lists to strings

    return polymorphic_codons, fixed_diff_codon;
    # polymorphic_codons:   A list of codons that incorporate all SNPs in the ingroup samples, one codon
    #                       per alternate allele per site. As is, this list will never have the reference codon in it,
    #                       and could be an empty list if there are no SNPs in the ingroup samples.
    # fixed_diff_codon:     A single codon string that incorporates all fixed differences in the outgroup species
    #                       relative to the reference codon. A fixed difference only occurs if all the alleles in
    #                       the outgroup samples 1) are not the reference allele and 2) never occur in the ingroup 
    #                       samples. As currently implemented, if there are no fixed differences this will return
    #                       the reference codon.

#############################################################################