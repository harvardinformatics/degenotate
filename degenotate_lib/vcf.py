#############################################################################
# Functions to handle vcf files
#############################################################################

import math
import os
import random
import sys
from collections import defaultdict
import degenotate_lib.core as CORE

#############################################################################

def read(globs):
# Reads a VCF file into a VariantFile object with pysam and stores object in globs

    try:
        from pysam import VariantFile
    except:
        print();
        CORE.errorOut("VCF1", "Missing pysam dependency. Please install and try again: https://anaconda.org/bioconda/pysam", globs);

    globs['vcf'] = VariantFile(globs['vcf-file']);
    # Read the VCF

    for header in globs['vcf'].header.contigs:
        if header not in globs['genome-seqs']:
            CORE.printWrite(globs['logfilename'], 3, "# WARNING: " + header + " is present in the VCF file but not the genome fasta file.");
            globs['warnings'] += 1;  
    # Add warnings for each header in the VCF file that isn't in the input fasta file

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
    globs['num-ingroup-chr'] = len(globs['vcf-ingroups']) * 2;
    if not globs['vcf-ingroups']:
        CORE.errorOut("VCF2", "No ingroup samples in VCF file -- all samples were specified as either outgroups (-u) or to be excluded (-e).", globs);
    # Check to make sure there are ingroup samples left after removing the excluded samples.

    if not globs['ingroup-maf-cutoff']:
        globs['ingroup-maf-cutoff'] = 1 / globs['num-ingroup-chr'];
    # Calculate the default MAF cutoff here if none is provided (1 / N)
    ## Ingroups

    return globs;

#############################################################################

def getVariants(globs, transcript, transcript_region, codons, extra_leading_nt, extra_trailing_nt):
# Gets variant codons in the in- and outgroups

    mk_codons = { c : { 'poly' : [], 'fixed' : list(codons[c]), 'fixed-flag' : False, 'AF' : [] } for c in range(len(codons)) };
    #mk_codons = { c : { 'ingroup-poly-samples' : [], 'ingroup-ref-samples' : [], 'poly' : defaultdict(int), 'fixed' : list(codons[c]), 'fixed-flag' : False } for c in range(len(codons)) };
    #mk_codons = { c : { 'ref' : { codons[c] : 0 }, 'poly' : {}, 'fixed' : list(codons[c]), 'fixed-flag' : False } for c in range(len(codons)) };
    # This dictionary will keep track of polymorphic codons and their frequencies (initialized as empty dict)
    # and fixed difference codons (initialized as reference codon) for each codon in the current transcript
    # The fixed flag is swapped to True if any fixed differences have actually been found, else
    # the codon is still the ref and shouldn't be counted
    # Also keeps track of which samples in the ingroup have at least 1 polymorphism in the codon (ingroup-poly-samples) and samples that are reference for each position in the codon (ingroup-ref-samples)

    strand = globs['annotation'][transcript]['strand'];
    # Need to get strand, to convert leading/trailing bases into start/end padding

    if strand == "+":
        start_pad = extra_leading_nt;
        end_pad = extra_trailing_nt;
    
    if strand == "-":
        start_pad = extra_trailing_nt;
        end_pad = extra_leading_nt;
    # Extra leading nt and trailing nt are defined (leading, trailing) based on transcript orientation
    # So in genomic coordinates, trailing nt on the minus strand is actually a shift of the start, and leading nt is a shift of the end
    # This is because for the transcript feature, start < end
    
    adj_ts_start = globs['annotation'][transcript]['start'] + start_pad;
    adj_ts_end = globs['annotation'][transcript]['end'] - end_pad;
    # Adjust the genomic start and end coordinates for this transcript based on 
    # the extra out of frame nts in the transcript

    transcript_records = globs['vcf'].fetch(transcript_region, adj_ts_start, adj_ts_end);
    # Look up all variant records in the range of the adjusted transcript start and end coordinates

    for rec in transcript_records:

        alt_nts = rec.alts;
        if not alt_nts:
            continue;
        # Look up the alleles at the current position and if there are no alternate alleles (invariant site), skip

        if not all(len(alt) == 1 for alt in alt_nts):
            CORE.printWrite(globs['logfilename'], 3, "# WARNING: invalid allele in VCF file at position " + str(rec.start + 1) + "....skipping");
            globs['warnings'] += 1;
            continue;
        # A warning for alternate alleles that are longer than one base (indels, svs)

        if strand == "-":
            alt_nts = [ globs['complement'][base] for base in alt_nts ];
        # For transcripts on the negative strand, the alternate alleles in the VCF
        # need to be complemented

        rec_pos = rec.start + 1
        # Adjust 0-based pysam coordinate to 1-based gff coordinate here

        if rec_pos not in globs['coords-rev'][transcript]:
            continue;
        # Skip any SNPs within the range of the transcript start and end,
        # but not in the CDS

        rec_transcript_pos = globs['coords-rev'][transcript][rec_pos];
        # Look up the position of the record relative to the start of
        # the transcript

        adj_rec_pos = rec_transcript_pos - start_pad;
        # Adjust the transcript position based on the number of extra leading nts

        if adj_rec_pos < 1:
            continue;
        #if the variant is in the partial starting codon, skip it

        rec_codon_pos = math.floor(adj_rec_pos / 3);
        # Look up the codon position of the adjusted record position
        # This is also the key for mk_codons

        codon_pos = (adj_rec_pos) % 3;
        # The position of the record within the codon, either 0, 1, or 2

        try:
            ref_codon = codons[rec_codon_pos];
        except IndexError:
            continue;
            # an IndexError here should mean that we are in a last partial codon and we should just ignore this variant
            # however we should probably add some code to formally confirm this before just skipping
        # Look up the codon at the record's codon position

        in_allele_counts = defaultdict(int);
        in_hom_alts = defaultdict(int);
        out_allele_counts = defaultdict(int);
        # Dicts for allele counts in each group

        ####################

        for sample in globs['vcf-ingroups']:

            for allele in rec.samples[sample]['GT']:
                if allele is None:
                    continue
                # Skip missing data
                in_allele_counts[allele] += 1;
                # Count alleles in the ingroups
                if (rec.samples[sample]['GT'][0] == rec.samples[sample]['GT'][1] and rec.samples[sample]['GT'] != (0,0)):
                    in_hom_alts[rec.samples[sample]['GT']] += 1;
                # Check if the sample is homozygous for an alternate allele and add that to the count of in_hom_alts here

        ## End ingroup allele counting block

        ## NOTE: Should we be comparing to ingroups with called genotypes here instead of number of all ingroups?

        if not globs['count-fixed-alt-ingroups'] and len(in_hom_alts) == 1 and sum(in_hom_alts.values()) == globs['num-ingroups']:
            pass;
        # If we don't want to count fixed ingroups, check to see if any alt alleles are fixed and skip if so
        else:
            AN = sum([in_allele_counts[allele] for allele in in_allele_counts]) # total allele number without missing data
            for allele in in_allele_counts:
                if allele == 0:
                    continue;
                # Construct a codon for each allele in the ingroups except the
                # reference allele (since that would match the reference codon)

                # if in_allele_counts[allele] / globs['num-ingroup-chr'] < globs['ingroup-maf-cutoff']:
                #     continue;
                # Make sure this allele is present at a frequency higher than the specified cutoff

                AC = in_allele_counts[allele] # derived allele count
                AA = rec.info['AA']
                derived = rec.alts[allele]
                if derived == AA:
                    AC = AN - AC
                # if ancestral allele is the same as derived, flip allele count
                AF =  AC / AN  # derived allele frequency
                print('AF = {}'.format(AF))
                
                polymorphic_codon = list(ref_codon);
                # Convert the ref_codon to a list so we can change nts by index

                polymorphic_codon[codon_pos] = alt_nts[int(allele)-1]
                mk_codons[rec_codon_pos]['poly'].append(polymorphic_codon)
                mk_codons[rec_codon_pos]['AF'].append(AF)
                # Add the polymorphic codon to the list of polymorphic codons
            ## End ingroup codon loop
        ## End ingroup codon block

        ####################

        for sample in globs['vcf-outgroups']:
            for allele in rec.samples[sample]['GT']:
                if allele is None:
                    continue
                out_allele_counts[allele] += 1;
        # Count alleles in the outgroup

        print('IN allele count = {}, OUT allele count = {}, IN hom alts = {}'.format(in_allele_counts, out_allele_counts, in_hom_alts))

        if 0 in out_allele_counts or not out_allele_counts:
            continue;
        # If there are any reference alleles (0) in the outgroup, then this cannot
        # be a fixed difference and it should be skipped
        # Likewise, in the case of no outgroup alleles (e.g. all missing data), this
        # should be skipped, else it would crash if there are also no ingroup alleles
        # (because of missing data or because all alt alleles are in the excluded samples)

        if all(alleles not in in_allele_counts for alleles in out_allele_counts):
        # Sites contain fixed differences only if all the alleles in the outgroup do not
        # exist in the ingroup

            mk_codons[rec_codon_pos]['fixed-flag'] = True;

            if len(out_allele_counts) == 1:
                mk_codons[rec_codon_pos]['fixed'][codon_pos] = alt_nts[list(out_allele_counts)[0]-1];
            # If there is only one allele in the outgroup (the site is not polymorphic in the outgroup)
            # add that allele to the fixed_diff_codon by looking it up in the alt_nts list

            else:
                max_allele = [ allele for allele, count in out_allele_counts.items() if count == max(out_allele_counts.values()) ];
                max_allele = random.choice(max_allele);
                mk_codons[rec_codon_pos]['fixed'][codon_pos] = alt_nts[max_allele-1];
            # If there is more than one alternate allele in the outgroup (the site is polymorphic in the outgroup)
            # use the allele at the highest frequency
            # If more than one allele exists at the highest frequency, pick one randomly
        ## End outgroup block
        ####################

    ## End transcript record loop
    ## End transcript record block

    ####################

    for codon_index in mk_codons:

        mk_codons[codon_index]['poly'] = [ "".join(poly_codon) for poly_codon in mk_codons[codon_index]['poly'] ];   
        mk_codons[codon_index]['fixed'] = "".join(mk_codons[codon_index]['fixed']);
        # Convert the fixed difference codons from lists to strings
    # End bookkeeping loop
    #sys.exit();

    return mk_codons, globs;

#############################################################################


