#############################################################################
# Functions to handle vcf files
#############################################################################

import sys
import os
import math
import random
from collections import defaultdict
import lib.core as CORE

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

    # print();

    # x = globs['vcf'].fetch("NC_046340.1", 71, 150);

    # print(x);
    # #print(list(x));
    # print(len(list(x)))

    # sys.exit();

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

    #sys.exit();
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

    genome_region = globs['annotation'][transcript]['header'];
    genome_start_pos = globs['coords'][transcript][transcript_position]-1;
    #genome_positions = list(range(genome_start_pos, genome_start_pos+3));
    # Convert the transcript coordinates to genomic ones for the VCF
    # Pysam coordinates are 0 based, so subtract 1
    ## TODO: Make sure this is consistent with our other data structures!!

    polymorphic_codons = [];
    fixed_diff_codon = ref_codon;
    # A list of polymorphic codons from the ingroups and a single codon string
    # for fixed differences in the outgroups
    
    records = globs['vcf'].fetch(genome_region, genome_start_pos, genome_start_pos+3)
    #get records for the current codon
    
    empty = next(records,None) is None
    #check if records returns anything
    
    if empty:
    #no record at this position, which means no variation at this codon_pos
        return polymorphic_codons, "".join(fixed_diff_codon);
    
    else:
    #one or more records at this position
    
        for rec in records:
            codon_pos = rec.start - genome_start_pos
            #print(rec.start, genome_start_pos, codon_pos) 
            #get codon position of this record
        
            ref_nt = rec.ref;
            alt_nts = rec.alts;
            # Look up the alleles at the current position

            in_allele_counts = defaultdict(int);
            out_allele_counts = defaultdict(int);
            # Dicts for allele counts in each group

            ####################

            for sample in globs['vcf-ingroups']:
                for allele in rec.samples[sample]['GT']:
                    if allele is None:
                        continue
                    #skip missing data
                    
                    in_allele_counts[allele] += 1;
            # Count alleles in the ingroups

            for allele in in_allele_counts:
                if allele == 0:
                    continue;
            
                # Construct a codon for each allele in the ingroups except the
                # reference allele (since that would match the reference codon)

                polymorphic_codon = ref_codon;
                # TODO: Make sure this doesn't copy just a pointer or something
                polymorphic_codon[codon_pos] = alt_nts[int(allele)-1];
                polymorphic_codons.append(polymorphic_codon);
                # Add the polymorphic codon to the list of polymorphic codons
            ## End ingroup allele loop

            ####################

            for sample in globs['vcf-outgroups']:
                for allele in rec.samples[sample]['GT']:
                    if allele is None:
                        continue
                    out_allele_counts[allele] += 1;
            # Count alleles in the outgroup

            if 0 in out_allele_counts:
                continue;
            # If all outgroup alleles are reference, add the reference allele to
            # the fixed_diff_codon

            if all(alleles not in in_allele_counts for alleles in out_allele_counts):
            # Sites contain fixed differences only if all the alleles in the outgroup do not
            # exist in the ingroup 
                if len(out_allele_counts) == 1:
                    fixed_diff_codon[codon_pos] = alt_nts[list(out_allele_counts)[0]-1];
                # If there is only one allele in the outgroup (the site is not polymorphic in the outgroup)
                # add that allele to the fixed_diff_codon by looking it up in the alt_nts list
                
                else:
                    max_allele = [ allele for allele, count in out_allele_counts.items() if count == max(out_allele_counts.values()) ];
                    max_allele = random.choice(max_allele);
                    fixed_diff_codon[codon_pos]= alt_nts[max_allele-1];
                # If there is more than one alternate allele in the outgroup (the site is polymorphic in the outgroup)
                # use the allele at the highest frequency
                # If more than one allele exists at the highest frequency, pick one randomly
            ## End outgroup block

            ####################
        ## End codon record block

    polymorphic_codons = [ "".join(poly_codon) for poly_codon in polymorphic_codons ];
    # Convert the polymorphic_codons from lists to strings

    return polymorphic_codons, "".join(fixed_diff_codon);
    # polymorphic_codons:   A list of codons that incorporate all SNPs in the ingroup samples, one codon
    #                       per alternate allele per site. As is, this list will never have the reference codon in it,
    #                       and could be an empty list if there are no SNPs in the ingroup samples.
    # fixed_diff_codon:     A single codon string that incorporates all fixed differences in the outgroup species
    #                       relative to the reference codon. A fixed difference only occurs if all the alleles in
    #                       the outgroup samples 1) are not the reference allele and 2) never occur in the ingroup 
    #                       samples. As currently implemented, if there are no fixed differences this will return
    #                       the reference codon.

#############################################################################

def getVariantsTranscript(globs, transcript, transcript_position, codons, extra_leading_nt, extra_trailing_nt):
# Gets variant codons in the in- and outgroups

    mk_codons = { c : { 'poly' : [], 'fixed' : codons[c], 'fixed-flag' : False } for c in range(len(codons)) };
    # This dictionary will keep track of polymorphic codons (initialized as empty list)
    # and fixed difference codons (initialized as reference codon) for each codon in the current
    # transcript
    # The fixed flag is swapped to True if any fixed differences have actually been found, else
    # the codon is still the ref and shouldn't be counted

    genome_region = globs['annotation'][transcript]['header'];
    # Look up the header for the transcript

    adj_ts_start = globs['annotation'][transcript]['start'] + extra_leading_nt - 1;
    adj_ts_end = globs['annotation'][transcript]['end'] - extra_trailing_nt;
    # Adjust the genomic start and end coordinates for this transcript based on 
    # the extra out of frame nts in the transcript

    transcript_records = globs['vcf'].fetch(genome_region, adj_ts_start, adj_ts_end);
    # Look up all variant records in the range of the adjusted transcript start and end coordinates

    empty = next(transcript_records,None) is None
    # Check if records returns anything
    if empty:
    # No record at this position, which means no variation at this codon_pos
        return mk_codons;
    else:
    # one or more records at this position    
        for rec in transcript_records:
            if rec.start not in globs['coords-rev'][transcript]:
                continue;
            # Skip any SNPs within the range of the transcript start and end,
            # but not in the CDS

            rec_transcript_pos = globs['coords-rev'][transcript][rec.start];
            # Look up the position of the record relative to the start of
            # the transcript

            adj_rec_pos = rec_transcript_pos - extra_leading_nt;
            # Adjust the transcript position based on the number of extra leading nts

            rec_codon_pos = math.floor(adj_rec_pos / 3);
            # Look up the codon position of the adjusted record position
            # This is also the key for mk_codons

            codon_pos = rec_codon_pos % 3;
            # The position of the record within the codon, either 0, 1, or 2

            ref_codon = codons[rec_codon_pos];
            # Look up the codon at the record's codon position
        
            ref_nt = rec.ref;
            alt_nts = rec.alts;
            # Look up the alleles at the current position

            in_allele_counts = defaultdict(int);
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

            for allele in in_allele_counts:
                if allele == 0:
                    continue;
                # Construct a codon for each allele in the ingroups except the
                # reference allele (since that would match the reference codon)         

                polymorphic_codon = list(ref_codon);
                # Convert the ref_codon to a list so we can change nts by index

                polymorphic_codon[codon_pos] = alt_nts[int(allele)-1];
                mk_codons[codon_pos]['poly'].append(polymorphic_codon);
                # Add the polymorphic codon to the list of polymorphic codons
            ## End ingroup allele loop

            ####################

            for sample in globs['vcf-outgroups']:
                for allele in rec.samples[sample]['GT']:
                    if allele is None:
                        continue
                    out_allele_counts[allele] += 1;
            # Count alleles in the outgroup

            if 0 in out_allele_counts:
                continue;
            # If all outgroup alleles are reference, add the reference allele to
            # the fixed_diff_codon

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
        ## End codon record block

    for codon_index in mk_codons:
        mk_codons[codon_index]['poly'] = [ "".join(poly_codon) for poly_codon in mk_codons[codon_index]['poly'] ];    
    # Convert the polymorphic_codons from lists to strings

    return mk_codons;

    # polymorphic_codons:   A list of codons that incorporate all SNPs in the ingroup samples, one codon
    #                       per alternate allele per site. As is, this list will never have the reference codon in it,
    #                       and could be an empty list if there are no SNPs in the ingroup samples.
    # fixed_diff_codon:     A single codon string that incorporates all fixed differences in the outgroup species
    #                       relative to the reference codon. A fixed difference only occurs if all the alleles in
    #                       the outgroup samples 1) are not the reference allele and 2) never occur in the ingroup 
    #                       samples. As currently implemented, if there are no fixed differences this will return
    #                       the reference codon.

#############################################################################