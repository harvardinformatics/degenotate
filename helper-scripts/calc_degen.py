#############################################################################
# Adds degeneracy info to a genetic code in .csv format
#############################################################################

import sys
import os

#############################################################################

usage_msg = """
# Given a .csv file with a genetic code (codon,aa), this adds a column of
# the degeneracy for each nucleotide in each codon.
# E.g. for the standard genetic code, the following entry in the input:
#       GGG,G
# will be converted to the following line in the output:
#       GGG,G,004
#
# Degeneracy codes:
# 0 = non-degenerate; any mutation will change the amino acid
# 2 = two nucleotides at the position code the same AA, so 1 of the three possible
#     mutations will be synonymous and 2 will be non-synonymous
# 3 = three nucleotides at the position code for the same AA, so 2 of the three possible
#     mutations will be synonymous and 1 will be non-synonymous
# 4 = four nucleotides at the position code for the same AA, so all 3 possible
#     mutations are synonymous
#

Usage: calc_degen.py [input file.csv] [output file.csv]
---> [input file.csv] must be a .csv file with two columns: codon,single letter amino acid code
---> [output file.csv] will be a .csv file with three columns, codon,single letter amino acid code, degeneracy of each nucleotide in the codon
"""
# A usage/help message

if len(sys.argv) != 3 or "-h" in sys.argv:
    print(usage_msg);
    sys.exit();
# Print the usage/help message if an incorrect number of inputs is specified, or if "-h" is specified

####################

input_file, output_file = sys.argv[1], sys.argv[2];
# Read the user inputs

if not os.path.isfile(input_file):
    print("\n\tError 1: Cannot find specified input file: " + input_file);
    sys.exit(1);
# Make sure the specified input file exists

if os.path.isfile(output_file):
    #print("\n\tError 2: Specified output file already exists. Move or delete it before running.");
    #sys.exit(1);

    overwrite = input("\n\tWarning: Specified output file already exists. Overwrite? (y/n): ");
    if overwrite.lower() not in ["y", "yes"]:
        print("\n\tNot overwriting output file. Exiting...");
        sys.exit(0);
    print();
# If the specified output file already exists, ask the user if they want to overwrite it

####################

nts = "ATCG";
aas = "*ACDEFGHIKLMNPQRSTVWY";
# Strings with characters for standard nucleotides and amino acids

print("# Reading codon table...");

codon_dict = {};
# The dictionary to read the codon info into

for line_str in open(input_file):
    line = line_str.strip().upper().split(",");
    line[0] = line[0].replace("U", "T");
    # Parse the line into a list, replacing U's in the codon with T's

    if len(line) > 2 or len(line[0]) != 3 or len(line[1]) != 1:
        print("\n\tError 3: Error reading the following line in the input file: " + line_str +"\n");
        print("\tLines must have 2 columns, with the first column being a codon string (3 nucleotide characters)");
        print("\tand the second being a single character amino acid code\n");
        sys.exit(1);
    # Check for several errors in each line

    if any(nt not in nts for nt in line[0]):
        print("\n\tError 4: The codon on the following line contains a non-nucleotide character: " + line_str +"\n");
        sys.exit(1);
    if line[1] not in aas:
        print("\n\tError 5: The character in the second column doesn't match a known amino acid: " + line_str +"\n");
        sys.exit(1);
    # Make sure only nucleotide and AA codes are given

    if line[0] in codon_dict:
        print("\n\tError 6: The codon on the following line appears more than once in the input file: " + line_str +"\n");
        sys.exit(1);
    # Make sure each codon only shows up once in the file

    codon_dict[line[0]] = line[1];
# Read the codon table into a dictionary as codon:aa and check for more errors

####################

print("# Calculating degeneracy per codon...");
with open(output_file, "w") as outfile:
    for codon in codon_dict:

        orig_aa = codon_dict[codon];
        degen_str = "";
        # Each codon has an AA and a degeneracy string

        for codon_pos in [0,1,2]:

            orig_nt = codon[codon_pos];
            muts_that_change_aa = 0;
            # Each position within a codon has a nucleotide and a count of how many mutations will change
            # the amino acid of that codon

            for new_nt in nts:
                if new_nt == orig_nt:
                    continue;
                # Check each other possible nucleotide to see if it changes
                # the current amino acid

                new_codon = list(codon);
                new_codon[codon_pos] = new_nt;
                new_codon = "".join(new_codon);
                new_aa = codon_dict[new_codon];
                # Replace the original nt with the new nt at the current codon
                # position and look up the new AA

                if new_aa != orig_aa:
                    muts_that_change_aa += 1;
                # If the new aa doesn't match the original one, increment the mut counter
            # End new_nt loop
            ##########

            if muts_that_change_aa == 0:
                degen_str += "4";
            elif muts_that_change_aa == 1:
                degen_str += "3";
            elif muts_that_change_aa == 2:
                degen_str += "2";
            elif muts_that_change_aa == 3:
                degen_str += "0"
            # Add the correct string to the degen_str for this codon position

        # End codon_pos loop
        ##########

        outline = ",".join([codon,orig_aa,degen_str]);
        outfile.write(outline + "\n");
        # Write to output

    # End codon loop
    ##########
# Close output file
##########
        