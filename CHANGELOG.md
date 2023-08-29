2023.08.29
- Fixed bug with processing CDS fasta files

2023.05.12
- Added warnings for and skip alternate alleles longer than 1bp in the VCF

2023.04.26
- Added the `-maf` option that allows the user to specify a minor allele frequency cutoff for ingroups in the MK tests; alleles below the specified frequency will be excluded. The default `-maf` cutoff is 1 / 2N where N is the number of ingroup samples.
- Added `--no-fixed-ingroup` option that will exclude sites from MK tests if all ingroup samples are fixed for an alternate allele, which is likely an error in the reference

2023.04.20
- Converted error for transcripts with exons on differing strands to a warning, and added warning for transcripts with no coding exons

2023.02.16
- Added path to python executable to runtime info to help track down version issues as they arise
- Now check for VCF indices with both .tbi and .csi extensions

2023.02.07
- Fixed bug that miscounted substitutions on transcripts on the negative strand since the variants in the VCF file are all reported on the positive strand with respect to the reference

2023.02.06
- Added `-m` option so user can specify minimum transcript length, with a default value (and global minimum) of 3bp

2023.02.02
- Added check to skip invariant sites in the VCF file
- Added a check to skip sites where there are no alleles in the outgroup, which could happen in the case of missing data or for variant sites in which all the alternate alleles are among excluded samples. In the latter case, the program would crash because it would try to select an allele from an empty list

2023.01.30
- We now skip transcripts with 0 length in the input annotation and print warnings about them in the log file
- Updated `environment.yml` to include the [scipy](https://anaconda.org/conda-forge/scipy) dependency

2023.01.26
- Fixed bug when reading feature info from a gxf file and the field splitter remained on the last entry, causing no exons to be read
- Added scipy dependency for Fisher's test
- Implemented MK test (outputs raw p-value and odds ratio, which is equivalent to the neutrality index) and DoS statistic (from https://doi.org/10.1093/molbev/msq249)
- Added `-ca` and `-la` to write amino acid sequences when specified, which necessitated adding the bioTranslator() and readCodonTable() functions
- Changed how `-c`, `-ca`, `-l`, and `-la` behaved, allowing users to provide a file name, but using a default if none is provided

2022.11.21
- Added `-l` to extract CDS sequences from longest transcripts and exit

2022.11.16
- Numerous bugfixes in McDonald-Kreitman test calculations
- Efficiency improvements, now can process about 1000 transcripts per minute (for a VCF with 20 individuals)
- Updated README
- Added `-x` option to extract sequences by degeneracy
- Change column headers of output to play nice with R

2022.10.07
- Added transcript summary output
- Added `-d` option for users to specify delimiters on which to split FASTA headers in the genome file
- Added error checking for dependencies for VCF functions (pysam and networkx)
- Added README
