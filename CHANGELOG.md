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