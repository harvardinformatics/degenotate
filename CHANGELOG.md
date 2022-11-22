2022.11.21
- Added `-l` to extract CDS sequences from longest transcripts and exit.

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