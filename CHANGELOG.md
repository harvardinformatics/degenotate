2026.11.03
- Added a strand column (column 6) to `degeneracy-all-sites.bed`, making it a BED6+ file, and a `strand` column to `transcript-counts.tsv` (#56)
- Clarified in the README that the reference nucleotides and the nucleotides in the mutation summary column are reported on the strand of the transcript, so they are the complement of the genome for transcripts on the `-` strand (#56)
- Added a pytest test suite (`tests/`) with small hand-checked inputs covering degeneracy annotation, sequence extraction, the command line options, and the MK tables and tests
- Added a GitHub Actions workflow that runs the tests on every push and pull request
- Fixed MK counts skipping variants at the first base of a transcript's first full codon
- Fixed MK counts for partial transcripts on the `-` strand with different numbers of leading and trailing out-of-frame bases, where variants were placed at the wrong codon position
- Fixed MK counts missing a variant at the first base of a transcript, which happens when a transcript has no UTR at that end
- Fixed pN and pS including alleles at or below the `-maf` cutoff when `-imp` is set below it
- `-maf 0` and `-imp 0` are no longer replaced by the default cutoffs
- The neutrality index (`odds_ni`, `imp.odds_ni`) is now calculated from the exact counts, which can be fractional. The p-values from Fisher's exact test use the counts rounded to the nearest integer. Previously, scipy truncated fractional counts for both
- Fixed `mrna_length` in `transcript-counts.tsv` being one base short. The minimum transcript length (`-m`) now uses the correct length, so transcripts exactly `-m` bases long are kept
- Blank lines in the annotation file are now skipped instead of crashing
- `-ca` and `-la` now translate codons with ambiguous bases as `X` instead of crashing, and now exit after writing sequences, like `-c` and `-l`
- Using `-v` with `-s` now gives an error (OP17) instead of crashing, and a record with no AA field in a polarized VCF now gives an error (VCF4) instead of a traceback
- Error codes used twice were renumbered: `-m` errors are now OP16 and missing outgroup samples are now VCF3
- `-imp` and `-sfs` without `-v` now give warnings, and the start of the log only lists the output files that will be written
- The log now shows the default `-maf` setting as `1 / 2N (default)` instead of `False`, and records the cutoff calculated from the VCF file. It also now shows the `-imp` setting
- Clarified the `-maf` and `-imp` documentation (alleles at or below the `-maf` cutoff are not counted), and fixed the OP2 error message and the `--version` help

2023.08.31
- Added a check for trailing semi-colons in GFF file info fields

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
