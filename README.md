
# degenotate

## Annotate degeneracy of sites in coding regions of a genome

# Authors

### Gregg Thomas and Timothy Sackton

## About

degenotate takes as input either a genome FASTA file and a corresponding annotation file (GFF or GTF) OR file or directory of files that contain coding sequences in FASTA format and outputs a bed-like file that contains the degeneracy score (0-, 2-, 3-, or 4-fold) of every coding site.

If given a corresponding VCF file with specified outgroup samples, degenotate can also count synonymous and non-synonymous polymorphisms and fixed differences for use in [MK tests](https://en.wikipedia.org/wiki/McDonald%E2%80%93Kreitman_test).

The program also offers coding sequence extraction from the input genome and (coming soon) extraction of sequences by degeneracy (e.g. extract only the 4-fold degenerate sites).

## Installation

### Install

Simply download the program by cloning this repo and run it as `python degenotate.py`. You may want to add the degenotate folder to your $PATH variable for ease of use.

### Dependencies

degenotate is a standalone program for its core function of annotating degeneracy on a site-by-site basis.

**The main dependency is Python 3+**

The VCF functionality (`-v`) specifically requires **Python version 3.10+** for the `itertools.pairwise()` function, as well as a couple of external packages.

1. [pysam](https://pysam.readthedocs.io/en/latest/api.html) is used for efficent reading of VCF files for MK test site counts. pysam can be easily [installed with conda](https://anaconda.org/bioconda/pysam), but if you don't use the VCF options the program should run fine without it.
2. [NetworkX](https://networkx.org/) is used to easily trace the effect of mutations on different codons. NetworkX can also be [installed with conda](https://anaconda.org/conda-forge/networkx). Again if you don't use the VCF options the program should run without it.

We facilitate the installation of these dependencies by providing a pre-made conda environment, `environment.yml`. To create this environment, run:

```bash
conda env create -f environment.yml
```

And then:

```bash
conda activate degenotate
```

to activate the environment.

## Usage

1. Annotate degeneracy of coding sites from a genome:

```bash
python degenotate.py -a [annotation file] -g [genome fasta file] -o [output directory]
```

2. Annotate degeneracy of coding sites from a directory of individual coding sequences in FASTA format:

```bash
python degenotate.py -s [directory containing fasta files] -o [output directory]
```

3. Annotate degeneracy of coding sites from a genome and output synonymous and non-synonomous polymorphisms and fixed differences for MK tests:

```
python degenotate.py -a [annotation file] -g [genome fasta file] -v [vcf file] -u [file containin outgroup samples] -o [output directory]
```

4. Extract coding sequences from genome:

```
python degenotate.py -a [annotation file] -g [genome fasta file] -c [output file for CDS sequences]
```
## Output 

### How degenotate classifies degeneracy

| Fold | Description |
| ---- | ----------- |
| 0    | non-degenerate; any mutation will change the amino acid |
| 2    | two nucleotides at the position code the same AA, so 1 of the three possible mutations will be synonymous and 2 will be non-synonymous |
| 3    | three nucleotides at the position code for the same AA, so 2 of the three possible mutations will be synonymous and 1 will be non-synonymous |
| 4    | four nucleotides at the position code for the same AA, so all 3 possible mutations are synonymous |

### Degeneracy per site (bed file)

Default name: `[output directory]/degeneracy-all-sites.bed`

This is the main output file for degenotate. It contains one line for every coding site in the input genome, formatted with the following columns:

| Scaffold | Start pos | End pos | Transcript ID | Degeneracy code | Reference nucleotide | Reference amino acid | Mutation summary |
| -------- | --------- | ------- | ------------- | --------------- | -------------------- | -------------------- | ---------------- |
| The assembly scaffold or chromosome | The start position of the site | The end position of the site | The transcript ID | [See above](#how-degenotate-classifies-degeneracy) | The nucleotide at this site as read from the genome | The amino acid translated from the codon in that this site is in in the current transcript | [See below](#mutation-summary-column) |

#### Mutation summary column

For non-degenerate sites (not 0-fold), the last column of the bed file contains information about how each mutation to non-degenerate nucleotides changes the amino acid. For example, if the final 4 columns of the bed file are:

```
2       A       E       T:D;C:D
```

This indicates that this site has 2-fold degeneracy, the nucleotide is A, and the codon that this nucleotide is part of in this transcript translates to E (Glutamic Acid). The final column shows the two nucleotides that change the amino acid at this position and the amino acid they change it to. Formatted as:

```
[nucleotide 1]:[corresponding amino acid 1];[nucleotide 2]:[corresponding amino acid 2]
```

For 3-fold sites, this would only have one `[nucleotide]:[amino acid]` entry and for 0-fold it would have three, each separated by a semi-colon.

### Transcript site counts (tab-delimited)

Default name: `[output directory]/transcript-counts.tsv`

In addition to the information for every coding site, degenotate also outputs summaries by transcript. The columns in this file are:

| transcript | gene | transcript length | 0-fold | 2-fold | 3-fold | 4-fold |
| ---------- | ---- | ----------------- | ------ | ------ | ------ | ------ |
| Transcript ID | Gene ID | Length of transcript | Count of 0-fold degenerate sites | Count of 2-fold degenerate sites | Count of 3-fold degenerate sites | Count of 4-fold degenerate sites | 

### MK site counts (tab delimited)

Default name: `[output directory]/mk.tsv`

When provided with a multi-sample VCF file and outgroup samples, degenotate counts polymorphic and fixed differences for MK tests. The output counts are put in a file with the following columns:

| transcript | pN | pS | dN | dS | 
| ---------- | -- | -- | -- | -- | 
| Transcript ID | Count of polymorphic non-synonymous sites | Count of polymorphic synonymous sites | Count of fixed non-synonymous sites | Count of fixed synonymous sites | 

## Options

| Option | Description | 
| ------ | ----------- |
| `-a` | A GFF or GTF file that contains the coordinates of transcripts in the provided genome file (`-g`). Only one of -`a`/`-g` OR `-s` is REQUIRED. |
| `-g` | A FASTA file containing a genome. `-a` must also be specified. Only one of `-a`/`-g` OR `-s` is REQUIRED. |
| `-s` | Either a directory containing individual, in-frame coding sequence files or a single file containing multipl in-frame coding sequences on which to calculate degeneracy. Only one of `-a`/`-g` OR `-s` is REQUIRED. |
| `-v` | Optional VCF file with in and outgroups to output polymorphic and fixed differences for MK tests. |
| `-u` | A comma separated list of sample IDs in the VCF file that make up the outgroup (e.g. 'sample1,sample2') or a file with one sample per line. |
| `-e` | A comma separated list of sample IDs in the VCF file to exclude (e.g. 'sample1,sample2') or a file with one sample per line. |
| `-o` |  Desired output directory. This will be created for you if it doesn't exist. Default: `degenotate-[date]-[time]` |
| `-d` | degenotate assumes the chromosome IDs in the GFF file exactly match the sequence headers in the FASTA file. If this is not the case, use this to specify a character at which the FASTA headers will be trimmed. |
| `-c` | If a file is provided, the program will extract CDS sequences from the genome and write them to the file and exit. |
| `--overwrite` | Set this to overwrite existing files. |
| `--appendlog` | Set this to keep the old log file even if `--overwrite` is specified. New log information will instead be appended to the previous log file. |
| `--info` |  Print some meta information about the program and exit. No other options required. |
| `--version` | Simply print the version and exit. Can also be called as `-version`, or `--v`. |
| `--quiet` | Set this flag to prevent degenotate from reporting detailed information about each step. |
