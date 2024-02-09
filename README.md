
# degenotate

## Annotate degeneracy of sites in coding regions of a genome

[![Install](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://bioconda.github.io/recipes/degenotate/README.html)
[![OS](https://anaconda.org/bioconda/degenotate/badges/platforms.svg)](#degenotate)
[![Version](https://img.shields.io/conda/vn/bioconda/degenotate?label=version)](https://bioconda.github.io/recipes/degenotate/README.html)
[![Release Date](https://anaconda.org/bioconda/degenotate/badges/latest_release_date.svg)](#degenotate)
[![Downloads](https://img.shields.io/conda/dn/bioconda/degenotate.svg?style=flat)](https://bioconda.github.io/recipes/degenotate/README.html)
[![License](https://anaconda.org/bioconda/degenotate/badges/license.svg)](https://github.com/harvardinformatics/degenotate/blob/develop/LICENSE)

# Authors

### Gregg Thomas and Timothy Sackton

# Table of Contents

- [About](#about)
- [Installation](#installation)
    - [Installing with bioconda](#installing-from-bioconda)
    - [Installing from source](#installing-from-source)
- [Usage](#usage)
- [Output](#output)
    - [How degenotate classifies degeneracy](#how-degenotate-classifies-degeneracy)
    - [Degeneracy per site (bed file)](#degeneracy-per-site-bed-file)
    - [Transcript site counts (tab delimited)](#transcript-site-counts-tab-delimited)
    - [MK site counts (tab delimited)](#mk-site-counts-and-tests-tab-delimited)
- [Options](#options)
- [Assumptions](#assumptions)

# About

degenotate takes as input either a genome FASTA file and a corresponding annotation file (GFF or GTF) OR file or directory of files that contain coding sequences in FASTA format and outputs a bed-like file that contains the degeneracy score (0-, 2-, 3-, or 4-fold) of every coding site.

If given a corresponding VCF file with specified outgroup samples, degenotate can also count synonymous and non-synonymous polymorphisms and fixed differences for use in [MK tests](https://en.wikipedia.org/wiki/McDonald%E2%80%93Kreitman_test).
By default tries to run both standard ([McDonald and Kreitman 1991](https://doi.org/10.1038/351652a0) and imputed ([Murga-Moreno et al 2022](https://doi.org/10.1093/g3journal/jkac206))) MK tests. imputed MK test provides increased power to detect positive selection by taking into account the distribution of fitness effects.

The program also offers coding sequence extraction from the input genome and extraction of sequences by degeneracy (e.g. extract only the 4-fold degenerate sites).

**Warning: This is an early release. While we have done extensive testing, we are not certain our tests have hit all possible edge cases, especially those involving partial transcripts. We welcome bug reports and feature suggestions and are actively working to do more validation and testing.**

# Installation

## Installing from bioconda

We recommend installing degenotate from [bioconda](https://bioconda.github.io/recipes/degenotate/README.html) with the package manager [`conda`](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) or [`mamba`](https://mamba.readthedocs.io/en/latest/installation.html):

```
conda install degenotate
```

After this, you can run the program as `degenotate.py`. All dependencies should be automatically installed along with degenotate.

## Installing from source

Alternatively, since degenotate is purely Python and does not need compilation, one could simply download the program by cloning this repo and run it as `python degenotate.py`. In this case, you may want to add the degenotate folder to your $PATH variable for ease of use.

degenotate is a standalone program for its core function of annotating degeneracy on a site-by-site basis.

**The main dependency is Python 3+**

The VCF functionality (`-v`) specifically requires **Python version 3.10+** for the `itertools.pairwise()` function, as well as a some external libraries:

1. [pysam](https://pysam.readthedocs.io/en/latest/api.html) is used for efficent reading of VCF files for MK test site counts. pysam can be easily [installed with conda](https://anaconda.org/bioconda/pysam), but if you don't use the VCF options the program should run fine without it.
2. [NetworkX](https://networkx.org/) is used to easily trace the effect of mutations on different codons. NetworkX can also be [installed with conda](https://anaconda.org/conda-forge/networkx). Again if you don't use the VCF options the program should run without it.
3. [scipy](https://scipy.org/)'s implementation of Fisher's exact test is used to do the MK test. scipy can also be [installed with conda](https://anaconda.org/conda-forge/scipy). Again if you don't use the VCF options the program should run without it.

These dependencies are all installed automatically when degenotate is installed from bioconda. But if you need to clone this repo to install degenotate, we facilitate the installation of these dependencies by providing a pre-made conda environment, `environment.yml`. To create this environment, run:

```bash
conda env create -f environment.yml
```

And then:

```bash
conda activate degenotate
```

to activate the environment.

# Usage

### 1. Annotate degeneracy of coding sites from a genome in FASTA format and its annotation (gtf or gff):

```
python degenotate.py -a [annotation file] -g [genome fasta file] -o [output directory]
```

### 2. Annotate degeneracy of coding sites from a directory of individual coding sequences in FASTA format:

```
python degenotate.py -s [directory containing fasta files] -o [output directory]
```

### 3. Annotate degeneracy of coding sites from a genome and perform the MK test on every coding transcript:

```
python degenotate.py -a [annotation file] -g [genome fasta file] -v [vcf file] -u [sample ID(s) of outgroup in VCF file] -o [output directory] -sfs
```

### 4. Extract all coding sequences from the genome as nucleotides:

```
python degenotate.py -a [annotation file] -g [genome fasta file] -c [output file for CDS sequences] -o [output directory]
```

### 5. Extract coding sequences from the longest transcript of each gene as nucleotides and amino acids:

```
python degenotate.py -a [annotation file] -g [genome fasta file] -l [output file for CDS sequences] -la [output file for amino acid sequences] -o [output directory]
```

### 6. Extract 4-fold degenerate sites from all coding sequences:

```
python degenotate.py -a [annotation file] -g [genome fasta file] -x 4 -o [output directory]
```

# Output 

## How degenotate classifies degeneracy

| Fold | Description |
| ---- | ----------- |
| 0    | non-degenerate; any mutation will change the amino acid |
| 2    | two nucleotides at the position code the same AA, so 1 of the three possible mutations will be synonymous and 2 will be non-synonymous |
| 3    | three nucleotides at the position code for the same AA, so 2 of the three possible mutations will be synonymous and 1 will be non-synonymous |
| 4    | four nucleotides at the position code for the same AA, so all 3 possible mutations are synonymous |

## Degeneracy per site (bed file)

Default name: `[output directory]/degeneracy-all-sites.bed`

This is the main output file for degenotate. It contains one line for every coding site in the input genome, formatted with the following columns:

| Scaffold | Start pos | End pos | Transcript ID | Degeneracy code | Reference nucleotide | Reference amino acid | Mutation summary |
| -------- | --------- | ------- | ------------- | --------------- | -------------------- | -------------------- | ---------------- |
| The assembly scaffold or chromosome | The start position of the site | The end position of the site | The transcript ID | [See above](#how-degenotate-classifies-degeneracy) | The nucleotide at this site as read from the genome | The amino acid translated from the codon in that this site is in in the current transcript | [See below](#mutation-summary-column) |

### Mutation summary column

For non-degenerate sites (not 0-fold), the last column of the bed file contains information about how each mutation to non-degenerate nucleotides changes the amino acid. For example, if the final 4 columns of the bed file are:

```
2       A       E       T:D;C:D
```

This indicates that this site has 2-fold degeneracy, the nucleotide is A, and the codon that this nucleotide is part of in this transcript translates to E (Glutamic Acid). The final column shows the two nucleotides that change the amino acid at this position and the amino acid they change it to. Formatted as:

```
[nucleotide 1]:[corresponding amino acid 1];[nucleotide 2]:[corresponding amino acid 2]
```

For 3-fold sites, this would only have one `[nucleotide]:[amino acid]` entry and for 0-fold it would have three, each separated by a semi-colon.

## Transcript site counts (tab delimited)

Default name: `[output directory]/transcript-counts.tsv`

In addition to the information for every coding site, degenotate also outputs summaries by transcript. The columns in this file are:

| transcript | gene | cds_length | mrna_length | is_longest | f0 | f2 | f3 | f4 |
| ---------- | ---- | ----------------- | ------ | ------ | ------ | ------ |------ | ------ |
| Transcript ID | Gene ID | Length of coding sequence | Length of transcript | Indicator of longest transcript per gene | Count of 0-fold degenerate sites | Count of 2-fold degenerate sites | Count of 3-fold degenerate sites | Count of 4-fold degenerate sites | 

## MK site counts and tests (tab delimited)

Default name: `[output directory]/mk.tsv`

When provided with a multi-sample VCF file that includes outgroup samples, degenotate counts polymorphic and fixed differences and performs one or two versions of MK test: 1) standard [MK test](https://en.wikipedia.org/wiki/McDonald%E2%80%93Kreitman_test)  ([McDonald and Kreitman 1991](https://doi.org/10.1038/351652a0)) and 2) imputed MK test ([Murga-Moreno et al 2022](https://doi.org/10.1093/g3journal/jkac206)), which is performed if provided VCF is polarized. It also calculates the direction of selection ([Stoletzki and Eyre-Walker 2010](https://doi.org/10.1093/molbev/msq249)) for both versions of MK test. The output counts and calculations  are put in a file with the following columns:

| transcript | pN | pS | dN | dS | pval	| odds_ni |	dos | imp.pval | imp.odds_ni | imp.dos | pn_af | ps_af |
| ---------- | -- | -- | -- | -- | ---- | ------- | --- | -------- | ----------- | ------- | ----- | ----- |
| Transcript ID | Count of polymorphic non-synonymous sites | Count of polymorphic synonymous sites | Count of fixed non-synonymous sites | Count of fixed synonymous sites | The raw p-value from the MK test | The odds-ratio from the MK test, which is equivalent to the neutrality index | The direction of selection | The raw p-value from the imputed MK test | The odds ratio (neutrality index) from the imputed MK test | allele frequencies of each non-synonymous polymorphism found in the transcript (if -sfs specified) | allele frequencies of each synonymous polymorphism found in the transcript (if -sfs specified)|

# Options

| Option | Description | 
| :-------------------- | -------- |
| `-h`, `--help` | Show this help message and exit |
| `-a` | A GFF or GTF file that contains the coordinates of transcripts in the provided genome file (`-g`). Only one of -`a`/`-g` OR `-s` is REQUIRED. |
| `-g` | A FASTA file containing a genome. `-a` must also be specified. Only one of `-a`/`-g` OR `-s` is REQUIRED. |
| `-s` | Either a directory containing individual, in-frame coding sequence files or a single file containing multiple in-frame coding sequences on which to calculate degeneracy. Only one of `-a`/`-g` OR `-s` is REQUIRED. |
| `-v` | Optional VCF file with in- and outgroups to output polymorphic and fixed differences for MK tests. The VCF should contain SNPs only (no indels or structural variants). |
| `-u` | A comma-separated list of sample IDs in the VCF file that make up the outgroup (e.g. 'sample1,sample2') or a file with one sample per line. |
| `-e` | A comma separated list of sample IDs in the VCF file to exclude (e.g. 'sample1,sample2') or a file with one sample per line. |
| `-o` |  Desired output directory. This will be created for you if it doesn't exist. Default: `degenotate-[date]-[time]` |
| `-sfs` | Set this flag to output raw allele frequencies in the mk table|
| `-d` | degenotate assumes the chromosome IDs in the GFF file exactly match the sequence headers in the FASTA file. If this is not the case, use this to specify a character at which the FASTA headers will be trimmed. |
| `-c` | If a file is provided, the program will extract CDS sequences from the genome and write them to the file and exit. If no file is given with the option, a file with the name of 'cds-nt.fa' will be written to the output directory. This option is equivalent to '-x 0234' except this stops the program before calculating degeneracy. |
| `-ca` |  The same as `-c`, but writes translated amino acid sequences instead. Both `-c` and `-ca` can be specified. Default file name is 'cds-aa.fa'. |
| `-l` | If a file is provided, the program will extract CDS sequences from the longest transcript for each gene and write them to the file and exit. If no file is given with the option, a file with the name of 'cds-nt-longest.fa' will be written to the output directory. Both `-c` and `-l` can be specified. |
| `-la` |  The same as `-l`, but writes translated amino acid sequences instead. Both `-l` and `-la` can be specified to write both files. Default file name is 'cds-aa-longest.fa'. |
| `-x` | Extract sites of a certain degeneracy. For instance, to extract 4-fold degenerate sites enter '4'. To extract 2- and 4-fold degenerate sites enter '24' and so on. | 
| `-m` | The minimum length of a transcript for it to be counted. Default (and global min): 3 | 
| `-maf` | The minor allele frequency cutoff for MK tests. Sites where alternate alleles in the ingroup are below this frequency will be excluded. Default: 1 / 2N, where N is the number of ingroup samples | 
| `-imp` | The minor allele frequency cutoff that distinguishes low and high allele frequencies for imputed MK test. Only used if provided VCF is polarized. Default: 0.15 |
| `--no-fixed-in` | Set this if you wish to exclude sites from the MK test in which all ingroup samples share the same alternate allele (only the reference differs). | 
| `--overwrite` | Set this to overwrite existing files. |
| `--appendlog` | Set this to keep the old log file even if `--overwrite` is specified. New log information will instead be appended to the previous log file. |
| `--info` |  Print some meta information about the program and exit. No other options required. |
| `--version` | Simply print the version and exit. Can also be called as `-version`, or `--v`. |
| `--quiet` | Set this flag to prevent degenotate from reporting detailed information about each step. |

## Assumptions

To compute MK tables of synonymous and nonsynonymous polymorphism and divergence requires making certain assumptions. Here are the assumptions built into degenotate. 

- Codons with multiple polymorphisms: If a codon has more than one variant segregating within a population (either because multiple positions at the codon have segregating sites, or because one position has a multi-allelic SNP), we treat each segregating variant as independent.
- Codons with multiple fixed differences: If there are multiple fixed differences in a single codon in the outgroup, we compute all possible mutational pathways between the ingroup codon and the outgroup codon, and take the average number of nonsynonymous and synonymous changes across these paths, weighted equally. This means we can have fractional numbers of synonymous and nonsynonymous divergence. 
- Defining fixed differences: We consider a site to be a fixed difference only if all the alleles in the outgroup do not exist in the ingroup. This means, for example, that a polymorphic site in the outgroup can still contain a fixed difference as long as both the alleles are different from any ingroup allele. In this case, we use the highest frequency outgroup allele.
- Defining polymorphic sites: We consider polymorphic sites to be any site where at least one ingroup individual has a non-reference allele. In most cases this is intuitive, however there are two edge cases worth pointing out. First, if all ingroup individuals are homozygous alternate, we will count that position as a polymorphism **unless the `--no-fixed-in` option is specified**. This has implications if the reference genome is from a different population or species than your ingroup sequence. Second, we consider sites with shared polymorphism between the ingroup and the outgroup to be polymorphic, since we do not consider outgroup sequence at all when defining polymorphisms.
- Site frequency: Use the `-maf` option to specify a minimum ingroup minor allele frequency to consider sites. For example, with `-maf 0.3`, sites will only be considered for polymorphism if the minor allele exists in 30% of haplotypes (2 times the number of individuals). For now, if you want to exclude singletons or other low-frequency **fixed differences** from MK calculations, you need to remove these sites from your VCF before running degenotate.
- Weakly deleterious polymorphisms (of non-synonymous class) are segregating at low allele frequencies and their fraction can be estimated from the polymorphisms of synonymous class (main assumption of imputed MK test).

**Any of these assumptions may change in future releases.** 
