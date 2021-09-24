### Initial Roadmap and Design Notes

The most obvious, pressing need is a program that will take a GFF file and a fasta file as input, and output annotations of degeneracy in some BED or BED-like format. While there are many possible extensions to other things (variant annotation, polarization of subs and counting nonsyn/syn subs, probably others), as a first pass the goal should be to do the simple thing well and robustly. Other extensions may come as part of a refined VCF->MK toolkit, though.

These are just loose notes, not even pseudocode yet.

#### Input

Pretty definitely want to keep this so that the only input is the GFF and the reference fasta. 

Very loosely, will need to:
- extract coordinates of coding regions of transcripts from GFF, keeping track of splicing and strand
- keep track of translation between within CDS coordinate and within reference genome coordinate
- extract nucleotide sequence (accounting for strand) for each transcript
- keep track of transcript <-> gene key
- keep track of frame for incomplete / partial transcripts
- properly filter / ignore transcripts marked as pseudogenes, perhaps have option for other filtering (ignore partial transcripts)
- will need to decide what to do with transcripts with premature stop codons (often arise from assembly or annotation errors)

#### Processing

For each transcript, need to iterate through codons, and for each position score degeneracy

Requires:
- A function to compute degeneracy for each position in a codon, which requires consideration of edge cases (e.g., stop codons, weird asymetric codons like ATx, multiple mutuations in one codon)

#### Output

Default output is probably a bed-like file that indicates degeneracy, transcript id(s), and gene id for each position in the reference genome. May want to consider more compact options though as for big genomes this could get large (although not reporting non-coding bases would help)

We might want to consider various kinds of summary information; e.g., it is useful to output also a summary of syn and nonsyn sites per gene (or have a script that will compute this from the BED-ish output, to allow filtering based on callable sites beds or the like before calculations). This requires not just summing per gene, but translating degeneracy into syn/nonsyn (which requires consideration of whether you assume a uniform symmetric mutation matrix or not). 

#### Extensions

While writing intial code, worth keeping in mind the possible extensions. In particular in the MK context, polarizing polymorphisms to ancestral/derived is very useful and I think not well supported by existing tool, as is counting syn/nonsyn substitutions accounting for codon paths in multi-mutation codons. 

Also, SnpEff and VEP and similar often feel a bit "heavyweight" for evolutionary analysis and the core functions could possibly be replicated easily in a simple tool like this.

Finally, could consider additional annotations, like CNEEs, UTRs, introns, etc. Could be useful for people who want to do something like compute pi/theta for different functional compartments.
