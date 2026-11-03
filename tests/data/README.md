# Test fixtures

Small, hand-checkable inputs for the test suite. All expected values were worked out from the standard genetic code ([genetic-code.csv](../../degenotate_lib/genetic-code.csv)), independently of degenotate's code. **Don't regenerate the files in `expected/` by copying degenotate's output.** If a test fails, first work out by hand which side is wrong.

| File | Contents |
| ---- | -------- |
| `tiny.fa` | One 287 bp chromosome, `chr1`. `N`s separate the transcripts. TX2's second CDS exon (84-94) is soft-masked (lowercase). |
| `tiny.gff3`, `tiny.gtf` | The same annotation in both formats. |
| `tiny-cds.fa` | TX1 and TX2 CDS sequences (transcript strand) for `-s` input, plus `OUT1`, which has 10 bp and is not a multiple of 3, so it is skipped. |
| `tiny.vcf` | 16 SNPs. Samples `in01`-`in10` are the ingroup (20 chromosomes, so the default MAF cutoff is 1/20) and `out1` is the outgroup. The tests bgzip and index it into a temporary directory. |
| `tiny-polarized.vcf` | The same SNPs with an `AA` (ancestral allele) INFO field. |
| `expected/` | Expected outputs for the degeneracy tests. |

Coordinates are 1-based, as in the GFF/VCF. The bed output is 0-based half-open, so genome position *g* is the bed row `g-1  g`.

## Transcripts

CDS sequences and codons are in transcript orientation. For `-` strand transcripts, the genome holds the reverse complement.

| ID | Gene | Strand | Transcript | CDS exons (transcript order) | Phase | CDS | What it tests |
| -- | ---- | ------ | ---------- | ---------------------------- | ----- | --- | ------------- |
| TX1 | G1 | + | 11-44 | 16-39 | 0 | `ATG GCT TTA ATT TGG CGA GAA TAA` | 0/2/3/4-fold sites; main MK cases |
| TX3 | G1 | + | 55-68 | 58-66 | 0 | `ATG GCT TAA` | shorter isoform of G1: `is_longest` = no |
| TX2 | G2 | - | 79-116 | 105-111, 84-94 | 0, 2 | `ATG AAA G\|GC CCC TTC TGA` | 2 exons split mid-codon (`\|`), reverse complement, soft-masking, − strand MK cases |
| TX4 | G3 | + | 127-138 | 127-138 | 2 | `GC AAA TGC GGT A` | partial: 2 leading and 1 trailing out-of-frame bases (`.` rows) |
| TX5 | G4 | - | 149-170 | 154-165 | 1 | `G CAT AAA TTT AC` | partial on − strand: 1 leading, 2 trailing bases, so variant positions must be offset by the leading bases |
| TX6 | G5 | + | 181-199 | 186-194 | 0 | `ATG CCC TAA` | variant at the first base of the first codon |
| TX7 | G6 | - | 210-223 | 210-218 | 0 | `ATG GGA TAG` | transcript starts at its CDS; variant at that base (the first base of the VCF lookup) |
| TX8 | G7 | + | 234-252 | 239-247 | 0 | `ATG CAT TAA` | two fixed differences in one codon |
| TXS | G8 | + | 263-264 | 263-264 | 0 | | shorter than the minimum length: warned about and skipped |
| TXN | G9 | + | 268-277 | none | | | no CDS: warned about and removed |

### Degeneracy per codon

For each position, count how many of the 3 possible mutations keep the amino acid: 0 → 0-fold, 1 → 2-fold, 2 → 3-fold, 3 → 4-fold. Stop codons count as an amino acid (`*`).

| Codon | AA | Fold | Reasoning |
| ----- | -- | ---- | --------- |
| ATG | M | 000 | the only Met codon |
| GCT | A | 004 | GCN are all Ala |
| TTA | L | 202 | CTA = Leu; TTG = Leu (TTT/TTC = Phe) |
| ATT | I | 003 | ATC/ATA = Ile, ATG = Met |
| TGG | W | 000 | the only Trp codon |
| CGA | R | 204 | AGA = Arg (TGA stop, GGA Gly); CGN all Arg |
| GAA | E | 002 | GAG = Glu (GAT/GAC = Asp) |
| TAA | * | 022 | TGA = stop; TAG = stop |
| AAA | K | 002 | AAG = Lys (AAT/AAC = Asn) |
| GGC, CCC, GGT, GGA | G, P, G, G | 004 | 4-fold third positions |
| TTC | F | 002 | TTT = Phe (TTA/TTG = Leu) |
| TGA | * | 020 | TAA = stop; TGT/TGC = Cys, TGG = Trp |
| TGC | C | 002 | TGT = Cys |
| CAT | H | 002 | CAC = His (CAA/CAG = Gln) |
| TTT | F | 002 | TTC = Phe |
| TAG | * | 002 | TAA = stop (TAT/TAC = Tyr) |

### Transcript counts

These are the sums of the folds above. `mrna_length` is end − start + 1. Out-of-frame bases (TX4, TX5) aren't counted.

| ID | cds_length | mrna_length | f0 | f2 | f3 | f4 |
| -- | ---------- | ----------- | -- | -- | -- | -- |
| TX1 | 24 | 34 | 15 | 6 | 1 | 2 |
| TX3 | 9 | 14 | 6 | 2 | 0 | 1 |
| TX2 | 18 | 38 | 13 | 3 | 0 | 2 |
| TX4 | 12 | 12 | 6 | 2 | 0 | 1 |
| TX5 | 12 | 22 | 6 | 3 | 0 | 0 |
| TX6 | 9 | 19 | 6 | 2 | 0 | 1 |
| TX7 | 9 | 14 | 7 | 1 | 0 | 1 |
| TX8 | 9 | 19 | 6 | 3 | 0 | 0 |

`test_degeneracy.py` also puts an `N` at chr1:62, making TX3's middle codon `GNT`. Those 3 sites get degeneracy `.`, and TX3's counts become f0 = 4, f2 = 2, f3 = 0, f4 = 0.

## VCF sites

"Poly" sites are polymorphic in the ingroup with the outgroup reference. "Fixed" sites have the outgroup `1/1` and the ingroup `0/0`. AF is the ALT count over the called ingroup alleles. For `-` strand transcripts, the VCF alleles are the complements of the transcript bases. Every TX2 site is at a 2-fold position, so leaving the ALT uncomplemented would flip synonymous ↔ nonsynonymous.

| Pos | REF>ALT | Transcript change | Type | AF | Counts as |
| --- | ------- | ----------------- | ---- | -- | --------- |
| 21 | T>C | TX1 GCT→GCC (A→A) | poly, in10 missing | 3/18 = 0.167 | pS |
| 24 | A>G | TX1 TTA→TTG (L→L) | fixed | | dS |
| 25 | A>G | TX1 ATT→GTT (I→V) | poly, singleton | 1/20 = 0.05 | nothing: AF is not above the 1/20 cutoff |
| 29 | G>C | TX1 TGG→TCG (W→S) | poly | 0.1 | pN |
| 32 | G>A | TX1 CGA→CAA (R→Q) | fixed | | dN |
| 33 | A>G | TX1 CGA→CGG (R→R) | all ingroups 1/1 | 1.0 | pS (nothing with `--no-fixed-in`) |
| 36 | A>G,T | TX1 GAA→GAG (E→E), GAA→GAT (E→D) | poly, multi-allelic | 0.1, 0.2 | pS + pN |
| 87 | G>A | TX2 TTC→TTT (F→F) | poly | 0.2 | pS |
| 89 | A>G | TX2 TTC→CTC (F→L) | fixed | | dN |
| 106 | T>C | TX2 AAA→AAG (K→K) | fixed | | dS |
| 108 | T>C | TX2 AAA→GAA (K→E) | poly; AA = C in the polarized file | 0.1 | pN |
| 161 | T>G | TX5 AAA→CAA (K→Q) | fixed | | dN |
| 186 | A>G | TX6 ATG→GTG (M→V) | poly | 0.1 | pN |
| 210 | C>T | TX7 TAG→TAA (*→*) | poly | 0.1 | pS |
| 243, 244 | A>G, T>A | TX8 CAT→CGA (H→R) | fixed, 2 in one codon | | dN 1.5, dS 0.5 |

TX8: the two shortest paths are CAT→CGT (N)→CGA (S) and CAT→CAA (N)→CGA (N). Averaged, that is dN = 3/2, dS = 1/2.

### Expected MK tables (default options)

| ID | pN | pS | dN | dS | DoS | odds_ni | pval |
| -- | -- | -- | -- | -- | --- | ------- | ---- |
| TX1 | 2 | 3 | 1 | 1 | 1/2 − 2/5 = 0.1 | (2·1)/(3·1) = 0.667 | 1.0 |
| TX2 | 1 | 1 | 1 | 1 | 0 | 1 | 1.0 |
| TX3, TX4 | 0 | 0 | 0 | 0 | | | |
| TX5 | 0 | 0 | 1 | 0 | | | |
| TX6 | 1 | 0 | 0 | 0 | | | |
| TX7 | 0 | 1 | 0 | 0 | | | |
| TX8 | 0 | 0 | 1.5 | 0.5 | NA | NA | NA |

- DoS = dN/(dN+dS) − pN/(pN+pS), and odds_ni = (pN·dS)/(pS·dN).
- TX1's p-value: in Fisher's exact test on [[2, 3], [1, 1]], every table with the same margins is at most as likely as the observed one, so p = 1.
- TX8's statistics are NA because pN + pS = 0.
- TX5: if variant positions were offset by the trailing bases instead of the leading ones, its fixed difference would be read as CAT→CAC, giving dN = 0, dS = 1.
- Other options:
  - `--no-fixed-in`: TX1 pS = 2.
  - `-maf 0.1`: TX1 pN = 1, pS = 2.
  - `-maf 0`: TX1 pN = 3, pS = 3 (the singleton is counted).
  - `-imp 0.01`: TX1 is unchanged (pN = 2, pS = 3). `-imp` only splits the polymorphisms above the MAF cutoff into low and high classes.
  - `-sfs` for TX1: pn_af = `0.05,0.1,0.2`, ps_af = `0.17,1.0,0.1` (codon order, including AFs below the cutoff).

### Imputed MK (`tiny-polarized.vcf`)

Here AA equals REF everywhere except at position 108, where it is the ALT, so that SNP's derived frequency flips from 0.1 to 0.9. Low frequency means 1/20 < AF ≤ 0.15; high means AF > 0.15. With p_wd = p_low − p_high/p0_high · p0_low:

- **TX1:**
  - pN: 1 low (0.1) and 1 high (0.2). pS: 1 low (0.1) and 2 high (0.167, 1.0).
  - p_wd = 1 − ½·1 = 0.5.
  - imp.dos = 1/2 − 1.5/(1.5+3) = 1/6.
  - imp.odds_ni = (1.5·1)/(3·1) = 0.5, from the exact counts. Fisher's test needs integers, so imp.pval comes from the rounded table [[2, 3], [1, 1]].
- **TX2:**
  - pN: 0 low and 1 high (0.9, after the flip). pS: 1 high (0.2).
  - p_wd = 0.
  - imp.dos = 1/2 − 1/2 = 0. Without the flip it would be 0.5.
- **`-imp 0`:** there is no low class, so p_wd = 0 and TX1's imp.dos equals its dos, 0.1.
