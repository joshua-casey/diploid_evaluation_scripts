Several diploid assembly evaluation scripts.

Requirements: python with matplotlib, BioConda, and pysam

### Plot Alignments

```plot_alignments.py```
This script will plot the alignments of reads in three measures, the direction of reads, clips, and coverage.

Accepts exactly one parameter <bam input file>

Adjust soft parameters (i.e. bin size) in lines 17 to 20 as needed.

Will write results to directories direction/, clips/, and coverage/ on current working directory

### Telomere signals plot 

```plot_telomere.py```
This script will plot the appearance of telomere signals.

Accepts exactly one parameter <fasta input file>

Adjust soft parameters (i.e. bin size and which telomere signal to find) in lines 17 to 20 as needed.

Will write results to directory telomere/ on current working directory

### Nanopore reads kmer phasing 

```haploid_kmers_ont.sh```
Need meryl in $PATH

Adjust file paths and K in lines 3-6

Outputs ont_haploid_kmers.txt

Use ```plot_haploid_kmers.py``` with ont_haploid_kmers.txt as input to plot.
