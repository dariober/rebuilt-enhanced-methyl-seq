Table of Contents
=================

  * [Table of Contents](#table-of-contents)
  * [Alignment and methylation calling](#alignment-and-methylation-calling)
  * [Alignment of reads from <em>Escherichia coli</em> libraries](#alignment-of-reads-from-escherichia-coli-libraries)
  * [Detecting methylation](#detecting-methylation)
  * [Identification of methylated blocks via hidden Markov model](#identification-of-methylated-blocks-via-hidden-markov-model)
  * [Observed and expected duplication rate](#observed-and-expected-duplication-rate)



<!---
TOC created with ./gh-md-toc from https://github.com/ekalinin/github-markdown-toc as
gh-md-toc ~/git_sblab/rebuilt-enhanced-methyl-seq/trunk/Method.md
-->

# Alignment and methylation calling

Prior to alignment, reads have been trimmed to remove adapter contamination and low quality 3’ ends using [trim_galore 0.3.7](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore)
with the following settings:

```
trim_galore --stringency 3 --paired -o fastq_trimmed $fq1 $fq2
```

Where `$fq1` and `$fq2` are paired raw fastq files.

Alignment was performed with [bwameth.py](https://github.com/brentp/bwa-meth) with underlying [bwa mem](https://github.com/lh3/bwa) version 0.7.10-
r789:

```
bwameth.py index $ref
bwameth.py -t 4 --reference $ref --prefix $bname.mm9_pb11.un $fq1 $fq2
```

Where `$ref` is a fasta file comprising both the Mus musculus genome version [mm9](http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz) and
the Plasmodium berghei genome version [11](http://plasmodb.org/common/downloads/release-11.0/PbergheiANKA/fasta/data/PlasmoDB-11.0_PbergheiANKA_Genome.fasta).
`$fq1` and `$fq2` are trimmed paired fastq files and `$bname` is a variable to identify each output file.

For methylation calling the aligned reads were first processed to remove reads showing
mismatch rate in excess of 10% using [resetHighMismatchReads.py](https://github.com/dariober/bioinformatics-cafe/blob/master/resetHighMismatchReads.py).
In order to avoid double counting coverage due to overlapping read pairs, one of the two reads was
soft clipped with `clipOverlap` in the [bamUtil](http://genome.sph.umich.edu/wiki/BamUtil) suite version 1.0.12:

```
resetHighMismatchReads.py -i $bname.mm9_pb11.un.bam \
| bam clipOverlap --in - --out $bname.mm9_pb11.bam --stats --storeOrig XC
```

Finally, the count of methylated, unmethylated and mismatching (i.e. not C or T) reads at each cytosine was obtained with the script
[bam2methylation.py](https://github.com/dariober/bioinformatics-cafe/blob/master/bam2methylation.py) version 0.3.
In addition, the nucleotide composition following each cytosine (i.e. the cytosine context)
was assigned with [addCGcontextToBdg.sh](https://github.com/dariober/bioinformatics-cafe/blob/master/addCGcontextToBdg.sh) version 0.2:

```
bam2methylation.py -mq 13 -mm -A --samargs '-q 15' -i $bam -r $ref -l $bed > $bname.met.bdg.tmp
addCGcontextToBdg.sh $bname.met.bdg.tmp $ref 7 > $bname.met.bdg
```

Here `$bam` indicates the alignment file generated above, `$ref` is the reference fasta file
for _P. berghei_, `$bed` is a bed file of chromosome sizes, and `$bname` a variable to identify
each output file. The options passed to `bam2methylation.py` specify that only read
bases with quality 13 or above (`-mq 13`) and with mapping quality 15 or above (`--samargs '-q 15'`)
are included; anomalous read pairs are included (`-A`). Option `-mm` is used to include in output the number of mismatches.

# Alignment of reads from _Escherichia coli_ libraries

Read trimming and alignment to the E. coli genome was performed as described above
for mouse and Plasmodium. Raw fastq files of the PCR-BS library were downsampled to
3534470 reads to be comparable with the ReBuilT library. The _Escherichia coli_ reference genome
with accession ID CP001509.3 was downloaded from GenBank.

# Detecting methylation

The alignment and methylation calling described above produces for each library a table
of count of methylated, unmethylated and mismatching reads at each cytosine in the
genome. To test whether a cytosine has more methylated reads than expected by chance
we applied the following Fisher test to each cytosine implemented in [R](https://cran.r-project.org/):

```
fisher.test(line, alternative = "greater", conf.int= FALSE, or= 0.5)
```

Where `line` is a 2x2 matrix with: count of methylated cytosine, count of unmethylated
cytosines, count of mismatches (not C or T), count of C and T.

To merge methylation calls within the two PCR-BS and three REBUiLT replicates we
combined p-values using [Stouffer’s method](http://en.wikipedia.org/wiki/Fisher%27s_method#Relation_to_Stouffer.27s_Z-score_method):

```
Stouffer.test <- function(p, w) {
    # p is a vector of p-values
    if (missing(w)) {
        w <- rep(1, length(p))/length(p)
    } else {
        if (length(w) != length(p))
            stop("Length of p and w must equal!")
    }
    Zi <- qnorm(1-p)
    Z <- sum(w*Zi)/sqrt(sum(w^2))
    p.val <- 1-pnorm(Z)
    return(c(Z = Z, p.value = p.val))
}
```

# Identification of methylated blocks via hidden Markov model

As described in the main text, runs of methylated cytosines were identified by applying a
hidden Markov model (HMM) to segment the genome in methylated and unmethylated
blocks. We used the p-values obtained from Stouffer’s method above as a response
variable to identify methylated blocks. In addition, p-values were converted to a discrete
variable. By visually inspecting the output (as for example in Supplementary Figure 11)
we found this strategy to be the most reliable in identifying runs of highly methylated
cytosines. In contrast, using as input variable the percentage methylation or the raw or log transformed p-values,
seemed to produce less consistent and credible results. Obviously, an experimental validation of the identified blocks would be necessary to
conclusively prove the reliability of the model(s). Nevertheless, the identified blocks can
serve as starting point in this direction and to eventually gain insights into their biological
relevance.

The HMM model was fitted with the R library [RHmm](https://cran.r-project.org/src/contrib/Archive/RHmm/) and data processing was
performed in R. First, p-values were converted to a discrete scale from 0 to 3 to reflect the likelihood of methylation:

```
recodePvals<- function(p){
    if (is.na(p)){
       o<- NA
    } else if (p > 0.1 ){
       o<- 0
    } else if (p <= 0.1 & p > 0.05) {
       o<- 1
    } else if (p <= 0.05 & p > 0.001) {
       o<- 2
    } else {
       o<- 3
    }
       return(o)
    }
```

The recoded p-values where then processed one chromosome at a time if the numbers
of cytosines in the chromosome does not exceed 40000. If a single chromosome
contains more than 40000 cytosines then the vector of p-values was split in the fewest
number of equal chunks each not exceeding 40000. This splitting was necessary to
contain the amount of computer memory required and to overcome a built-in limit in the RHmm library.
The HMM model was fitted with function HMMFit :

```
hmmfit<- HMMFit(obs, nStates= 2, dis= 'DISCRETE')
vit<- viterbi(hmmfit, obs)
```

Where obs is a vector of recoded p-values belonging to the same chromosome or
chunks from the same chromosome. As shown, two states were assumed and they were
decoded via Viterbi algorithm.

# Observed and expected duplication rate

Duplication rates in ReBuilT and PCR-BS libraries were measured by extracting read mate 1 only and downsampling
to the size of the smallest library. Downsampling and duplicate detetction were performed with the [Picard package](http://broadinstitute.github.io/picard/)

```
REFCOUNT=`samtools view -c -F 2944 grm038_BS_plasb2_AD16.clean.pb11.bam` # Smallest library: 8933900 reads

java -XX:-UseGCOverheadLimit -Xmx9g -jar DownsampleSam.jar I=$inbam O=$bname.ds.bam P=\$p R=1234 VALIDATION_STRINGENCY=SILENT
java -jar -Xmx2g MarkDuplicates.jar I=$bname.ds.bam O=/dev/null M=$bname.md.txt AS=true VALIDATION_STRINGENCY=SILENT
```

Where `$inbam` is a bam file containing only read mate 1 mapping to the _Plasmodium_ genome.

The expected duplication rate was calculated under the assumption that reads map randomly to the genome and
therefore the number of duplicates follows a Poisson distribution (see also the Illumina technical note [Estimating Sequencing Coverage](http://www.illumina.com/documents/products/technotes/technote_coverage_calculation.pdf)
and Lander ES, Waterman MS., 1988, Genomic mapping by fingerprinting random clones: a mathematical analysis, Genomics). Given the genome size of _P. berghei_, approximately 18 Mb, and
8930000 single end reads mapped the expetced duplication rate was calculated in R as follows:

```R
G<- 18000000
N<- 8930000 / 2 ## /2 because F and R strand are different 

depth<- 0:10
expect<- data.frame(
    pct_pos= dpois(depth, N/G),
    depth= depth    
)
expect$npos<- expect$pct_pos * G
expect$nreads<- expect$npos * expect$depth 

pct_usable<- sum(expect$npos[expect$depth > 0]) / N
pct_dups<- 1 - pct_usable
```

