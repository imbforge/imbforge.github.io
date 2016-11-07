# Merge bigWig track files 

Generating a single track summarizing NGS coverage from several replicates is a common request among scientists.

Usually, the scientist will be happy if the bioinformatician just merges together the *coverage per million* **normalized** tracks. This is a trivial task handled with the tools from the guys of the (https://genome.ucsc.edu/util.html UCSC Genome Bioinformatics Group) \(aka the kent utilities\).

All the magic is done with the `bigWigMerge` tool, which puts together the signal from the several bigwig tracks, (bedGraph output), to be eventually converted back to bigwig with the `bedGraphToBigWig` tool. The *trick* here is to divide to amount of signal in each position by the number of replicates. Something like this will do the job:

```bash
bigWigMerge sample1_rep1.bw sample1_rep2.bw sample1_rep3.bw stdout \
    | awk -v NORM=1/3 '$4=(NORM)*$4' \
    | LC_COLLATE=C sort -k1,1 -k2,2n > sample1.bed
bedGraphToBigWig sample1.bed chr.sizes sample1.bw && rm sample1.bed
```

Chromosome sizes can be retrieved from the UCSC Genome Browser tools with the tool `fetchChromSizes`, also from the kent utilities. Alternatively, `samtools idxstats sample1.bam | cut -f1-2 > ./chr.sizes` will work.

### Requirements

* `bedGraphToBigWig` from the UCSC Genome Browser tools.
* the cromosome sizes, which can be retrieved from the UCSC Genome Browser tools with the fetchChromSizes. Alternatively, `samtools idxstats sample.bam | cut -f1-2 > ./chr.sizes`.

### Source

The following script loops over serveral samples defined in a array, and identifies the replicates assuming they follow a pattern (matched by the ?? wildcards):

```bash
#!/bin/bash
UCSC=/opt/ucsc/latest/
INPUTDIR=/project/tracks
SAMPLES=(sample1_?? sample2_?? sample3_??)
CHRSIZES=/project/annotation/chr.sizes  # got with: samtools idxstats sample1.bam | cut -f1-2 > ./chr.sizes

cd $INPUTDIR

for SAMPLE in ${SAMPLES[@]}; do
    echo $SAMPLE

    REP=$(ls ${SAMPLE}*.bw) # get list of replicates to merge
    N=$(echo $REP | wc -w)  # how many?
    NORM=$(echo "1 / $N" | bc -l)   # average each signal by the number of replicates
    OUT=$(echo $SAMPLE | sed 's/??_//')

    ${UCSC}/bigWigMerge $REP stdout \                    # merge replicates and output to stdout
        | awk -v NORM=$NORM '$4=NORM*$4' \               # normalize (4th field contains the coverage)
        | LC_COLLATE=C sort -k1,1 -k2,2n > ${OUT}.bed    # sort by chr-pos
    ${UCSC}/bedGraphToBigWig ${OUT}.bed $CHRSIZES ${OUT}.bw && rm ${OUT}.bed
done
```
