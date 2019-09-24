# Merge track files

## bigWig tracks

Generating a single track summarizing NGS coverage from several replicates is a common request among scientists.

Usually, the scientist will be happy if the bioinformatician just merges together the *coverage per million* **normalized** tracks. This is a trivial task handled with the tools from the guys of the UCSC Genome Bioinformatics Group (aka [*the kent utilities*](https://genome.ucsc.edu/util.html)).

All the magic is done with the `bigWigMerge` tool, which puts together the signal from the several bigwig tracks, (bedGraph output), to be eventually converted back to bigwig with the `bedGraphToBigWig` tool. The *trick* here is to divide to amount of signal in each position by the number of replicates. Something like this will do the job:

```bash
bigWigMerge sample1_rep1.bw sample1_rep2.bw sample1_rep3.bw stdout \
    | awk -v NORM=1/3 '$4=(NORM)*$4' \
    | LC_COLLATE=C sort -k1,1 -k2,2n > sample1.bed
bedGraphToBigWig sample1.bed chr.sizes sample1.bw && rm sample1.bed
```

Chromosome sizes can be retrieved from the UCSC Genome Browser tools with the tool `fetchChromSizes`, also from the kent utilities. Alternatively, `samtools idxstats sample1.bam | cut -f1-2 > ./chr.sizes` will work.

### Requirements

* `bedGraphToBigWig` from the [kent utilities](https://genome.ucsc.edu/util.html).
* the cromosome sizes, which can be retrieved from the UCSC Genome Browser tools with the fetchChromSizes. Alternatively, `samtools idxstats sample1.bam | cut -f1-2 > ./chr.sizes`.

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

## BAM files

Sometimes the researcher might be interested in keeping information at read level, e.g. to explore splicing. In this case, usual formats to represent continuous data in tracks are of no use, and the bioinformatician has to fall back to BAM.

The main basic thing to consider is avoid overrepresentation of any library over the rest. Thus, we have to first calculate the library size of the replicates, and subsample from the original bam files to the smallest fraction before merging.

### Steps

* get the smallest library size of all replicates of the sample
* downsample the rest of the replicates to the smallest library
* merge the downsampled replicates together

### Requirements

* [samtools](http://www.htslib.org/)

### Source

The following script loops over serveral samples defined in a array, and identifies the replicates assuming they follow a pattern (matched by the ?? wildcards):

```bash
#!/bin/bash
SAMTOOLS=/opt/samtools/1.3/samtools
INPUTDIR=/project/mapped
SAMPLES=(sample1_?? sample2_?? sample3??)
SEED=666

cd $INPUTDIR

for SAMPLE in ${SAMPLES[@]}; do
    OUT=$(echo $SAMPLE | sed 's/??_//')
    echo "== $OUT =="

    # get minimum number of reads mapped
    REPS=() # will fill out later
    MAPPEDREADS=()  # will fill out later
    for REP in ${SAMPLE}*.bam; do
        M=$($SAMTOOLS flagstat $REP | grep mapped | cut -f1 -d" " | head -n1)
        MAPPEDREADS+=($M)
        REPS+=($REP)
        echo "$REP has $M reads"
    done
    MIN=$(echo ${MAPPEDREADS[@]} | tr " " "\n" | sort -n | head -n1)

    # downsample to the minimum
    REPSDOWNSAMPLED=()
    for i in $(seq 1 ${#REPS[@]}); do
        X=${REPS[$i-1]}
        REPSDOWNSAMPLED+=(${X%.bam}.subsampled.bam)
        if [ "$MIN" -eq "${MAPPEDREADS[$i-1]}" ]; then
            echo "not downsampling $X"
            $SAMTOOLS view -F4 -bh $X > ${REPSDOWNSAMPLED[$i-1]}
        else
            PROB=$(echo "$MIN / ${MAPPEDREADS[$i-1]}" | bc -l | cut -c2-3)
            echo "downsampling $X to ${PROB}% of reads"
            $SAMTOOLS view -s ${SEED}.${PROB} -F4 -bh $X > ${REPSDOWNSAMPLED[$i-1]}
        fi
    done

    # merge the downsampled bams
    echo "merging ${REPSDOWNSAMPLED[@]} into ${OUT}.bam"
    $SAMTOOLS merge ${OUT}.bam ${REPSDOWNSAMPLED[@]} && rm ${REPSDOWNSAMPLED[@]}
done
```

And happy sashimi ;-)

![happy sashimi!](figures/2016-11-04-Merge-track-files_fig1.png)
