#!/bin/bash

for file in `ls -1 ./rnaseq`
do
	hisat2 -q -p 8 -x mm_index -U ./rnaseq/$file --no-softclip --known-splicesite-infile mss.txt \
		| samtools sort -@ 8 | samtools view -@ 8 -b > ./bams/${file::-5}sorted.bam
done
