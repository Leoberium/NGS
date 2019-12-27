#!/bin/bash

for file in `ls -1 ./bams/*.bam`
do
	file=$(basename $file)
	echo $file
	stringtie ./bams/$file -o ./gtfs/${file::-10}gtf -G GRCm38.98.filtered.gtf -p 8
done
