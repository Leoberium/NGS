#!/bin/bash

for file in `ls -1 ./rebams/*.bam`
do
	file=$(basename $file)
	echo $file
	stringtie ./rebams/$file -B -e -G merged.gtf -p 8 > ./covs/${file::-10}cov.gtf
done
