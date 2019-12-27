#!/bin/bash
  
for file in `ls -1 ./bedGraph`
do
	./bedGraphToBigWig ./bedGraph/$file hg38.chrom.sizes ./bigWig/${file::-8}bigWig
done

