#!/bin/bash
  
for file in `ls -1 ./bed`
do
	awk '$1 ~ /^chr[0-9YX][0-9]?$/ {print $1"\t"$2"\t"$3"\t"$5}' ./bed/$file \
		| sort -k1,1 -k2,2n > ./bedGraph/${file::-3}bedGraph
done

