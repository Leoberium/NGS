#!/bin/bash

for file in `ls -1 ./covs/*.gtf`
do
	echo "$(basename ${file::-8}) $file" >> sample_list.txt
done
