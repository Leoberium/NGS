#!/bin/bash

for file in `ls -1 ./rebams`
do
        echo $file
        samtools index ./rebams/$file
done
