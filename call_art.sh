#!/bin/bash

genome=$1
read_length=$2
cover=$3
frag_size=$4
stdev=$5
output=$6
pp=$7
ql=$8
qu=$9

if [ $pp == "p" ]
then
	./art_illumina -na -q -p -i $genome -l $read_length -f $cover -m $frag_size -s $stdev -o ${output}. -qL $ql -qU $qu
else
	./art_illumina -na -q -i $genome -l $read_length -f $cover -o ${output} -qL $ql -qU $qu
fi