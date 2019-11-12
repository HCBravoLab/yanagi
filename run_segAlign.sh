#!/bin/bash -x

SRC=/PATH/TO/rapmap
SEGSFASTA=/PATH/TO/segs.fa.meta
INDEX=/PATH/TO/segs_quasiindex

SAMPLES='sample1 sample2'
SAMPLEPATH=/PATH/TO/readsFastaFiles

OUTPUT=output/

mkdir $OUTPUT

for s in $SAMPLES
do
	echo $s
	# Single-end Use
	#python yanagi.py align -ref $SEGSFASTA -o $OUTPUT${s}_segcounts.txt -cmd1 "$SRC quasimap -i $INDEX -r ${SAMPLEPATH}${s}.fa"

	# Paired-end Use
	python yanagi.py align -ref $SEGSFASTA -o $OUTPUT${s}_segcounts.txt -cmd1 "$SRC quasimap -i $INDEX -r ${SAMPLEPATH}${s}_1.fa" -cmd2 "$SRC quasimap -i $INDEX -r ${SAMPLEPATH}${s}_2.fa"
done
