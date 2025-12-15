#!/bin/bash

$BDNMRHOME/salem/salem-tools/probe-finder/probe-finder 1_${projectName}-min.pdb ${id1} ${id2} ${id3} ${id4}

rm -f zmatrix.zmt

nf=$(find . -name "probe_NH_*" |wc -l)
if [ "$nf" -gt 0 ]
then
	for i in probe_NH_*; do
		mv $i 5_${i}
	done
fi

nf=$(find . -name "probe_CH_*" |wc -l)
if [ "$nf" -gt 0 ]
then
	for i in probe_CH_*; do
		mv $i 5_${i}
	done
fi

nf=$(find . -name "probe_CH2_*" |wc -l)
if [ "$nf" -gt 0 ]
then
	for i in probe_CH2_*; do
		mv $i 5_${i}
	done
fi

