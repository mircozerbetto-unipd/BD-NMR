#!/bin/bash

# Read trj parameters if in "auto" mode
if [ "$nSteps" -lt 0 ] || [ "$nDump" -lt 0 ] || [ "$dt" -lt 0 ]
then
	while IFS=' ' read -r key v1; do
	        if [ "$key" == "nSteps" ]
		then
                	nSteps=$v1
	        elif [ "$key" == "nDump" ]
	        then
        	        nDump=$v1
	        elif [ "$key" == "dt" ]
	        then
        	        dt=$v1
		elif [ "$key" == "nChunks" ]
		then
			nRep=$v1
		fi
	done < 3_trjEstimatedInput.dat

	let ntrj=${nSteps}/${nDump}
fi

# Read the number of coordinates
read -r nz < 1_${projectName}-min.hes
let nz=$nz-6

# Prepare input
cat > ./4_zToCart.inp << EOF
${nz}
${ntrj}
./1_${projectName}-min.pdb
${id1} ${id2} ${id3} ${id4}
./2_T.dat
./2_sqKmat.dat
./3_trj_rc.dat
EOF

# Run
$BDNMRHOME/z-to-Cartesian/z2cart < 4_zToCart.inp

mv zmatrix.zmt 4_zmatrix.zmt
mv q.dat 4_q.dat
mv trj.xyz 4_trj.xyz


