#!/bin/bash

####################################################
# SETUP THE ENVIRONMENT AS NECESSARY TO USE GCC>=9 #
####################################################
GCChome=/mnt/data02/software/GCC/gcc-9
export PATH=$GCChome/bin:$PATH
export INCLUDE=$GCChome/include:$INCLUDE
export LD_LIBRARY_PATH=$GCChome/lib:$GCChome/lib64:$LD_LIBRARY_PATH
export MANPATH=$GCChome/man:$MANPATH
####################################################

###echo "Time stretching constant = "
###read timeAlpha
timeAlpha=1.0

fstring=$(echo "frequency ${f1} ${f2} ${f3} ${f4}" | xargs) 
cstring=$(echo "calculate ${R1} ${R2} ${R3} ${R4}" | xargs) 

cat > ./6_nmr_relax_data.inp << EOF
nucleus ${nucleus}
bondLength ${rXH}
deltaCSA ${CSA}
${fstring}
${cstring}
EOF

# Read trj parameters if in "auto" mode
if [ "$nSteps" -lt 0 ] || [ "$nDump" -lt 0 ] || [ "$(echo "$dt < 0.0" | bc -l)" -eq 1 ]
then
	while IFS=' ' read -r key v1; do
	        if [ "$key" == "nSteps" ]
		then
                	nSteps=$v1
	        elif [ "$key" == "nDump" ]
	        then
        	        nDump=$v1
		elif [ "$key" == "nChunks" ]
		then
			nRep=$v1
		fi
	done < 3_trjEstimatedInput.dat

	let ntrj=${nSteps}/${nDump}
fi

# Determine number of atoms
read -r nAtoms < 1_${projectName}-min.hes
let nAtoms=$nAtoms/3

# Which type of probe?
if [ "$nucleus" == "15N1H" ]
then
	probe=1
	baseName="5_probe_NH_*"
elif [ "$nucleus" == "13C1H" ]
then
	probe=2
	baseName="5_probe_CH_*"
else
	probe=3
	baseName="5_probe_CH2_*"
fi

# Run
nf=$(find . -name "$baseName" |wc -l)

cat > ./6_nmr_trj_data.inp <<-EOF
${nAtoms}
$ntrj
$nRep
./4_trj.xyz
$nExp
$timeAlpha
${probe}
${nf}
EOF

if [ "$nf" -gt 0 ]
then
	for i in $baseName; do
		{
		    read
		    IFS=' ' read key a1 a2 a3 a4
		} < ${i}

		echo "${a1} ${a2} ${a3} ${betaCSA}" >> ./6_nmr_trj_data.inp
	done
fi

date

echo ""
echo "Calculating D200 functions"
echo ""
$BDNMRHOME/nmr_relaxation/getOmega 6_nmr_relax_data.inp < 6_nmr_trj_data.inp >& 6_D200.out
date

echo ""
echo "Calculating NMR relaxation data"
echo ""
$BDNMRHOME/nmr_relaxation/relax 6_nmr_relax_data.inp < 6_nmr_trj_data.inp >& 6_nmr.out
date
