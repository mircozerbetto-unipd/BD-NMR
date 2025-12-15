#!/bin/bash

read -r nz < 1_${projectName}-min.hes
let nz=$nz-6

cat > ./3_sfb.inp << EOF
nSteps ${nSteps}
nDump ${nDump}
nIntCoor ${nz}
RCDiff ./2_rcdiff_z.dat
dt ${dt}
seed ${seed}
EOF

##$BDNMRHOME/SFB/sfb_RK 3_sfb.inp  # IMPLEMENTATION OF RK SCHEME HAS TO BE CHECKED
$BDNMRHOME/SFB/sfb 3_sfb.inp  # USE EULER FOR CAUTION

mv trj.dat 3_trj_rc.dat

if [ "$nSteps" -lt 0 ] || [ "$nDump" -lt 0 ] || [ "$dt" -lt 0 ]; then
	mv trjEstimatedInput.dat 3_trjEstimatedInput.dat
fi

