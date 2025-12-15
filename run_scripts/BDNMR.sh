#!/bin/bash

echo ''
echo ''
echo '__/\\\\\\\\\\\\\____/\\\\\\\\\\\\_______________/\\\\\_____/\\\__/\\\\____________/\\\\____/\\\\\\\\\_____        '
echo ' _\/\\\/////////\\\_\/\\\////////\\\____________\/\\\\\\___\/\\\_\/\\\\\\________/\\\\\\__/\\\///////\\\___       '
echo '  _\/\\\_______\/\\\_\/\\\______\//\\\___________\/\\\/\\\__\/\\\_\/\\\//\\\____/\\\//\\\_\/\\\_____\/\\\___      '
echo '   _\/\\\\\\\\\\\\\\__\/\\\_______\/\\\___________\/\\\//\\\_\/\\\_\/\\\\///\\\/\\\/_\/\\\_\/\\\\\\\\\\\/____     '
echo '    _\/\\\/////////\\\_\/\\\_______\/\\\___________\/\\\\//\\\\/\\\_\/\\\__\///\\\/___\/\\\_\/\\\//////\\\____    '
echo '     _\/\\\_______\/\\\_\/\\\_______\/\\\___________\/\\\_\//\\\/\\\_\/\\\____\///_____\/\\\_\/\\\____\//\\\___   '
echo '      _\/\\\_______\/\\\_\/\\\_______/\\\____________\/\\\__\//\\\\\\_\/\\\_____________\/\\\_\/\\\_____\//\\\__  '
echo '       _\/\\\\\\\\\\\\\/__\/\\\\\\\\\\\\/_____________\/\\\___\//\\\\\_\/\\\_____________\/\\\_\/\\\______\//\\\_ '
echo '        _\/////////////____\////////////_______________\///_____\/////__\///______________\///__\///________\///__'
echo '                                                                                                                  '
echo ''
echo ''

##############
# ADJUST ENV #
##############

export PATH=${BDNMRHOME}/external_packages/tinker/bin/:$PATH

#################
# SOME DEFAULTS #
#################

gradTol=0.01
seed=-1
nucleus=15N1H
nSteps=-1
nDump=-1
dt=-1
nRep=1
nExp=4
betaCSA=0
tStretch=1.0
ignoreNegEigs=0

####################
# PARSE INPUT FILE #
####################

inpFile=$1
while IFS=' ' read -r key v1 v2 v3 v4
do
	if [ "$key" == "project" ]
	then
		projectName=$v1
	elif [ "$key" == "ff" ]
	then
		ff=$v1
	elif [ "$key" == "refAtoms" ]
	then
		id1=$v1	
		id2=$v2	
		id3=$v3	
		id4=$v4	
	elif [ "$key" == "gradtol" ]
	then
		gradTol=$v1
	elif [ "$key" == "Reff" ]
	then
		Reff=$v1
	elif [ "$key" == "C" ]
	then
		C=$v1
	elif [ "$key" == "viscosity" ]
	then
		eta=$v1
	elif [ "$key" == "temperature" ]
	then
		T=$v1
	elif [ "$key" == "nSteps" ]
	then
		nSteps=$v1
	elif [ "$key" == "nDump" ]
	then
		nDump=$v1
	elif [ "$key" == "dt" ]
	then
		dt=$v1
	elif [ "$key" == "seed" ]
	then
		seed=$v1
	elif [ "$key" == "nucleus" ]
	then
		nucleus=$v1
	elif [ "$key" == "bondLength" ]
	then
		rXH=$v1
	elif [ "$key" == "deltaCSA" ]
	then
		CSA=$v1
	elif [ "$key" == "betaDipCSA" ]
	then
		betaCSA=$v1
	elif [ "$key" == "frequency" ]
	then
		f1=$v1
		f2=$v2
		f3=$v3
		f4=$v4
	elif [ "$key" == "calculate" ]
	then
		R1=$v1
		R2=$v2
		R3=$v3
		R4=$v4
	elif [ "$key" == "nChunks" ]
	then
		nRep=$v1
	elif [ "$key" == "nExp" ]
	then
		nExp=$v1
	elif [ "$key" == "tStretch" ]
	then
		tStretch=$v1
	elif [ "$key" == "ignore_negative_hess_eigs" ]
	then
		ignoreNegEigs=$v1
	fi
done < "$inpFile"

let ntrj=${nSteps}/${nDump}

#######################
# SELECT STEPS TO RUN #
#######################

echo "Please, select which step(s) to run"
OPTIONS="all minimization diffusion-tensor BD-trajectory z-to-Cartesian probe-finder NMR select-range exit"
select opt in $OPTIONS; do
	if [ "$opt" = "all" ]; then
		firstStep=1
		lastStep=6
		break
	elif [ "$opt" = "minimization" ]; then
	firstStep=1
		lastStep=1
		break
	elif [ "$opt" = "diffusion-tensor" ]; then
		firstStep=2
		lastStep=2
		break
	elif [ "$opt" = "BD-trajectory" ]; then
		firstStep=3
		lastStep=3
		break
	elif [ "$opt" = "z-to-Cartesian" ]; then
		firstStep=4
		lastStep=4
		break
	elif [ "$opt" = "probe-finder" ]; then
		firstStep=5
		lastStep=5
		break
	elif [ "$opt" = "NMR" ]; then
		firstStep=6
		lastStep=6
		break
	elif [ "$opt" = "select-range" ]; then
		echo "Select the numbers of first and last steps with 1 = minimization, 2 = diffusion, 3 = trajetory, 4 = convesion to Cartesian, 5 = find probes, and 6 = nmr"
		echo "First step: "
		read firstStep
		echo "Last step: "
		read lastStep
		break
	elif [ "$opt" == "exit" ]; then
		exit
	fi
done

if [ "$firstStep" == 1 ]; then

	echo ""
	echo "*************************"
	echo "** ENERGY MINIMIZATION **"
	echo "*************************"
	echo ""

	source "${BDNMRHOME}/run_scripts/minimize.sh"

fi

if [ "$firstStep" -le 2 ] && [ "$lastStep" -ge 2 ]; then

	echo ""
	echo "*****************************"
	echo "** SCALED DIFFUSION TENROR **"
	echo "*****************************"
	echo ""

	source "${BDNMRHOME}/run_scripts/salem.sh"

fi

if [ "$firstStep" -le 3 ] && [ "$lastStep" -ge 3 ]; then

	echo ""
	echo "*************************"
	echo "** BROWNIAN TRAJECTORY **"
	echo "*************************"
	echo ""

	source "${BDNMRHOME}/run_scripts/sfb.sh"

fi

if [ "$firstStep" -le 4 ] && [ "$lastStep" -ge 4 ]; then

	echo ""
	echo "**************************************"
	echo "** CONVERT TO CARTESIAN COORDINATES **"
	echo "**************************************"
	echo ""

	source "${BDNMRHOME}/run_scripts/rc2xyz.sh"

fi

if [ "$firstStep" -le 5 ] && [ "$lastStep" -ge 5 ]; then

	echo ""
	echo "*********************"
	echo "** FIND NMR PROBES **"
	echo "*********************"
	echo ""

	source "${BDNMRHOME}/run_scripts/probes.sh"

fi

if [ "$firstStep" -le 6 ] && [ "$lastStep" -ge 6 ]; then

	echo ""
	echo "************************"
	echo "** CALCULATE NMR DATA **"
	echo "************************"
	echo ""

	source "${BDNMRHOME}/run_scripts/nmr.sh"

fi

