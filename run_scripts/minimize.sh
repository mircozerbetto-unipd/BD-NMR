#!/bin/bash

# Commented out to use implicit solvation
#cat > ./tinker.key << EOF
#verbose
#parameters ${BDNMRHOME}/external_packages/tinker/params/${ff}
#digits 8
#EOF

# set filename
sys=$projectName

# origin dir
wd=`pwd`

# Adjust atom types based on first conversion
cat > ./edit.inp << EOF
1_${sys}-min.xyz
10
${sys}.xyz


EOF

# generate .pdb, .gra and .hes files in AF

echo "Hessian will be rotated to AF frame defined from atoms: $id1, $id2, $id3, and $id4"

# minimize, rotate, get gradient and hessian
echo ""
echo "Running energy minimization and Hessian calculation"
echo ""
echo "minimize ${sys}.xyz  $gradTol -k ./tinker.key"
echo ""
minimize ${sys}.xyz $gradTol -k ./tinker.key

# clean
mv ${sys}.xyz_2 1_${sys}-min.xyz
# Adjust atom types based on first conversion
xyzedit < edit.inp
mv 1_${sys}-min.xyz_2 1_${sys}-min.xyz

# Convert to SALEM
echo ""
echo "Convert results to SALEM format\n"
echo ""

# Obtain PDB from XYZ
echo "Running: xyzpdb 1_${sys}-min.xyz -key ./tinker.key"
xyzpdb -key ./tinker.key 1_${sys}-min.xyz

# Guess PDB elements
pdbguesselement 1_${sys}-min.pdb guessed-${sys}-min.pdb
mv guessed-${sys}-min.pdb 1_${sys}-min_LF.pdb

# Rotation to AF frame
pdbaf 1_${sys}-min_LF.pdb 1_${sys}-min.pdb $id1 $id2 $id3 $id4
rm 1_${sys}-min_LF.pdb

# Convert to XYZ
pdbxyz -key ./tinker.key 1_${sys}-min.pdb
mv 1_${sys}-min.xyz_2 1_${sys}-min.xyz
# Adjust atom types based on first conversion
xyzedit < edit.inp
mv 1_${sys}-min.xyz_2 1_${sys}-min.xyz

# Recalculate Hessian in AF
echo ""
echo "Test gradient"
echo ""
printf 'y\nn\nn' | testgrad 1_${sys}-min.xyz > 1_${sys}-min.tgra
echo ""
echo "Test Hessian"
echo ""
printf 'y\nn\nn' | testhess 1_${sys}-min.xyz
mv 1_${sys}-min.hes 1_${sys}-min.thes
echo "\nConvert to SALEM\n"
tinker2salemgra 1_${sys}-min.pdb 1_${sys}-min.tgra 1_${sys}-min.gra
tinker2salemhes 1_${sys}-min.pdb 1_${sys}-min.thes 1_${sys}-min.hes

# set occupancies to 1.0
pdboccupy 1_${sys}-min.pdb ${sys}-occ.pdb
mv ${sys}-occ.pdb 1_${sys}-min.pdb

# adjust the name of the .seq file
mv ${sys}.seq 1_${sys}.seq
