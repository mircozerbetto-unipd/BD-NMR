#!/bin/bash

cat > ./2_salem.inp << EOF
pdb ./1_${projectName}-min.pdb
refAtoms ${id1} ${id2} ${id3} ${id4}
Reff ${Reff}
C ${C}
viscosity ${eta}
temperature ${T}
hydroint 1
kmatrix hessian ./1_${projectName}-min.hes
activeatoms 0
ignore_negative_hess_eigs ${ignoreNegEigs}
EOF

$BDNMRHOME/salem/src/salem-srb-diff_GPU 2_salem.inp

mv zmatrix.zmt 2_zmatrix.zmt
mv rebuild.xyz 2_rebuild.xyz
mv sqKmat.dat 2_sqKmat.dat
mv rcdiff_z.dat 2_rcdiff_z.dat
mv E.dat 2_E.dat
mv T.dat 2_T.dat
mv times.dat 2_times.dat

