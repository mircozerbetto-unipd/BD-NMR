#!/bin/sh
################################################
##                   ______                   ##
#                   |      |                   #
#                   |  GO  |                   #
#                   |______|                   #
##                                            ##
# Change the following variables as needed     #
##                                            ##
################################################
############
# HARDWARE #
############
# Change list according to hardware ##############
GPU="Hopper"
##################################################
###############
# ENVIRONMENT #
###############
# Change according to environment ###############
export CUDADIR=/cm/shared/apps/cuda/12.6/
export LAPACKDIR=/cm/shared/apps/hpc_sdk/Linux_x86_64/24.9/compilers/lib/
export OPENBLASDIR=/cm/shared/apps/openblas/0.3.18/
####################
# GFORTRAN LIBRARY #
####################
GFORTRANLIB=/cm/local/apps/gcc/13.1.0/lib64/libgfortran.so
################################################
##                   ______                   ##
#                   |      |                   #
#                   | STOP |                   #
#                   |______|                   #
##                                            ##
# Do not touch the following, unless necessary #
##                                            ##
################################################
################
# WHERE WE ARE #
################
BDNMRHOME=$(pwd)
LIBHOME=$(pwd)/lib/
######################
# Update environment #
######################
export PATH=${CUDADIR}:$PATH
export INCLUDE=${CUDADIR}/include:$INCLUDE
export LD_LIBRARY_PATH=${CUDADIR}:$LD_LIBRARY_PATH
# Determine CUDA architecture(s)
CUDA_ARCH_=""
# Check for Kepler
case "$GPU" in
  *Kepler*)
    CUDA_ARCH_="$CUDA_ARCH_ -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35"
    ;;
esac

# Check for Maxwell
case "$GPU" in
  *Maxwell*)
    CUDA_ARCH_="$CUDA_ARCH_ -gencode arch=compute_50,code=sm_50"
    ;;
esac

# Check for Pascal
case "$GPU" in
  *Pascal*)
    CUDA_ARCH_="$CUDA_ARCH_ -gencode arch=compute_60,code=sm_60"
    ;;
esac

# Check for Volta
case "$GPU" in
  *Volta*)
    CUDA_ARCH_="$CUDA_ARCH_ -gencode arch=compute_70,code=sm_70"
    ;;
esac

# Check for Turing
case "$GPU" in
  *Turing*)
    CUDA_ARCH_="$CUDA_ARCH_ -gencode arch=compute_75,code=sm_75"
    ;;
esac

# Check for Ampere
case "$GPU" in
  *Ampere*)
    CUDA_ARCH_="$CUDA_ARCH_ -gencode arch=compute_80,code=sm_80"
    ;;
esac

# Check for Ada
case "$GPU" in
  *Ada*)
    CUDA_ARCH_="$CUDA_ARCH_ -gencode arch=compute_89,code=sm_89"
    ;;
esac

# Check for Hopper
case "$GPU" in
  *Hopper*)
    CUDA_ARCH_="$CUDA_ARCH_ -gencode arch=compute_90,code=sm_90 -gencode arch=compute_90a,code=sm_90a"
    ;;
esac

if [ -z "$CUDA_ARCH_" ]; then
	echo "Cannot determine CUDA architecture for GPU ${GPU}, please fix the GPU variable to be one (or more) of: Kepler, Maxwell, Pascal, Volta, Ampere, Ada, Hopper"
	exit 1
else
	echo "Compiling CUDA architecture(s): \"$CUDA_ARCH_\""
fi

###########################################################################
# make.in CONTAINS INFORMATION FOR LIB AND INCLUDE, AND CUDA ARCHITECTURE #
###########################################################################
cat > ./conf/make.in << EOF
LL=-L$BDNMRHOME/lib/CBLAS -L$BDNMRHOME/lib/cqp -L$BDNMRHOME/lib/EulerQuaternion/ -L$BDNMRHOME/lib/Wigner/ -L$BDNMRHOME/lib/magma-2.9.0/lib/ -L$BDNMRHOME/lib/zmatlib/ -L$BDNMRHOME/lib/fftw-3.3/.libs/ -L$BDNMRHOME/lib/cminpack-1.3.11/ -L${LAPACKDIR} -L$OPENBLASDIR/lib
II=-I$BDNMRHOME/lib/CBLAS/src -I$BDNMRHOME/lib/cqp -I$BDNMRHOME/include/ -I$BDNMRHOME/lib/EulerQuaternion/ -I$BDNMRHOME/lib/magma-2.9.0/include/ -I$BDNMRHOME/lib/zmatlib/include/ -I$BDNMRHOME/lib/fftw-3.3/api/ -I$BDNMRHOME/lib/cminpack-1.3.11/
CUDA_ARCH=$CUDA_ARCH_
EOF
#############################
# INSTALL CPU-ONLY PACKAGES #
#############################
cd $BDNMRHOME
source ./install_cpu_packages.sh
############################
# INSTALL CPU-GPU PACKAGES #
############################
cd $BDNMRHOME
source ./install_gpu_packages.sh
#############################
# OUTPUT ENVIRONMENT SCRIPT #
#############################
cd $BDNMRHOME
cat > ./run_scripts/set_BD-NMR_env.sh << EOF
#!/bin/bash
###########################################################
# Use this script to set the environment variables.       #
#                                                         #
# Please, note that if this build directory is moved to   #
# to a different path after the program has been          #
# compiled, the paths in the following commands must      #
# be changed accordingly to the new location of the files #
###########################################################
export PATH=${CUDADIR}:${BDNMRHOME}/salem/salem-tools/pdbtk/bin:$PATH
export INCLUDE=${CUDADIR}/include:$INCLUDE
export LD_LIBRARY_PATH=${CUDADIR}:$LAPACKDIR:$OPENBLASDIR/lib:${LIBHOME}/magma-2.9.0/lib:$LD_LIBRARY_PATH
export ZMATLIB_HOME=${LIBHOME}/zmatlib/
export BDNMRHOME=${BDNMRHOME}
EOF
chmod u+x ./run_scripts/set_BD-NMR_env.sh
echo 'Welcome to...'
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
echo 'To run the calculations, source first the ./run_scripts/set_BD-NMR.env.sh script.'
echo ''
echo 'Next, as the input geometry and configuration files are ready (see the examples)'
echo 'run the ./run_scripts/BDNMR.sh script to handle the calculation workflow.'
echo ''
