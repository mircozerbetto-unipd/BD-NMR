#!/bin/bash
###########################
# BUILD THE CBLAS LIBRARY #
###########################
cd ${LIBHOME}
tar xfz CBLAS.tgz
cd CBLAS/
cat ./Makefile.LINUX ../../conf/make.COMPILERS.CPU ../../conf/make.in > ./Makefile.in
make alllib
cp ./src/cblas_LINUX.a ./libcblas_LINUX.a
# Check if CBLAS has been built
log=`find -name "libcblas_LINUX.a"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the CBLAS library (lib/CBLAS/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
###########################
# BUILD THE FFTW3 LIBRARY #
###########################
cd ${LIBHOME}/
tar zxf fftw-3.3.tar.gz
cd fftw-3.3
./configure
make -j4
# Check if FFTW3 has been built
cd ./.libs/
log=`find -name "libfftw3.a"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the FFTW3 library (lib/fftw3.3/.libs/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
#####################################
# BUILD THE EULERQUATERNION LIBRARY #
#####################################
cd ${LIBHOME}/
tar zxf EulerQuaternion.tgz
cd EulerQuaternion/
make
# Check if EulerQuaternion has been built
log=`find -name "libEulerQuaternion.a"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the EulerQuaternion library (lib/EulerQuaternion/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
############################
# BUILD THE WIGNER LIBRARY #
############################
cd ${LIBHOME}/
tar zxf Wigner.tgz
cd Wigner
make
# Check if Wigner has been built
log=`find -name "libWigner.a"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the Wigner library (lib/Wigner/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
##############################
# BUILD THE CMINPACK LIBRARY #
##############################
cd ${LIBHOME}/
tar zxf cminpack-1.3.11.tar.gz
cd cminpack-1.3.11/
make -j4
# Check if Wigner has been built
log=`find -name "libcminpack.a"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the CMINPACK library (lib/cminpack-1.3.11/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
################################
# BUILD NMR_RELAXATION PACKAGE #
################################
cd ${BDNMRHOME}/nmr_relaxation/
make -j4 all
# Check if nmr has been built
log=`find -name "nmr"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the NMR package (nmr_relaxation/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
########################
# BUILD TINKER PACKAGE #
########################
cd ${BDNMRHOME}/external_packages
tar zxf tinker-8.11.3.tgz
cd tinker/source/
oldline="LINKFLAGS = "
modline="LINKFLAGS = \$(OPTFLAGS) $GFORTRANLIB"
sed -i "s@.*$oldline.*@$modline@" "Makefile"
make -j4 all
# Check if Tinker has been built
log=`find -name "minimize.x"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the TINKER package (external_packages/tinker). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
make install
