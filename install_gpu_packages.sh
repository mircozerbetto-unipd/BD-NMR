#!/bin/bash
###########################
# BUILD THE MAGMA LIBRARY #
###########################
cd ${LIBHOME}
tar zxf magma-2.9.0.tgz
cd magma-2.9.0/
cp ./make.inc-examples/make.inc.openblas ./make.inc
echo GPU_TARGET=${GPU} >> make.inc
make -j4 lib
# Check if MAGMA has been built
cd lib
log=`find -name "libmagma.a"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the MAGMA library (lib/magma-2.9.0/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
########################
# BUILD THE SESTO TOOL #
########################
cd ${LIBHOME}
tar zxf zmatlib.tgz
cd zmatlib/sesto_1.0/
make
# Check if SESTO has been built
log=`find -name "sesto"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build SESTO (lib/zmatlib/sesto_1.0/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
#############################
# BUILD THE ZMATLIB LIBRARY #
#############################
cd ${LIBHOME}/zmatlib/
make
# Check if ZMAT has been built
log=`find -name "libzmat.a"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the ZMAT library (lib/zmatlib/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
#########################
# BUILD THE CQP LIBRARY #
#########################
cd ${LIBHOME}
tar zxf cqp.tgz
cd cqp/
make
# Check if CQP has been built
log=`find -name "libcqp++.a"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the CQP library (lib/cqp/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
###############
# BUILD SALEM #
###############
cd ${BDNMRHOME}/salem
SALEMEXE="salem-srb-diff_GPU"
###########################################
# PREPARATION OF THE make.COMPILERS FILE ##
###########################################
cat ../conf/make.COMPILERS.CPU-GPU ../conf/make.in > ./Makefile.in
cd src/
make -j4 ${SALEMEXE}
# Check if SALEM has been built
log=`find -name "salem-srb-diff_GPU"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the SALEM executable (src/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
cd ../
cat > ./salem << EOF
#!/bin/bash
###########################################################
# Use this script to run SALEM-SRB                        #
#                                                         #
# Please, note that if the SALEM build directory is       #
# moved to a different path after the program has been    #
# compiled, the paths in the following commands must      #
# be changed accordingly to the new location of the files #
###########################################################
export PATH=${CUDADIR}:$PATH
export INCLUDE=${CUDADIR}/include:$INCLUDE
export LD_LIBRARY_PATH=${CUDADIR}:$LD_LIBRARY_PATH
export ZMATLIB_HOME=${ZMATPATH}
${BDNMRHOME}/salem/src/${SALEMEXE} \$1
EOF
chmod u+x ./salem
#############
# END SALEM #
#############
echo ""
echo "====================================================================================================================="
echo "SALEM succesfully compiled. The program can be run as:"
echo ""
echo "./salem input.dat"
echo ""
echo "Please note that 'salem' is a bash script pointing to the actual position of this directory:"
echo ""
echo $(pwd)
echo ""
echo "If this directory is moved to another location, the 'salem' script should be modified accordingly to the new path."
echo "====================================================================================================================="
echo ""
#####################
# BUILD SALEM TOOLS #
#####################
cd ${BDNMRHOME}/salem/salem-tools/pdbtk/
cat > ./config.mk << EOF
PDBTKPATH=${BDNMRHOME}/salem/salem-tools/pdbtk/
INSTALLPATH=/usr/local
FC=gfortran
FFLAGS= -Wall
EOF
make -j4
# Check if pdbtk has been built
cd bin
log=`find -name "pdbaf"`
if [ -z "$log" ]; then
        echo ""
        echo "ERROR: it was not possible to build the pdbtk executables (salem/salem-tools/pdbtk). Please, check the Makefile therein."
        echo ""
        echo "Build script stopped"
        exit
fi
log=`find -name "tinker2salemhes"`
if [ -z "$log" ]; then
        echo ""
        echo "ERROR: it was not possible to build the pdbtk executables (salem/salem-tools/pdbtk). Please, check the Makefile therein."
        echo ""
        echo "Build script stopped"
        exit
fi

cd ${BDNMRHOME}/salem/salem-tools/probe-finder/
make
# Check if probe-finder has been built
log=`find -name "probe-finder"`
if [ -z "$log" ]; then
        echo ""
        echo "ERROR: it was not possible to build the probe-finder executable (salem/salem-tools/probe-finder/). Please, check the Makefile therein."
        echo ""
        echo "Build script stopped"
        exit
fi
###################
# BUILD Z-To-Cart #
###################
cd ${BDNMRHOME}/z-to-Cartesian
make
# Check if z-to-Cartesian has been built
log=`find -name "z2cart"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the z2cart executable (z-to-Cartesian/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
#############
# BUILD SFB #
#############
cd ${BDNMRHOME}/SFB
make
# Check if SFB has been built
log=`find -name "sfb"`
if [ -z "$log" ]; then
	echo ""
	echo "ERROR: it was not possible to build the SFB executable (SFB/). Please, check the Makefile therein."
	echo ""
	echo "Build script stopped"
	exit
fi
