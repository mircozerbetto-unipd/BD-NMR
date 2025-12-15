#!/bin/bash
###########################################################
# Use this script to run Z2CART                           #
#                                                         #
# Please, note that if the build directory is             #
# moved to a different path after the program has been    #
# compiled, the paths in the following commands must      #
# be changed accordingly to the new location of the files #
###########################################################
export PATH=/usr/local/cuda-7.5:/usr/local/cuda-7.5:/usr/local/cuda-7.5/bin:/home/mirco/work/software/amber22/bin:/home/mirco/work/software/amber22/bin:/usr/share/Modules/bin:/usr/local/bin:/usr/local/sbin:/usr/bin:/usr/sbin:/bin:/sbin:/home/mirco/.local/bin:/home/mirco/bin
export INCLUDE=/usr/local/cuda-7.5/include:/usr/local/cuda-7.5/include:/usr/local/cuda-7.5/include
export LD_LIBRARY_PATH=/usr/local/cuda-7.5:/usr/local/cuda-7.5:/usr/local/cuda-7.5/lib64:/home/mirco/work/software/amber22/lib:/home/mirco/work/software/amber22/lib
export ZMATLIB_HOME=/home/mirco/work/development/BD-NMR/lib/zmatlib
/home/mirco/work/development/BD-NMR/z-to-Cartesian/z2cart < $1
