#!/bin/sh
export ZMATLIB_HOME=../../../lib/zmatlib

probe-finder:
	g++ main.cpp rotomat.cpp searchNH.cpp searchCH.cpp searchCH2.cpp -I${ZMATLIB_HOME}/include -L${ZMATLIB_HOME} -lzmat -o probe-finder

clean:
	rm -f ./probe-finder


