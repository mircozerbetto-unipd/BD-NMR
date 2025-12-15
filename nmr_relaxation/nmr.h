/**********************************************/
/*               NMR  class header            */
/*                                            */
/* This object hndles the calculation of NMR  */
/* parameters given nuclear properties and    */
/* the spetral densities.                     */
/*                                            */
/* Author: Mirco Zerbetto                     */
/*         mirco.zerbetto@unip.it             */
/**********************************************/

#ifndef NMR_H_
#define NMR_H_

#undef DEBUG

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>

#include "spline.h"

#define gyroH 2.675198e8    /* A m N^-1 s^-1 */
#define hbar 1.054494e-34   /* N m  s        */
#define MU0_OVER_4PI 1.0e-7 /* N A^-2        */
#define TWO_OVER_FIFTHEEN 2.0/15.0
#define ONE_OVER_SIX 1.0/6.0

#define T1  1
#define T2  2
#define NOE 3

typedef std::vector<int> vectorOfIntegers;
typedef std::vector<double> vectorOfDoubles;

class nmr
{
	public:

		nmr();
		nmr(std::string, std::string);

		~nmr();

		void setDataFile(std::string);
		void setSpdFile(std::string);

		void readData();
		void readSpd();

		void setNMRparameters(vectorOfIntegers);
		void setBondLength(double);
		void setDeltaCsa(double);
		void setNucleus(std::string);
		void setFrequency(vectorOfDoubles);

		void larmor(double,std::string,vectorOfDoubles*);

		void calcNMR();

		std::string toString();

	private:

		vectorOfIntegers dataToCalculate;

		double rXH;
		double dCSA;
		double gyromag;
		double wH, wX;
		vectorOfDoubles freq;
		vectorOfDoubles nmrData;

		vectorOfDoubles w;
		vectorOfDoubles Jd1d1, Jd2d2, Jd3d3, Jcc, Jqq;
		vectorOfDoubles Jd1d2, Jd1d3, Jd1c;
		vectorOfDoubles Jd2d3, Jd2c;
		vectorOfDoubles Jd3c;

		std::string nucleus;
		std::string dataFileName;
		std::string spdFileName;

		void setMagnetogyricRatio(std::string);
		void clearJ();
		void readLine(std::string,int);
};

#endif
