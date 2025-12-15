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

#include "nmr.h"

/* Initiate object */

nmr::nmr()
{
}

nmr::nmr(std::string f1, std::string f2)
{
	dataFileName = f1;
	spdFileName  = f2;
}

/* Destroy obect */

nmr::~nmr()
{
}

/* Set input files */

void nmr::setDataFile(std::string f1)
{
	dataFileName = f1;
	return;
}

void nmr::setSpdFile(std::string f2)
{
	spdFileName = f2;
}

/* Read in the files */
void nmr::readData()
{
	bool nucleusFound    = false;
	bool bondLengthFound = false;
	bool deltaCSAFound   = false;
	bool frequencyFound  = false;
	bool calculateFound  = false;

	int position, position2;
	vectorOfIntegers nmrData;

	vectorOfDoubles f;

	char fileLine[2048];
	std::string str, keyword, value, value2;

	std::fstream dataFile;

	dataFile.open(dataFileName.c_str(), std::ios::in);

	while (dataFile.getline(fileLine,2048))
	{
		if (fileLine[0] != '#' & fileLine[0] != ' ')
		{
			str.assign(fileLine);
                        for (unsigned int i=0; i<str.length(); i++) str[i] = tolower(str[i]);

                        position = str.find_first_of(" ",0);
                        keyword  = str.substr(0,position);
                        value    = str.substr(position+1);

			if (!keyword.compare("nucleus"))
			{
				setNucleus(value);
				setMagnetogyricRatio(value);
				nucleusFound = true;
			}
			else if (!keyword.compare("bondlength"))
			{
				setBondLength(atof(value.c_str()) * 1.0e-10);
				bondLengthFound = true;
			}
			else if (!keyword.compare("deltacsa"))
			{
				setDeltaCsa(atof(value.c_str()));
				deltaCSAFound = true;
			}
			else if (!keyword.compare("frequency"))
			{
				position2 = value.find_first_of(" ",0);
				while (position2 > 0)
				{
					value2 = value.substr(0,position2);
					f.push_back(atof(value2.c_str()) * 1.0e6);
					value  = value.substr(position2+1);
					position2 = value.find_first_of(" ",0);
				}
				f.push_back(atof(value.c_str()) * 1.0e6);

				frequencyFound = true;
			}
			else if (!keyword.compare("calculate"))
			{
				position2 = value.find_first_of(" ",0);
				while (position2 > 0)
				{
					value2 = value.substr(0,position2);
					if (!value2.compare("t1"))
						nmrData.push_back(T1);
					else if (!value2.compare("t2"))
						nmrData.push_back(T2);
					else if (!value2.compare("noe"))
						nmrData.push_back(NOE);
					else
					{
						std::cout << std::endl << "ERROR : NMR parameter " << value2 << " not implemented" << std::endl;
						exit(1);
					}
					value  = value.substr(position2+1);
					position2 = value.find_first_of(" ",0);
				}

				if (!value.compare("t1"))
					nmrData.push_back(T1);
				else if (!value.compare("t2"))
					nmrData.push_back(T2);
				else if (!value.compare("noe"))
					nmrData.push_back(NOE);
				else
				{
					std::cout << std::endl << "ERROR : NMR parameter " << value << " not implemented. Therefore it will not be calculated" << std::endl;
					exit(1);
				}

				calculateFound = true;
			}
		}
	}

	dataFile.close();

	/* Check input integrity */
	
	if (!nucleusFound)
	{
		std::cout << std::endl << "ERROR: nucleus not found in " << dataFileName << std::endl;
		exit(1);
	}
	if (!bondLengthFound)
	{
		std::cout << std::endl << "ERROR: bondLength not found in " << dataFileName << std::endl;
		exit(1);
	}
	if (!deltaCSAFound)
	{
		std::cout << std::endl << "ERROR: deltaCsa not found in " << dataFileName << std::endl;
		exit(1);
	}
	if (!frequencyFound)
	{
		std::cout << std::endl << "ERROR: frequency not found in " << dataFileName << std::endl;
		exit(1);
	}
	if (!calculateFound)
	{
		std::cout << std::endl << "ERROR: no NMR to calculate defined in " << dataFileName << std::endl;
		exit(1);
	}
	
	/* Conclude setup */

	setNMRparameters(nmrData);
	setFrequency(f);

#ifdef DEBUG
	std::cout << "------------------------------------------------------" << std::endl;
	std::cout << "Echoing " << dataFileName << " file content:" << std::endl << std::endl;
	std::cout << nucleus << std::endl;
	std::cout << rXH << std::endl;
	std::cout << dCSA << std::endl;
	std::cout << gyromag << std::endl;
	for (unsigned int i = 0; i < nmrData.size(); i++)
		std::cout << nmrData[i] << std::endl;
	for (unsigned int i = 0; i < f.size(); i++)
		std::cout << f[i] << std::endl;
	std::cout << "------------------------------------------------------" << std::endl;
#endif

	return;
}

void nmr::readSpd()
{
	// Reset spectral densities arrays
	clearJ();

	// Read file

	int ncol = 1;
	int position;

	char fileLine[2048];
	std::string str, str0;

	std::fstream spdFile;
	spdFile.open(spdFileName.c_str(), std::ios::in);

	// n of columns	

	spdFile.getline(fileLine,2048);
	str.assign(fileLine);
	str0.assign(fileLine);
	position = str.find_first_of(" ",0);

	while (position > 0)
	{
		ncol ++;
		str = str.substr(position+1);
		position = str.find_first_of(" ",0);
	}

	std::cout << std::endl << "Number of columns in " << spdFileName << " = " << ncol << std::endl;

	readLine(str0,ncol);

	while (spdFile.getline(fileLine,2048))
	{
		str.assign(fileLine);
		readLine(str,ncol);
	}
	spdFile.close();

#ifdef DEBUG
	std::cout << "------------------------------------------------------" << std::endl;
	std::cout << "Echoing " << spdFileName << " file content: " << std::endl;
	for (unsigned int i = 0; i < w.size(); i++)
	{
		std::cout << w[i];
		if (Jd1d1.size() > 0) std::cout << "\t" << Jd1d1[i];
		if (Jd2d2.size() > 0) std::cout << "\t" << Jd1d1[i];
		if (Jd3d3.size() > 0) std::cout << "\t" << Jd1d1[i];
		if (Jcc.size()   > 0) std::cout << "\t" << Jd1d1[i];
		if (Jqq.size()   > 0) std::cout << "\t" << Jd1d1[i];
		if (Jd1d2.size() > 0) std::cout << "\t" << Jd1d1[i];
		if (Jd1d3.size() > 0) std::cout << "\t" << Jd1d1[i];
		if (Jd1c.size()  > 0) std::cout << "\t" << Jd1d1[i];
		if (Jd2d3.size() > 0) std::cout << "\t" << Jd1d1[i];
		if (Jd2c.size()  > 0) std::cout << "\t" << Jd1d1[i];
		if (Jd3c.size()  > 0) std::cout << "\t" << Jd1d1[i];
		std::cout << std::endl;
	}
	std::cout << "------------------------------------------------------" << std::endl;
#endif
	return;
}

/* Set which NMR parameters should be calculated */

void nmr::setNMRparameters(vectorOfIntegers v)
{
	dataToCalculate.clear();
	dataToCalculate = v;
	return;
}

/* Set various physical data */

void nmr::setBondLength(double r)
{
	rXH = r;
	return;
}

void nmr::setDeltaCsa(double d)
{
	dCSA = d;
	return;
}

void nmr::setNucleus(std::string n)
{
	nucleus = n;
	return;
}

void nmr::setFrequency(vectorOfDoubles f)
{
	freq.clear();
	freq = f;
	return;
}

void nmr::setMagnetogyricRatio(std::string n)
{
	for (unsigned int i = 0; i < n.size(); i++) n[i] = tolower(n[i]);
	if (!n.compare("15n1h"))
		gyromag = -2.7116e7;
	else if (!n.compare("13c1h"))
		gyromag = 6.7283e7;
	else
	{
		std::cout << std::endl << "ERROR : probe " << n << " not implemented in nmr::setMagnetogyricRatio " << std::endl;
		exit(1);
	}
}
 
/* Calculate NMR parameters */

void nmr::calcNMR()
{

	// Some initializations

	nmrData = vectorOfDoubles(freq.size()*dataToCalculate.size(),0.0);
	vectorOfDoubles larmorFreq = vectorOfDoubles(5,0.0);
	vectorOfDoubles jj;
	spline spl = spline();

	// Cycle over fields
	for (unsigned int ifreq = 0; ifreq < freq.size(); ifreq++)
	{

		// Calculate needed frequencies

		for (unsigned int i = 0; i < 5; i++) larmorFreq[i] = 0.0;
		larmor(freq[ifreq],nucleus,&larmorFreq);

		// Retrive J(w) at the frequencies via cubic spline
		
		if (!nucleus.compare("15n1h") | !nucleus.compare("13c1h"))
		{
			if (Jd1c.size() > 0)
				jj = vectorOfDoubles(15,0.0);
			else
				jj = vectorOfDoubles(10,0.0);

			// Dip - Dip

			spl.setExperimentalPoints(w,Jd1d1);
			spl.makeSpline();
			jj[0] = spl.splint(larmorFreq[0]);
			jj[1] = spl.splint(larmorFreq[1]);
			jj[2] = spl.splint(larmorFreq[2]);
			jj[3] = spl.splint(larmorFreq[3]);
			jj[4] = spl.splint(larmorFreq[4]);

			// CSA - CSA

			spl.setExperimentalPoints(w,Jcc);
			spl.makeSpline();
			jj[5] = spl.splint(larmorFreq[0]);
			jj[6] = spl.splint(larmorFreq[1]);
			jj[7] = spl.splint(larmorFreq[2]);
			jj[8] = spl.splint(larmorFreq[3]);
			jj[9] = spl.splint(larmorFreq[4]);
	
			// Dip - CSA

			if (Jd1c.size() > 0)
			{
				spl.setExperimentalPoints(w,Jd1c);
				spl.makeSpline();
				jj[10] = spl.splint(larmorFreq[0]);
				jj[11] = spl.splint(larmorFreq[1]);
				jj[12] = spl.splint(larmorFreq[2]);
				jj[13] = spl.splint(larmorFreq[3]);
				jj[14] = spl.splint(larmorFreq[4]);
			}
			
		}
		else
		{
			std::cout << std::endl << "ERROR : probe " << nucleus << " not implemented in nmr::calcNMR - J(w) spline calculation" << std::endl << std::endl;
			exit(1);
		}

		// Calculate NMR parameters

		if (!nucleus.compare("15n1h") | !nucleus.compare("13c1h"))
		{
			double dipFactor, csaFactor;

			dipFactor =  MU0_OVER_4PI*gyroH*gyromag*hbar/(rXH*rXH*rXH);
			dipFactor *= dipFactor;
			dipFactor *= 0.1;

			csaFactor = wX*wX*dCSA*dCSA*TWO_OVER_FIFTHEEN*1.0e-12;

			for (unsigned int iNMR = 0; iNMR < dataToCalculate.size(); iNMR++)
			{
				switch (dataToCalculate[iNMR])
				{
					case (T1):
					{
						double R1 = dipFactor*(jj[3]+3.0*jj[2]+6.0*jj[4]) + csaFactor*jj[7];
						nmrData[ifreq*dataToCalculate.size()+iNMR] = 1.0 / R1;
						break;
					}
					case (T2):
					{
						double R2 = 0.5*dipFactor*(4.0*jj[0]+jj[3]+3.0*jj[2]+6.0*jj[1]+6.0*jj[4])+ONE_OVER_SIX*csaFactor*(3.0*jj[7]+4.0*jj[5]);
						nmrData[ifreq*dataToCalculate.size()+iNMR] = 1.0 / R2;
						break;
					}
					case (NOE):
					{
						double R1 = dipFactor*(jj[3]+3.0*jj[2]+6.0*jj[4]) + csaFactor*jj[7];
						nmrData[ifreq*dataToCalculate.size()+iNMR] = 1.0 + (gyroH/gyromag) * dipFactor * (6.0*jj[4] - jj[3]) / R1;
						break;
					}
					default:
					{
						std::cout << std::endl << "ERROR : NMR data code " << dataToCalculate[iNMR] << " not implemented in nmr::calcNMR" << std::endl << std::endl;
						exit(1);
					}
				}
			} // iNMR
		}
		else
		{
			std::cout << std::endl << "ERROR : probe " << nucleus << " not implemented in nmr::calcNMR - NMR data calculation" << std::endl << std::endl;
			exit(1);
		}

	} // ifreq
	return;
}

/* Output results */

std::string nmr::toString()
{
	std::ostringstream outStr;
	for (unsigned int i = 0; i < freq.size(); i++)
	{
		outStr << std::endl << "Frequency: " << freq[i]*1.0e-6 << " Mhz" << std::endl;
		for (unsigned int j = 0; j < dataToCalculate.size(); j++)
		{
			switch (dataToCalculate[j])
			{
				case (T1):
				{
					outStr << "T1 = " << nmrData[i*dataToCalculate.size()+j] << " s" << std::endl;
					break;
				}
				case (T2):
				{
					outStr << "T2 = " << nmrData[i*dataToCalculate.size()+j] << " s" << std::endl;
					break;
				}
				case (NOE):
				{
					outStr << "NOE = " << nmrData[i*dataToCalculate.size()+j] << std::endl;
					break;
				}
				default:
				{
					std::cout << std::endl << "ERROR : NMR data code " << j << " not implemented in nmr::toString" << std::endl << std::endl;
					exit(1);
				}
			} // switch
		}
	}

	return outStr.str();
}

/* Other routines */

void nmr::clearJ()
{
	w.clear();

	Jd1d1.clear();
	Jd2d2.clear();
	Jd3d3.clear();
	Jcc.clear();
	Jqq.clear();

	Jd1d2.clear();
	Jd1d3.clear();
	Jd1c.clear();
	Jd2d3.clear();
	Jd2c.clear();
	Jd3c.clear();
	
	return;
}

void nmr::readLine(std::string str, int ncol)
{

	if (!nucleus.compare("15n1h"))
	{
		double tmpW, tmpDD, tmpCC, tmpDC;
		if (ncol < 3)
		{
			std::cout << std::endl << "ERROR : at least 3 columns are required in " << spdFileName << " for 15N-1H probe" << std::endl << std::endl;
			exit(1);
		}
		else if (ncol == 3)
		{
			sscanf(str.c_str(),"%le %le %le",&tmpW,&tmpDD,&tmpCC);
			w.push_back(tmpW);
			Jd1d1.push_back(tmpDD);
			Jcc.push_back(tmpCC);
		}
		else
		{
			sscanf(str.c_str(),"%le %le %le %le",&tmpW,&tmpDD,&tmpCC,&tmpDC);
			w.push_back(tmpW);
			Jd1d1.push_back(tmpDD);
			Jcc.push_back(tmpCC);
			Jd1c.push_back(tmpDC);
		}
	} // NH

	else if (!nucleus.compare("13c1h"))
	{
		double tmpW, tmpDD, tmpCC, tmpDC;
		if (ncol < 3)
		{
			std::cout << std::endl << "ERROR : a minimum of 3 columns are required in " << spdFileName << " for 13C-1H probe" << std::endl << std::endl;
			exit(1);
		}
		else if (ncol == 3)
		{
			sscanf(str.c_str(),"%le %le %le",&tmpW,&tmpDD,&tmpCC);
			w.push_back(tmpW);
			Jd1d1.push_back(tmpDD);
			Jcc.push_back(tmpCC);
		}
		else
		{
			sscanf(str.c_str(),"%le %le %le %le",&tmpW,&tmpDD,&tmpCC,&tmpDC);
			w.push_back(tmpW);
			Jd1d1.push_back(tmpDD);
			Jcc.push_back(tmpCC);
			Jd1c.push_back(tmpDC);
		}
	} // CH
	else
	{
		std::cout << std::endl << "ERROR : probe " << nucleus << " not implemented in nmr::readLine " << std::endl << std::endl;
		exit(1);
	}

	return;
}

void nmr::larmor(double w0, std::string probe, vectorOfDoubles* larmfreq)
{
	for (unsigned int i = 0; i < probe.size(); i++) probe[i] = tolower(probe[i]);

	if (!probe.compare("15n1h") | !probe.compare("13c1h"))
	{
		wH = w0 * 2.0*M_PI;
		wX = gyromag * wH / gyroH;
		larmfreq->at(0) = 0.0;
		larmfreq->at(1) = wH;
		larmfreq->at(2) = (!probe.compare("15n1h") ? -wX : wX);
		larmfreq->at(3) = wH - wX;
		larmfreq->at(4) = wH + wX;
#ifdef DEBUG
		std::cout << wH << std::endl << wX << std::endl << wH+wX << std::endl;
#endif
	}
	else
	{
		std::cout << std::endl << "ERROR : probe " << probe << " not implemented in nmr::larmor" << std::endl << std::endl;
		exit(1);
	}

	return;
}
