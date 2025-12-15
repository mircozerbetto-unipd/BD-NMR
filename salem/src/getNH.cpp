#include "getNH.h"

/***************/
/* Constructor */
/***************/

NHsites::NHsites()
{
	return;
}

/**************/
/* Destructor */
/**************/

NHsites::~NHsites()
{
	return;
}

/***********/
/* Methods */
/***********/

void NHsites::locateNHsites(molecule *mol)
{
	hydrogens.clear();
	nitrogens.clear();
	carbons.clear();

	int nA = mol->getNAtoms();
	atom a, a1;

	int resNum;
	std::string aNamePDB, aResidue;

	for (int i = 1; i <= nA; ++i)
	{
		a = mol->getAtomOld(i);

		aNamePDB = a.getAtomTypePDB();
		aNamePDB.erase(std::remove_if(aNamePDB.begin(), aNamePDB.end(), isspace), aNamePDB.end());

		aResidue = a.getResidueName();
		aResidue.erase(std::remove_if(aResidue.begin(), aResidue.end(), isspace), aResidue.end());

		if (!aNamePDB.compare("HN") && a.getResidueNumber() > 1)
			hydrogens.push_back(a.id0);
		else if (!aNamePDB.compare("N") && aResidue.compare("PRO") && a.getResidueNumber() > 1)
			nitrogens.push_back(a.id0);
		else if (!aNamePDB.compare("C"))
		{
			resNum = a.getResidueNumber();
			if (resNum < mol->getNResidues())
			{
				for (int j = i + 1; j <= nA; ++j)
				{
					a1 = mol->getAtomOld(j);
					aResidue = a1.getResidueName();
					if (resNum != a1.getResidueNumber() && aResidue.compare("PRO"))
					{
						carbons.push_back(a.id0);
						break;
					}
				}
			}
		}
	}

	nSites = hydrogens.size();

	if (nitrogens.size() != nSites)
	{
		std::cout << "ERROR in file " << __FILE__ << ", line " << __LINE__ << ": number of N atoms (" << nitrogens.size() << ") is different from number of H atoms (" << nSites << ")." << std::endl;
		exit(1);
	}
	else if (carbons.size() != nSites)
	{
		std::cout << "ERROR in file " << __FILE__ << ", line " << __LINE__ << ": number of C atoms (" << carbons.size() << ") is different from number of H atoms (" << nSites << ")." << std::endl;
		exit(1);
	}
	else
	{
		std::cout << "* Found " << nSites << " N-H sites in PDB file:" << std::endl;
		for (int i = 0; i < nSites; ++i)
			std::cout << "  - Site " << i + 1 << " H(" << hydrogens[i] << "), N(" << nitrogens[i] << "), C(" << carbons[i] << ")" << std::endl;
	}
	return;
}

int NHsites::getNumberOfSites(void)
{
	return nSites;
}

// NOTE:
// i is the site number, 0-based;
// adID will be the 1-based ID in the **OLD** numeration (i.e., same as input PDB)
void NHsites::getAtomsOfSite(int i, int *atID)
{
	atID[0] = hydrogens[i];
	atID[1] = nitrogens[i];
	atID[2] = carbons[i];
	return;
} 

// NOTE: 1-based indexes in **OLD** numeration (i.e., same as input PDB)
vectorOfIntegers NHsites::getHydrogens(void)
{
	return hydrogens;
}

vectorOfIntegers NHsites::getNitrogens(void)
{
	return nitrogens;
}

vectorOfIntegers NHsites::getCarbons(void)
{
	return carbons;
}

