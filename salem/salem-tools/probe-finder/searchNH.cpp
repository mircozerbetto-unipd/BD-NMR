#include "searchNH.h"

NHprobeVec searchNH(molecule mol)
{
	std::ofstream ofs;
	std::stringstream sts;
	NHprobeVec nhp;
	NHprobe p;
	int nprobes = 0;
	int nH, nC;
	vectorOfIntegers aBonds;
	vectorOfDoubles r;
	atom a1, a2, a3, a4, a0;
	zmatline za;

	vectorOfIntegers atomsConnected, atomsConnected3;

	int n = mol.getNAtoms();
	for (int i = 0; i < n; ++i)
	{
		a2 = mol.getAtom(i+1);
		if (!a2.getAtomType().compare("N"))
		{
			nH = 0;
			nC = 0;
			atomsConnected = mol.getConnectivity(a2.getID());
			for (int nb = 0; nb < atomsConnected.size(); ++nb)
			{
				a0 = mol.getAtom(atomsConnected[nb]);
				if (!a0.getAtomType().compare("H"))
				{
					nH++;
					a1 = a0;
				}
				else if (!a0.getAtomType().compare("C") && !nC)
				{
					a3 = a0;
					atomsConnected3 = mol.getConnectivity(a3.getID());
					for (int nb3 = 0; nb3 < atomsConnected3.size(); ++nb3)
					{
						a0 = mol.getAtom(atomsConnected3[nb3]);
						if (!a0.getAtomType().compare("O"))
						{
							a4 = a0;
							nC++;
						}
					}
				}
			}
			if (nH == 1 && nC)
			{
				p.Hid = a1.getID();
				p.Nid = a2.getID();
				p.Cid = a3.getID();
				p.Oid = a4.getID();
				p.T   = a2.getPosition();
				p.R   = rotomat(a1, a2, a3);
				p.resNumber = a2.getResidueNumber();
				p.resName = a2.getResidueName();
				nhp.push_back(p);

				sts.str("");
				sts << "probe_NH_" << nhp.size() << ".dat";
				ofs.open(sts.str().c_str(), std::ofstream::out);
				ofs << "probe NH" << std::endl;
				ofs << "refAtoms " << p.Hid << " "<< p.Nid << " " << p.Cid << " " << p.Oid <<std::endl;
				ofs << "Residue " << p.resNumber << " " << p.resName << std::endl;
				ofs << "Translation " << std::setprecision(16) << std::scientific << p.T[0] << " " << p.T[1] << "  " << p.T[2] << std::endl;
				ofs << "Rotation " << std::setprecision(16) << std::scientific << p.R[0] << " " << p.R[1] << "  " << p.R[2] << " " << p.R[3] << " " << p.R[4] << " " << p.R[5] << " " << p.R[6] <<  " " << p.R[7] << " " << p.R[8] << std::endl;
				ofs.close();
			}
		}
	}

	return nhp;
}

