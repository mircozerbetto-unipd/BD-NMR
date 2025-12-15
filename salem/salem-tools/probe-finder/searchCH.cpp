#include "searchCH.h"

CHprobeVec searchCH(molecule mol)
{
	std::ofstream ofs;
	std::stringstream sts;
	CHprobeVec chp;
	CHprobe p;
	int nprobes = 0;
	int nH, nX;
	vectorOfIntegers aBonds;
	vectorOfDoubles r;
	atom a1, a2, a3, a4, a0;

	int n = mol.getNAtoms();
	for (int i = 0; i < n; ++i)
	{
		a2 = mol.getAtom(i+1);
		if (!a2.getAtomType().compare("C"))
		{
			nH = 0;
			nX = 0;
			aBonds = mol.getConnectivity(a2.getID());
			for (int nb = 0; nb < aBonds.size(); ++nb)
			{
				a0 = mol.getAtom(aBonds[nb]);
				if (!a0.getAtomType().compare("H"))
				{
					nH++;
					a1 = a0;
				}
				else if (!nX)
				{
					a3 = a0;
					nX++;
				}
			}
			aBonds = mol.getConnectivity(a3.getID());
			for (int nb = 0; nb < aBonds.size(); ++nb)
			{
				if ((aBonds[nb] != a2.getID()) && a3.getAtomType().compare("H"))
				{
					a4 = mol.getAtom(aBonds[nb]);
					break;
				}
			}

			if (nH == 1)
			{
				p.Hid = a1.getID();
				p.Cid = a2.getID();
				p.Xid = a3.getID();
				p.Yid = a4.getID();
				p.T   = a2.getPosition();
				p.R   = rotomat(a1, a2, a3);
				p.resNumber = a2.getResidueNumber();
				p.resName = a2.getResidueName();
				chp.push_back(p);

				sts.str("");
				sts << "probe_CH_" << chp.size() << ".dat";
				ofs.open(sts.str().c_str(), std::ofstream::out);
				ofs << "probe CH" << std::endl;
				ofs << "refAtoms " << p.Hid << " "<< p.Cid << " " << p.Xid << " " << p.Yid <<std::endl;
				ofs << "Residue " << p.resNumber << " " << p.resName << std::endl;
				ofs << "Translation " << std::setprecision(16) << std::scientific << p.T[0] << " " << p.T[1] << "  " << p.T[2] << std::endl;
				ofs << "Rotation " << std::setprecision(16) << std::scientific << p.R[0] << " " << p.R[1] << "  " << p.R[2] << " " << p.R[3] << " " << p.R[4] << " " << p.R[5] << " " << p.R[6] <<  " " << p.R[7] << " " << p.R[8] << std::endl;
				ofs.close();
			}
		}
	}

	return chp;
}

