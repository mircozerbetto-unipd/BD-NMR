===================================================

 NMR - calculate NMR data from spectral densities

===================================================

[A] USAGE

    nmr <data file> <spectral density file>

[B] DATA FILE

    Input data is introduced with the following keywords:

    nucleus    X           [ X  (string) = 15N1H, 13C1H, 13C1H2, 13C1H3, 2H                ]
    bondLength r           [ r  (floating point) is the X-H bond length in Angstroms       ]
    deltaCsa   d           [ d  (floating point) is the deltaCSA for nucleus X in ppm      ]
    frequency  f1 f2 ...   [ fn (floating point) is the n-th spectrometer frequency in MHz ]
    calculate  P1 P2 ...   [ Pn (string) = T1, T2, NOE                                     ]
    
    NOTES:

    - Only 15N1H and 13C1H are implementend at present.
    - Order of keywords is not relevant.
    - Input is not case sensitive.
    - Lines starting with a pound (#) or a space are considered as comments.
    - For CHn groups the C-H bond length is assumed to be the same for all the n bonds.
    - Keywords and values must be separated by one space.

[C] SPECTRAL DENSITIES FILE

    The file with spectral densities is expected to contain N+1 columns. The first column contains frequencies in Hz.
    The successive M <= N columns contain auto-relaxation spectral densities (in s). Following P = min{M*(M-1)/2 , N}
    columns contain cross-relaxation spectral densities (in s).

    Based on the probe, the columns are interpreted as:

    1. 15N1H

       Auto-relaxation  (M = 2) : Dip-Dip, Csa-Csa
       Cross-relaxation (P = 1) : Dip-Csa

    2. 13C1H(n)

       Auto-relaxation (M = n+1) : Dip1-Dip1, ..., Dipn-Dipn, Csa-Csa
       Cross-relaxation          : Dip1-Dip2, ..., Dip1-Dipn, Dip1-Csa, Dip2-Dip3, Dip2-Csa, Dip3-Csa

    3. 2H

       Auto-relaxation (M = 1) : Quad-Quad
       Cross-relaxation        : -

