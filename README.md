
_/\\\\\\\\\\\\\____/\\\\\\\\\\\\_______________/\\\\\_____/\\\__/\\\\____________/\\\\____/\\\\\\\\\_____        
_\/\\\/////////\\\_\/\\\////////\\\____________\/\\\\\\___\/\\\_\/\\\\\\________/\\\\\\__/\\\///////\\\___       
 _\/\\\_______\/\\\_\/\\\______\//\\\___________\/\\\/\\\__\/\\\_\/\\\//\\\____/\\\//\\\_\/\\\_____\/\\\___      
  _\/\\\\\\\\\\\\\\__\/\\\_______\/\\\___________\/\\\//\\\_\/\\\_\/\\\\///\\\/\\\/_\/\\\_\/\\\\\\\\\\\/____     
   _\/\\\/////////\\\_\/\\\_______\/\\\___________\/\\\\//\\\\/\\\_\/\\\__\///\\\/___\/\\\_\/\\\//////\\\____    
    _\/\\\_______\/\\\_\/\\\_______\/\\\___________\/\\\_\//\\\/\\\_\/\\\____\///_____\/\\\_\/\\\____\//\\\___   
     _\/\\\_______\/\\\_\/\\\_______/\\\____________\/\\\__\//\\\\\\_\/\\\_____________\/\\\_\/\\\_____\//\\\__  
      _\/\\\\\\\\\\\\\/__\/\\\\\\\\\\\\/_____________\/\\\___\//\\\\\_\/\\\_____________\/\\\_\/\\\______\//\\\_
       _\/////////////____\////////////_______________\///_____\/////__\///______________\///__\///________\///__


***************
* BD-NMR v1.0 *
***************

Brownian Dynamics NMR (BD-NMR) is a suite of programs developed by Mirco Zerbetto, Sergio Rampino, and Antonino Polimeno
at the Theoretical Chemistry Group of the Department of Chemical Sciences of the University of Padova (Italy).

The suite is tailored for the calculation of NMR relaxation data (T1, T2, and NOE) of 15N-1H and 13C-1H probes in
flexible molecules. The calculation is based on a harmonic modeling of the internal dynamics of the molecule and on a
Brownian dynamics simulation of the roto-conformational motions.

Therefore the basic approximations are
- the internal dynamics of the molecule is modeled as harmonic: this is called the semi-flexible body (SFB) model
- the solvent is treated implicitly as generator of fluctuation-dissipation to the roto-conformational degrees of freedom.

Thanks to these approximations, the BD simulation is 10-100 times faster than the all-atom moleculard dynamics simulations,
allowing one to calculate trajectories long enough to catch the coupled effect of internal motions and global tumbling.

BD-NMR merges the interpretative power of stochastic modeling with the simplicity of treatment of the dynamics within a 
framework that is similar to that of standard molecular dynamics simulations.

****************
* Installation *
****************

1. Edit the install.sh file and make the necessary changes to these environment variables:
   - GPU: the kind of GPU where BD-NMR is being compiled
   - CUDADIR: the location of the CUDA development kit
   - LAPACKDIR: place where to locate LAPACK
   - OPENBLASDIR: place where to locate OPENBLAS
   - GFORTRANLIB: the libgfortran.so library to be used

2. Run:
   
     sh ./install.sh
 
   from the BD-BMR_v1.0/ package directory. The script will handle the installation of all the libraries and packages
   of the suite. Also, it prepares a script to set the correct environment to run BD-NMR.

   If any problem is encountered during the installation, it is possible to manually compile the library or package,
   comment out the compilation instructions for that library/package into instal_X_packages.sh (X = cpu, gpu), and
   then run again the install.sh script.


*********
* Usage *
*********

The calculation is divided in 6 steps.

1. Energy minimization and Hessian calculation
  
   By default this is hanlded using the Tinker package. However, the user can provide a file with the Hessian calculated
   in another way. We suggest to run the default protocol to see the files that are neeeded and their format. Files
   produced at this stage have the names starting with 1_.

2. Calculation of the diffusion tensor
 
   The Hessian of the energy and the hydrodynamic information are used to calculate the diffusion tensor. All the
   files produced in this step have names starting with 2_.

3. Calculation of the Brownian trajectory
 
   At this stage, the Brownian trajectory is calculated in internal coordinates. The 3_trj_rc.dat file contains
   the trajectory in the rotational (Euler angles) and internal (z coordinates). This file can be used to check
   convergence by plotting the histograms of the Euler angles and those of the internal coordinates.
   Hisograms of alpha, cos(beta), and gamma must resemble a uniform distribution. Histograms of each of the z
   coordinates must resemble a Gaussian distribution with 0 mean and variace equal to 1.
   All the files names start as 3_.

4. Conversion of the trajectory

   The trajectory is converted into Cartesian coordinates. The output file is 4_trj.xyz.

5. Probes recognition

   Based on the type of probes that the user is studying, a script automatically recognizes the possible probes
   in the molecule. For example, if M N-H groups are found, then the script writes out files named
   5_probe_NH_j.dat, for j = 1,...,M. Each file contains information about the ID's of the atoms and roto-translation
   from the molecular frame to the frame that diagonalizes the dipolar interaction.

6. Calculation of NMR relaxation

   This is the last stage where the XYZ trajectory and the information on the probes are used to calculate the
   dipolar-dipolar and CSA-CSA autocorrelation function to be used to access the spectral densities, and then
   the T1, T2, and NOE data. The name of the files produced in this step start with 6_.

   WARNING: at the moment, the script takes the 5_probe_* files with an order that depends on the file name.
            For example, if it finds 11 probes, than the probe numbering in the finel 6_nmr.out file will be:

            Probe | File
            ------------
            1     | 5_probe_NH_10.dat
            2     | 5_probe_NH_11.dat
            3     | 5_probe_NH_1.dat
            4     | 5_probe_NH_2.dat
            5     | 5_probe_NH_3.dat
            6     | 5_probe_NH_4.dat
            7     | 5_probe_NH_5.dat
            8     | 5_probe_NH_6.dat
            9     | 5_probe_NH_7.dat
            10    | 5_probe_NH_8.dat
            11    | 5_probe_NH_9.dat

            This is an issue that will be solved in the next release of BD-NMR.

The BDNMR.sh script in the run_scripts/ subfolder allows the user to select whether to run the full workflow
described above, or parts of the protocol.

IMPORTANT
=====================================================================================================================
The full calculation using default choices in the BD-NMR suite is a good starting point. However, it is important
to always check the convergence of the trajectory before relying on the results.
=====================================================================================================================
