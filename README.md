------------------------------------------------------------------
                         FMO-MS-RMD
                   Author: Adrian W. Lange
				  Modified by Olaseni Sode
				  
------------------------------------------------------------------

Fragment Molecular Orbital Multi-state Reactive Molecular Dynamics


This is my interface to Q-Chem for running FMO molecular dynamics.
However, since traditional FMO is not reactive, I added a cool
hook to the theory that uses multiple FMO fragmentations and a
model Hamiltonian, which is diagonalized to compute the energy
and forces via the Hellman-Feynman theorem. This allows chemical
reactions to take place between fragments, like proton transfers. 


This code is only set up to work for protonated water clusters 
for the moment. 

Also, this code works with a development version of Q-Chem, 
NWChem6.3 with python module, and Gamess.


Sample input script:

Run_Type        moldyn
QChem           rungms-xt
Scratch         /lustre/scratch/oosode/scr/
Gamess_version  asis
Gamess_ncores   12
Atoms_File      3w.xyz
PrintLevel      1
Max_Hops        1
vOption         rand
Temperature     10.0
nTimeSteps      1000
timestep        1.0
vseed           999
Thermostat      berendsen
tau             10.0
Correlation     mp2
Exchange        hf
Basis_Set       sto-3g
SCF_algorithm   diis_gdm
cellA           8.5
cellB           8.5
cellC           8.5
Periodic        0 0 0
Field           0 0 0

Sample atoms file:

 10 3 1
  H   0.67903311   0.73601772  -0.06096849 1
  O   1.76689739   1.83821361  -0.21398170 0
  H   1.72704940   2.75750807   0.02702567 0
  H   2.62005625   1.67860366  -0.60423304 0
  O   0.00000000   0.00000000   0.00000000 1
  H  -0.26044680  -0.26560957   0.93153659 1
  H  -0.75957202   0.15196230  -0.56290433 1
  O  -0.70307821  -0.76977643   2.33717800 2
  H  -1.07741374  -0.27497044   3.05814772 2
  H  -0.55395385  -1.65832267   2.64356545 2




