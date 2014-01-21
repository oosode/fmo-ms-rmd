
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

Also, this code only works with a development version of Q-Chem.

