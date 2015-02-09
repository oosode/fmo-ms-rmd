/* AWGL */
#ifndef FMR_ATOM_H
#define FMR_ATOM_H

#include "pointers.h"

#define MAX_STATES 30

namespace FMR_NS {

class Atom : protected Pointers {
 public:

   // ** Constructor/destructor ** //
   Atom(FMR *);
   ~Atom();

   // ** Variables ** //
   int natoms;           // Number of atoms (static)
   int nfragments;       // Number of fragments (static)
   int nstates;          // Number of fragmentation states (dynamic)
   int prev_nstates;     // Number of fragmentation states in previous step 
   int ireactive;        // Index of initial reactive fragment
   double *coord;        // Cartesian coordinates of atoms in Angstroms
   double *force;        // Forces on atoms (IMPORTANT: actually the gradient)
   double *veloc;        // Velocities of atoms
   double *mass;         // Mass of atoms
   double totalMass;     // Sum of all atom masses
   char *symbol;         // Atomic symbol of atoms
   int *fragment;        // Fragment index of atom for each fragmentation state
   int *reactive;        // Index specifying if atom belongs to a reactive fragment in each state
   int *available;       // Availability of atom in state search
   int *hop;             // How many hops this atom is from the pivot state reactive fragment
   int *environment;     // Do I belong to an environment fragment? 1=yes, 0=no

   // periodic boundary conditions      
   double cellA;
   double cellB;
   double cellC; 
   double xprd;
   double yprd;
   double zprd;
   double xprd_half;
   double yprd_half;
   double zprd_half;
   int na;
   int nb;
   int nc;
   int afield;
   int bfield;
   int cfield;

   // MM charge parameters
   double qO_SPCE;        // SPC/E charge oxygen 
   double qH_SPCE;        // SPC/E charge hydrogen
   double qO_hydronium;   // Hydronium charge oxygen
   double qH_hydronium;   // Hydronium charge hydrogen
   double qCl;            // Chloiride charge

   // ** Functions ** //
   bool AtomInFragment(int iatom, int ifrag, int istate){ 
     return fragment[istate*natoms + iatom] == ifrag;
   }
   bool AtomInFragment(int iatom, int ifrag, int istate, int icella, int icellb, int icellc) {

     int ra = 2*na+1;
     int rb = 2*nb+1;
     int rc = 2*nc+1;

     // seems unnecessary
     int statesize = ra*rb*rc*natoms;

     int ar =   (iatom % statesize) / (rb*rc*natoms); 			       	ar = ar - na; 
     int br =  ((iatom % statesize) % (rb*rc*natoms)) / (rc*natoms); 		br = br - nb;
     int cr = (((iatom % statesize) % (rb*rc*natoms)) % (rc*natoms)) / natoms; 	cr = cr - nc;
     //printf("%d %d %d %d %d %d %d %d\n",ar,br,cr,icella,icellb,icellc,iatom,natoms);
     //if (ar> na) exit(0);

     if (fragment[istate*natoms + iatom%natoms]==ifrag && ar==icella && br==icellb && cr==icellc) return true;

     return false;
   }
   bool AtomInFragment(int iatom, int ifrag, int istate, int icella, int icellb, int icellc, int afield, int bfield, int cfield) {

     int ra = 2*afield + 1;
     int rb = 2*bfield + 1;
     int rc = 2*cfield + 1;

     //seems unnecessary
     int statesize = ra*rb*rc*natoms;

     int ar =   (iatom % statesize) / (rb*rc*natoms); 				ar = ar - afield;
     int br =  ((iatom % statesize) % (rb*rc*natoms)) / (rc*natoms); 		br = br - bfield;
     int cr = (((iatom % statesize) % (rb*rc*natoms)) % (rc*natoms)) / natoms; 	cr = cr - cfield;
     //
     if (fragment[istate*natoms + iatom%natoms]==ifrag && ar==icella && br==icellb && cr==icellc) return true;
     
     return false;
   }
   bool AtomInCell(int iatom, int istate, int icella, int icellb, int icellc) {

     int ra = 2*na+1;
     int rb = 2*nb+1;
     int rc = 2*nc+1;

     int ar =   iatom / (rb*rc*natoms); ar = ar - na;
     int br =  (iatom % (rb*rc*natoms)) / (rc*natoms); br = br - nb;
     int cr = ((iatom % (rb*rc*natoms)) % (rc*natoms)) / natoms; cr = cr - nc;
     
     if (ar==icella && br==icellb && cr==icellc) return true;

     return false;
   }
   bool AtomInCell(int iatom, int istate, int icella, int icellb, int icellc, int afield, int bfield, int cfield) {

     int ra = 2*afield+1;
     int rb = 2*bfield+1;
     int rc = 2*cfield+1;

     int ar =   iatom / (rb*rc*natoms); 				ar = ar - afield;
     int br =  (iatom % (rb*rc*natoms)) / (rc*natoms); 			br = br - bfield;
     int cr = ((iatom % (rb*rc*natoms)) % (rc*natoms)) / natoms; 	cr = cr - cfield;

     if (ar==icella && br==icellb && cr==icellc) return true;

     return false;
   }
   int getAtomPosition(int istate, int icella, int icellb, int icellc, int iatom) {

     int ra = 2*na+1;
     int rb = 2*nb+1;
     int rc = 2*nc+1;

     int ar = icella + na;
     int br = icellb + nb;
     int cr = icellc + nc;

     int pos = istate*ra*rb*rc*natoms + ar*rb*rc*natoms + br*rc*natoms + cr*natoms + iatom;
	
     return pos;	
   }
   int getMonomerIndex(int istate, int icella, int icellb, int icellc, int ifrag) {

     int ra = 2*na+1;
     int rb = 2*nb+1;
     int rc = 2*nc+1;

     int ar = icella + na;
     int br = icellb + nb;
     int cr = icellc + nc;

     int imon=istate*ra*rb*rc*nfragments + ar*rb*rc*nfragments + br*rc*nfragments + cr*nfragments + ifrag;
     //printf("imon:%5d a:%2d b:%2d c:%2d nfrag:%d state:%d\n", imon, icella, icellb, icellc, ifrag, istate);     
     return istate*ra*rb*rc*nfragments + ar*rb*rc*nfragments + br*rc*nfragments + cr*nfragments + ifrag;
   }
   int getDimerIndex(int istate, int icella, int icellb, int icellc, int ifrag, int jfrag) {

     int ra = 2*na+1;
     int rb = 2*nb+1;
     int rc = 2*nc+1;

     int nf2 = nfragments * nfragments;

     int ar = icella + na;
     int br = icellb + nb;
     int cr = icellc + nc;

     return istate*ra*rb*rc*nf2 + ar*rb*rc*nf2 + br*rc*nf2 + cr*nf2 + ifrag*nfragments + jfrag;

   }
   void getMonomerIndices(int imon, int &istate, int &icella, int &icellb, int &icellc, int &ifrag) {

     int ra = 2*na+1;
     int rb = 2*nb+1;
     int rc = 2*nc+1;
     //printf("na:%3d\n",na);

     istate =     imon / (ra*rb*rc*nfragments);
     icella =    (imon % (ra*rb*rc*nfragments)) / (rb*rc*nfragments) - na;
     icellb =   ((imon % (ra*rb*rc*nfragments)) % (rb*rc*nfragments)) / (rc*nfragments) - nb;
     icellc =  (((imon % (ra*rb*rc*nfragments)) % (rb*rc*nfragments)) % (rc*nfragments)) / (nfragments) - nc;
     ifrag  =  (((imon % (ra*rb*rc*nfragments)) % (rb*rc*nfragments)) % (rc*nfragments)) % (nfragments);

     //printf("imon:%5d a:%2d b:%2d c:%2d nfrag:%d\n",imon,icella,icellb,icellc,ifrag); 
     return;
   }
   void getDimerIndices(int idim, int &istate, int &icella, int &icellb, int &icellc, int &ifrag, int &jfrag) {

     int ra = 2*na+1;
     int rb = 2*nb+1;
     int rc = 2*nc+1;

     int nf2 = nfragments * nfragments;
     istate =      idim / (ra*rb*rc*nf2);
     icella =     (idim % (ra*rb*rc*nf2)) / (rb*rc*nf2) - na;
     icellb =    ((idim % (ra*rb*rc*nf2)) % (rb*rc*nf2)) / (rc*nf2) - nb;
     icellc =   (((idim % (ra*rb*rc*nf2)) % (rb*rc*nf2)) % (rc*nf2)) / (nf2) - nc;
     ifrag  =  ((((idim % (ra*rb*rc*nf2)) % (rb*rc*nf2)) % (rc*nf2)) % (nf2)) / nfragments;
     jfrag  = (((((idim % (ra*rb*rc*nf2)) % (rb*rc*nf2)) % (rc*nf2)) % (nf2)) % nfragments);

     //printf("dimer indices\n");
     return;
   }

   double getCharge(int iatom, int istate) {
     double mmq = 0.0;
     if (symbol[iatom] == 'H') {
       if (reactive[istate*natoms + iatom]) mmq = qH_hydronium; 
       else                                 mmq = qH_SPCE;
     } else if (symbol[iatom] == 'O') {
       if (reactive[istate*natoms + iatom]) mmq = qO_hydronium; 
       else                                 mmq = qO_SPCE;
     } else if (symbol[iatom] == 'L') {
       mmq = qCl; 
     } else {
       printf("Charge not found for atom %d of state %d\n", iatom, istate);
     }

     return mmq;
   }

   double getAtomMass(int iatom) {
     // Return atom mass in AMU. Must convert later to a.u. 
     double mass = 0.0;
     if      (symbol[iatom] == 'H') mass = 1.00783;
     else if (symbol[iatom] == 'O') mass = 15.99491;
     else if (symbol[iatom] == 'L') mass = 35.4527;
     return mass;     
   }

   double getAtomicNumber(int iatom) {
     // Return atomic number.
     double number = 0.0;
     if      (symbol[iatom] == 'H') number = 1.0;
     else if (symbol[iatom] == 'O') number = 8.0;
     else if (symbol[iatom] == 'L') number = 17.0;
     return number;
   }

   void setGlobalBox();
   void minimum_image(double &, double &, double &);
   void setAtomMasses(); // in atom.cpp

};

}

#endif
