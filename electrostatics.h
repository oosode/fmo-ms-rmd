/* AWGL */
#ifndef FMR_ELECTROSTATICS_H
#define FMR_ELECTROSTATICS_H

#include "pointers.h"
#include <complex>

namespace FMR_NS {


class Electrostatics : protected Pointers {
 public:

   // ** Constructor/destructor ** //
   Electrostatics(FMR *);
   ~Electrostatics();

   // ** Variables ** // 
   double eps_ewald;
   double ratio_ewald;
   double alpha_ewald;
   double rcut_ewald;
//   double r_cec[3];
//   double **r_coc;
//   double *qsum_coc;
//   int natom_coc[MAX_STATES];
   int *nbox_ewald;
   int *lmax_ewald;

   complex<double> **eigax;
   complex<double> **eigay;
   complex<double> **eigaz;
   complex<double> **eigbx;
   complex<double> **eigby;
   complex<double> **eigbz;
   complex<double> **eigcx;
   complex<double> **eigcy;
   complex<double> **eigcz;   


   // ** Functions ** //
   void coulomb(); 

   void ewald();
   void setup();
   void rs_ewald();
   void fs_ewald();
   void charge_ewald();   
   void self_ewald();
 
};

}

#endif
