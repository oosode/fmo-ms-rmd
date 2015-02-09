/* AWGL */
#ifndef FMR_ELECTROSTATICS_H
#define FMR_ELECTROSTATICS_H

#include "pointers.h"

namespace FMR_NS {


class Electrostatics : protected Pointers {
 public:

   // ** Constructor/destructor ** //
   Electrostatics(FMR *);
   ~Electrostatics();

   // ** Variables ** //
//   double r_cec[3];
//   double **r_coc;
//   double *qsum_coc;
//   int natom_coc[MAX_STATES];


   // ** Functions ** //
   void coulomb();
   void ewald();

};

}

#endif
