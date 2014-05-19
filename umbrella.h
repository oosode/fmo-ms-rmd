/* AWGL */
#ifndef FMR_UMBRELLA_H
#define FMR_UMBRELLA_H

#include "pointers.h"

namespace FMR_NS {

#define GTP_COORD  0
#define GTP_CEC    1
#define GTP_ATOM   2
#define GTP_CM     3
#define GTP_CG     4
#define GTP_CECV2  5

#define COORD_CARTESIAN    0
#define COORD_SPHERICAL    1
#define COORD_CYLINDRICAL  2
#define COORD_PT           3

class Umbrella : protected Pointers {
 public:

   // ** Constructor/destructor ** //
   Umbrella(FMR *);
   ~Umbrella();

   // ** Variables ** //
   double *deltaE;
   double **koppa;    // Eq. 24 double summation
   double **upsilon;  // Eq. 25 double summation
   double **beta;     // Eq. 21 triple summation 
   double **omega;    // Eq. 22 triple summation

   /* Potential setting */
   int    do_umbrella_sampling;
   char   umb_typ[256];     // Umbrella potential coordinate name
   int    umb_coord;        // Umbrella potential coordinate index
   int    di[3];            // Are potential on x, y, z - direction
   double k[3];             // force constants
   double ref[3];           // Reference vector

   double center[3],f[2][3],dx[3],dx2[3];
   double energy, ff[3];
   double virial[6];

   double diff;             // displacement along umb vector

   /* output */
   char   umbFile[256];
     
   // ** Functions ** //
   void setup();
   void compute();
   void decompose_force(double *);
   void decompose_energy(double);
   void partial_C_N2(double *);
   void partial_C_N3(double *);
   void partial_C_N4(double *);

   void writeStepUmb(int);
};

}

#endif
