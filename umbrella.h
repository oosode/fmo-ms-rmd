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

#define COORD_CART      0
#define COORD_SPHERICAL 1
#define COORD_CYLINDER  2
#define COORD_PT        3

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
   int umb_typ;             // Umbrella-potential coordinate
   int di[3];               // Are potential on x, y, z - direction
   double k[3];             // force constants
   double ref[3];           // Reference distance

   double center[3],f[2][3],dx[3],dx2[3];
   double energy, ff[3];
   double virial[6];

  /* output */
  int next_out,freq_out;
     
   //int next_pivot_state; // The state index of the next step's pivot state
   //int max_hops;         // Maximum number of hops in search
   //double cut_OH;        // Distance cutoff in state search between O and H atoms
   //int flag_read_MOs;    // Flag to indicate if it is safe to read MO coefficients from file from previous step
   //int flag_state_number_change; // Flag to indicate a change in the number of states b/w steps

   // ** Functions ** //
   void compute();
   void decompose_force(double *);
   void decompose_energy(double);
   void partial_C_N2(double *);

   void write_log();
};

}

#endif