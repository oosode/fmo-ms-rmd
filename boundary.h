#ifndef FMR_BOUNDARY_H
#define FMR_BOUNDARY_H

#include "pointers.h"

namespace FMR_NS {
    
#define COORD_CARTESIAN    0
#define COORD_SPHERICAL    1
#define COORD_CYLINDRICAL  2
#define COORD_PT           3

class Boundary : protected Pointers {
 public:

   // ** Constructor/destructor ** //
   Boundary(FMR *);
   ~Boundary();

   // ** Variables ** //
   int    do_boundary_conditions; // Do boundary conditions
   int    bound_coord;        // Boundary potential coordinate index
   double di[3];              // Are potential on x, y, z - direction
   double k[3];               // force constants
   double radius[3];             // Boundary radius in x, y, z - direction
    
   double center[3];
   double f[2][3],dx[3],dx2[3];
   double energy, ff[3];
   double virial[6];
    
   double diff;             // displacement along umb vector
   double diff2;

   // ** Functions ** //
   void setup();
   void compute();
   void decompose_force(double *);
   void decompose_energy(double);

};

}

#endif
