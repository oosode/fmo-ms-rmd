/* AWGL */
#ifndef FMR_CEC_H
#define FMR_CEC_H

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

class Cec : protected Pointers {
 public:

   // ** Constructor/destructor ** //
   Cec(FMR *);
   ~Cec();

   // ** Variables ** //
   double r_cec[3];
   //double r_coc[MAX_STATES][3];        

   double **r_coc;
   double *qsum_coc;
   int natom_coc[MAX_STATES];
   double *deltaE;
   double **koppa;    // Eq. 24 double summation
   double **upsilon;  // Eq. 25 double summation
   double **beta;     // Eq. 21 triple summation 
   double **omega;    // Eq. 22 triple summation
   double **gamma;    //  

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
   void compute_coc();
   void compute_cec();
   void compute2();
   void decompose_force2(double *);
   void decompose_energy2(double);
   void partial_C_N22(double *);
   //void state_search();           // The breadth-first search for fragmentation states
   //void write_qchem_inputs(int);  // Writes the inputs for Q-Chem
   //void write_nwchem_inputs(int); // Writes the inputs for NWChem
   //void write_nwchem_inputs_cutoff(int);
   //void write_gamess_inputs(int); // Writes the inputs for Gamess
   
   //void updatePivotState();       // Update pivot state information *after* matrix diagonalization
   //void updateCoordinates();      // Update geometry coordinates of atoms in unit cell.
};

}

#define VECTOR_ZERO(a) a[0]=a[1]=a[2]=0.0
#define VECTOR_PBC(a) domain->minimum_image(a[0],a[1],a[2])
#define VECTOR_R2(b,a) b=a[0]*a[0]+a[1]*a[1]+a[2]*a[2]
#define VECTOR_R(b,a) b=sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])
#define VECTOR_SELF_ADD(b,a) b[0]+=a[0];b[1]+=a[1];b[2]+=a[2]
#define VECTOR_SELF_SUB(b,a) b[0]-=a[0];b[1]-=a[1];b[2]-=a[2]
#define VECTOR_SCALE(c,b) c[0]*=b;c[1]*=b;c[2]*=b
#define VECTOR_SUB(c,a,b) c[0]=a[0]-b[0];c[1]=a[1]-b[1];c[2]=a[2]-b[2]
#define VECTOR_SCALE_SUB(c,a,b) c[0]-=(a[0]*b);c[1]-=(a[1]*b);c[2]-=(a[2]*b)
#define VECTOR_ADD(c,a,b) c[0]=a[0]+b[0];c[1]=a[1]+b[1];c[2]=a[2]+b[2]
#define VECTOR_SCALE_ADD(c,a,b) c[0]+=(a[0]*b);c[1]+=(a[1]*b);c[2]+=(a[2]*b)
#define VECTOR_COPY(a,b) a[0]=b[0];a[1]=b[1];a[2]=b[2];

#endif
