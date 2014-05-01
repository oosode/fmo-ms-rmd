/* AWGL */
#ifndef FMR_CEC_H
#define FMR_CEC_H

#include "pointers.h"

namespace FMR_NS {

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
     
   //int next_pivot_state; // The state index of the next step's pivot state
   //int max_hops;         // Maximum number of hops in search
   //double cut_OH;        // Distance cutoff in state search between O and H atoms
   //int flag_read_MOs;    // Flag to indicate if it is safe to read MO coefficients from file from previous step
   //int flag_state_number_change; // Flag to indicate a change in the number of states b/w steps

   // ** Functions ** //
   void compute_coc();
   void compute();
   void decompose_force(double *);
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
