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
   double **r_coc;
   double *qsum_coc;
   int natom_coc[MAX_STATES];



   // ** Functions ** //
   void compute_coc();
   void compute_cec();

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
