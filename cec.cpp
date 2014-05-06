/* AWGL */
#include "fmr.h"
#include "fmr_math.h"
#include "atom.h"
#include "state.h"
#include "cec.h"
#include "run.h"
#include "matrix.h"
#include <vector>
#include <string>

using namespace FMR_NS;

/*-----------------------------------------------------------------
  Constructor
-----------------------------------------------------------------*/
Cec::Cec(FMR *fmr) : Pointers(fmr)
{
    // Basic initializations
    qsum_coc   = NULL;
    r_coc      = NULL;

}

/*-----------------------------------------------------------------
  Destructor
-----------------------------------------------------------------*/
Cec::~Cec()
{
    int nstates = fmr->atom->nstates;
    
    if (r_coc != NULL) {
        for (int i=0; i<nstates; ++i)
            delete [] r_coc[i];
        delete [] r_coc;
    }

    if (qsum_coc != NULL) delete [] qsum_coc;

}

/*-----------------------------------------------------------------
  State search algorithm 
-----------------------------------------------------------------*/
void Cec::compute_coc()
{
    
    Atom *atom	 = fmr->atom;
    int natoms       = atom->natoms;
    int nstates      = atom->nstates;
    int prev_nstates = atom->prev_nstates;
    
    // ** Allocate qsum_coc array ** //
    // If qsum_coc is already allocated from previous step, de-allocate
    if (qsum_coc != NULL) delete [] qsum_coc;
    
	// Allocate qsum_coc array based on number of states
    qsum_coc = new double [nstates];
    for (int i=0; i<nstates; ++i) qsum_coc[i] = 0.0;
    
    // ** Allocate r_coc array ** //
    // If r_coc is already allocated from previous step, de-allocate
	if (r_coc != NULL) {
	    for (int i=0; i<prev_nstates; ++i) {
            delete [] r_coc[i];
        }
        delete [] r_coc;
    }
    
    // Allocate r_coc array based on number of states for this step
    r_coc = new double*[nstates];
    for (int i=0; i<nstates; ++i) {
        r_coc[i] = new double[3];
        r_coc[i][0] = 0.0;
        r_coc[i][1] = 0.0;
        r_coc[i][2] = 0.0;
    }
    
    if (fmr->master_rank) {
        
        double ref[3];
        /********************************************/
        /*** Assign COC information  ****************/
        /********************************************/
        for (int istate=0; istate<nstates; ++istate) {
            int atoms_per_ireactive = 0;
            for (int iatom=0; iatom<natoms; ++iatom) {
                if (atom->reactive[istate*natoms + iatom]) {
                    
                    if (atom->symbol[iatom] == 'O') {
                        ref[0] = atom->coord[3*iatom];
                        ref[1] = atom->coord[3*iatom + 1];
                        ref[2] = atom->coord[3*iatom + 2];
                    }
                    atoms_per_ireactive++;
                    qsum_coc[istate] += fabs(atom->getCharge(iatom,istate));
                    
                }
            }
            
            //printf("reference state %d: %f %f %f\n",istate,ref[0],ref[1],ref[2]);
            //natom_coc[istate] = atoms_per_ireactive;
            for (int iatom=0; iatom<natoms; ++iatom) {
                if (atom->reactive[istate*natoms + iatom]) {
                    
                    double rr[3];
                    double dr[3];
                    
                    rr[0] = atom->coord[3*iatom];
                    rr[1] = atom->coord[3*iatom + 1];
                    rr[2] = atom->coord[3*iatom + 2];
                    
                    //printf("current position: %f %f %f\n",rr[0],rr[1],rr[2]);
                    
                    VECTOR_SUB(dr,rr,ref);
                    //VECTOR_PBC(dr);
                    
                    //printf("charge: %f/%f = %f\n",atom->getCharge(iatom,istate),qsum_coc[istate],atom->getCharge(iatom,istate)/qsum_coc[istate]);
                    VECTOR_SCALE(dr,fabs(atom->getCharge(iatom,istate))/qsum_coc[istate]);
                    VECTOR_ADD(r_coc[istate],r_coc[istate],dr);
                    //printf("coc: %f %f %f\n",r_coc[istate][0],r_coc[istate][1],r_coc[istate][2]);
                    
                }
            }
            VECTOR_ADD(r_coc[istate],r_coc[istate],ref);
            printf("COC position for state %2d: %15.10f %15.10f %15.10f\n",istate,r_coc[istate][0],r_coc[istate][1],r_coc[istate][2]);
            
        }
    }
}

void Cec::compute_cec()
{

  double *GSCoeffs = fmr->matrix->GSCoeffs;
  int natoms       = atom->natoms;
  int nstates      = atom->nstates;

  if (fmr->master_rank) {

    r_cec[0]=r_cec[1]=r_cec[2]=0.0;
    
    double ref[3];
    ref[0] = r_coc[0][0];
    ref[1] = r_coc[0][1];
    ref[2] = r_coc[0][2];

    for (int istate=1; istate<nstates; istate++) {

      double dr[3];
      double C2 = GSCoeffs[istate]*GSCoeffs[istate];

      VECTOR_SUB(dr,r_coc[istate],ref);
      //VECTOR_PBC(dr);
      VECTOR_SCALE_ADD(r_cec,dr,C2);

    }
    VECTOR_ADD(r_cec,r_cec,ref);

    printf("CEC position:              %15.10f %15.10f %15.10f\n",r_cec[0],r_cec[1],r_cec[2]);

  }
}

