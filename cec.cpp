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
  //next_pivot_state = 0;
  //cut_OH   = 2.5; // Angstroms
  //cut_OH   = 2.4; // Angstroms
  //max_hops = 2;
  //flag_read_MOs = 0; // default to not okay to read
  //flag_state_number_change = 0;
}

/*-----------------------------------------------------------------
  Destructor
-----------------------------------------------------------------*/
Cec::~Cec()
{

}

/*-----------------------------------------------------------------
  State search algorithm 
-----------------------------------------------------------------*/
void Cec::compute_coc()
{
  
  Atom *atom	= fmr->atom;;
  int natoms     = atom->natoms;
  int nstates    = atom->nstates;
  
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
        qsum_coc[iatom] += atom->getCharge(iatom,istate);      

      }
    }

    //natom_coc[istate] = atoms_per_ireactive;

    for (int iatom=0; iatom<natoms; ++iatom) {
      if (atom->reactive[istate*natoms + iatom]) {

        double rr[3];
        double dr[3];

        rr[0] = atom->coord[3*iatom]; 
        rr[1] = atom->coord[3*iatom+1]; 
        rr[2] = atom->coord[3*iatom+2];

        VECTOR_SUB(dr,rr,ref);
        //VECTOR_PBC(dr);
        VECTOR_SCALE(dr,atom->getCharge(iatom,istate)/qsum_coc[istate]);
        VECTOR_ADD(r_coc[iatom],r_coc[iatom],dr);

      }
    }
    VECTOR_ADD(r_coc[istate],r_coc[istate],ref);

  }
}

void Cec::compute()
{

  double *C2     = fmr->matrix->GSCoeffs;
  int natoms     = atom->natoms;
  int nstates    = atom->nstates;

  if (fmr->master_rank) {

    r_cec[0]=r_cec[1]=r_cec[2]=0.0;
    
    double ref[3];
    ref[0] = r_coc[0][0];
    ref[1] = r_coc[0][1];
    ref[2] = r_coc[0][2];

    for (int istate=1; istate<nstates; istate++) {

      double dr[3];
      VECTOR_SUB(dr,r_coc[istate],ref);
      //VECTOR_PBC(dr);
      VECTOR_SCALE_ADD(r_cec,dr,C2[istate]);

    }
    VECTOR_ADD(r_cec,r_cec,ref);

    printf("CEC position: %f %f %f\n",r_cec[0],r_cec[1],r_cec[2]);

  }
}


