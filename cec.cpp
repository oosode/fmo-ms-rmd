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
  
  Atom *atom;

  /********************************************/
  /*** Assign COC information  ****************/
  /********************************************/
  for (int i=0; i<nstates; ++i) {
  for (int i=0; i<natoms; ++i) {
    if (atom->reactive[istate*natoms + i]) {
      chgfrag = atom->fragment[istate*natoms + i];
      break;
    }
  }

  natom_coc[i] = 
  qsum_coc[i]  = 0.0;
  
  // Get atom ID and charge
  for(int j=0; j<natom_coc[i]; j++)
  {
    id_coc[i][j] = molecule_map[mol_id][index[j]];
    qi_coc[i][j] = fabs(q[id_coc[i][j]]);
    qsum_coc[i]+=qi_coc[i][j];
  }

  // Reset the qi
  for(int j=0; j<natom_coc[i]; j++) qi_coc[i][j]/=qsum_coc[i];

  //Calculate COC
  for(int j=1; j<natom_coc[i]; j++)
  {
      double dr[3];
      VECTOR_SUB(dr,x[id_coc[i][j]],ref);
      VECTOR_PBC(dr);
      VECTOR_SCALE(dr,qi_coc[i][j]);
      VECTOR_ADD(r_coc[i],r_coc[i],dr);
  }
  VECTOR_ADD(r_coc[i],r_coc[i],ref); 
}


